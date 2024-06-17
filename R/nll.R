##' Get univariate densities and uniform order statistics
##'
##' @param x vector of observations
##' @param eta,phi linear component and dispersion parameters
##' @param df degrees of freedom (only for t-distribution)
##' @param family numeric indicator of family
##' @param link link function
##'
##' @details `fam` follows the usual numeric pattern: 1=normal,
##' 2=t-distribution and 3=Gamma with a log-link.
##'
##' @return A list with entries being the numeric vectors `u` (the
##' quantiles of the input values) and `ld` (the log density of each
##' observation).
##'
##'
glm_dens <- function (x, eta, phi, other_pars, family=1, link) {

  ## deal with case that family is a 'causl_family'
  if (is(family, "causl_family")) {
    if (!missing(link)) warning("Using link from 'causl_family' object")
    link <- family$link

    if (!(family$name %in% names(links_list))) {
      # if (!(link %in% family$links_list)) stop(paste0(link, " is not a valid link function for ", family$name, " family"))
      stop(paste0("Family ", family$name, " is not a valid and registered family"))
    }
    if (!(link %in% links_list[[family$name]])) stop(paste0(link, " is not a valid link function for ", family$name, " family"))

    ## get link function
    mu <- link_apply(eta, link, family$name)

    pars <- list(mu=mu)
    if ("phi" %in% family$pars) pars <- c(pars, list(phi=phi))
    if ("par2" %in% family$pars) pars <- c(pars, list(par2=other_pars$par2))

    lp <- do.call(family$ddist, c(list(x=x, log=TRUE), pars))
    u <- do.call(family$pdist, c(list(x=x), pars))
  }
  else if (is.numeric(family)) {
    if (missing(link)) link <- links_list[[family_vals[family_vals$val==family,"family"]]][1]

    ## get the densities for x
    if (family == 1) {
      if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else if (link=="log") mu <- exp(eta)
      else stop("Not a valid link function for Gaussian distribution")

      lp <- dnorm(x, mu, sd=sqrt(phi), log=TRUE)
      u <- pnorm(x, mu, sd=sqrt(phi))
    }
    else if (family == 2) {
      if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else if (link=="log") mu <- exp(eta)
      else stop("Not a valid link function for t-distribution")
      df <- other_pars$par2

      lp <- dt((x - mu)/sqrt(phi), df=df, log=TRUE) - log(sqrt(phi))
      u <- pt((x - mu)/sqrt(phi), df=df)
    }
    else if (family == 3) {
      if (link=="log") mu <- exp(eta)
      else if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else stop("Not a valid link function for gamma distribution")

      lp <- dgamma(x, shape=1/phi, scale=phi*mu, log=TRUE)
      u <- pgamma(x, shape=1/phi, scale=phi*mu)
    }
    else if (family == 5) {
      if (link=="logit") mu <- expit(eta)
      else if (link=="probit") mu <- pnorm(eta)
      else stop("Not a valid link function for Bernoulli distribution")

      lp <- x*log(mu) + (1-x)*log(1-mu)
      lp[is.nan(lp)] <- 0
      u <- x
    }
    else stop("Only Gaussian, t, gamma and Bernoulli distributions are allowed")
  }

  return(list(u=u, ld=lp))
}

##' @describeIn glm_dens old name
univarDens <- function (x, eta, phi, other_pars, family=1, link) {
  deprecate_soft("0.8.8", "univarDens()", "glm_dens()")
  glm_dens(x=x, eta=eta, phi=phi, other_pars=other_pars, family=family, link=link)
}

##' Negative log-likelihood
##'
##' @param theta concatenated vector of parameters (`beta` followed by `phi`)
##' @param dat matrix of data
##' @param mm model matrix for use with `beta`
##' @param beta (sparse) matrix of regression parameters for each variable and copula
##' @param phi vector of dispersion parameters
##' @param inCop vector of integers giving variables in `dat` to be included in copula
##' @param fam_cop,family integer and integer vector for copula and distribution families respectively
##' @param link vector of link functions
##' @param par2 degrees of freedom for t-distribution
##' @param useC logical: should Rcpp functions be used?
##'
##' @details The number of columns of `beta` should be the number of columns
##' in `dat` plus the number required to parameterize the copula.  The first
##' few columns and the entries in `phi` are assumed to be in the order of
##' those in `dat`.  If the \eqn{i}th
##' family for a variable does not require a dispersion parameter then the value of
##' `phi[i]` is ignored.
##'
## @importFrom Matrix Matrix
##'
nll2 <- function(theta, dat, mm, beta, phi, inCop, fam_cop=1,
                 family, link, par2=NULL, useC=TRUE) {
  np <- sum(beta > 0)

  beta[beta > 0] <- theta[seq_len(np)]
  phi[phi > 0] <- theta[-seq_len(np)]

  -sum(ll(dat, mm=mm, beta=beta, phi=phi, inCop=inCop, fam_cop=fam_cop,
      family=family, link=link, par2=par2, useC=useC))
}


ll <- function(dat, mm, beta, phi, inCop, fam_cop=1,
                 family=rep(1,nc), link, par2=NULL, useC=TRUE,
                exclude_Z = FALSE, outcome = "y") {

  if (missing(inCop)) inCop <- seq_along(dat)

  if (missing(link)) {
    fams <- family_vals[match(family, family_vals$val),2]
    link <- sapply(fams, function(x) links_list[[x]][1])
  }

  if (any(phi < 0)) return(-Inf)

  nv <- length(phi)
  nc <- length(inCop)
  if (length(family) != nc) stop(paste0("family should have length ", nc))
  else if (nv != nc) stop(paste0("phi should have length ", nv))

  if (fam_cop == 2 && any(par2 <= 0)) stop("par2 must be positive for t-copula")

  ## number of discrete variables
  ndisc <- sum(family %in% c(5,0))
  ## compute etas for each variable
  eta <- mm %*% beta

  ## get the densities for each observation
  log_den <- dat_u <- matrix(NA, nrow(dat), nc)
  # if (length(family) < nc) family <- rep_len(family, nc)

  ## get univariate densities
  for (i in which(family != 5)) {
    tmp <- univarDens(dat[,i], eta[,i], phi=phi[i], family=family[i])
    log_den[,i] <- tmp$ld
    dat_u[,i] <- pmax(pmin(tmp$u,1-1e-10),1e-10)
  }
  # wh_trunc = 0
  ## deal with discrete variables separately
  for (i in which(family == 5)) {
    # wh_trunc <- wh_trunc + 1
    tmp <- univarDens(dat[,i], eta[,i], family=family[i])
    log_den[,i] <- tmp$ld
    # log_den[,i] <- 0 #### CHANGED HERE XI
    dat_u[,i] <- tmp$u
  }

  ## now set up copula parameters
  ncv <- length(inCop)
  par <- vector(mode="list", length=choose(ncv,2))

  for (i in seq_len(choose(ncv,2))) {
    if (fam_cop <= 2 || fam_cop == 11) {
      par[[i]] <- pmin(pmax(2*expit(eta[,i+nv])-1, -1+1e-10), 1-1e-10)
    }
    else if (fam_cop == 3) {
      par[[i]] <- exp(eta[,i+nv])
    }
    else if (fam_cop == 4 || fam_cop == 6) {
      par[[i]] <- exp(eta[,i+nv])+1
    }
    else if (fam_cop == 5) {
      par[[i]] <- eta[,i+nv]
    }
  }

  ## define appropriate matrices
  if (fam_cop != 11) {
    if (ncv > 2 || fam_cop <= 2) {
      Sigma <- rep(diag(ncv), length(par[[1]]))
      dim(Sigma) <- c(ncv,ncv,length(par[[1]]))
      k = 1
      for (j in seq_len(ncv)[-1]) for (i in seq_len(j-1)) {
        Sigma[i,j,] <- Sigma[j,i,] <- par[[k]]
        k <- k+1
      }

      ## deal with Gaussian and t-copulas
      if (fam_cop == 1) {
        if (any(family == 5 | family == 0)) {
          # new_ord <- order(family[inCop])
          dat_u2 <- dat_u[,inCop,drop=FALSE]#[,new_ord,drop=FALSE]
          # Sigma <- Sigma[new_ord,new_ord,,drop=FALSE]
          eta2 <- eta
          # columns of eta that correspond to discrete variables
          eta_disc <- eta[,nc-ndisc + seq_len(ndisc),drop=FALSE]

          # link functions that correspond to discrete variables
          link_disc <- link[nc-ndisc + seq_len(ndisc)]
          # which columns to convert
          eta_disc[,which(link_disc == "logit")] <- qnorm(expit(eta_disc[,which(link_disc == "logit")]))

          eta2[,nc-ndisc + seq_len(ndisc)] <- eta_disc

          # conversion from logit to probit scale
          # eta2[,(nc - ndisc+1):nc] <- qnorm(expit(eta[,(nc - ndisc+1):nc]))

          cop <- dGaussDiscCop(dat_u2, m = ndisc, Sigma=Sigma, eta=eta2[,inCop,drop=FALSE], log=TRUE, useC=useC)
        }
        else cop <- dGaussCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, log=TRUE, useC=useC)

        # cop <- dGaussCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, log=TRUE, useC=useC)
      }
      else if (fam_cop == 2) {
        #q_dat <- qt(as.matrix(dat_u), df=par2)
        cop <- dtCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, df=par2, log=TRUE)
      }
      else stop("Only Gaussian and t-copulas implemented for more than two dimensions")
    }
    else if (nc == 2) {
      ### MODIFY THIS TO USE copula PACKAGE!
      cop <- log(VineCopula::BiCopPDF(dat_u[,1], dat_u[,2], family=fam_cop, par=par[[1]], par2=par2))
    }
    else stop("should have that nc is an integer >= 2")
  }
  else {
    cop <- log(causl::dfgmCopula(dat_u[,1], dat_u[,2], alpha=par[[1]]))
  }

  if (exclude_Z == FALSE) {
    out <- cop + rowSums(log_den)

  } else{
    wh_y <- which(colnames(dat) == outcome)
    out <- cop + log_den[,wh_y]
  }

  out
}

##' Extract parameter estimates and standard errors
##'
##' @param fit output of `optim`
##' @param beta output of `initializeParams2`
##' @param merged_formula formula with all variables on RHS
##' @param kwd keyword for copula variable
##'
ests_ses <- function(fit, beta, merged_formula, kwd) {

  if (is.null(fit$par)) {
    warning("No fitted values found")
    return(fit)
  }

  ## get number of parameters and varible
  np <- sum(beta$beta_m > 0)
  nv <- sum(regexpr(kwd, colnames(beta$beta_m)) != 1L)
  if (sum(regexpr(kwd, colnames(beta$beta_m)[seq_len(nv)]) != 1L) != nv) stop("beta_m poorly formatted")

  pars <- vector(mode="list", length=nv+1)
  if (!missing(merged_formula)) names(pars) <- c(colnames(beta$beta_m)[seq_len(nv)], kwd)

  nphi <- length(beta$phi_m)
  if (nphi != nv) stop("Something has gone wrong here!")

  beta_out <- beta$beta_m
  beta_out[beta_out > 0] <- fit$par[seq_len(np)]
  for (i in seq_len(nv)) {
    pars[[i]]$beta <- beta_out[beta$beta_m[,i] > 0, i]
    if (i <= nphi && beta$phi_m[i] > 0) pars[[i]]$phi <- fit$par[np+i]
  }
  pars[[nv+1]]$beta <- beta_out[beta$beta_m[,nv+1] > 0, nv+seq_len(ncol(beta$beta_m)-nv), drop=FALSE]

  ## now deal with standard errors
  if (!is.null(fit$se)) {
    beta_out <- beta$beta_m
    beta_out[beta_out > 0] <- fit$se[seq_len(np)]

    for (i in seq_len(nv)) {
      pars[[i]]$beta_se <- beta_out[beta$beta_m[,i] > 0, i]
      if (beta$phi_m[i] > 0) pars[[i]]$phi_se <- fit$se[np+i]
    }
    pars[[nv+1]]$beta_se <- beta_out[beta$beta_m[,nv+1] > 0, nv+seq_len(ncol(beta$beta_m)-nv), drop=FALSE]
    # if (ncol(pars[[nv+1]]$beta_se) == 1) pars[[nv+1]]$beta_se <- c(pars[[nv+1]]$beta_se)
  }

  ## now deal with sandwich standard errors
  if (!is.null(fit$sandwich_se)) {
    beta_out <- beta$beta_m
    beta_out[beta_out > 0] <- fit$sandwich_se[seq_len(np)]

    for (i in seq_len(nv)) {
      pars[[i]]$beta_sandwich <- beta_out[beta$beta_m[,i] > 0, i]
      if (beta$phi_m[i] > 0) pars[[i]]$phi_sandwich <- fit$sandwich_se[np+i]
    }
    pars[[nv+1]]$beta_sandwich <- beta_out[beta$beta_m[,nv+1] > 0, nv+seq_len(ncol(beta$beta_m)-nv), drop=FALSE]
  }

  for (i in seq_len(nv)) {
    rnms <- if(attr(merged_formula$old_forms[[i]],"intercept") == 1) "(intercept)" else character(0)
    rnms <- c(rnms, attr(merged_formula$old_forms[[i]],"term.labels"))

    ## reorder tables according to original formulas
    pars[[i]]$beta <- pars[[i]]$beta[rank(merged_formula$wh[[i]])]
    names(pars[[i]]$beta) <- rnms
    if (!is.null(fit$se)) {
      pars[[i]]$beta_se <- pars[[i]]$beta_se[rank(merged_formula$wh[[i]])]
      names(pars[[i]]$beta_se) <- rnms
    }
    if (!is.null(fit$sandwich_se)) {
      pars[[i]]$beta_sandwich <- pars[[i]]$beta_sandwich[rank(merged_formula$wh[[i]])]
      names(pars[[i]]$beta_sandwich) <- rnms
    }

    ### CHECK IF THIS SHOULD BE order()!!!
  }

  {
    ## do the above for copula params
    # dimnames(pars[[nv+1]]) <- list(rownames(beta_m[beta_m[,nv+1] > 0,,drop=FALSE]),
    #                                colnames(beta_m[,nv+seq_len(ncol(beta_m)-nv)]))
  }

  fit$pars <- pars

  return(fit)
}
