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
##' @param cop_pars other parameters for copula
##' @param use_cpp logical: should Rcpp functions be used?
##' @param other_pars other parameters to pass to `glm_dens`
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
                 family, link, cop_pars=NULL, use_cpp=TRUE, other_pars=list()) {
  np <- sum(beta > 0)

  beta[beta > 0] <- theta[seq_len(np)]
  phi[phi > 0] <- theta[-seq_len(np)]

  -sum(ll(dat, mm=mm, beta=beta, phi=phi, inCop=inCop, fam_cop=fam_cop,
      family=family, link=link, cop_pars=cop_pars, use_cpp=use_cpp, other_pars=other_pars))
}


ll <- function(dat, mm, beta, phi, inCop, fam_cop=1,
                 family=rep(1,nc), link, cop_pars=NULL, use_cpp=TRUE,
                exclude_Z = FALSE, other_pars=list(), outcome = "y",
               tol=c(sd=sqrt(.Machine$double.eps), quan=sqrt(.Machine$double.eps))) {

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

  if (fam_cop == 2 && any(cop_pars <= 0)) stop("degrees of freedom must be positive for t-copula")

  ## number of discrete variables
  family[family == 0] <- 5
  ndisc <- sum(family == 5)
  ## compute etas for each variable
  eta <- mm %*% beta

  ## get the densities for each observation
  log_den <- dat_u <- matrix(NA, nrow(dat), nc)
  # if (length(family) < nc) family <- rep_len(family, nc)

  ## get univariate densities
  for (i in which(family != 5)) {
    if (family[i] == 2) {
      vr_nm <- names(dat)[i]
      tmp <- glm_dens(dat[,i], eta[,i], phi=phi[i], other_pars = other_pars[[vr_nm]], family=family[i])
    }
    else tmp <- glm_dens(dat[,i], eta[,i], phi=phi[i], family=family[i])

    log_den[,i] <- tmp$ld
    dat_u[,i] <- pmax(pmin(tmp$u,1-tol[["quan"]]),tol[["quan"]])
  }
  # wh_trunc = 0
  ## deal with discrete variables separately
  for (i in which(family == 5)) {
    # wh_trunc <- wh_trunc + 1
    tmp <- glm_dens(dat[,i], eta[,i], family=family[i])
    log_den[,i] <- tmp$ld
    # log_den[,i] <- 0 #### CHANGED HERE XI
    dat_u[,i] <- tmp$u
  }

  ## now set up copula parameters
  ncv <- length(inCop)
  par <- vector(mode="list", length=choose(ncv,2))

  for (i in seq_len(choose(ncv,2))) {
    if (nrow(eta) > 1 && sd(eta[,i+nv]) < tol[["sd"]]) sel <- 1
    else sel <- TRUE

    if (fam_cop <= 2 || fam_cop == 11) {
      par[[i]] <- pmin(pmax(2*expit(eta[sel,i+nv])-1, -1+tol[["quan"]]), 1-tol[["quan"]])
    }
    else if (fam_cop == 3) {
      par[[i]] <- exp(eta[sel,i+nv])
    }
    else if (fam_cop == 4 || fam_cop == 6) {
      par[[i]] <- exp(eta[sel,i+nv])+1
    }
    else if (fam_cop == 5) {
      par[[i]] <- eta[sel,i+nv]
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

          cop <- dGaussDiscCop(dat_u2, m = ndisc, Sigma=Sigma, eta=eta2[,inCop,drop=FALSE], log=TRUE, use_cpp=use_cpp)
        }
        else cop <- dGaussCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, log=TRUE, use_cpp=use_cpp)

        # cop <- dGaussCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, log=TRUE, use_cpp=use_cpp)
      }
      else if (fam_cop == 2) {
        #q_dat <- qt(as.matrix(dat_u), df=cop_pars)
        cop <- dtCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, df=cop_pars, log=TRUE)
      }
      else stop("Only Gaussian and t-copulas implemented for more than two dimensions")
    }
    else if (nc == 2) {
      ### MODIFY THIS TO USE copula PACKAGE!
      cop <- log(VineCopula::BiCopPDF(dat_u[,1], dat_u[,2], family=fam_cop, par=par[[1]], par2=cop_pars))
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
