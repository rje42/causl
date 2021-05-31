##' Get univariate densities and uniform order statistics
##'
##' @param x vector of observations
##' @param mu,phi mean and dispersion parameters
##' @param df degrees of freedom (only for t-distribution)
##' @param fam numeric indicator of family
##'
##' @details \code{fam} follows the usual numeric pattern: 1=normal,
##' 2=t-distribution and 3=Gamma with a log-link.
##'
##' @return A list with entries being the numeric vectors \code{u} (the
##' quantiles of the input values) and \code{ld} (the log density of each
##' observation).
##'
univarDens <- function (x, mu, phi, df, fam=1) {

  ## get the densities for x
  if (fam == 1) {
    lp <- dnorm(x, mu, sd=sqrt(phi), log=TRUE)
    u <- pnorm(x, mu, sd=sqrt(phi))
  }
  else if (fam == 2) {
    lp <- dt((x - mu)/sqrt(phi), df=df, log=TRUE) - log(sqrt(phi))
    u <- pt((x - mu)/sqrt(phi), df=df)
  }
  else if (fam == 3) {
    lp <- dgamma(x, shape=1/phi, scale=phi*exp(mu), log=TRUE)
    u <- pgamma(x, shape=1/phi, scale=phi*exp(mu))
  }
  else if (fam == 5) {
    lp <- log(mu[x+1])
    u <- x
  }
  else stop("Only t, normal and gamma distributions are allowed")

  return(list(u=u, ld=lp))
}



##' Negative log-likelihood of multivariate conditional copula
##'
##' @param beta vector of parameters
##' @param dat output variables
##' @param mms list of model matrices and indices created by \code{fitCausal}
##' @param fam_cop copula family
##' @param fam univariate families for distributions
##' @param par2 second parameter for t-distributions
##' @param useC logical: should C++ routine be used?
##'
##' @details
##' \code{fam_cop} has: 1=Gaussian, 2=t, 3:6 Archimidean, 11 FGM, and
##'  \code{fam} has 1=Gaussian, 2=t, 3=Gamma (with log-link), 5=truncated normal.
##'
##' @return Numeric value of the negative log-likelihood.
##'
nll <- function(beta, dat, mms, fam_cop=1, fam=rep(1,nc), par2=NULL, useC=TRUE) {

  nc <- ncol(dat)
  if (length(mms) != nc + 3) stop(paste0("mms should have length ", nc+2))
  if (length(fam) != nc) stop(paste0("fam should have length ", nc))

  if (fam_cop == 2 && any(par2 <= 0)) stop("par2 must be positive")

  eta <- vector(mode="list", length=nc)
  eta_cor <- vector(mode="list", length=choose(nc,2))
  ## compute etas for each variable
  for (i in seq_len(nc)) {
    eta[[i]] <- mms[[i]] %*% beta[mms$wh[[2*i - 1]]]
  }
  ## now compute copula eta parameters
  for (i in seq_len(choose(nc,2))) {
    if (ncol(mms[[length(mms)-2]]) == 1) eta_cor[[i]] <- mms[[length(mms)-2]][1,] %*% beta[mms$wh[[2*nc+i]]]
    else eta_cor[[i]] = mms[[length(mms)-2]] %*% beta[mms$wh[[2*nc+i]]]
  }

  ## get the densities for each observation
  log_den <- dat_u <- matrix(NA, nrow(dat), ncol(dat))
  # if (length(fam) < nc) fam <- rep_len(fam, nc)

  for (i in which(fam != 5)) {
    tmp <- univarDens(dat[[i]], eta[[i]], phi=beta[mms$wh[[2*i]]]^2, fam=fam[i])
    log_den[,i] <- tmp$ld
    dat_u[,i] <- tmp$u
  }
  wh_trunc = 0
  for (i in which(fam == 5)) {
    wh_trunc <- wh_trunc + 1
    tmp <- univarDens(dat[[i]], mms$trunc[[wh_trunc]], fam=fam[i])
    log_den[,i] <- tmp$ld
    dat_u[,i] <- tmp$u
  }

  ## now set up copula parameters
  par <- vector(mode="list", length=choose(nc,2))

  for (i in seq_len(choose(nc,2))) {
    if (fam_cop <= 2 || fam_cop == 11) {
      par[[i]] <- pmin(pmax(2*expit(eta_cor[[i]])-1, -1+1e-10), 1-1e-10)
    }
    else if (fam_cop == 3) {
      par[[i]] <- exp(eta_cor[[i]])
    }
    else if (fam_cop == 4 || fam_cop == 6) {
      par[[i]] <- exp(eta_cor[[i]])+1
    }
    else if (fam_cop == 5) {
      par[[i]] <- eta_cor[[i]]
    }
  }

  ## define appropriate matrices
  if (fam_cop != 11) {
    if (nc > 2 || fam_cop <= 2) {
      Sigma <- rep(diag(nc), length(par[[1]]))
      dim(Sigma) <- c(nc,nc,length(par[[1]]))
      k = 1
      for (j in seq_len(nc)[-1]) for (i in seq_len(j-1)) {
        Sigma[i,j,] <- Sigma[j,i,] <- par[[k]]
        k <- k+1
      }

      ## deal with Gaussian and t-copulas
      if (fam_cop == 1) {
        if (any(fam == 5)) {
          new_ord <- c(1,order(fam)+1)
          dat_u <- dat_u[,new_ord,drop=FALSE]
          Sigma <- Sigma[new_ord,new_ord,,drop=FALSE]
          cop <- dGaussDiscCop(dat_u, Sigma=Sigma, trunc=mms$trunc, log=TRUE, useC=useC)
        }
        else cop <- dGaussCop(dat_u, Sigma=Sigma, log=TRUE, useC=useC)
      }
      else if (fam_cop == 2) {
        #q_dat <- qt(as.matrix(dat_u), df=par2)
        cop <- dtCop(dat_u, Sigma=Sigma, df=par2, log=TRUE)
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

  out <- -sum(cop) - sum(log_den)

  out
}


##' Fit multivariate copula regression model
##'
##' @param dat data frame of observations
##' @param formulas list of model formulae, for Y, for the Z variables, and
##' finally for the copula
##' @param family families for the Y and Z distributions, and the copula. Should
##' be the same length as \code{formulas}
##' @param par2 additional parameters if required
##' @param sandwich logical: should sandwich standard errors be returned?
##' @param useC logical: should C++ routines be used?
##' @param init should linear models be used to initialize starting point?
##' @param control list of parameters to be passed to \code{optim}
##'
##' @details \code{forms} is list of three or more formulae giving
##' predictors of y-margin, z-margin(s) and interaction
##' parameters.  Fit is by maximum likelihood.
##'
##' \code{control} has the same arguments as the argument in \code{optim}, as well
##' as \code{newton}, a logical which controls whether Newton iterates should be performed at the end.
##' Useful for altering are \code{trace} (1 shows steps of optimization) and
##' \code{maxit} for the number of steps, as well as \code{fact} which controls
##' the proportion of correlation put into a Gaussian or t-Copula's starting
##' values.
##'
##' @return Returns a list of class \code{cop_fit}.
##'
##' @export
fitCausal <- function(dat, formulas=list(y~x, z~1, ~x),
                      family=rep(1,length(formulas)), par2,
                      sandwich=TRUE, useC=TRUE, init=FALSE,
                      control=list()) {

  # get control parameters for optim or use defaults
  con <- list(newton = FALSE, trace = 0, fnscale = 1, maxit = 10000L, abstol = -Inf,
              reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5,
              gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1,
              lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10)
  matches = pmatch(names(control), names(con))
  con[matches[!is.na(matches)]] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ", paste(names(control[is.na(matches)]), sep = ", "))

  newton <- con$newton
  con <- con[names(con) != "newton"]

  d <- length(formulas) - 1
  if (length(family) != length(formulas)) {
    stop("Should be a family variable for each formula")
  }

  forms <- formulas
  fam_cop <- last(family)

  ## make sure all the univariate formulae have left-hand sides
  if (any(lengths(forms[-length(forms)]) < 3)) {
    wh <- which(lengths(forms[-length(forms)]) < 3)
    for (i in seq_along(wh)) {
      nm <- paste0("V", i)
      forms[[wh[i]]] <- update.formula(forms[[wh[i]]], nm ~ .)
    }
  }

  LHS <- lhs(forms[-length(forms)])
  mms = lapply(forms, model.matrix, data=dat)

  ## set secondary parameter to 4 if in a t-Copula model
  if (missing(par2)) {
    if (fam_cop == 2) {
      par2 <- 4
      message("par2 set to 4\n")
    }
    else par2 = 0
  }

  ## parameters, in order:
  ## regression for Y, standard deviation for Y,
  ## regression 1, st dev 1,
  ## regression 2, st dev 2,
  ## ...
  ## regression d, st dev d
  ## correlation.
  wh = list()
  wh[[1]] = seq_len(ncol(mms[[1]]))
  wh[[2]] = last(wh[[1]])+1

  for (i in seq_len(d-1)) {
    wh[[2*i+1]] = seq_len(ncol(mms[[i+1]]))+wh[[2*i]]
    wh[[2*i+2]] = last(wh[[2*i+1]])+1
  }
  for (i in seq_len(choose(d,2))) {
    wh[[2*d+i]] = seq_len(ncol(mms[[length(mms)]]))+last(wh[[2*d+i-1]])
  }
  mms = c(mms, list(wh=wh))

  ## for discrete variables, plug in empirical probabilities
  trunc <- list()
  wh_disc <- which(family[-length(family)] == 5)

  for (i in seq_along(wh_disc)) {
    trunc[[i]] <- tabulate(dat[[LHS[wh_disc[i]+1]]] + 1)
    if (sum(trunc[[i]]) == 0) stop("tabulation of values failed")
    trunc[[i]] <- trunc[[i]]/sum(trunc[[i]])
  }
  mms <- c(mms, list(trunc=trunc))

  beta_start <- initializeParams(dat, formulas=forms, family=family, init=init,
                                 LHS=LHS, wh=mms$wh)


  other_args2 <- list(dat=dat[, LHS, drop=FALSE], mms=mms,
                     fam_cop=fam_cop, fam=family[-length(family)], par2=par2,
                     useC=useC)

#  ll(beta_start, y=dat$y, z=dat[,grep("z", names(dat))], mms=mms, fam=fam_cop, fam_y=fam_y, fam_z=fam_z, par2=par2)
#  do.call(ll, c(list(beta=beta_start), other_args))

  # out <- do.call(optim, c(list(fn=ll, par=beta_start), other_args, list(control=con)))
  out <- do.call(optim, c(list(fn=nll, par=beta_start), other_args2, list(control=con)))
  curr_val = out$value
  if (out$convergence != 0) {
    cat("Error in convergence of optim\n")
    # return(out)
  }
  if (sandwich || newton) gr <- do.call(grad, c(list(nll, x=out$par), other_args2))
  else gr <- NULL

  ## finish with Newton's method if so required
  it = 0

  if (newton) {
    if (con$trace > 0) {
      cat("Newton's method, iteration: ")
    }

    while (it == 0 || (max(abs(gr)) > 1e-05 && it < 100)) {
      if (con$trace > 0) printCount(it+1)
      Hess <- do.call(hessian, c(list(nll, x = out$par),
                                 other_args2))
      if (rcond(Hess) > 1e-16) iHess <- solve.default(Hess)
      else iHess <- MASS::ginv(Hess)

      new = out$par - c(iHess %*% gr)
      new_val <- do.call(nll, c(list(beta = new), other_args2))
      if (is.na(new_val)) {
        cat("Newton's algorithm gives NA value, exiting\n")
        break
      }
      else if (curr_val < new_val) {
        cat("Newton's algorithm failed, exiting\n")
        break
      }
      else {
        out$par = new
        curr_val = new_val
      }
      gr = do.call(grad, c(list(nll, x = out$par), other_args2))
      it = it + 1
    }
  }

  ## if sandwich == TRUE then compute sandwich standard errors
  if (sandwich) {
    other_args2a <- other_args2
    gr2 <- matrix(0, nrow=length(gr), ncol=length(gr))

    for (i in seq_len(nrow(dat))) {
      other_args2a$dat <- other_args2$dat[i,]
      for (j in seq_along(other_args2$mms[-length(mms)+0:1])) {
        other_args2a$mms[[j]] <- other_args2$mms[[j]][i,,drop=FALSE]
      }
      tmp <- do.call(grad, c(list(nll, x=out$par), other_args2a))
      gr2 <- gr2 + outer(tmp, tmp)
    }
  }

  # ## construct output
  out$counts = c(out$counts, newton_its = it)

  out$value = do.call(nll, c(list(beta=out$par), other_args2))
  out$ll = -out$value
  out$grad = gr
  out$FI = do.call(hessian, c(list(nll, x=out$par), other_args2))

  ## get rid of standard error parameters for non-Gaussian models
  # rm <- integer(0)
  # if (fam_z != 1) {
  #   rm <- wh[[4]]
  #   out$grad <- out$grad[-rm]
  #   out$FI <- out$FI[-rm, -rm, drop=FALSE]
  # }
  # if (fam_y != 1) {
  #   rm <- wh[[2]]
  #   out$grad <- out$grad[-rm]
  #   out$FI <- out$FI[-rm, -rm, drop=FALSE]
  # }

  if (!any(is.na(out$FI)) && rcond(out$FI) > 1e-16) {
    invFI <- solve.default(out$FI)
  }
  else {
    invFI <- tryCatch(MASS::ginv(out$FI), error = function(e) NA)
  }
  if (is.numeric(invFI)) {
    out$se = sqrt(pmax(0, diag(invFI)))
    if (out$se[2] < 1e-3) stop("Error here")
  }
  else out$se <- NULL

  ## regroup estimates by variables
  out$pars <- vector(length=length(formulas), "list")
  names(out$pars) <- c(lhs(formulas[-length(formulas)]), "cop")
  # if (length(wh) != 2*length(formulas) - 1) message("problem with relisting")

  out$pars <- vector(length=length(formulas)-1, mode="list")
  for (i in seq_along(out$pars)) out$pars[[i]] <- list(beta=numeric(ncol(mms[[i]])), phi=numeric(1))
  names(out$pars) <- lhs(formulas[-length(formulas)])

  ln <- length(unlist(out$pars))
  out$pars <- relist(out$par[seq_len(ln)], out$pars)
  for (i in seq_along(formulas)[-length(formulas)]) {
    names(out$pars[[i]]$beta) <- colnames(mms[[i]])
  }
  out$pars$cop <- matrix(out$par[-seq_len(ln)], nrow=ncol(mms[[length(formulas)]]))
  rownames(out$pars$cop) <- colnames(mms[[length(formulas)]])
  sapply(combn(LHS, 2, simplify = FALSE),
         function(x) paste(x,collapse=""))
  colnames(out$pars$cop) <- sapply(combn(LHS, 2, simplify = FALSE),
                                   function(x) paste(x,collapse=""))
  if (ncol(out$pars$cop) > 3) {
    M <- matrix(NA_character_, length(formulas)-1, length(formulas)-1)
    M[upper.tri(M)] <- colnames(out$pars$cop)
    colnames(out$pars$cop) <- t(M)[lower.tri(M)]
  }

  ## construct Sandwich estimates:
  if (sandwich) {
    if(!is.null(out$se)) {
      out$sandwich <- invFI %*% gr2 %*% invFI
      out$sandwich_se <- sqrt(diag(out$sandwich))
      out$sandwich_se <- relist(out$sandwich_se, out$pars)
    }
    else {
      out$sandwich <- out$sandwich_se <- NULL
    }

  }
  else out$sandwich <- out$sandwich_se <- NULL

  if (!is.null(out$se)) out$se <- relist(out$se, out$pars)

  out$mms = mms
  out$formulas = forms
  class(out) = "cop_fit"

  out
}

## @param sandwich logical indicating whether to print sandwich standard errors (if present)
##' @export
print.cop_fit <- function(x, sandwich = TRUE, digits=3, ...) {
  # wh = x$mms$wh
  # d = length(x$pars)
  # nv <- c(1,v,choose(v+1,2))


  sandwich = sandwich && !is.null(x$sandwich_se)

  cat("log-likelihood: ", x$ll, "\n")

  # nform <- length(x$forms)

  ## print parameters for univariate regressions
  for (i in seq_along(x$formulas[-1])) {
    print(x$formulas[[i]])

    if (sandwich) tab = cbind(par=x$pars[[i]]$beta, se=x$se[[i]]$beta, sandwich=x$sandwich_se[[i]]$beta)
    else tab = cbind(par=x$pars[[i]]$beta, se=x$se[[i]]$beta)
    rownames(tab) = names(x$pars[[i]]$beta)
    print(tab, digits=digits)

    if (!is.null(x$pars[[i]]$phi)) {
      if (sandwich) cat("  residual s.e.: ", signif(c(x$pars[[i]]$phi, x$se[[i]]$phi, x$sandwich_se[[i]]$phi),digits), "\n\n")
      else cat("  residual s.e.: ", signif(c(x$pars[[i]]$phi, x$se[[i]]$phi),digits), "\n\n")
    }
  }

  ## print parameters for copula
  if (!is.null(x$pars$cop)) {
    cat("copula parameters:\n")
    print(x$formulas[[length(x$formulas)]])
  }
  for (i in seq_len(ncol(x$pars$cop))) {
    cat(colnames(x$pars$cop)[i], ":\n", sep="")
    if (sandwich) tab = cbind(par=x$pars$cop[,i], se=x$se$cop[,i], sandwich=c(x$sandwich_se$cop[,i]))
    else tab = cbind(par=x$pars$cop[,i], se=x$se$cop[,i])

    rownames(tab) = rownames(x$pars$cop)
    print(tab, digits=digits)

    cat("\n")
  }

  invisible(x)
}
