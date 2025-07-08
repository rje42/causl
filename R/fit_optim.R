##' Fit using maximum likelihood estimation
##'
##' @inheritParams fit_causl
fit_optim <- function (dat, full_form, family, fam_cop, link, mm, cop_pars, LHS,
                       other_pars, control, use_cpp) {

  ## ensure copula keyword is not a variable name
  kwd <- control$cop
  if (kwd %in% names(dat)) stop(paste("Must not have a variable named '", kwd, "'", sep=""))
  sandwich <- control$sandwich
  newton <- control$newton

  control <- control[-c(1:4)]

  ## set secondary parameter to 4 if any copula is t-Copula
  if (missing(cop_pars)) {
    if (any(fam_cop == 2)) {
      cop_pars <- 4
      message("degrees of freedom set to 4\n")
    }
    else cop_pars <- 0
  }

  ## obtain LHSs if not supplied
  if (missing(LHS)) {
    LHS <- lhs(full_form$reforms)
  }

  ## fitting code
  ## get some intitial parameter values
  beta_start2 <- initializeParams2(dat[, LHS, drop=FALSE], full_form=full_form,
                                   family=family, link=link, kwd=kwd)
  theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0], beta_start2$phi[beta_start2$phi_m > 0])

  # beta_start2 <- initializeParams2(dat, formulas=formulas, family=family, link=link,
  #                                  full_form=full_form, kwd=kwd)
  # theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0], beta_start2$phi[beta_start2$phi_m > 0])

  ## other arguments to nll2()
  other_args2 <- list(dat=dat[, LHS, drop=FALSE], mm=mm,
                      beta = beta_start2$beta_m, phi = beta_start2$phi_m,
                      inCop = seq_along(beta_start2$phi_m),
                      fam_cop=fam_cop, fam=family[-length(family)], cop_pars=cop_pars,
                      use_cpp=use_cpp, link = link, other_pars = other_pars)

  ## parameters to
  maxit <- control$maxit
  conv <- FALSE
  # if (!is.null(control$start)) out2 <- list(par = start)
  # else
  out2 <- list(par = theta_st)
  control <- control[names(control) != 'start']

  while (!conv) {
    control$maxit <- min(maxit, 5e3)
    out <- do.call(optim, c(list(fn=nll2, par=out2$par), other_args2, list(method="Nelder-Mead", control=control)))
    control$maxit <- min(max(maxit - 5e3, 1e3), maxit)
    out2 <- tryCatch(do.call(optim, c(list(fn=nll2, par=out$par), other_args2, list(method="BFGS", control=control))),
                     warning=function(e) NA, error=function(e) NA)
    if (!isTRUE(is.na(out2))) {
      out <- out2
      conv  <- TRUE
    }
    else out2 <- list(par = out$par)
  }
  curr_val = out$value
  if (out$convergence != 0) {
    cat("Error in convergence of optim\n")
    # return(out)
  }
  if (sandwich || newton) gr <- do.call(grad, c(list(nll2, x=out$par), other_args2))
  else gr <- NULL

  ## finish with Newton's method if so required
  it <- 0

  if (newton) {
    if (control$trace > 0) {
      cat("Newton's method, iteration: ")
    }

    while (it == 0 || (max(abs(gr)) > 1e-05 && it < 100)) {
      if (control$trace > 0) printCount(it+1)
      Hess <- do.call(hessian, c(list(nll2, x = out$par),
                                 other_args2))
      if (rcond(Hess) > 1e-16) iHess <- solve.default(Hess)
      else iHess <- MASS::ginv(Hess)

      new = out$par - c(iHess %*% gr)
      new_val <- do.call(nll2, c(list(theta = new), other_args2))
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
      gr = do.call(grad, c(list(nll2, x = out$par), other_args2))
      it = it + 1
    }
    if (control$trace > 0) {
      cat(" - done\n")
    }
  }

  # ## construct output
  out$counts = c(out$counts, newton_its = it)

  out$value = do.call(nll2, c(list(theta=out$par), other_args2))
  out$ll = -out$value
  out$grad = gr
  out$OI = do.call(hessian, c(list(nll2, x=out$par), other_args2))

  ## check that observed information is well conditioned
  if (!any(is.na(out$OI)) && rcond(out$OI) > 1e-16) {
    invOI <- solve.default(out$OI)
  }
  else {
    invOI <- tryCatch(MASS::ginv(out$OI), error = function(e) NA)
  }
  if (is.numeric(invOI)) {
    out$se = sqrt(pmax(0, diag(invOI)))
    # if (out$se[2] < 1e-3) stop("Error here")
  }
  else out$se <- NULL

  ## construct Sandwich estimates:
  if (sandwich) {
    ## if sandwich == TRUE then compute sandwich standard errors
    sand_method <- "simple"
    gr2 <- sandwich_errors(func=nll2, pars=out$par, method=sand_method,
                           trace=control$trace, other_args2)

    if(!is.null(out$se)) {
      out$sandwich <- invOI %*% gr2 %*% invOI
      out$sandwich_se <- sqrt(diag(out$sandwich))
      # out$sandwich_se <- relist(out$sandwich_se, out$pars)
    }
    else {
      out$sandwich <- out$sandwich_se <- NULL
    }
  }
  else out$sandwich <- out$sandwich_se <- NULL

  ## record values
  out <- ests_ses(out, beta_start2, full_form, kwd=kwd)

  ## add on additional objects used
  out$mm <- mm
  out$formulas <- full_form
  out$method <- "optim"
  class(out) <- "cop_fit"

  return(out)
}

##' Compute robust standard errors
##'
##' @param func function to be optimized
##' @param pars point at which function should be optimized
##' @param method taken from `grad()`; defaults to `"simple"` for computational reasons
##' @param trace provides more details if > 0.
##' @param other_args other arguments
##'
sandwich_errors <- function(func, pars, method="simple",
                       trace=0L, other_args) {
  if (trace > 0) cat("Computing sandwhich errors...")

  ## get matrix for results
  gr <- matrix(0, nrow=length(pars), ncol=length(pars))
  other_args2 <- other_args

  for (i in seq_len(nrow(other_args$dat))) {
    ## for each observation, compute the gradient
    other_args2$dat <- other_args$dat[i,]
    other_args2$mm <- other_args$mm[i,,drop=FALSE]

    tmp <- do.call(grad, c(list(func, x=pars), other_args2))
    gr <- gr + outer(tmp, tmp)
  }

  if (trace > 0) cat("done\n")

  return(gr)
}

