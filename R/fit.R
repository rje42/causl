##' Fit multivariate copula regression model
##'
##' @param dat data frame of observations
##' @param formulas list of model formulae, for Y, for the Z variables, and
##' finally for the copula
##' @param family families for the Y and Z distributions, and the copula. Should
##' be the same length as `formulas`
##' @param link link functions for each variable
##' @param cop_pars additional parameters for copula if required
##' @param use_cpp logical: should C++ routines be used?
## @param init should linear models be used to initialize starting point?
##' @param control list of parameters to be passed to `optim`
##'
##' @details `forms` is list of three or more formulae giving
##' predictors of y-margin, z-margin(s) and interaction
##' parameters.  Fit is by maximum likelihood.
##'
##' `control` has the same arguments as the argument in `optim`, as well

##' as `sandwich`, a logical indicating if sandwich estimates of standard errors
##' should be computed, `newton`, a logical which controls whether Newton iterates should be
##' performed at the end, and `cop` which can edit the restricted variable name
##' for the left-hand side of formulae.
##' Useful for altering are `trace` (1 shows steps of optimization) and
##' `maxit` for the number of steps.

##'
##' **Warning** By default, none of the variables should be called `cop`, as
##' this is reserved for the copula.  The reserved word can be changed using
##' the argument `cop` within control.
##'
##' @return Returns a list of class `cop_fit`.
##'
##' @export
fit_causl <- function(dat, formulas=list(y~x, z~1, ~x),
                      family=rep(1,length(formulas)), link, cop_pars, use_cpp=TRUE,
                      control=list()) {

  # get control parameters for optim or use defaults
  con <- list(sandwich = TRUE, method = "BFGS", newton = FALSE, cop="cop", trace = 0, fnscale = 1, maxit = 10000L,
              abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5,
              gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1,
              lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10, start = NULL)
  matches = pmatch(names(control), names(con))
  con[matches[!is.na(matches)]] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ", paste(names(control[is.na(matches)]), sep = ", "))

  ## ensure copula keyword is not a variable name
  sandwich <- con$sandwich
  kwd <- con$cop
  if (kwd %in% names(dat)) stop(paste("Must not have a variable named '", kwd, "'", sep=""))
  newton <- con$newton
  method <- con$method
  con <- con[-c(1:4)]

  ## may need to fix this to allow more flexibility in copula
  d <- length(formulas) - 1
  if (length(family) != length(formulas)) {
    stop("Should be a family variable for each formula")
  }

  ## tidy up the formulae
  forms <- tidy_formulas(formulas, kwd=kwd)
  fam_cop <- last(family)
  link <- link_setup(link, family = family[-length(family)])
  LHS <- lhs(forms[-length(forms)])

  ## put discrete variables at the end
  if (any(family[-length(forms)] %in% c(0,5))) {
    tmp <- process_discrete_dens(dat = dat, family = family[-length(forms)],
                                 LHSs=LHS)
    trunc <- tmp$trunc
    forms <- forms[c(tmp$order, length(forms))]
    LHS <- LHS[tmp$order]
    family <- c(tmp$family, last(family))
    link <- link[tmp$order]
  }
  else trunc <- list()

  ## get a joint formula
  full_form <- merge_formulas(forms)
  mm <- model.matrix(full_form$formula, data=dat)
  # wh <- full_form$wh
  # dat[full_form$formula]

  ## handle missingness cleanly
  if (nrow(mm) < nrow(dat)) {
    nlost <- nrow(dat) - nrow(mm)
    message(paste0(nlost, " rows deleted due to missing covariates"))
    mm_vars <- attr(terms(full_form$formula), "variables")
    dat <- dat[complete.cases(with(dat, eval(mm_vars))),]
  }
  # mms = lapply(forms, model.matrix, data=dat)
  ## attach truncation values as an attribute of the model matrix
  attr(mm, "trunc") <- trunc

  ## set secondary parameter to 4 if in a t-Copula model
  if (missing(cop_pars)) {
    if (fam_cop == 2) {
      cop_pars <- 4
      message("degrees of freedom set to 4\n")
    }
    else cop_pars = 0
  }


  ## fitting code
  ## get some intitial parameter values
  beta_start2 <- initializeParams2(dat, formulas=forms, family=family, link=link,
                                   full_form=full_form, kwd=kwd)
  theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0], beta_start2$phi[beta_start2$phi_m > 0])
  # theta_st <- c(0,0,-0.4,0.3,0.5,0.5,1,1,1,1)

  ## other arguments to nll2()
  other_args2 <- list(dat=dat[, LHS, drop=FALSE], mm=mm,
                      beta = beta_start2$beta_m, phi = beta_start2$phi_m,
                      inCop = seq_along(beta_start2$phi_m),
                      fam_cop=fam_cop, fam=family[-length(family)], cop_pars=cop_pars,
                      use_cpp=use_cpp,
                      link = link)

  ## get some intitial parameter values
  beta_start2 <- initializeParams2(dat[, LHS, drop=FALSE], formulas=forms, family=family, link=link,
                                   full_form=full_form, kwd=kwd)
  theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0], beta_start2$phi[beta_start2$phi_m > 0])


  ## parameters to
  maxit <- con$maxit
  conv <- FALSE
  # if (!is.null(con$start)) out2 <- list(par = start)
  # else
    out2 <- list(par = theta_st)
  con <- con[names(con) != 'start']

  while (!conv) {
    con$maxit <- min(maxit, 5e3)
    out <- do.call(optim, c(list(fn=nll2, par=out2$par), other_args2, list(method="Nelder-Mead", control=con)))
    con$maxit <- min(max(maxit - 5e3, 1e3), maxit)
    out2 <- tryCatch(do.call(optim, c(list(fn=nll2, par=out$par), other_args2, list(method="BFGS", control=con))),
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
  it = 0

  if (newton) {
    if (con$trace > 0) {
      cat("Newton's method, iteration: ")
    }

    while (it == 0 || (max(abs(gr)) > 1e-05 && it < 100)) {
      if (con$trace > 0) printCount(it+1)
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
    if (con$trace > 0) {
      cat(" - done\n")
    }
  }

  ## if sandwich == TRUE then compute sandwich standard errors
  if (sandwich) {
    if (con$trace > 0) cat("Computing gradients for sandwich estimates...")
    other_args2a <- other_args2
    gr2 <- matrix(0, nrow=length(gr), ncol=length(gr))

    for (i in seq_len(nrow(dat))) {
      other_args2a$dat <- other_args2$dat[i,]
      # for (j in seq_along(other_args2$mms[-length(mms)+0:1])) {
      other_args2a$mm <- other_args2$mm[i,,drop=FALSE]
      # }
      tmp <- do.call(grad, c(list(nll2, x=out$par), other_args2a))
      gr2 <- gr2 + outer(tmp, tmp)
    }

    if (con$trace > 0) cat("done\n")
  }
  # ## construct output
  out$counts = c(out$counts, newton_its = it)

  out$value = do.call(nll2, c(list(theta=out$par), other_args2))
  out$ll = -out$value
  out$grad = gr
  out$FI = do.call(hessian, c(list(nll2, x=out$par), other_args2))

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
    # if (out$se[2] < 1e-3) stop("Error here")
  }
  else out$se <- NULL

  ## construct Sandwich estimates:
  if (sandwich) {
    if(!is.null(out$se)) {
      out$sandwich <- invFI %*% gr2 %*% invFI
      out$sandwich_se <- sqrt(diag(out$sandwich))
      # out$sandwich_se <- relist(out$sandwich_se, out$pars)
    }
    else {
      out$sandwich <- out$sandwich_se <- NULL
    }
    if (con$trace > 0) cat("done\n")
  }
  else out$sandwich <- out$sandwich_se <- NULL

  ## record values
  out <- ests_ses(out, beta_start2, full_form, kwd=kwd)

  # ## record parameter values
  # beta_out <- ses <- beta_start2$beta_m
  # np <- sum(beta_out > 0)
  # beta_out[beta_out > 0] <- out$par[seq_len(np)]
  # ses[ses > 0] <- out$se[seq_len(np)]
  # phi_out <- beta_start2$phi_m
  # phi_out[phi_out > 0] <- out$par[-seq_len(np)]
  # ses_phi[ses_phi > 0] <- out$se[-seq_len(np)]
  #
  # ## regroup estimates by variables
  # out$pars <- vector(length=length(formulas), mode="list")
  # for (i in seq_len(nf-1)) {
  #   if (beta_start2$phi_m[i] == 1) out$pars[[i]] <- list(beta=beta_out[beta_start2$beta_m[,i] > 0,i], phi=phi_out[i])
  #   else out$pars[[i]] <- list(beta=beta_out[beta_start2$beta_m[,i] > 0,i])
  # }
  # out$pars[[nf]] <- list(beta=beta_out[beta_start2$beta_m[,nf] > 0,nf])
  # names(out$pars) <- lhs(formulas)

  # ln <- length(unlist(out$pars))
  # out$pars <- relist(out$par[seq_len(ln)], out$pars)
  # for (i in seq_along(formulas)[-length(formulas)]) {
  #   names(out$pars[[i]]$beta) <- colnames(mms[[i]])
  # }
  # out$pars$cop <- matrix(out$par[-seq_len(ln)], nrow=ncol(mms[[length(formulas)]]))
  # rownames(out$pars$cop) <- colnames(mms[[length(formulas)]])
  # sapply(combn(LHS, 2, simplify = FALSE),
  #        function(x) paste(x,collapse=""))
  # colnames(out$pars$cop) <- sapply(combn(LHS, 2, simplify = FALSE),
  #                                  function(x) paste(x,collapse=""))
  # if (ncol(out$pars$cop) > 3) {
  #   M <- matrix(NA_character_, length(formulas)-1, length(formulas)-1)
  #   M[upper.tri(M)] <- colnames(out$pars$cop)
  #   colnames(out$pars$cop) <- t(M)[lower.tri(M)]
  # }

  # if (!is.null(out$se)) out$se <- relist(out$se, out$pars)

  out$mm = mm
  out$formulas = full_form
  class(out) = "cop_fit"

  out
}

##' @describeIn fit_causl old name
##'
##' @param par2 former name for `cop_pars` argument
##' @param sandwich logical: should sandwich standard errors be returned?
##'
##' @export
fitCausal <- function(dat, formulas=list(y~x, z~1, ~x),
                      family=rep(1,length(formulas)), link, par2,
                      sandwich=TRUE, use_cpp=TRUE, control=list()) {
  deprecate_soft("0.8.8", "fitCausal()", with="fit_causl()")

  ## put sandwich option into control list
  control = c(list(sandwich=sandwich), control)

  fit_causl(dat=dat, formulas=formulas, family=family, link=link, cop_pars=par2,
            use_cpp=use_cpp, control=control)
}


## @param sandwich logical indicating whether to print sandwich standard errors (if present)
##' @export
print.cop_fit <- function(x, sandwich = TRUE, digits=3, ...) {
  # wh = x$mms$wh
  # d = length(x$pars)
  # nv <- c(1,v,choose(v+1,2))

  formulas <- x$formulas$old_forms
  formulas <- lapply(formulas, function(x) {attributes(x) <- NULL; x})
  sandwich <- sandwich && !is.null(x$sandwich_se)

  cat("log-likelihood: ", x$ll, "\n")

  # nform <- length(x$forms)

  ## print parameters for univariate regressions
  for (i in seq_along(formulas[-1])) {
    print(formulas[[i]])

    if (sandwich) {
      tab = cbind(par=x$pars[[i]]$beta, se=x$pars[[i]]$beta_se, sandwich=x$pars[[i]]$beta_sandwich)
      colnames(tab) = c("est.", "s.e.", "sandwich") #names(x$pars[[i]]$beta)
    }
    else {
      tab = cbind(par=x$pars[[i]]$beta, se=x$pars[[i]]$beta_se)
      colnames(tab) = c("est.", "s.e.")
    }

    print(tab, digits=digits)

    if (!is.null(x$pars[[i]]$phi)) {
      if (sandwich) cat("  residual s.e.: ", signif(c(x$pars[[i]]$phi, x$pars[[i]]$phi_se, x$pars[[i]]$phi_sandwich),digits), "\n\n")
      else cat("  residual s.e.: ", signif(c(x$pars[[i]]$phi, x$pars[[i]]$phi_se),digits), "\n\n")
    }
  }

  ## print parameters for copula
  if (!is.null(x$pars$cop)) {
    cat("copula parameters:\n")
    print(formulas[[length(formulas)]])

    if (nrow(x$pars[[i+1]]$beta) == 1) {
      if (sandwich) {
        tab = cbind(est=c(x$pars[[i+1]]$beta), se=c(x$pars[[i+1]]$beta_se), sandwich=c(x$pars[[i+1]]$beta_sandwich))
        if (ncol(x$pars[[i+1]]$beta) == 1) {
          rownames(tab) <- rownames(x$pars[[i+1]]$beta)
        }
        else rownames(tab) <- colnames(x$pars[[i+1]]$beta)
        colnames(tab) <- c("est.", "s.e.", "sandwich")
      }
      else {
        tab = cbind(est=c(x$pars[[i+1]]$beta), se=c(x$pars[[i+1]]$beta_se))
        rownames(tab) <- colnames(x$pars[[i+1]]$beta)
        colnames(tab) <- c("est.", "s.e.")
      }
      print(tab, digits=digits)
    }
    else {
      for (j in seq_len(ncol(x$pars[[i+1]]$beta))) {
        if (sandwich) {
          tab = cbind(par=x$pars[[i+1]]$beta[,j], se=x$pars[[i+1]]$beta_se[,j], sandwich=x$pars[[i+1]]$beta_sandwich[,j])
          colnames(tab) = c("est.", "s.e.", "sandwich") #names(x$pars[[i]]$beta)
        }
        else {
          tab = cbind(par=x$pars[[i+1]]$beta[,j], se=x$pars[[i+1]]$beta_se[,j])
          colnames(tab) = c("est.", "s.e.")
        }
        cat(paste0(colnames(x$pars[[i+1]]$beta)[j],":\n"))
        print(tab, digits=digits)
      }
    }
  }
  # for (i in seq_len(ncol(x$pars$cop))) {
  #   cat(colnames(x$pars$cop)[i], ":\n", sep="")
  #   # if (sandwich) tab = cbind(par=x$pars$cop[,i], se=x$se$cop[,i], sandwich=c(x$sandwich_se$cop[,i]))
  #   # else tab = cbind(par=x$pars$cop[,i], se=x$se$cop[,i])
  #   #
  #   # rownames(tab) = rownames(x$pars$cop)
  #   # print(tab, digits=digits)
  #
  #   cat("\n")
  # }

  invisible(x)
}
