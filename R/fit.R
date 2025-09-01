##' Fit multivariate copula regression model
##'
##' @param dat data frame of observations
##' @param formulas list of model formulae, for outcome, covariates, and
##' finally for the copula
##' @param family families for variables and copula as above. Should be the
##' same length as `formulas`, and in same order
##' @param link link functions for each variable
##' @param cop_pars additional parameters for copula if required
##' @param use_cpp logical: should C++ routines be used?
## @param init should linear models be used to initialize starting point?
##' @param control list of parameters to be passed to `optim`
##' @param other_pars list of other parameters to use (e.g. degrees of freedom for a t-distribution)
##' @param method currently only `"optim"` is supported
##'
##' @details `formulas` is list of three or more formulae giving predictors of
##' y-margin, z-margin(s) and interaction parameters.  Fit is by maximum
##' likelihood.
##'
##' `family` takes numeric values of the following forms:
##'
##' | val|family      |
##' |---:|:-----------|
##'   |   0|binomial    |
##'   |   1|gaussian    |
##'   |   2|t           |
##'   |   3|Gamma       |
##'   |   4|beta        |
##'   |   5|binomial    |
##'   |   6|lognormal   |
##'   |  10|categorical |
##'   |  11|ordinal     |
##'
##' `control` has the same arguments as the argument in `optim`, as well
##' as `sandwich`, a logical indicating if sandwich estimates of standard errors
##' should be computed, `newton`, a logical which controls whether Newton iterates should be
##' performed at the end, and `cop` which can edit the restricted variable name
##' for the left-hand side of formulae.
##' Useful for altering are `trace` (1 shows steps of optimization) and
##' `maxit` for the number of steps.
##'
##' `cop_pars` is currently just a numeric scalar, generally representing the
##' degrees of freedom to use in a t-copula.
##'
##' The list `other_pars` should be named with the relevant variables, and
##' each entry should be a named list containing the relevant parameters.  As
##' an example, for a t-copula, the required parameter is the degrees of
##' freedom, `df`.
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
                      control=list(), other_pars=list(), method = "optim") {

  # get control parameters for optim or use defaults
  con <- list(sandwich = TRUE, method = "BFGS", newton = FALSE, cop="cop", trace = 0, fnscale = 1, maxit = 10000L,
              abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5,
              gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1,
              lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10, start = NULL)
  matches = pmatch(names(control), names(con))
  con[matches[!is.na(matches)]] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ", paste(names(control[is.na(matches)]), sep = ", "))

  # ## get more new versions of family variables
  # if (is.atomic(family)) {
  #   fam <- list(family[-length(family)], last(family))
  # }
  # dim_fm <- lengths(fam[-length(fam)])
  # fam <- process_family(fam[-length(fam)], dims = dim_fm)
  # fam <- unlist(fam, recursive = FALSE)

  ## may need to fix this to allow more flexibility in copula
  d <- length(formulas) - 1
  # if (sum(dim_fm) != d) {
  #   stop("Should be a family variable for each formula")
  # }

  fam <- fam_chk(family[seq_len(d)], d)
  new_fams <- is.list(fam)
  fam <- insert_lev(fam, target_class="causl_family")
  fam[[length(fam)+1]] <- family[[d+1]]

  ## tidy up the formulae
  forms <- tidy_formulas(formulas, kwd=con$cop)
  fam_cop <- last(family)
  if (!new_fams) link <- link_setup(link, family = fam[-length(fam)])
  else {
    if (!missing(link)) warning("Links should not be specified separately from 'causl_family' objects. 'link' argument will be ignored")
    fam <- rmv_lev(fam, target_class="causl_family")
    link <- sapply(fam[-length(fam)], `[[`, "link")
  }
  LHS <- lhs(forms[-length(forms)])

  if (fam_cop == 0) method <- "dissmann"

  if (method == "optim") {
    ## put discrete variables at the end
    if (any(sapply(family[-length(forms)], is_discrete))) {
      tmp <- process_discrete_dens(dat = dat, family = family[-length(forms)],
                                   LHSs=LHS)
      trunc <- tmp$trunc
      forms <- forms[c(tmp$order, length(forms))]
      link <- link[tmp$order]
      LHS <- LHS[tmp$order]
      fam <- fam[c(tmp$order, length(forms))]
      family <- c(tmp$family, last(family))
      if (!new_fams) link <- link[tmp$order]
    }
    else trunc <- list()
  }

  ## get a joint formula
  full_form <- merge_formulas(forms)
  mm <- model.matrix(full_form$formula, data=dat)
  # wh <- full_form$wh

  # ## get list of interactions already parameterized
  # frm_excl <- form_excl(formulas, LHS = LHS)

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

  if (method == "optim") {
    out <- fit_optim(dat=dat, full_form=full_form, family=fam, fam_cop=fam_cop,
                     link=link, mm=mm, cop_pars=cop_pars, LHSs=LHS,
                     other_pars=other_pars, control=con, use_cpp=use_cpp)
  }
  else if (method == "dissmann") {
    stop("Dissmann method not implemented yet")
  }

  return(out)
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
