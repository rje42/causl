##' Log-likelihood for frugal parameterization
##' 
##' @param pars parameter values
##' @param dat 
##' @inheritParams rfrugalParam
##' 
##' @export
ll_frugal <- function (pars, dat, formulas, family) {
  
  nf <- length(formulas) - 1
  if (nf != 4) stop("Function only works for original frugal parameterization")
  
  ## get formulae for copula
  forms_mod <- unlist(formulas[c(1,3)])
  
  ## tidy up the formulae
  forms <- tidy_formulas(forms_mod, kwd=kwd)
  fam_cop <- last(family)
  link <- link_setup(link, family = family[-length(family)])
  
  ## reorder variables so that discrete ones come last
  ## for discrete variables, plug in empirical probabilities
  disc <- family[-length(family)] == 5 | family[-length(family)] == 0
  LHS <- lhs(forms[-length(forms)])
  inCop <- unlist(LHS)
  trunc <- list()
  
  if (any(disc)) {
    wh_disc <- which(disc)
    
    # ninCop <- setdiff(names(dat), inCop)
    # dat <- dat[,c(inCop, ninCop)]
    
    ## tabulate discrete variables
    for (i in seq_along(wh_disc)) {
      trunc[[i]] <- tabulate(dat[[LHS[wh_disc[i]]]] + 1)
      if (sum(trunc[[i]]) == 0) stop("tabulation of values failed")
      trunc[[i]] <- trunc[[i]]/sum(trunc[[i]])
    }
    
    ## then move discrete variables to the end
    wh_cnt <- which(!disc)
    new_ord <- c(wh_cnt, wh_disc, length(forms))  # adjust for multiple copula formulae
    new_ord0 <- new_ord[-length(new_ord)]
    
    LHS <- LHS[new_ord0]
    forms <- forms[new_ord]
    family <- family[new_ord]
    link <- link[new_ord0]
  }
  
  full_form <- merge_formulas(forms)
  # wh <- full_form$wh
  # dat[full_form$formula]
  
  mm <- model.matrix(full_form$formula, data=dat)
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
  if (missing(par2)) {
    if (fam_cop == 2) {
      par2 <- 4
      message("par2 set to 4\n")
    }
    else par2 = 0
  }
  
  msk <- masks(forms, family = unlist(family[c(1,3)]), wh = full_form$wh, LHS=LHS)
  prs <- msk
  names(prs) <- c("beta", "phi")
  for (i in seq_along(forms)) {
    ### consolidate this with pars2mask()
    var <- LHS[i]
    be <- pars[[var]]$beta
    prs$beta[full_form$wh[[var]],i] <- be
    if (msk$phi_m[i] > 0) prs$phi[i] <- pars[[var]]$phi
  }
  
  ## other arguments to nll2()
  other_args2 <- list(dat=dat[, LHS, drop=FALSE], mm=mm,
                      beta = beta_start2$beta_m, phi = beta_start2$phi_m,
                      inCop = seq_along(inCop),
                      fam_cop=fam_cop, fam=family[-length(family)], par2=par2,
                      useC=useC,
                      link = link)
  
  ll(dat=dat, mm=mm, beta)
  
  ## parameters to
  maxit <- con$maxit
  conv <- FALSE
  # if (!is.null(con$start)) out2 <- list(par = start)
  # else 
  out2 <- list(par = theta_st)
  
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
  
}