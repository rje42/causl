##' Function to generate random copula parameters for simulation
##'
##' @param formulas formulas as specified in `rfrugalParam`
##' @param data dataset to obtain parameterization for
##' @param range range of parameters to target
##'
##' @return A list suitable for the `cop` entry of the `pars` argument of
##' `rfrugalParam`
##'
##' @description Attempts to ensure that values after passing through the
##' standard link function used for Gaussian copulas will have the specified
##' value.  For other copulas this will not target the correct range, but it
##' can still be used by considering how the relevant link functions work for
##' the Gaussian and other copula.
##'
##' @export
gen_cop_pars <- function (formulas, data, range=c(-1,1)) {
  if (missing(data)) {
    stop("Data set should be supplied")
  }


  if (is(formulas[[4]], "formula")) {
    form_cop <- formulas[[4]]
  }
  else if (length(formulas[[4]]) > 1) {
    stop("At most one formula please!")
  }
  else if (length(formulas[[4]]) == 1) {
    form_cop <- formulas[[4]][[1]]
  }

  samp <- model.matrix(form_cop, data=data)
  np <- ncol(samp)

  mns <- colMeans(samp)
  mns[1] <- 0
  sds <- apply(samp, 2, sd)

  mid_pt <- (range[1]+range[2])/2
  sprd <- (range[2]-range[1])/2
  lin_mid_pt <- logit((mid_pt+1)/2)

  ## get means and standard deviations for random parameters
  par_mns <- par_sds <- numeric(np)
  par_mns[1] <- lin_mid_pt - sum(mns)
  par_sds[1] <- sprd/np
  par_sds[-1] <- pmin(sprd/(2*np), 1/sds[-1] * sprd/np)

  ## get variables to put in copula
  vars <- c(lhs(formulas[[2]]), lhs(formulas[[3]]))
  wh_q <- setdiff(c(lhs(formulas[[1]]), unlist(causl:::rhs_vars(formulas[[2]]))),
                  c(unlist(causl:::rhs_vars(formulas[[3]])), vars))
  nv <- length(wh_q)

  ## simulate values and put in appropriate structure
  beta_m <- matrix(rnorm(np*nv, mean=par_mns, sd=par_sds), np, nv)
  # which(samp %*% beta_m
  beta_pars <- apply(beta_m, 2, c, simplify = FALSE)
  beta_pars <- lapply(beta_pars, function(x) list(beta=x))
  names(beta_pars) <- wh_q

  ## now add
  beta_pars <- list(Y=beta_pars)

  return(beta_pars)
}




