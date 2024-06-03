##' Simulate a single variable using the inversion method
##'
##' @param n sample size
##' @param formulas list consisting of a formula for the output variables and a list of formulae for the pair-copula
##' @param family list containing family variable
##' @param pars list with two entries, first a list of parameters for response, and second a further list of parameters for pair-copula
##' @param link list of same form as `family`
##' @param dat data frame of current variables
##' @param quantiles data frame of quantiles
##'
##' @return The data frame `dat` with an additional column given by the left-hand side of `formula[[1]]`.
##'
##' @description
##' Each entry `formulas`, `family`, `pars`, `link` is a list
##' with two entries, the first referring to the variable being simulated and the
##' second to the pair-copulas being used.
##'
##' @export
sim_variable <- function (n, formulas, family, pars, link, dat, quantiles) {
  qY <- runif(n)

  ## get variables
  vnm <- lhs(formulas[[1]])
  if (length(vnm) != 1L) stop("Unable to extract variable name")

  quantiles[[vnm]] <- qY

  LHS_cop <- lhs(formulas[[2]])

  for (i in rev(seq_along(formulas[[2]]))) {
    X <- model.matrix(formulas[[2]][[i]], data=dat)
    # eta <- X %*% pars[[2]][[i]]$beta

    ## rescale quantiles for pair-copula
    qs <- cbind(quantiles[[LHS_cop[[i]]]], qY)
    qY <- rescale_cop(qs, X=X, beta=pars[[2]][[i]]$beta, family=family[[2]][i],
                     par2=pars[[2]][[i]]$par2) #, link=link[[2]][j])
  }

  ## now rescale to correct margin
  X <- model.matrix(delete.response(terms(formulas[[1]])), data=dat)

  #Y <- rescale_var(qY, X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  Y <- rescale_var(runif(n), X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  dat[[vnm]] <- Y
  # quantiles[[vnm]] <- qY
  attr(dat, "quantiles") <- quantiles

  return(dat)
}
