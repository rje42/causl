##' Simulate for single time-step
##'
##' @param out data frame for output
##' @param proc_inputs output of `process_inputs()`
## @param control list of control parameters
##'
##' @details `sim_inversion` and `sim_rejection` correspond to
##' performing the sampling by inversion or using rejection sampling.
##'
##' @export
sim_inversion <- function (out, proc_inputs) {
  ## unpack proc_inputs
  formulas <- proc_inputs$formulas
  pars <- proc_inputs$pars
  family <- proc_inputs$family
  link <- proc_inputs$link
  dZ <- proc_inputs$dim[1]; dX <- proc_inputs$dim[2]; dY <- proc_inputs$dim[3]
  LHS_Z <- proc_inputs$LHSs$LHS_Z; LHS_X <- proc_inputs$LHSs$LHS_X; LHS_Y <- proc_inputs$LHSs$LHS_Y; kwd <- proc_inputs$kwd
  famZ <- proc_inputs$family[[1]]; famX <- proc_inputs$family[[2]]; famY <- proc_inputs$family[[3]]; famCop <- proc_inputs$family[[4]]

  vars <- proc_inputs$vars

  ## get quantiles, if any available
  if (is.null(proc_inputs$quantiles)) {
    quantiles <- out
  }
  else {
    quantiles <- cbind(proc_inputs$quantiles, out[vars])
  }


  ## code to get causal order
  order <- proc_inputs$order

  ## sample size
  n <- nrow(out)

  for (i in seq_along(order)) {
    vnm <- vars[order[i]]

    ## code to get Y quantiles conditional on different Zs
    if (order[i] > dZ+dX) {
      ## simulate Y variable
      # qY <- runif(n)
      wh <- order[i] - dZ - dX
      # print(wh)

      ## code to use sim_variable
      forms <- list(formulas[[3]][[wh]], formulas[[4]][[wh]])
      fams <- list(family[[3]][[wh]], family[[4]][[wh]])
      prs <- list(pars[[vnm]], pars[[kwd]][[vnm]])
      lnk <- list(link[[3]][wh], list()) # link[[4]][[wh]])

      if (any(is.na(lhs(forms[[2]])))) {
        forms[[2]] <- `lhs<-`(forms[[2]], c(LHS_Z[rank(order[seq_len(dZ)])],
                                            LHS_Y[rank(order[dZ+dX+seq_len(i-1 - dZ-dX)])]))
      }

      out <- sim_variable(n=nrow(out), formulas=forms, family=fams, pars=prs,
                          link=lnk, dat=out, quantiles=quantiles)
      quantiles <- attr(out, "quantiles")
      attr(out, "quantiles") <- NULL
    }
    else {
      ## code to simulate Z and X variables in causal order
      curr_link <- unlist(link)[order[i]]

      if (vnm %in% LHS_X) {
        curr_form <- formulas[[2]][[order[i]-dZ]]
        curr_fam <- famX[[order[i]-dZ]]
      }
      else {
        curr_form <- formulas[[1]][[order[i]]]
        curr_fam <- famZ[[order[i]]]
      }
      trm <- terms(curr_form)
      # curr_form2 <- delete.response(terms(curr_form))
      MM <- model.matrix(delete.response(trm), data=out)
      if (nrow(MM) != nrow(out)) {
        if (length(attr(trm, "factors")) == 0) {
          if (attr(trm, "intercept") == 1) MM <- matrix(1, nrow=nrow(out), ncol=1)
          else MM <- matrix(0, nrow=nrow(out), ncol=0)
        }
        else warning(paste0("Missing entries for ", vnm))
      }
      eta <- MM %*% pars[[vnm]]$beta
      curr_phi <- pars[[vnm]]$phi
      oth_pars <- pars[[vnm]][!(names(pars[[vnm]]) %in% c("beta", "phi"))]
      tmp <- glm_sim(family=curr_fam, eta=eta, phi=curr_phi, other_pars=oth_pars, link=curr_link)
      if (vnm %in% LHS_Z) quantiles[[vnm]] <- attr(tmp, "quantile")
      attr(tmp, "quantile") <- NULL
      out[[vnm]] <- tmp
    }
  }

  attr(out, "qZ") <- quantiles

  return(out)
}


##' Simulate a single variable using the inversion method
##'
##' @param n sample size
##' @param formulas list consisting of a formula for the output variables and a list of formulae for the pair-copula
##' @param family list containing family variable
##' @param pars list with two entries, first a list of parameters for response, and second a further list of parameters for pair-copula
##' @param link list of same form as `family`
##' @param dat data frame of current variables
##' @param quantiles data frame of quantiles
##' @param tol tolerance for quantile closeness to 0 or 1
##'
##' @return The data frame `dat` with an additional column given by the left-hand side of `formula[[1]]`.
##'
##' @description
##' Each entry `formulas`, `family`, `pars`, `link` is a list
##' with two entries, the first referring to the variable being simulated and the
##' second to the pair-copulas being used.
##'
##' If the quantile is smaller than `tol` or bigger than 1 - `tol` the
##' associated value is moved to `tol` and 1 - `tol` respectively.
##'
##' @export
sim_variable <- function (n, formulas, family, pars, link, dat, quantiles,
                          tol=10*.Machine$double.eps) {
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
    qY <- rescale_cop(qs, X=X, beta=pars[[2]][[i]]$beta, family=family[[2]][[i]],
                      df=pars[[2]][[i]]$df) #, link=link[[2]][j])
    ##
    if (max(qY) > 1 - tol) {
      warning(paste0("Quantiles numerically 1, reducing by ", tol))
      qY <- pmin(qY, 1 - tol)
    }
    if (min(qY) < tol) {
      qY <- pmax(qY, tol)
      warning(paste0("Quantiles numerically 0, increasing by", tol))
    }
  }

  ## now rescale to correct margin
  X <- model.matrix(delete.response(terms(formulas[[1]])), data=dat)

  Y <- rescale_var(qY, X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  # Y <- rescale_var(runif(n), X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  dat[[vnm]] <- Y
  # quantiles[[vnm]] <- qY
  attr(dat, "quantiles") <- quantiles

  return(dat)
}
