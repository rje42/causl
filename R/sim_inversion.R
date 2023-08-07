##' Simulate for single time-step
##'
##' @param out data frame for output
##' @param proc_inputs output of \code{process_inputs()}
##' @param control list of control parameters
##'
##' @details \code{sim_inversion} and \code{sim_rejection} correspond to
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

  ## code to get causal order
  order <- proc_inputs$order
  qZs <- out[LHS_Z]

  ## sample size
  n <- nrow(out)

  for (i in seq_along(order)) {
    ## code to get Y quantiles conditional on different Zs
    if (order[i] > dZ+dX) {
      ## simulate Y variable
      qY <- runif(n)
      wh <- order[i] - dZ - dX
      # print(wh)

      for (j in seq_len(dZ)) {
        curr_qZ <- qZs[[vars[j]]]
        X <- model.matrix(formulas[[4]][[wh]][[j]], data=out)
        curr_fam <- family[[4]][wh,j]
        curr_par <- pars[[kwd]]$beta[[wh]][[j]]
        # eta <- X %*% curr_par
        qY <- rescaleCop(cbind(curr_qZ,qY), X=X, pars=curr_par, family=curr_fam) #, link=link[[4]][i,j])
      }
      ##
      X <- model.matrix(formulas[[3]][[wh]], data=out)
      qY <- rescaleVar(qY, X=X, family=famY[wh], pars=pars[[LHS_Y[wh]]],
                       link=link[[3]][wh])

      out[[vars[order[i]]]] <- qY
    }
    else {
      ## code to simulate Z and X variables in causal order
      vnm <- vars[order[i]]
      curr_link <- unlist(link)[order[i]]

      if (vnm %in% LHS_X) {
        curr_form <- formulas[[2]][[order[i]-dZ]]
        curr_fam <- famX[order[i]-dZ]
      }
      else {
        curr_form <- formulas[[1]][[order[i]]]
        curr_fam <- famZ[order[i]]
      }
      MM <- model.matrix(curr_form, data=out)
      if (nrow(MM) != nrow(out)) {
        trm <- terms(curr_form)
        if (length(attr(trm, "factors")) == 0) {
          if (attr(trm, "intercept") == 1) MM <- matrix(1, nrow=nrow(out), ncol=1)
          else MM <- matrix(0, nrow=nrow(out), ncol=0)
        }
        else warning(paste0("Missing entries for ", vnm))
      }
      eta <- MM %*% pars[[vnm]]$beta
      curr_phi <- pars[[vnm]]$phi
      tmp <- glm_sim(fam=curr_fam, eta=eta, phi=curr_phi, link=curr_link,
                     par2=pars[[vnm]]$par2)
      out[[vnm]] <- tmp
      if (vnm %in% LHS_Z) qZs[[vnm]] <- attr(tmp, "quantile")
    }
  }

  return(out)
}
