##' Simulate for single time-step
##'
##' @param out data frame for output
##' @param proc_inputs output of \code{process_inputs()}
##' @param control list of control parameters
##'
##' @details `sim_inversion` and `sim_rejection` correspond to
##' performing the sampling by inversion or using rejection sampling.`sim_multi` first 
##' simulates from the copula then transforms to the correct margins in the correct causal ordering
##'
##' @export
sim_multi <- function (out, proc_inputs) {
  ## unpack proc_inputs
  formulas <- proc_inputs$formulas
  pars <- proc_inputs$pars
  family <- proc_inputs$family
  link <- proc_inputs$link
  dZ <- proc_inputs$dim[1]; dX <- proc_inputs$dim[2]; dY <- proc_inputs$dim[3]
  LHS_Z <- proc_inputs$LHSs$LHS_Z; LHS_X <- proc_inputs$LHSs$LHS_X; LHS_Y <- proc_inputs$LHSs$LHS_Y; kwd <- proc_inputs$kwd
  famZ <- proc_inputs$family[[1]]; famX <- proc_inputs$family[[2]]; famY <- proc_inputs$family[[3]]; famCop <- proc_inputs$family[[4]]
  # dC <-
  
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
  
  ## Simulate Copula
  #TODO: add in the other COPULA types, besides normal, T, FGM
  
  beta_vector <- pars$cop$Y[[1]]$beta
  beta_vector <- 2 * expit(beta_vector) - 1
  
  # infer dimension (only do marginal uncoditional on X) for now
  SIGMA <- diag(dZ + 1)
  SIGMA[upper.tri(SIGMA)] <- beta_vector
  SIGMA[lower.tri(SIGMA)] <- t(SIGMA[upper.tri(SIGMA)])
 

  copulafams <- family[[4]]
  # TODO: allow for more than one copula fam only do one copula fam for now
  if(length(unique(copulafams)) != 1) stop()
  # TODO: allow for different copula's for Y, Z_i
  # 
  # us <- matrix(nrow = n)
  # for(copulafam in unlist(copulafams)){
  #   if(copulafam == 1){
  #     us_i <- rGaussCop(n, SIGMA)
  #   }
  #   else if(copulafam == 2){
  #     df <- pars$cop$Y[[1]]$df
  #     us_i <- rtCop(n, SIGMA, df)
  #   }
  #   else if(copulafam == 11){
  #     alpha <- pars$cop$Y[[1]]$alpha
  #     us_i <- rfgmCopula(n, alpha)
  #   }
  #   us <- cbind(us, us_i)
  # }
  
  copulafam <- unlist(copulafams)[1]
  if(copulafam == 1){
    us <- rGaussCop(n, SIGMA)
  }
  else if(copulafam == 2){
    df <- pars$cop$Y[[1]]$df
    us <- rtCop(n, SIGMA, df)
  }
  else if(copulafam == 11){
    alpha <- pars$cop$Y[[1]]$alpha
    us <- rfgmCopula(n, alpha)
  }
  
  
  for (i in seq_along(order)) {
    vnm <- vars[order[i]]
    
    ## code to get Y quantiles conditional on different Zs
    if (order[i] > dZ+dX) {
      ## simulate Y variable
      # qY <- runif(n)
      wh <- order[i] - dZ - dX
      wh_u <- order[i] - dZ
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
      
      
      # now rescale to correct margin

      X <- model.matrix(delete.response(terms(forms[[1]])), data=out)
      Y <- rescale_var(us[,wh_u], X=X, family=fams[[1]], pars=prs[[1]], link=lnk[[1]])
      out[[vnm]] <- Y
    
      
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
      if(vnm %in% LHS_X){
        out[[vnm]] <- rescale_var(runif(n), X=MM, family=curr_fam, pars=pars[[vnm]], 
                                  link=curr_link)
      }
      else{
        out[[vnm]] <- rescale_var(us[,i], X=MM, family=curr_fam, pars=pars[[vnm]], 
                                  link=curr_link)
      }

      
    }
  }
  
  attr(out, "qZ") <- quantiles
  
  return(out)
}
