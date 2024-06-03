##' Simulate for single time-step
##'
##' @param out data frame for output
##' @param proc_inputs output of \code{process_inputs()}
##' @param control list of control parameters
##'
##' @details `sim_inversion` and `sim_rejection` correspond to
##' performing the sampling by inversion or using rejection sampling.
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
  #TODO: fix for more than just normal copula
  
  beta_vector <- pars$cop$Y[[1]]$beta
  beta_vector <- 2 * expit(beta_vector) - 1
  # infer dimension
  cop <- normalCopula(beta_vector, dim = 2, dispstr = 'un')
  
  # simulate 
  us <- rCopula(n,cop)

  
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
      
      # out[[vars[order[i]]]] <- sim_Y(n, formulas=formulas[[4]][[wh]],
      #                                family=family[[4]][[wh]],
      #                                pars=pars[[kwd]][[wh]],
      #                                formY = formulas[[3]][[wh]],
      #                                famY=family[[3]][wh],
      #                                parsY=pars[[LHS_Y[wh]]],
      #                                linkY=link[[3]][wh], qZ=quantiles, vars=vars,
      #                                dat=out)
      
      # for (j in seq_len(dZ)) {
      #   curr_qZ <- qZs[[vars[j]]]
      #   X <- model.matrix(formulas[[4]][[wh]][[j]], data=out)
      #   curr_fam <- family[[4]][wh,j]
      #   curr_par <- pars[[kwd]]$beta[[wh]][[j]]
      #   # eta <- X %*% curr_par
      #   qY <- rescale_cop(cbind(curr_qZ,qY), X=X, pars=curr_par, family=curr_fam) #, link=link[[4]][i,j])
      # }
      # ##
      # X <- model.matrix(formulas[[3]][[wh]], data=out)
      # qY <- rescale_var(qY, X=X, family=famY[[wh]], pars=pars[[LHS_Y[wh]]],
      #                  link=link[[3]][wh])
      #
      # out[[vars[order[i]]]] <- qY
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
      if(order[i] > dZ){
        out[[vnm]] <- rescale_var(runif(n), X=MM, family=curr_fam, pars=pars[[vnm]], 
                                  link=curr_link)
      }
      else{
        out[[vnm]] <- rescale_var(us[,i], X=MM, family=curr_fam, pars=pars[[vnm]], 
                                  link=curr_link)
      }

      
    }
  }
  
  
  # Z1 <- qnorm(us[,1], mean = 0, sd = 1)
  # Z2 <- qnorm(us[,2], mean = 0, sd = 1)  # changed to be standard normal
  # 
  # X <- rbinom(n = n, size = 1, prob = expit(Z1 + Z2))
  # # X | Z_1, Z_2 ~ bernoulli(p = expit(Z1 + Z2))
  # 
  # Y <- qnorm(us[,3], mean = 1 + 2*X , sd = 1)
  # # Y | do(X) ~ N(1 + 2*X, 1)
  # # combine XYZ into a new dataset
  # dat1 <- data.frame(Z1=Z1,Z2=Z2,X=X,Y=Y)
  
  attr(out, "qZ") <- quantiles
  
  return(out)
}

# ##' Generate outcome data by inversion
# ##'
# ##' @param n number of samples
# ##' @param formulas list of formulae for copula parameters
# ##' @param family vector of integers for distributions of variables
# ##' @param pars vector of parameters for copula
# ##' @param formY formulae for response variable
# ##' @param famY family for response variable
# ##' @param parsY regression parameters for response variable
# ##' @param linkY link function for response variable
# ##' @param qZ quantiles of covariate distribution
# ##' @param vars character vector of variable names
# ##' @param dat data frame containing variables in `formulas` and `formY`
# ##'
# sim_Y <- function(n, formulas, family, pars, formY, famY, parsY, linkY, qZ, vars, dat) {
#
#   qY <- runif(nrow(qZ))  # uniform quantiles
#
#   for (j in seq_along(formulas)) {
#     ## for each Z variable
#     quZ <- qZ[[vars[j]]]
#     X <- model.matrix(delete.response(terms(formulas[[j]])), data=dat)  ## covariates matrix
#     # curr_fam <- family[j]
#     # curr_par <- pars[[j]]
#     # eta <- X %*% curr_par
#     qY <- rescale_cop(cbind(quZ,qY), X=X, beta=pars[[j]]$beta, family=family[j],
#                      par2=pars[[j]]$par2) #, link=link[j])
#   }
#   ##
#   X <- model.matrix(formY, data=dat)
#   qY <- rescale_var(qY, X=X, family=famY, pars=parsY, link=linkY)
#
#   qY
# }
