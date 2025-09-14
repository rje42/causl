##' Estimate the power under a causl_model and estimator
##'
##' @param causl_model the specified causal model
##' @param n number of observations drawn from causal model
##' @param estimator method to estimate the coefficient and standard errors
##' @param monte_carlo_sims the number of simulations. Defaults to 500
##' @param alpha_level significance level one reject

##' @export
estimate_power <- function(causl_model, n, estimator,monte_carlo_sims = 500, 
                           alpha_level = 0.05){
  
  pvals <- numeric(length = monte_carlo_sims)
  for(i in seq(monte_carlo_sims)){
    dat <- rfrugal(n, causl_model  = causl_model)
    causal_estimate <- get_causal_estimate(dat, estimator)
    pvals[i] <- causal_estimate$pval
  }

  # Create a summary data.frame with short column names
  res_table <- data.frame(
    Power = round(mean(pvals < alpha_level), 3),
    Alpha = alpha_level,
    N = n,
    Sims = monte_carlo_sims,
    stringsAsFactors = FALSE
  )
  return(res_table)

}

##' A helper function to calculate estimates from data. Make better to allow
##' for more estimators.
##'
##' @param data the data
##' @param estimator method to estimate the coefficient and standard errors
get_causal_estimate <- function(data, estimator){
  outcome_info <- estimator$outcome_model
  ipw_info <- estimator$ipw_model
  
  outcome_model <- outcome_info$model
  outcome_se <- outcome_info$se
  outcome_formula <- outcome_info$formula
  outcome_family <- outcome_info$family

  
  ipw_weights <- get_ipw_estimate(data, ipw_info = ipw_info)
  
  design <- survey::svydesign(~1, probs = ipw_weights,  data=data)
  mod <- outcome_model(outcome_formula, design = design)
  confint_X <- confint(mod)[2,]
  return(list(pval = summary(mod)$coefficients[2,"Pr(>|t|)"],
              confint = confint_X))
}

##' A helper function to calculate ipw weights from data. 
##'
##' @param data the data
##' @param ipw_info information about the propensity score methods
get_ipw_estimate <- function(data, ipw_info){
  
  ipw_formula <- ipw_info$formula
  ipw_model <- ipw_info$model
  ipw_family <- ipw_info$family
  
  ipw_weights <- fitted(ipw_model(ipw_formula, data = data,family = ipw_family))
  return(ipw_weights)
}

