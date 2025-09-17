##' Estimate statistical power for a causal effect estimator using Monte Carlo simulation
##'
##' This function simulates data from a specified causal model multiple times, applies a given estimator to each simulated dataset,
##' and calculates the proportion of times the null hypothesis is rejected at a specified significance level (i.e., the estimated power).
##'
##' @param causl_model An object specifying the data-generating causal model from which to simulate datasets.
##' @param n Integer. The number of observations to generate in each simulated dataset.
##' @param estimator A list or function specifying the estimation procedure, including the outcome model and standard error computation.
##' @param monte_carlo_sims Integer. The number of Monte Carlo simulations to perform. Defaults to 500.
##' @param alpha_level Numeric. The significance level at which to test the null hypothesis (e.g., 0.05).
##'
##' @return A data.frame summarizing the estimated power, significance level, sample size, and number of simulations.
##' @export
estimate_power <- function(causl_model, n, estimator, monte_carlo_sims = 500, 
                          alpha_level = 0.05) {
    pvals <- numeric(length = monte_carlo_sims)
    for (i in seq(monte_carlo_sims)) {
        dat <- rfrugal(n, causl_model = causl_model)
        causal_estimate <- est_causal(dat, estimator)
        pvals[i] <- causal_estimate$pval
    }

    # Create a summary data.frame with short column names
    res_table <- data.frame(
        Power = round(mean(pvals < alpha_level), 3),
        alpha = alpha_level,
        N = n,
        Sims = monte_carlo_sims,
        stringsAsFactors = FALSE
    )

    return(res_table)
}

##' Helper function to compute causal effect estimates from data using a specified estimator.
##'
##' This function applies the provided estimator to the input data to estimate the causal effect
##' and its associated statistics (e.g., p-value, confidence interval). The estimator should
##' specify the outcome model and the inverse probability weighting (IPW) model, allowing for
##' flexible estimation strategies.
##'
##' @param data A data.frame containing the observed data for estimation.
##' @param estimator A list specifying the estimation procedure, including the outcome model,
##'   standard error computation, model formula, family, and IPW model details.
##'
##' @return A list containing the estimated p-value and confidence interval for the causal effect.
est_causal <- function(data, estimator) {
    outcome_info <- estimator$outcome_model
    ipw_info <- estimator$ipw_model

    outcome_model <- outcome_info$model
    outcome_se <- outcome_info$se
    outcome_formula <- outcome_info$formula
    outcome_family <- outcome_info$family

    ipw_weights <- ipw_weights_from(data, ipw_info = ipw_info)

    design <- survey::svydesign(~1, probs = ipw_weights, data = data)
    mod <- outcome_model(outcome_formula, design = design)
    confint_X <- confint(mod)[2, ]
    return(list(
        pval = summary(mod)$coefficients[2, "Pr(>|t|)"],
        confint = confint_X
    ))
}

##' Computes inverse probability weights (IPW) for each observation in the data.
##'
##' This function fits the specified propensity score model to the data and returns
##' the predicted probabilities (propensity scores) as weights for use in IPW estimation.
##'
##' @param data A data.frame containing the variables required by the propensity score model.
##' @param ipw_info A list containing the propensity score model specification, including:
##'   \itemize{
##'     \item \code{formula}: a formula specifying the model for the treatment assignment.
##'     \item \code{model}: a function to fit the propensity score model (e.g., \code{glm}).
##'     \item \code{family}: the error distribution and link function to be used in the model.
##'   }
##' @return A numeric vector of fitted propensity scores (IPW weights) for each observation.
ipw_weights_from <- function(data, ipw_info) {
    ipw_formula <- ipw_info$formula
    ipw_model <- ipw_info$model
    ipw_family <- ipw_info$family

    ipw_weights <- fitted(ipw_model(ipw_formula, data = data, family = ipw_family))
    return(ipw_weights)
}
