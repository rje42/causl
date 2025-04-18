##' Sample from a causal model
##'
##' Obtain samples from a causal model parameterized as in Evans and Didelez (2024).
##'
##' @param n number of samples required
##' @param causl_model object of class `causl_model`
##' @param formulas list of lists of formulas
##' @param family families for variables and copula
##' @param pars list of lists of parameters
##' @param link list of link functions
##' @param dat optional data frame of covariates
##' @param estimand quantity to control, default is `"ATE"`
##' @param method either `"inversion"` (the default), `"inversion_mv"`, or `"rejection"`
##' @param control list of options for the algorithm
##' @param ... other arguments, such as custom families
## @param seed random seed used for replication
##'
##' @details Samples from a given causal model under the frugal
##' parameterization.
##'
##' The entries for  `formula` and `family` should each be a
##' list with four entries, corresponding to the \eqn{Z}, \eqn{X}, \eqn{Y} and
##' the copula.  `formula` determines the model, so it is crucial that
##' every variable to be simulated is represented there exactly once.  Each
##' entry of that list can either be a single formula, or a list of formulae.
##' Each corresponding entry in `family` should be the same length as the
##' list in `formula` or of length 1 (in which case it will be repeated
##' for all the variables therein).
##'
##' We use the following codes for different families of distributions:
##'   | val | family      |
##'   |----:|:------------|
##'   |   0 | binomial    |
##'   |   1 | gaussian    |
##'   |   2 | t           |
##'   |   3 | Gamma       |
##'   |   4 | beta        |
##'   |   5 | binomial    |
##'   |   6 | lognormal   |
##'   |  11 | ordinal     |
##'   |  10 | categorical |
##'
##' The family variables for the copula are also numeric and taken from
##' `VineCopula`. Use, for example, 1 for Gaussian, 2 for t, 3 for Clayton,
##' 4 for Gumbel, 5 for Frank, 6 for Joe and 11 for FGM copulas.
##'
##' `pars` should be a named list containing variable names that correspond
##' to the LHS of
##' formulae in `formulas`.  Each of these should themselves be a list
##' containing `beta` (a vector of regression parameters) and (possibly)
##' `phi`, a dispersion parameter.  For any discrete variable that is a
##' treatment, you can also specify `p`, an initial proportion to simulate
##' from (otherwise this defaults to 0.5).
##'
##' Link functions for the Gaussian, t and Gamma distributions can be the
##' identity, inverse or log functions.  Gaussian and t-distributions default to
##' the identity, and Gamma to the log link.  For the Bernoulli the logit,
##' probit, and log links are available.
##'
##' A variety of sampling methods are implemented.  The
##' inversion method with pair-copulas is the default (`method="inversion"`),
##' but we cam also use a multivariate copula (`method="inversion_mv"`) or even
##' rejection sampling (`method="rejection"`).  Note that the `inveresion_mv`
##' method simulates the entire copula, so it cannot depend upon intermediate
##' variables.
##'
##' The only control parameters are `cop`: which gives a keyword for the
##' copula that defaults to `"cop"`; `quiet` which defaults to `FALSE` but will
##' reduce output if set to `TRUE`; and (if rejection sampling is selected)
##' `careful`: this logical enables one to implement the full rejection
##' sampling method, which means we do get exact samples (note this method
##' is generally very slow, especially if we have an outlying value, so the
##' default is `FALSE`).
##'
##' @examples
##' pars <- list(z=list(beta=0, phi=1),
##'              x=list(beta=c(0,0.5), phi=1),
##'              y=list(beta=c(0,0.5), phi=0.5),
##'              cop=list(beta=1))
##' rfrugalParam(100, pars = pars)
##'
##'
##' @return A data frame containing the simulated data.
##' @export
rfrugal <- function (n, causl_model, control=list()) {

  # get control parameters or use defaults
  con <- list(careful = FALSE, cop="cop", quiet = FALSE, quan_tol = 1e3*.Machine$double.eps,
              pm_cond = TRUE, pm_cor_thresh = 0.25, pm_nlevs = 5)
  matches <- match(names(control), names(con))
  con[matches] <- control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ",
                                   paste(names(control[is.na(matches)]),
                                         sep = ", "))
  ## get keyword for copula formula
  kwd <- con$cop

  ## set up output
  if (missing(n)) n <- nrow(causl_model$dat)
  out <- data.frame(matrix(0, ncol=length(causl_model$vars), nrow=n))
  names(out) <- causl_model$vars

  if (!is.null(causl_model$dat)) {
    out <- cbind(causl_model$dat, out)
    # proc_inputs$vars <- c(names(dat), proc_inputs$vars)
  }

  ## obtain method
  method <- causl_model$method

  ## choose appropriate method
  if (method == "inversion") {  ## inversion method
    if (!con$quiet) message("Inversion method selected: using pair-copula parameterization")
    out <- sim_inversion(out, causl_model)
  }
  else if (method == "rejection") { ## rejection sampling method
    if (!con$quiet) message("Rejection method selected: using multivariate copula parameterization")
    out <- sim_rejection(out, causl_model, careful=con$careful)
  }
  else if (method == "inversion_mv") { # multivariate copula
    if (!con$quiet) message("Inversion with multivariate-copula selected")
    out <- sim_multi(out, causl_model)
  }
  else stop("Not a recognised method")

  rownames(out) <- NULL

  return(out)
}


##' @describeIn rfrugal old function for simulation
##' @export
rfrugalParam <- function(n, formulas = list(list(z ~ 1), list(x ~ z), list(y ~ x), list( ~ 1)),
                       family = c(1,1,1,1), pars, link=NULL, dat=NULL,
                       method="inversion", control=list(), ...) {

  # get control parameters or use defaults
  con <- list(careful = FALSE, quiet = FALSE, cop="cop", quan_tol = 1e3*.Machine$double.eps,
              pm_cond = TRUE, pm_cor_thresh = 0.25, pm_nlevs = 5)
  matches <- match(names(control), names(con))
  con[matches] <- control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ",
                                   paste(names(control[is.na(matches)]),
                                         sep = ", "))
  ## get keyword for copula formula
  kwd <- con$cop

  ## infer sample size from dat if possible
  datNULL <- is.null(dat)
  if (datNULL && missing(n)) stop("Must specify a sample size or supply a data frame")
  else if (missing(n)) n <- nrow(dat)
  else if (!datNULL && n != nrow(dat)) {
    warning("Dimension of supplied data does not match n, using size of supplied data frame")
    n <- nrow(dat)
  }

  ## process the four main arguments (change to use causl_model())
  proc_inputs <- process_inputs(formulas=formulas, family=family, pars=pars,
                                link=link, dat=dat, kwd=kwd, method=method, control=con)

  out <- rfrugal(n=n, causl_model=proc_inputs, control=con)

  return(out)
}

