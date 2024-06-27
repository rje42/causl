##' Sample from a causal model
##'
##' Obtain samples from a causal model using the rejection sampling approach of
##' Evans and Didelez (2024).
##'
##' @param n number of samples required
##' @param formulas list of lists of formulas
##' @param family families for variables and copula
##' @param pars list of lists of parameters
##' @param link list of link functions
##' @param dat optional data frame of covariates
##' @param careful logical: should full rejection sampling method be used with
##' correctly computed weight?
##' @param method either `"inversion"` (the default) or `"rejection"`
##' @param control list of options for the algorithm
##' @param ... other arguments, such as custom families
##' @param seed random seed used for replication
##'
##' @details Samples from a given causal model using rejection sampling (or,
##' if everything is discrete, direct sampling).
##'
##' The logical `careful` enables one to implement the full rejection
##' sampling method, which means we do get exact samples.  However, this method
##' is generally very slow, and in particular if we have an outlying value it
##' may be interminably so.
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
##' 0 or 5 = binary;
##' 1 = normal;
##' 2 = t-distribution;
##' 3 = gamma;
##' 4 = beta;
##' 6 = log-normal.
##'
##' The family variables for the copula are also numeric and taken from
##' `VineCopula`.
##' Use, for example, 1 for Gaussian, 2 for t, 3 for Clayton, 4 for Gumbel,
##' 5 for Frank, 6 for Joe and 11 for FGM copulas.
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
##' the identity, and Gamma to the log link.  For the Bernoulli the logit and
##' probit links are available.
##'
##' Control parameters are `trace` (default value 0, increasing to 1
##' increases verbosity of output), `warn` (which currently does nothing),
##' `max_wt` which is set to 1, and increases each time the function
##' is recalled.
## (if weights empirically appear not to have an upper bound, this warns if set
## to 1 (the default) and stops if set to 2), ...
##' Control parameters also include `cop`, which gives a keyword for the
##' copula that defaults to `"cop"`.
##'
##' The logical `careful` enables one to implement the full rejection
##' sampling method, which means we do get exact samples.  However, this method
##' may be slow, and in particular if we have an outlying value it may run very
##' slowly indeed.
##'
##' @examples
##' pars <- list(z=list(beta=0, phi=1),
##'              x=list(beta=c(0,0.5), phi=1),
##'              y=list(beta=c(0,0.5), phi=0.5),
##'              cop=list(beta=1))
##' rfrugalParam(100, pars = pars)
##'
## @importFrom frugalSim sim_chain
##'
##' @return A data frame containing the simulated data.
##'
##' @export
rfrugalParam <- function(n, formulas = list(list(z ~ 1), list(x ~ z), list(y ~ x), list( ~ 1)),
                       family = c(1,1,1,1), pars, link=NULL, dat=NULL, careful=FALSE,
                       method="inversion", control=list(), ..., seed) {

  # get control parameters or use defaults
  con = list(max_wt = 1, warn = 1, cop="cop", trace=0)
  matches = match(names(control), names(con))
  con[matches] = control[!is.na(matches)]
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

  ## inform user of parameterization being used
  if (method == "inversion") message("Inversion method selected: using pair-copula parameterization")
  else if (method == "rejection") message("Rejection method selected: using multivariate copula parameterization")

  ## process the four main arguments
  proc_inputs <- process_inputs(formulas=formulas, pars=pars, family=family, link=link,
                        dat=dat, kwd=kwd, method=method)

  # if (method == "particle") {
  #   forms <- tidy_formulas(unlist(formulas), kwd)
  #   form_all <- merge_formulas(forms)
  #   msks <- masks(forms, family=family, wh=form_all$wh, cp=length(output))
  #
  #   theta <- pars2mask(pars, msks)
  #
  #   ## construct suitable log-likelihood function
  #   llhd <- function(x) {
  #     if (!is.matrix(x)) x <- matrix(x, nrow=1)
  #     colnames(x) <- vars2
  #     mm <- model.matrix.default(form_all$formula, data=as.data.frame(x))
  #     ll(x, mm, beta=theta$beta, phi=theta$phi,
  #        inCop = match(c(LHS_Z,LHS_Y), vars2), fam_cop = unlist(last(family)),
  #        family = unlist(family)[seq_along(vars)])
  #   }
  #
  #   vars2 <- vars
  #
  #   dat <- as.data.frame(matrix(NA, n, d))
  #   names(dat) <- vars2
  #
  #   ## get parameters for particle simulator
  #   d <- sum(lengths(family[1:3]))
  #   dc <- d - sum(family[[1]] == 0 | family[[1]] == 5)
  #   n_state <- rep(2, d-dc)
  #   prob <- rep(list(c(0.5,0.5)), d-dc)
  #   attrib <- list(dc=dc, n_state=n_state, prob=prob)
  #
  #   for (i in seq_len(n)) {
  #     dat[i,] <- sim_chain(n=n, d=d, ltd=llhd, ns=1, sd1=2, attrib,
  #                          sim="importance", control=list(B=10))$x
  #     if (con$trace > 0) rje::printPercentage(i,n)
  #   }
  #
  #   return(dat)
  # }

  ## set up output
  out <- data.frame(matrix(0, ncol=length(proc_inputs$vars), nrow=n))
  names(out) <- proc_inputs$vars

  if (!datNULL) {
    out <- cbind(dat, out)
    # proc_inputs$vars <- c(names(dat), proc_inputs$vars)
  }
  if (anyDuplicated(na.omit(proc_inputs$vars))) stop("duplicated variable names")

  ## choose appropriate method
  if (method == "inversion") {  ## inversion method
    out <- sim_inversion(out, proc_inputs)
  }
  else if (method == "rejection") { ## rejection sampling method
    out <- sim_rejection(out, proc_inputs, careful, control)
  }
  else if (method == "multi_copula"){ # multivariate copula
    out <- sim_multi(out, proc_inputs)
  }
  else stop("Not a recognised method")

  rownames(out) <- NULL

  return(out)
}

