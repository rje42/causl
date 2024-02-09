# ##' Process inputs given in linear form
# process_inputs2(formulas, pars, family, link, kwd, ordering=FALSE) {
#   if (length(formulas))
# }

##' Sample from a causal model
##'
##' Obtain samples from a causal model using the rejection sampling approach of
##' Evans and Didelez (2024).
##'
##' @param n number of samples required
##' @param formulas list of lists of formulas
##' @param pars list of lists of parameters
## @param data optional data frame of covariates
##' @param family families for Z,X,Y and copula
##' @param link list of link functions
##' @param dat data frame of covariates
##' @param method only \code{"rejection"} is valid
##' @param control list of options for the algorithm
##' @param seed random seed used for replication
##'
##' @details Samples from a given causal model using rejection sampling (or,
##' if everything is discrete, direct sampling).
##'
##' The entries for  \code{formula} and \code{family} should each be a
##' list with four entries, corresponding to the \eqn{Z}, \eqn{X}, \eqn{Y} and
##' the copula.  \code{formula} determines the model, so it is crucial that
##' every variable to be simulated is represented there exactly once.  Each
##' entry of that list can either be a single formula, or a list of formulae.
##' Each corresponding entry in \code{family} should be the same length as the
##' list in \code{formula} or of length 1 (in which case it will be repeated
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
##' \code{VineCopula}.
##' Use, for example, 1 for Gaussian, 2 for t, 3 for Clayton, 4 for Gumbel,
##' 5 for Frank, 6 for Joe and 11 for FGM copulas.
##'
##' \code{pars} should be a named list containing: either entries \code{z},
##' \code{x}, \code{y} and \code{cop}, or variable names that correspond to the
##' LHS of formulae in \code{formulas}.  Each of these should themselves be a list
##' containing \code{beta} (a vector of regression parameters) and (possibly)
##' \code{phi}, a dispersion parameter.  For any discrete variable that is a
##' treatment, you can also specify \code{p}, an initial proportion to simulate
##' from (otherwise this defaults to 0.5).
##'
##' Link functions for the Gaussian, t and Gamma distributions can be the
##' identity, inverse or log functions.  Gaussian and t-distributions default to
##' the identity, and Gamma to the log link.  For the Bernoulli the logit and
##' probit links are available.
##'
##' Control parameters are \code{oversamp} (default value 10), \code{trace} (default
##' value 0, increasing to 1 increases verbosity of output),
##' \code{max_oversamp} (default value 1000), \code{warn} (which currently does
##' nothing), \code{max_wt} which is set to 1, and increases each time the function
##' is recalled.
## (if weights empirically appear not to have an upper bound, this warns if set
## to 1 (the default) and stops if set to 2), ...
##' Control parameters also include \code{cop}, which gives a keyword for the
##' copula that defaults to \code{"cop"}.
##'
##' This function is kept largely for the replication of simulations from
##' Evans and Didelez (2024).
##'
##' @examples
##' pars <- list(z=list(beta=0, phi=1),
##'              x=list(beta=c(0,0.5), phi=1),
##'              y=list(beta=c(0,0.5), phi=0.5),
##'              cop=list(beta=1))
##' causalSamp(100, pars = pars)
##'
## @importFrom frugalSim sim_chain
##'
##' @return A data frame containing the simulated data.
##'
##' @references Evans, R.J. and Didelez, V. Parameterizing and simulating from
##' causal models (with discussion). _Journal of the Royal Statistical Society, Series B_,
##' 2024.
##'
##' @export
causalSamp <- function(n, formulas = list(list(z ~ 1), list(x ~ z), list(y ~ x), list( ~ 1)),
                        pars, family, link=NULL, dat=NULL, method="rejection",
                        control=list(), seed) {

  # get control parameters or use defaults
  con = list(oversamp = 10, max_oversamp=1000, max_wt = 1, warn = 1, cop="cop", trace=0)
  matches = match(names(control), names(con))
  con[matches] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ",
                                   paste(names(control[is.na(matches)]),
                                         sep = ", "))
  if (round(con$oversamp) != con$oversamp) {
    con$oversamp <- ceiling(con$oversamp)
    message("oversamp not an integer, rounding up")
  }
  if (round(con$max_oversamp) != con$max_oversamp) {
    con$max_oversamp <- ceiling(con$max_oversamp)
    message("max_oversamp not an integer, rounding up")
  }
  ## get keyword for copula formula
  kwd <- con$cop

  ## process the four main arguments
  tmp <- process_inputs(formulas=formulas, pars=pars, family=family, link=link,
                        kwd=kwd)
  formulas <- tmp$formulas
  pars <- tmp$pars
  family <- tmp$family
  link <- tmp$link
  dZ <- tmp$dim[1]; dX <- tmp$dim[2]; dY <- tmp$dim[3]
  LHS_Z <- tmp$LHSs$LHS_Z; LHS_X <- tmp$LHSs$LHS_X; LHS_Y <- tmp$LHSs$LHS_Y
  famZ <- tmp$family[[1]]; famX <- tmp$family[[2]]; famY <- tmp$family[[3]]; famCop <- tmp$family[[4]]

  datNULL <- is.null(dat)

  output <- c(LHS_Z, LHS_Y)
  vars <- c(LHS_Z, LHS_X, LHS_Y)
  if (!datNULL) {
    vars <- c(names(dat), vars)
  }
  if (anyDuplicated(na.omit(vars))) stop("duplicated variable names")
  d <- length(vars)

  # ## start to simulate
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


  ## set seed to 'seed'
  if (missing(seed)) {
    seed <- round(1e9*runif(1))
  }
  # else .Random.seed <- seed
  set.seed(seed)

  ## set up data frame for output
  oversamp <- con$oversamp

  if (!datNULL) {
    dat2 <- dat[rep(seq_len(nrow(dat)), each=oversamp), ]
  }

  out <- data.frame(matrix(0, ncol=dZ+dX+dY, nrow=n*oversamp))
  names(out) <- c(LHS_Z, LHS_X, LHS_Y)

  if (!datNULL) {
    out <- cbind(dat2, out)
  }

  ## obtain treatment values
  tmp <- gen_X_values(n*oversamp, famX=famX, pars=pars, LHS_X=LHS_X, dX=dX)
  out[LHS_X] <- tmp$datX[LHS_X]
  qden <- tmp$qden

  # ## give default coefficients
  # if (is.null(pars2$z$beta)) pars2$z$beta = 0
  # if (is.null(pars2$z$phi)) pars2$z$phi = 1

  ## get linear predictors for outcomes
  mms <- vector(mode = "list", length=3)
  mms[[3]] = lapply(formulas[[3]], model.matrix, data=out)
  for (i in seq_along(mms[[3]])) {
    if (ncol(mms[[3]][[i]]) != length(pars[[LHS_Y[i]]]$beta)) stop(paste0("dimension of model matrix for ", LHS_Y[i], " does not match number of coefficients provided"))
  }

  # etas <- vector(mode="list", length=3)
  # for (i in c(1,3)) {
  #   etas[[i]] <- mapply(function(x, y) x %*% pars[[y]]$beta, mms[[i]], lhs(formulas[[i]]), SIMPLIFY = FALSE)
  # }
  copMM <- model.matrix(formulas[[4]][[1]], out)
  if (is.matrix(pars$cop$beta)) {
    if (nrow(pars$cop$beta) != ncol(copMM)) stop(paste0("dimension of model matrix for copula (", ncol(copMM), ") does not match number of coefficients provided (", nrow(pars$cop$beta),")"))
  }
  else if (is.atomic(pars$cop$beta)) {
    if (length(pars$cop$beta) != ncol(copMM)) stop(paste0("dimension of model matrix for copula (", ncol(copMM), ") does not match number of coefficients provided (", length(pars$cop$beta),")"))
  }

  # eta <- list()
  # eta$z <- mms[[1]] %*% pars2$z$beta
  # eta$y <- mms[[2]] %*% pars2$y$beta
  # mms[[3]] <- model.matrix(update.formula(formulas[[4]], NULL ~ . ), out)

  ## note that code will be slow if continuous covariates used in vine copula
  if (length(famCop) > 1) {
    if (nrow(unique(copMM)) > 25) warning("using vine copulas with continuous covariates may be very slow")
  }
  ## get copula data and then modify distributions of Y and Z
  out[,output] <- sim_copula(out[,output], family=famCop,
                             par = pars$cop, par2=pars$cop$par2, model_matrix=copMM)
  for (i in seq_along(LHS_Z)) {
    mms[[1]][[i]] <- model.matrix(formulas[[1]][[i]], data=out)
    out[[LHS_Z[i]]] <- rescale_var(out[[LHS_Z[i]]], X=mms[[1]][[i]],
                                                            family=famZ[[i]], pars=pars[[LHS_Z[i]]],
                                                            link=link[[1]][i])
  }
  for (i in seq_along(LHS_Y)) out[[LHS_Y[i]]] <- rescale_var(out[[LHS_Y[i]]], X=mms[[3]][[i]],
                                                            family=famY[[i]], pars=pars[[LHS_Y[i]]],
                                                            link=link[[3]][i])

  mms[[2]] = lapply(formulas[[2]], model.matrix, data=out)
  for (i in seq_along(mms[[2]])) {
    if (ncol(mms[[2]][[i]]) != length(pars[[LHS_X[i]]]$beta)) stop(paste0("dimension of model matrix for ", LHS_X[i], " does not match number of coefficients provided"))
  }

  ## perform rejection sampling
  wts <- rejectionWeights(out[LHS_X], mms[[2]], family=famX, pars=pars[LHS_X], qden = qden, link=link[[2]])
  con$max_wt <- max(max(wts), con$max_wt)
  wts <- wts/con$max_wt
  # if (mean(wts > 0.2) < 0.01) {
  #   if (con$warn == 1) warning("weights may be unbounded")
  #   else if (con$warn == 2) stop("weights may be unbounded")
  # }

  ## different behaviour depending upon whether covariates were supplied
  if (!datNULL) {
    done <- matrix(runif(nrow(out)) < wts, nrow=oversamp)
    wh2 <- apply(done, 2, function(x) which(x)[1])
    nr <- sum(colSums(done) > 0)

    while (any(!done)) {
      ## try again...
    }
  }
  else {
    kp <- runif(nrow(out)) < wts
    out2 <- out[kp, ]
    nr2 <- nrow(out2)

    if (nr2 < n) {
      ## not enough rows, so consider increasing oversampling rate
      if (con$oversamp == con$max_oversamp) {
        ## already at maximum oversample, so just return what we have
        warning(paste0("Only sampled ", nr2, " observations"))
      }
      else {
        con$oversamp <- min(con$max_oversamp, ceiling(con$oversamp*(n/nr2*1.1)))
        out2 <- Recall(n, formulas, pars=pars, family=family, link=link, dat=dat, control=con, seed=seed)
      }
    }
    else if (nr2 > n) {
      ## too many rows, so take a subsample
      out2 <- out2[seq_len(n), ]
    }
  }

  rownames(out2) <- NULL

  return(out2)
}
