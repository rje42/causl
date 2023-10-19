##' Process formulas, families and parameters
##'
##' @param formulas list of lists of formulas
##' @param pars list of lists of parameters
##' @param family families for Z,X,Y and copula
##' @param link list of link functions
##' @param kwd keyword for copula
##' @param ordering logical: should an ordering of variables be computed?
##'
process_inputs <- function (formulas, pars, family, link, kwd, ordering=FALSE) {

  # ## check list of formula
  # if (ordering && !is.numeric(family[[1]]) && ("formula" %in% class(formulas[[1]]))) {
  #   process_inputs2(formulas, pars, family, link, kwd)
  # }

  ## check we have four groups of formulas
  if (length(formulas) != 4) stop("formulas must have length 4")
  if (missing(pars)) stop("Must supply parameter values")

  ## ensure all formulas are given as lists
  if (any(sapply(formulas, class) == "formula")) {
    wh <- which(sapply(formulas, class) == "formula")
    for (i in wh) formulas[[i]] <- list(formulas[[i]])
  }

  dim <- lengths(formulas[1:3])

  if (missing(family)) {
    ## assume everything is Gaussian
    family = lapply(lengths(formulas), rep.int, x=1)
  }
  else if (is.list(family)) {
    if (!all(lengths(family[1:3]) == lengths(formulas[1:3]))) stop("mismatch in family and formulae specifications")
  }
  else if (length(family) == 4) family <- as.list(family)
  else stop("family should be a list, or vector of length 4")

  ## check that supplied parameters are sufficient
  nm_pars <- names(pars)
  wh_cop <- which(nm_pars == kwd)
  if (is.na(wh_cop)) stop("No parameters specified for copula")
  nms <- lapply(pars[-wh_cop], names)
  bpres <- sapply(nms, function(x) "beta" %in% x)
  if (!all(bpres)) {
    plur <- sum(!bpres) > 1
    stop(paste(ifelse(plur, "Variables", "Variable"), paste(names(pars)[!bpres], collapse=", "), ifelse(plur, "lack", "lacks"), "a beta parameter vector"))
  }

  # if (all(unlist(family) == 0)) {
  #   message("Perhaps better to simulate this using the MLLPs package")
  # }

  dZ <- dim[1]
  dX <- dim[2]
  dY <- dim[3]

  famZ <- family[[1]]
  famX <- family[[2]]
  famY <- family[[3]]
  famCop <- family[[4]]

  ## check variable names
  LHS_Z <- lhs(formulas[[1]])
  LHS_X <- lhs(formulas[[2]])
  LHS_Y <- lhs(formulas[[3]])
  rep_pars <- match(c(LHS_Z,LHS_X,LHS_Y), nm_pars, nomatch = 0L)
  if (any(rep_pars == 0)) {
    wh_nrep <- c(LHS_Z,LHS_X,LHS_Y)[rep_pars==0]
    stop(paste0("Variable ", paste(wh_nrep, collapse=", "), "not represented in the parameters list"))
  }

  formsZ <- lapply(formulas[[1]], terms)
  formsX <- lapply(formulas[[2]], terms)
  formsY <- lapply(formulas[[3]], terms)

  for (i in seq_along(formulas[[1]])) {
    npar <- length(attr(formsZ[[i]], "term.labels")) + attr(formsZ[[i]], "intercept")
    if (length(pars[[LHS_Z[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_Z[i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[2]])) {
    npar <- length(attr(formsX[[i]], "term.labels")) + attr(formsX[[i]], "intercept")
    if (length(pars[[LHS_X[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_X[i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[3]])) {
    npar <- length(attr(formsY[[i]], "term.labels")) + attr(formsY[[i]], "intercept")
    if (length(pars[[LHS_Y[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_Y[i], " does not match number of coefficients provided"))
  }

  ## check which variables are already in the formula for another, and set
  ## corresponding copula correlation to zero
  # formsZ <- merge_formulas(formulas[[1]])$old_forms
  # setmatch(LHS_Z, rhs_vars(formsZ))
  # formsY <- merge_formulas(formulas[[3]])$old_forms

  if (ordering) {
    ord_mat <- matrix(0, nrow=sum(dim), ncol=sum(dim))

    ## get adjacency matrix of dependencies
    vars <- c(LHS_Z, LHS_X, LHS_Y)
    ord_mat[seq_len(dZ), ] <- 1*t(sapply(formsZ, function(x) vars %in% attr(x, "term.labels")))
    ord_mat[dZ + seq_len(dX), ] <- 1*t(sapply(formsX, function(x) vars %in% attr(x, "term.labels")))
    ord_mat[dZ + dX + seq_len(dY), ] <- 1*t(sapply(formsY, function(x) vars %in% attr(x, "term.labels")))

    order <- topOrd(ord_mat)
    if (any(is.na(order))) stop("Formulae contain cyclic dependencies")
  }
  else {
    if (any(unlist(lapply(rhs_vars(formulas[[3]]),
                          function(x) any(x %in% LHS_Z))))) stop("Covariates cannot be direct predictors for outcomes")
    if (any(unlist(lapply(rhs_vars(formulas[[1]]),
                          function(x) any(x %in% LHS_Y))))) stop("Outcomes cannot be direct predictors for covariates")
    if (any(unlist(lapply(rhs_vars(formulas[[3]]),  # could relax this
                          function(x) any(x %in% LHS_Y))))) stop("Outcomes cannot directly predict one another")
    order <- NULL
  }

  #
  # if (any(sapply(formsZ,
  #        function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Y)))) stop("Outcomes cannot be direct predictors for covariates")
  # if (any(sapply(formsY,
  #                function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Z)))) stop("Covariates cannot be direct predictors for outcomes")
  # if (any(sapply(formsY,  # could weaken this?
  #                function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Y)))) stop("Outcomes cannot directly predict one another")

  ## check if covariates predict each other
  if (length(formulas[[1]]) > 1) {
    idx <- lapply(rhs_vars(formulas[[1]]), function(x) na.omit(match(x,LHS_Z)))
    # idx <- lapply(formsZ,
    #               function(x) match(intersect(LHS_Z, attr(x, "term.labels")[attr(x, "order") == 1]), LHS_Z))
    ignr <- do.call(cbind, mapply(function(x,y) {
      if (length(x) > 0) rbind(x,y)
      else matrix(NA, 2, 0)
    }, idx, seq_along(LHS_Z)))
    if (length(ignr) > 0) {
      wh_swt <- ignr[1,] > ignr[2,]
      ignr[,wh_swt] <- ignr[2:1, wh_swt]
      in_form <- setmatch(apply(ignr, 2, c, simplify = FALSE),
                          combn(length(LHS_Z), 2, simplify = FALSE))
    }
    else in_form <- numeric(0)
  }
  else in_form <- numeric(0)

  ## check that variables in families 1,2,3 have a dispersion parameter
  for (i in seq_along(famZ))
    if (famZ[i] <= 3 && is.null(pars[[LHS_Z[i]]]$phi))
      stop(paste(LHS_Z[i], "needs a phi parameter"))
  for (i in seq_along(famX))
    if (famX[i] <= 3 && is.null(pars[[LHS_X[i]]]$phi))
      stop(paste(LHS_X[i], "needs a phi parameter"))
  for (i in seq_along(famY))
    if (famY[i] <= 3 && is.null(pars[[LHS_Y[i]]]$phi))
      stop(paste(LHS_Y[i], "needs a phi parameter"))


  if (ordering) {

    ## put some validation code in here for inversion method

    ## get formulas in right format
    if (!is.list(formulas[[4]])) {
      formulas[[4]] <- rep(list(formulas[[4]]), dY)
    }
    if (!is.list(formulas[[4]][[1]])) {
      formulas[[4]] <- lapply(formulas[[4]], function(x) rep(list(x), dZ))
    }

    ## get parameters in right format
    if (!setequal(names(pars[[kwd]]), LHS_Y)) {
      if ("beta" %in% names(pars[[kwd]])) {
        pars[[kwd]] <- list(list(list(beta=pars[[kwd]]$beta)))
        names(pars[[kwd]]) <- LHS_Y
      }
      else {
        nrep <- setdiff(LHS_Y, names(pars[[kwd]]))
        if (length(nrep) > 0) stop(paste0("Variable ", paste(nrep, collapse=", "), "not represented in the copula parameters list"))
        rep <- setdiff(names(pars[[kwd]]), LHS_Y)
        if (length(rep) > 0) stop(paste0("Variable ", paste(nrep, collapse=", "), "represented in copula parameters list but not a response variable"))
        stop("Shouldn't get here")
      }
    }

    # ## TRY TO CODE UP CHECK THAT EVERY COMBINATION HAS A beta PARAMETER
    # if (any(!sapply(lapply(pars[[kwd]], names), "beta" %in% x))) {
    #   stop("")
    # }

    ## get family in right format
    if (!is.matrix(famCop)) {
      if (dY == 1) {
        family[[4]] <- list(unlist(family[[4]]))
        names(family[[4]]) <- LHS_Y
      }
      else if (is.matrix(family[[4]]) && nrow(family[[4]] == dY)) {
        family[[4]] <- apply(family[[4]], 1, c, simplify = FALSE)
      }
      else if (is.list(family[[4]]) && length(family[[4]] == dY)) {
        if (any(lengths(family[[4]]) != dZ + seq_len(dY) - 1)) stop("Incorrect format of family parameters for copula with inversion method")
      }
    }
  }
  else {
    ## ensure that copula parameters are in correct matrix format
    if (!is.matrix(pars[[kwd]]$beta)) {
      pars[[kwd]]$beta <- matrix(pars[[kwd]]$beta, ncol=choose(dZ+dY,2))
      form_cop <- terms(formulas[[4]][[1]])
      if (nrow(pars[[kwd]]$beta) != length(attr(form_cop, "term.labels"))+attr(form_cop, "intercept")) stop("Wrong number of regression parameters for copula")
    }
    if (length(in_form) > 0) {
      if (any(pars[[kwd]]$beta[,in_form] != 0)) {
        warning("Parameters in copula for objects already related by regression; setting copula parameters to zero")
        pars[[kwd]]$beta[,in_form] <- 0
      }
    }

    ## code to check if Y or Z is included in copula formula
    if (any(c(LHS_Z,LHS_Y) %in% rhs_vars(formulas[[4]])[[1]])) stop("copula cannot depend upon Z or Y variables")
  }

  ## set up link functions
  link <- link_setup(link, family[1:3], vars=list(LHS_Z,LHS_X,LHS_Y))

  ## other useful variables
  output <- c(LHS_Z, LHS_Y)
  vars <- c(LHS_Z, LHS_X, LHS_Y)

  return(list(formulas=formulas, pars=pars, family=family, link=link,
              LHSs=list(LHS_Z=LHS_Z, LHS_X=LHS_X, LHS_Y=LHS_Y), kwd=kwd,
              dim=dim, in_form=in_form, vars=vars, output=output,
              order=order))
}

# ##' Process inputs given in linear form
# process_inputs2(formulas, pars, family, link, kwd, ordering=FALSE) {
#   if (length(formulas))
# }

##' Sample from a causal model
##'
##' Obtain samples from a causal model using the rejection sampling approach of
##' Evans and Didelez (2023).
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
  out[,output] <- sim_CopVal(out[,output], family=famCop,
                             par = pars$cop, par2=pars$cop$par2, model_matrix=copMM)
  for (i in seq_along(LHS_Z)) {
    mms[[1]][[i]] <- model.matrix(formulas[[1]][[i]], data=out)
    out[[LHS_Z[i]]] <- rescaleVar(out[[LHS_Z[i]]], X=mms[[1]][[i]],
                                                            family=famZ[i], pars=pars[[LHS_Z[i]]],
                                                            link=link[[1]][i])
  }
  for (i in seq_along(LHS_Y)) out[[LHS_Y[i]]] <- rescaleVar(out[[LHS_Y[i]]], X=mms[[3]][[i]],
                                                            family=famY[i], pars=pars[[LHS_Y[i]]],
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
