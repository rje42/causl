##' Get maximum weight for each segment of a distribution
##'
##' @param pars list with all regression parameters
##' @param forms_X formulae for treatments
##' @param fam_X,fam_Z vector of families for treatments and covariates
##' @param LHS_Z variables in covariates
##' @param qden density of proposals
##' @param ranges range of segments
##' @param link link functions for treatments
##'
get_max_weights <- function (pars, forms_X, fam_X, qden, fam_Z, LHS_Z, ranges, link, ...) {
  # if (!missing(link)) warning("Link argument is ignored")

  ## get unique list of segments
  # m_wts <- unique(Z0q)
  # m_wts$max_wt <- NA

  LHS_X <- lhs(forms_X)
  dX <- length(LHS_X)
  full_form <- merge_formulas(forms_X)

  ## get states for discrete treatments
  disX <- sum(fam_X == 5)
  # disZ <- sum(fam_Z == 5)
  combs <- expand.grid(rep(list(0:1), disX))
  names(combs) <- names(LHS_X[fam_X == 5])

  ## build dummy data frame
  Zvars <- setdiff(unique(unlist(rhs_vars(forms_X))), LHS_X)
  fam_Z <- fam_Z[match(LHS_Z, Zvars, nomatch = 0L)]  # get only relevant Z variables
  nms <- c(LHS_X, Zvars)
  args <- rep(list(NA), length(nms))
  names(args) <- nms
  out <- do.call(data.frame, args)


  ## get matrix form for parameters
  msks <- initializeParams2(formulas=forms_X, family=fam_X, full_form=full_form,
                            only_masks=TRUE, nc=0, inc_cop=FALSE)
  beta <- beta_m <- msks$beta_m
  phi <- phi_m <- msks$phi_m

  for (i in seq_len(ncol(beta_m))) {
    beta[full_form$wh[[i]],i] <- pars[[LHS_X[i]]]$beta
    if (fam_X[i] %in% c(1:3,6)) phi[i] <- pars[[LHS_X[i]]]$phi
  }
  if (any(beta != 0 & beta_m == 0) ||
      any(phi != 0 & phi_m == 0)) stop("Mask problem")

  disc <- c(fam_X==5, fam_Z==5)

  func <- function (x, dis_val=integer(0)) {
    ## put variables in data frame
    ZX <- out
    ZX[1, !disc] <- x
    ZX[1, disc] <- dis_val
    ZX2 <- ZX[,seq_along(fam_X),drop=FALSE]
    # X <- x[seq_len(fam_X)]
    # Z <- x[-seq_len(fam_X)]
    # beta

    mm <- model.matrix(full_form$formula, ZX)
    eta <- mm %*% beta

    out <- get_X_density(dat = ZX2, eta = eta, phi = phi, par2=0*phi,
                         qden = qden, family=fam_X, link = link)

    -out
  }

  if (length(ranges) >= 1) rg <- lapply(purrr::transpose(ranges), unlist)
  else if (length(ranges) == 1) rg <- list(ranges[[1]][1],ranges[[1]][2])
  else rg <- list()

  ## need to make move in several directions
  best <- list(value=func(rep(0.1,dX+length(fam_Z)-sum(disc)),rep(0,sum(disc))))
  n_strts <- 10
  for (i in seq_len(n_strts)) {
    strt <- runif(length(ranges), min = rg[[1]], max = rg[[2]])
    strt <- c(rnorm(dX-disX)/10, strt)
    if (nrow(combs) > 0) {
      for (dx in seq_len(nrow(combs))) {
        opt <- tryCatch(optim(par=strt, fn=func, dis_val = combs[dx,], method = "L-BFGS-B",
                              lower=c(rep(-Inf, dX-disX),rg[[1]]), upper=c(rep(Inf, dX-disX),rg[[2]])),
                        error=function(e) return(list(value=Inf)))
        if (opt$value < best$value) best <- opt
      }
    }
    else {
      opt <- tryCatch(optim(par=strt, fn=func, dis_val = integer(0), method = "L-BFGS-B",
                            lower=c(rep(-Inf, dX),rg[[1]]), upper=c(rep(Inf, dX),rg[[2]])),
                      error=function(e) return(list(value=Inf)))
      if (opt$value < best$value) best <- opt
    }

  }
  # print(-best$value)
  # print(ranges)
  return(-best$value)

  # ### get maximum weights in each segment
  # for (i in m_wts) {
  #   rg
  #   opt <- optim(par)
  #
  # }
}

# ##' Sample from a causal model
# ##'
# ##' Obtain samples from a causal model using the rejection sampling approach of
# ##' Evans and Didelez (2021).
# ##'
# ##' @param n number of samples required
# ##' @param formulas list of lists of formulas
# ##' @param pars list of lists of parameters
# ##' @param data optional data frame of covariates
# ##' @param family families for Z,X,Y and copula
# ##' @param link list of link functions
# ##' @param dat data frame of covariates
# ##' @param method only \code{"rejection"} is valid
# ##' @param control list of options for the algorithm
# ##' @param seed random seed used for replication
# ##'
# ##' @details Samples from a given causal model using rejection sampling (or,
# ##' if everything is discrete, direct sampling).
# ##'
# ##' The entries for  \code{formula} and \code{family} should each be a
# ##' list with four entries, corresponding to the \eqn{Z}, \eqn{X}, \eqn{Y} and
# ##' the copula.  \code{formula} determines the model, so it is crucial that
# ##' every variable to be simulated is represented there exactly once.  Each
# ##' entry of that list can either be a single formula, or a list of formulae.
# ##' Each corresponding entry in \code{family} should be the same length as the
# ##' list in \code{formula} or of length 1 (in which case it will be repeated
# ##' for all the variables therein).
# ##'
# ##' We use the following codes for different families of distributions:
# ##' 0 or 5 = binary;
# ##' 1 = normal;
# ##' 2 = t-distribution;
# ##' 3 = gamma;
# ##' 4 = beta;
# ##' 6 = log-normal.
# ##'
# ##' The family variables for the copula are also numeric and taken from
# ##' \code{VineCopula}.
# ##' Use, for example, 1 for Gaussian, 2 for t, 3 for Clayton, 4 for Gumbel,
# ##' 5 for Frank, 6 for Joe and 11 for FGM copulas.
# ##'
# ##' \code{pars} should be a named list containing: either entries \code{z},
# ##' \code{x}, \code{y} and \code{cop}, or variable names that correspond to the LHS of
# ##' formulae in \code{formulas}.  Each of these should themselves be a list
# ##' containing \code{beta} (a vector of regression parameters) and (possibly)
# ##' \code{phi}, a dispersion parameter.  For any discrete variable that is a
# ##' treatment, you can also specify \code{p}, an initial proportion to simulate
# ##' from (otherwise this defaults to 0.5).
# ##'
# ##' Link functions for the Gaussian, t and Gamma distributions can be the
# ##' identity, inverse or log functions.  Gaussian and t-distributions default to
# ##' the identity, and Gamma to the log link.  For the Bernoulli the logit and
# ##' probit links are available.
# ##'
# ##' Control parameters are \code{trace} (default value 0, increasing to 1
# ##' increases verbosity of output), \code{warn} (which currently does nothing),
# ##' \code{max_wt} which is set to 1, and increases each time the function
# ##' is recalled.
# ## (if weights empirically appear not to have an upper bound, this warns if set
# ## to 1 (the default) and stops if set to 2), ...
# ##' Control parameters also include \code{cop}, which gives a keyword for the
# ##' copula that defaults to \code{"cop"}.
# ##'
# ##' @examples
# ##' pars <- list(z=list(beta=0, phi=1),
# ##'              x=list(beta=c(0,0.5), phi=1),
# ##'              y=list(beta=c(0,0.5), phi=0.5),
# ##'              cop=list(beta=1))
# ##' causalSamp(100, pars = pars)
# ##'
# ## @importFrom frugalSim sim_chain
# ##'
# ##' @return A data frame containing the simulated data.
# ##'
# ##' @importFrom dplyr ntile
# ##' @export
# rfrugalParam2 <- function(n, formulas = list(list(z ~ 1), list(x ~ z), list(y ~ x), list( ~ 1)),
#                           pars, family, link=NULL, dat=NULL, method="rejection",
#                           control=list(), seed) {
#
#   # get control parameters or use defaults
#   con = list(max_wt = 1, warn = 1, cop="cop", trace=0)
#   matches = match(names(control), names(con))
#   con[matches] = control[!is.na(matches)]
#   if (any(is.na(matches))) warning("Some names in control not matched: ",
#                                    paste(names(control[is.na(matches)]),
#                                          sep = ", "))
#   ## get keyword for copula formula
#   kwd <- con$cop
#
#   ## infer sample size from dat if possible
#   datNULL <- is.null(dat)
#   if (datNULL && missing(n)) stop("Must specify a sample size or supply a data frame")
#   else if (missing(n)) n <- nrow(dat)
#   else if (!datNULL && n != nrow(dat)) {
#     warning("Dimension of supplied data does not match n, using size of supplied data frame")
#     n <- nrow(dat)
#   }
#
#   ## process the four main arguments
#   tmp <- process_inputs(formulas=formulas, pars=pars, family=family, link=link,
#                         kwd=kwd)
#   formulas <- tmp$formulas
#   pars <- tmp$pars
#   family <- tmp$family
#   link <- tmp$link
#   dZ <- tmp$dim[1]; dX <- tmp$dim[2]; dY <- tmp$dim[3]
#   LHS_Z <- tmp$LHSs$LHS_Z; LHS_X <- tmp$LHSs$LHS_X; LHS_Y <- tmp$LHSs$LHS_Y
#   famZ <- tmp$family[[1]]; famX <- tmp$family[[2]]; famY <- tmp$family[[3]]; famCop <- tmp$family[[4]]
#
#   output <- c(LHS_Z, LHS_Y)
#   vars <- c(LHS_Z, LHS_X, LHS_Y)
#   if (!datNULL) {
#     vars <- c(names(dat), vars)
#   }
#   if (anyDuplicated(na.omit(vars))) stop("duplicated variable names")
#   d <- length(vars)
#
#
#   ## set up output
#   out <- data.frame(matrix(0, ncol=dZ+dX+dY, nrow=n))
#   names(out) <- c(LHS_Z, LHS_X, LHS_Y)
#
#   if (!datNULL) {
#     out <- cbind(dat, out)
#   }
#
#   ## get sample Z values
#   Z0s <- gen_X_values(n, famX=famZ, pars=pars, LHS_X=LHS_Z, dX=dZ)$datX
#   unb2_cts <- famZ %in% c(1,2)
#   unb_cts <- famZ %in% c(3,6)
#   b01 <- famZ == 4
#   rg <- list()
#
#   ## get range of bins for unbounded continuous variables
#   for (i in which(unb2_cts)) {
#     rg[[i]] <- range(Z0s[,i])
#     rg[[i]][1] <- floor(rg[[i]][1])
#     rg[[i]][2] <- ceiling(rg[[i]][2])
#   }
#   for (i in which(unb_cts)) {
#     rg[[i]] <- c(0, ceiling(max(Z0s[,i])))
#   }
#   for (i in which(b01)) {
#     rg[[i]] <- c(0,1)
#   }
#
#   ## then find constant needed over this space
#   tmp <- gen_X_values(n, famX=famX, pars=pars, LHS_X=LHS_X, dX=dX, sim=FALSE)
#   M <- get_max_weights(pars=pars, forms_X=formulas[[2]], fam_X = famX,
#                        qden = tmp$qden, fam_Z=famZ, LHS_Z=LHS_Z, ranges=rg,
#                        link=link[[2]])
#
#   ## ready for loop
#   OK <- rep(FALSE, n)
#   nr <- n
#
#   while (nr > 0) {
#     ## obtain treatment values
#     tmp <- gen_X_values(nr, famX=famX, pars=pars, LHS_X=LHS_X, dX=dX)
#     out[!OK,][LHS_X] <- tmp$datX[LHS_X]
#     qden <- tmp$qden
#
#     # ## give default coefficients
#     # if (is.null(pars2$z$beta)) pars2$z$beta = 0
#     # if (is.null(pars2$z$phi)) pars2$z$phi = 1
#
#     ## get linear predictors
#     mms <- vector(mode = "list", length=3)
#     mms[c(1,3)] = rapply(formulas[c(1,3)], model.matrix, data=out[!OK,,drop=FALSE], how = "list")
#     for (i in seq_along(mms[[1]])) {
#       if (ncol(mms[[1]][[i]]) != length(pars[[LHS_Z[i]]]$beta)) stop(paste0("dimension of model matrix for ", LHS_Z[i], " does not match number of coefficients provided"))
#     }
#     for (i in seq_along(mms[[3]])) {
#       if (ncol(mms[[3]][[i]]) != length(pars[[LHS_Y[i]]]$beta)) stop(paste0("dimension of model matrix for ", LHS_Y[i], " does not match number of coefficients provided"))
#     }
#
#     # etas <- vector(mode="list", length=3)
#     # for (i in c(1,3)) {
#     #   etas[[i]] <- mapply(function(x, y) x %*% pars[[y]]$beta, mms[[i]], lhs(formulas[[i]]), SIMPLIFY = FALSE)
#     # }
#     copMM <- model.matrix(formulas[[4]][[1]], out[!OK,,drop=FALSE])
#     if (is.matrix(pars$cop$beta)) {
#       if (nrow(pars$cop$beta) != ncol(copMM)) stop(paste0("dimension of model matrix for copula (", ncol(copMM), ") does not match number of coefficients provided (", nrow(pars$cop$beta),")"))
#     }
#     else if (is.atomic(pars$cop$beta)) {
#       if (length(pars$cop$beta) != ncol(copMM)) stop(paste0("dimension of model matrix for copula (", ncol(copMM), ") does not match number of coefficients provided (", length(pars$cop$beta),")"))
#     }
#
#     # eta <- list()
#     # eta$z <- mms[[1]] %*% pars2$z$beta
#     # eta$y <- mms[[2]] %*% pars2$y$beta
#     # mms[[3]] <- model.matrix(update.formula(formulas[[4]], NULL ~ . ), out)
#
#     ## note that code will be slow if continuous covariates used in vine copula
#     if (length(famCop) > 1) {
#       if (nrow(unique(copMM)) > 25) warning("using vine copulas with continuous covariates may be very slow")
#     }
#     ## get copula data and then modify distributions of Y and Z
#     out[!OK,output] <- sim_CopVal(out[!OK,output,drop=FALSE], family=famCop,
#                                   par = pars$cop, par2=pars$cop$par2, model_matrix=copMM)
#     OB <- rep(FALSE, nr)
#     for (i in seq_along(LHS_Z)) {
#       out[[LHS_Z[i]]][!OK] <- rescaleVar(out[[LHS_Z[i]]][!OK], X=mms[[1]][[i]],
#                                          family=famZ[i], pars=pars[[LHS_Z[i]]],
#                                          link=link[[1]][i])
#       tmpZ <- out[[LHS_Z[i]]][!OK]
#       if (famZ[i] <= 2) OB <- OB | (tmpZ < rg[[i]][1]) | (tmpZ > rg[[i]][2])
#       else if (famZ[i] %in% c(3,6)) OB <- OB | (tmpZ > rg[[i]][2])
#
#     }
#     for (i in seq_along(LHS_Y)) out[[LHS_Y[i]]][!OK] <- rescaleVar(out[[LHS_Y[i]]][!OK], X=mms[[3]][[i]],
#                                                                    family=famY[i], pars=pars[[LHS_Y[i]]],
#                                                                    link=link[[3]][i])
#
#
#     mms[[2]] = lapply(formulas[[2]], model.matrix, data=out[!OK,,drop=FALSE])
#     for (i in seq_along(mms[[2]])) {
#       if (ncol(mms[[2]][[i]]) != length(pars[[LHS_X[i]]]$beta)) stop(paste0("dimension of model matrix for ", LHS_X[i], " does not match number of coefficients provided"))
#     }
#
#     ## perform rejection sampling
#     wts <- rejectionWeights(out[LHS_X][!OK,,drop=FALSE], mms[[2]], family=famX, pars=pars[LHS_X], qden = qden, link=link[[2]])
#     wts[OB] <- 0  # for out of bounds values of Z
#     # con$max_wt <- max(max(wts), con$max_wt)
#     wts <- wts/M
#     if (any(wts > 1)) stop(paste("Weights ", paste(which(wts > 1), collapse=", "), " are > 1", sep=""))
#
#     OK[!OK] <- runif(nr) < wts
#     nr <- sum(!OK)
#   }
#
#   rownames(out) <- NULL
#
#   return(out)
# }

