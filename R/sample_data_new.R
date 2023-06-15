##' Sample from a causal model
##'
##' Obtain samples from a causal model using the rejection sampling approach of
##' Evans and Didelez (2021).
##'
##' @param n number of samples required
##' @param formulas list of lists of formulas
##' @param pars list of lists of parameters
##' @param family families for Z,X,Y and copula
##' @param link list of link functions
##' @param dat optional data frame of covariates
##' @param careful logical: should full rejection sampling method be used with
##' correctly computed weight?
##' @param method either \code{"rejection"} (the default) or \code{"inversion"}
##' @param control list of options for the algorithm
##' @param seed random seed used for replication
##'
##' @details Samples from a given causal model using rejection sampling (or,
##' if everything is discrete, direct sampling).
##'
##' The logical \code{careful} enables one to implement the full rejection
##' sampling method, which means we do get exact samples.  However, this method
##' may be slow, and in particular if we have an outlying value it may run very
##' slowly indeed.
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
##' \code{pars} should be a named list containing variable names that correspond
##' to the LHS of
##' formulae in \code{formulas}.  Each of these should themselves be a list
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
##' Control parameters are \code{trace} (default value 0, increasing to 1
##' increases verbosity of output), \code{warn} (which currently does nothing),
##' \code{max_wt} which is set to 1, and increases each time the function
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
rfrugalParam <- function(n, formulas = list(list(z ~ 1), list(x ~ z), list(y ~ x), list( ~ 1)),
                       pars, family, link=NULL, dat=NULL, careful=FALSE, method="rejection",
                       control=list(), seed) {

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

  ## process the four main arguments
  tmp <- process_inputs(formulas=formulas, pars=pars, family=family, link=link,
                        kwd=kwd, ordering = (method == "inversion"))
  formulas <- tmp$formulas
  pars <- tmp$pars
  family <- tmp$family
  link <- tmp$link
  dZ <- tmp$dim[1]; dX <- tmp$dim[2]; dY <- tmp$dim[3]
  LHS_Z <- tmp$LHSs$LHS_Z; LHS_X <- tmp$LHSs$LHS_X; LHS_Y <- tmp$LHSs$LHS_Y
  famZ <- tmp$family[[1]]; famX <- tmp$family[[2]]; famY <- tmp$family[[3]]; famCop <- tmp$family[[4]]

  output <- c(LHS_Z, LHS_Y)
  vars <- c(LHS_Z, LHS_X, LHS_Y)
  if (!datNULL) {
    vars <- c(names(dat), vars)
  }
  if (anyDuplicated(na.omit(vars))) stop("duplicated variable names")
  d <- length(vars)


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
  out <- data.frame(matrix(0, ncol=dZ+dX+dY, nrow=n))
  names(out) <- vars

  if (!datNULL) {
    out <- cbind(dat, out)
  }

  ## choose appropriate method
  if (method == "inversion") {  ## inversion method
    ## code to get causal order
    order <- tmp$order
    qZs <- out[LHS_Z]

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
        eta <- MM %*% pars[[vnm]]$beta
        curr_phi <- pars[[vnm]]$phi
        tmp <- sim_glm(fam=curr_fam, eta=eta, phi=curr_phi, link=curr_link)
        out[[vnm]] <- tmp
        if (vnm %in% LHS_Z) qZs[[vnm]] <- attr(tmp, "quantile")
      }
    }
  }
  else if (method == "rejection") { ## rejection sampling method
    if (careful) {
      ## get sample Z values
      Z0s <- gen_X_values(n, famX=famZ, pars=pars, LHS_X=LHS_Z, dX=dZ)$datX
      unb2_cts <- famZ %in% c(1,2)
      unb_cts <- famZ %in% c(3,6)
      b01 <- famZ == 4
      rg <- list()

      ## get range of bins for unbounded continuous variables
      for (i in which(unb2_cts)) {
        rg[[i]] <- range(Z0s[,i])
        rg[[i]][1] <- floor(rg[[i]][1])
        rg[[i]][2] <- ceiling(rg[[i]][2])
      }
      for (i in which(unb_cts)) {
        rg[[i]] <- c(0, ceiling(max(Z0s[,i])))
      }
      for (i in which(b01)) {
        rg[[i]] <- c(0,1)
      }

      ## then find constant needed over this space
      tmp <- gen_X_values(n, famX=famX, pars=pars, LHS_X=LHS_X, dX=dX, sim=FALSE)
      M <- get_max_weights(pars=pars, forms_X=formulas[[2]], fam_X = famX,
                           qden = tmp$qden, fam_Z=famZ, LHS_Z=LHS_Z, ranges=rg,
                           link=link[[2]])
    }

    ## ready for loop
    OK <- rep(FALSE, n)
    nr <- n

    while (nr > 0) {
      ## obtain treatment values
      tmp <- gen_X_values(nr, famX=famX, pars=pars, LHS_X=LHS_X, dX=dX)
      out[!OK,][LHS_X] <- tmp$datX[LHS_X]
      qden <- tmp$qden

      # ## give default coefficients
      # if (is.null(pars2$z$beta)) pars2$z$beta = 0
      # if (is.null(pars2$z$phi)) pars2$z$phi = 1

      ## get linear predictors
      mms <- vector(mode = "list", length=3)
      mms[c(1,3)] = rapply(formulas[c(1,3)], model.matrix, data=out[!OK,,drop=FALSE], how = "list")
      for (i in seq_along(mms[[1]])) {
        if (ncol(mms[[1]][[i]]) != length(pars[[LHS_Z[i]]]$beta)) stop(paste0("dimension of model matrix for ", LHS_Z[i], " does not match number of coefficients provided"))
      }
      for (i in seq_along(mms[[3]])) {
        if (ncol(mms[[3]][[i]]) != length(pars[[LHS_Y[i]]]$beta)) stop(paste0("dimension of model matrix for ", LHS_Y[i], " does not match number of coefficients provided"))
      }

      # etas <- vector(mode="list", length=3)
      # for (i in c(1,3)) {
      #   etas[[i]] <- mapply(function(x, y) x %*% pars[[y]]$beta, mms[[i]], lhs(formulas[[i]]), SIMPLIFY = FALSE)
      # }
      copMM <- model.matrix(formulas[[4]][[1]], out[!OK,,drop=FALSE])
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
      out[!OK,output] <- sim_CopVal(out[!OK,output,drop=FALSE], family=famCop,
                                    par = pars$cop, par2=pars$cop$par2, model_matrix=copMM)
      if (careful) OB <- rep(FALSE, nr)

      for (i in seq_along(LHS_Z)) {
        out[[LHS_Z[i]]][!OK] <- rescaleVar(out[[LHS_Z[i]]][!OK], X=mms[[1]][[i]],
                                           family=famZ[i], pars=pars[[LHS_Z[i]]],
                                           link=link[[1]][i])
        if (careful) {
          tmpZ <- out[[LHS_Z[i]]][!OK]
          if (famZ[i] <= 2 || famZ[i] == 4) OB <- OB | (tmpZ < rg[[i]][1]) | (tmpZ > rg[[i]][2])
          else if (famZ[i] %in% c(3,6)) OB <- OB | (tmpZ > rg[[i]][2])
        }
      }
      for (i in seq_along(LHS_Y)) out[[LHS_Y[i]]][!OK] <- rescaleVar(out[[LHS_Y[i]]][!OK], X=mms[[3]][[i]],
                                                                     family=famY[i], pars=pars[[LHS_Y[i]]],
                                                                     link=link[[3]][i])

      mms[[2]] = lapply(formulas[[2]], model.matrix, data=out[!OK,,drop=FALSE])
      for (i in seq_along(mms[[2]])) {
        if (ncol(mms[[2]][[i]]) != length(pars[[LHS_X[i]]]$beta)) stop(paste0("dimension of model matrix for ", LHS_X[i], " does not match number of coefficients provided"))
      }

      ## perform rejection sampling
      wts <- rejectionWeights(out[LHS_X][!OK,,drop=FALSE], mms[[2]], family=famX, pars=pars[LHS_X], qden = qden, link=link[[2]])
      if (careful) {
        wts[OB] <- 0  # for out of bounds values of Z
        wts <- wts/M
        if (any(wts > 1)) stop(paste("Weights ", paste(which(wts > 1), collapse=", "), " are > 1", sep=""))
      }
      else {
        con$max_wt <- max(max(wts), con$max_wt)
        wts <- wts/con$max_wt
      }

      OK[!OK] <- runif(nr) < wts
      nr <- sum(!OK)
    }
  }

  rownames(out) <- NULL

  return(out)
}

