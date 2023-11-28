##' Get weights for rejection sampling
##'
##' @param dat data frame of variables to change conditional distribution of
##' @param mms list of model matrices
##' @param family vector of distribution families
##' @param pars parameters for new distributions
##' @param qden functions for densities used to simulate variables
##' @param link link functions for GLMs
##'
##' @return a numeric vector of weights
##'
##' @export
rejectionWeights <- function (dat, mms,# formula,
                              family, pars, qden, link) {

  d <- ncol(dat)

  if (d != length(mms)) stop("Inconsistent length of dat and mms")
  if (d != length(family)) stop("Inconsistent length of dat and family")
  if (d != length(pars)) stop("Inconsistent length of dat and pars")
  if (d != length(qden)) stop("Inconsistent length of dat and family")

  betas <- lapply(pars, function(x) x$beta)
  eta <- mapply(function(X,y) X %*% y, mms, betas, SIMPLIFY = FALSE)

  nms <- names(dat)

  ## collect phi and par2 parameters
  phi <- par2 <- numeric(d)
  for (i in which(family %in% c(1:3,6))) {
    phi[i] <- pars[[nms[i]]]$phi
    if (family[i] == 2) par2[i] <- pars[[nms[i]]]$par2
  }

  wts <- get_X_density(dat, eta=eta, phi=phi, qden=qden, family=family,
                       link=link, par2=par2)

  if (any(is.na(wts))) stop("Problem with weights")

  wts
}


##' @param careful should full, slower method be used?
##' @describeIn sim_inversion Rejection sampling code
sim_rejection <- function (out, proc_inputs, careful, control) {

  n <- nrow(out)

  ## unpack proc_inputs
  formulas <- proc_inputs$formulas
  pars <- proc_inputs$pars
  family <- proc_inputs$family
  link <- proc_inputs$link
  dZ <- proc_inputs$dim[1]; dX <- proc_inputs$dim[2]; dY <- proc_inputs$dim[3]
  LHS_Z <- proc_inputs$LHSs$LHS_Z; LHS_X <- proc_inputs$LHSs$LHS_X; LHS_Y <- proc_inputs$LHSs$LHS_Y
  famZ <- proc_inputs$family[[1]]; famX <- proc_inputs$family[[2]]; famY <- proc_inputs$family[[3]]; famCop <- proc_inputs$family[[4]]

  output <- c(LHS_Z, LHS_Y)
  vars <- names(out)

  if (careful) {
    if (length(famZ) == 0) stop("Should have at least one Z variable")

    ## get sample Z values
    Z0s <- gen_X_values(n, famX=famZ, pars=pars, LHS_X=LHS_Z, dX=dZ)$datX
    if (is(famZ[[1]], "causl_family")) {
      famChr <- sapply(famZ, function(x) x$name)
      unb2_cts <- famChr %in% c("gaussian","t")
      unb_cts <- famChr %in% c("Gamma","lnormal")
      b01 <- famChr == "beta"
    }
    else {
      unb2_cts <- famZ %in% c(1,2)
      unb_cts <- famZ %in% c(3,6)
      b01 <- famZ == 4
    }
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
      mms[[1]][[i]] <- model.matrix(formulas[[1]][[i]], data=out[!OK,,drop=FALSE])
      out[[LHS_Z[i]]][!OK] <- rescaleVar(out[[LHS_Z[i]]][!OK], X=mms[[1]][[i]],
                                         family=famZ[[i]], pars=pars[[LHS_Z[i]]],
                                         link=link[[1]][i])
      if (careful) {
        tmpZ <- out[[LHS_Z[i]]][!OK]
        if (famZ[i] <= 2 || famZ[i] == 4) OB <- OB | (tmpZ < rg[[i]][1]) | (tmpZ > rg[[i]][2])
        else if (famZ[i] %in% c(3,6)) OB <- OB | (tmpZ > rg[[i]][2])
      }
    }
    for (i in seq_along(LHS_Y)) {
      mms[[3]][[i]] <- model.matrix(formulas[[3]][[i]], data=out[!OK,,drop=FALSE])
      out[[LHS_Y[i]]][!OK] <- rescaleVar(out[[LHS_Y[i]]][!OK], X=mms[[3]][[i]],
                                         family=famY[[i]], pars=pars[[LHS_Y[i]]],
                                         link=link[[3]][i])
    }

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
      control$max_wt <- max(max(wts), control$max_wt)
      wts <- wts/control$max_wt
    }

    OK[!OK] <- runif(nr) < wts
    nr <- sum(!OK)
  }

  return(out)
}
