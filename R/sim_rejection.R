##' Simulate initial X values
##'
##' @param n number of observations
##' @param fam_x number for distribution family
##' @param theta parameters for model
##' @param offset optional mean correction
##' @param sim should variables be simulated?
##'
##' @details Returns a list that includes a data frame containing a column
##' `x`, as well as the density that was used to generate it.  Possible
##' families are Gaussian (=1), t (=2), Exponential (=3), beta (=4)
##' Bernoulli/categorical (=5) and log-normal (=6).
##'
##' For the exponential distribution, `theta` is the mean.
##' Beta can take one or two parameters, and if there is only
##' one it is just repeated.
##'
##' The `offset` parameter alters the median for the normal and t-distributions,
##' or the median of the logarithm in the case of a log-normal.
##'
##' @return A list with two entries: `x` a vector of the simulated
##' values, and `qden`, which contains a function that evaluates to the
##' density of the distribution used to generate those values.
##'
##' @export
sim_X <- function(n, fam_x, theta, offset, sim=TRUE) {

  if (missing(offset)) offset <- 0

  if (is.numeric(fam_x)) {
    if (fam_x == 1) {
      if (sim) X <- rnorm(n, mean=offset, sd=sqrt(theta))
      qden <- function(x) dnorm(x, mean=offset, sd=sqrt(theta))
    }
    else if (fam_x == 6) {
      if (sim) X <- exp(rnorm(n, mean=offset, sd=sqrt(theta)))
      qden <- function(x) dnorm(log(x), mean=offset, sd=sqrt(theta))/x
    }
    else if (fam_x == 2) {
      if (sim) X <- sqrt(theta[1])*(rt(n, df=theta[2]) + offset)
      qden <- function(x) dt((x-offset)/sqrt(theta[1]), df=theta[2])/sqrt(theta[1])
    }
    else if (fam_x == 3) {
      if (sim) X <- rgamma(n, shape=1, rate=1/theta)
      qden <- function(x) dgamma(x, shape=1, rate=1/theta)
    }
    else if (fam_x == 4) {
      if (length(theta) == 1) theta <- c(theta, theta)
      if (sim) X <- rbeta(n, theta[1], theta[2])
      qden <- function(x) dbeta(x, theta[1], theta[2])
    }
    else if (fam_x == 5) {
      if (sum(theta) > 1) stop("Probabilities must be < 1")
      if (length(theta) == 1) theta <- c(theta, 1-theta)
      if (any(theta < 0)) stop("Negative parameters not allowed")
      if (!isTRUE(all.equal(sum(theta), 1))) {
        warning("Parameters do not sum to 1, rescaling")
        theta <- theta/sum(theta)
      }
      if (sim) X <- sample(length(theta), size=n, replace=TRUE, prob=theta)-1
      qden <- function(x) theta[x+1]
    }
    else if (fam_x == 11) {
      p <- c(1,theta)
      for (i in seq_along(theta)[-1]) p[i+1] <- p[i+1]/p[i]
      p <- p/sum(p)

      if (sim) X <- sample(length(p), size=n, replace=TRUE, prob=p)
      qden <- function(x) p[x]
    }
    else stop("X distribution must be normal (1), t (2), Gamma (3), Beta (4) or categorical (5)")
  }
  else if (is(fam_x, "causl_family")) {

    pars <- fam_x$default(theta)[fam_x$pars]
    # if (fam_x$name == "ordinal") fam_x

    if (sim) {
      X <- do.call(fam_x$rdist, c(list(n=n), pars))
    }
    qden <- function (x) do.call(fam_x$ddist, c(list(x=x), pars))
  }

  if (sim) out <- list(x=X, qden=qden)
  else out <- list(qden=qden)

  return(out)
}


gen_X_values <- function (n, famX, pars, LHS_X, dX, sim=TRUE) {
  ## list for densities used to simulate X's
  if (missing(dX)) dX <- length(LHS_X)
  qden <- vector(mode="list", length=dX)
  out <- data.frame(matrix(0, ncol=dX, nrow=n))
  names(out) <- LHS_X

  if (is.numeric(famX)) {
    ## use old style approach
    for (i in seq_len(dX)) {
      ## get parameters for X
      if (famX[i] == 0 || famX[i] == 5) {
        if (!is.null(pars[[LHS_X[i]]][["p"]])) theta <- pars[[LHS_X[i]]]$p
        else theta <- 0.5
        famX[i] <- 5
      }
      else if (famX[i] == 1 || famX[i] == 2) {
        theta <- 2*pars[[LHS_X[i]]]$phi
      }
      else if (famX[i] == 6) {
        theta <- 1.5*pars[[LHS_X[i]]]$phi
      }
      else if (famX[i] == 3) {
        theta <- 2*pars[[LHS_X[i]]]$phi
      }
      else if (famX[i] == 4) {
        theta = c(1,1)
      }
      else if (famX[i] == 11) {
        nl <- pars[[LHS_X[i]]]$nlevel
        theta <- rep(1, nl-1)
      }
      else stop(paste0("Family ", famX[i], " not recognized"))

      ## obtain data for X's
      tmp <- sim_X(n, fam_x = famX[i], theta=theta, sim=sim)
      out[LHS_X[[i]]] <- tmp$x
      qden[[i]] <- tmp$qden
    }
  }
  else if (is.list(famX) && length(famX) > 0 && is(famX[[1]], "causl_family")) {
    for (i in seq_len(dX)) {
      fam <- famX[[i]]$name
      ## get parameters for X
      if (fam == "binomial") {
        if (!is.null(pars[[LHS_X[i]]][["p"]])) theta <- pars[[LHS_X[i]]]$p
        else theta <- 0.5
      }
      else if (fam == "gaussian" || fam == "t") {
        theta <- 2*pars[[LHS_X[i]]]$phi
      }
      else if (fam == "lognormal") {
        theta <- 1.5*pars[[LHS_X[i]]]$phi
      }
      else if (fam == "Gamma") {
        theta <- 2*pars[[LHS_X[i]]]$phi
      }
      else if (fam == "beta") {
        theta = c(1,1)
      }
      else if (fam == "ordinal") {
        nl <- pars[[LHS_X[i]]]$nlevel
        theta <- rep(1, nl-1)
      }
      else stop(paste0("Family for entry ", i, " not recognized"))

      ## obtain data for X's
      tmp <- sim_X(n, fam_x = famX[[i]], theta=theta, sim=sim)
      out[LHS_X[[i]]] <- tmp$x
      qden[[i]] <- tmp$qden
    }
  }
  else stop("famX should be a causl_family list or integer vector")

  return(list(datX = out, qden = qden))
}

##' Get density of treatments
##'
##' @inherit rejectionWeights
##' @param eta list (or matrix) of linear forms
##' @param phi vector of dispersion coefficients
## @param par2 vector of degrees of freedom
##' @param other_par other parameters for family
##' @param log logical: should log-density be returned?
##'
get_X_density <- function (dat, eta, phi, qden, family, link, other_par,
                           log = FALSE) {

  wts <- rep(1, nrow(dat))

  if (is.matrix(eta)) eta <- apply(eta, 2, function(x) x, simplify = FALSE)

  if (length(family) > 0 && is(family[[1]], "causl_family")) {
    for (i in seq_len(ncol(dat))) {
      mu <- link_apply(eta[[i]], link[i], family_nm = family[[i]]$name)

      pars <- list(mu = mu)
      if ("phi" %in% family[[i]]$pars) pars <- c(pars, list(phi=phi))
      pars <- c(pars, other_par[[i]])
      # if ("df" %in% family[[i]]$pars) pars <- c(pars, list(df=other_par[[i]]$df))

      wts <- wts*do.call(family[[i]]$ddist, c(list(x=dat[,i]), pars))
      if (is.numeric(qden[[i]])) wts <- wts/qden[[i]]
      else wts <- wts/qden[[i]](dat[,i])
    }
  }
  else if (is.numeric(family)) {
    ## take the old-style approach
    for (i in seq_len(ncol(dat))) {
      if (family[i] == 1) {
        if (link[i] == "identity") mu <- eta[[i]]
        else if (link[i] == "log") mu <- exp(eta[[i]])
        else if (link[i] == "inverse") mu <- 1/eta[[i]]
        else stop("not a valid link function for Gaussian distribution")

        fam <- get_family(i)
        dens <- fam()$ddist

        if (is.numeric(qden[[i]])) wts <- wts*dens(dat[,i], mu=mu, phi=phi[i])/qden[[i]]
        else wts <- wts*dens(dat[,i], mu=mu, phi=phi[i])/qden[[i]](dat[,i])

        # if (is.numeric(qden[[i]])) wts <- wts*dnorm(dat[,i], mean=mu, sd=sqrt(phi[i]))/qden[[i]]
        # else wts <- wts*dnorm(dat[,i], mean=mu, sd=sqrt(phi[i]))/qden[[i]](dat[,i])
      }
      else if (family[i] == 2) {
        if (link[i] == "identity") mu <- eta[[i]]
        else if (link[i] == "log") mu <- exp(eta[[i]])
        else if (link[i] == "inverse") mu <- 1/eta[[i]]
        else stop("not a valid link function for t-distribution")

        fam <- get_family(i)
        dens <- fam()$ddist

        if (is.numeric(qden[[i]])) wts <- wts*dens(dat[,i], mu=mu, phi=phi[i], df=pars[[i]]$df)/qden[[i]]
        else wts <- wts*dens(dat[,i], mu=mu, phi=phi[i], df=pars[[i]]$df)/qden[[i]](dat[,i])

        # if (is.numeric(qden[[i]])) wts <- wts*dt((dat[,i] - mu)/sqrt(phi[i]), df=pars[[i]]$df)/(sqrt(phi[i])*qden[[i]])
        # else wts <- wts*dt((dat[,i] - mu)/sqrt(phi[i]), df=pars[[i]]$df)/(sqrt(phi[i])*qden[[i]](dat[,i]))
      }
      else if (family[i] == 3) {
        if (link[i] == "identity") mu <- eta[[i]]
        else if (link[i] == "log") mu <- exp(eta[[i]])
        else if (link[i] == "inverse") mu <- 1/eta[[i]]
        else stop("not a valid link function for Gamma distribution")

        fam <- get_family(i)
        dens <- fam()$ddist

        if (is.numeric(qden[[i]])) wts <- wts*dens(dat[,i], mu=mu, phi=phi[i])/qden[[i]]
        else wts <- wts*dens(dat[,i], mu=mu, phi=phi[i])/qden[[i]](dat[,i])
        # if (is.numeric(qden[[i]])) wts <- wts*dgamma(dat[,i], rate=1/(mu*phi[i]), shape=1/phi[i])/qden[[i]]
        # else wts <- wts*dgamma(dat[,i], rate=1/(mu*phi[i]), shape=1/phi[i])/qden[[i]](dat[,i])
      }
      else if (family[i] == 4) {
        mu <- expit(eta[[i]])
        if (is.numeric(qden[[i]])) wts <- wts*dbeta(dat[,i], shape1=1+phi[i]*mu, shape2=1+phi[i]*(1-mu))/qden[[i]]
        else wts <- wts*dbeta(dat[,i], shape1=1+phi[i]*mu, shape2=1+phi[i]*(1-mu))/qden[[i]](dat[,i])
      }
      else if (family[i] == 5) {
        if (link[i] == "probit") mu <- qnorm(eta[[i]])
        else if (link[i] == "logit") mu <- expit(eta[[i]])
        else stop("not a valid link function for binomial distribution")

        if (is.numeric(qden[[i]])) wts <- wts*dbinom(dat[,i], prob=mu, size=1)/qden[[i]]
        else wts <- wts*dbinom(dat[,i], prob=mu, size=1)/qden[[i]](dat[,i])
      }
      else if (family[i] == 6) {
        if (link[i] == "identity") lmu <- eta[[i]]
        else if (link[i] == "identity") lmu <- log(eta[[i]] - phi/2)
        else stop("not a valid link function for log-normal distribution")
        if (is.numeric(qden[[i]])) wts <- wts*dnorm(log(dat[,i]), mean=lmu, sd=sqrt(phi[i]))/(dat[,i]*qden[[i]])
        else wts <- wts*dnorm(log(dat[,i]), mean=lmu, sd=sqrt(phi[i]))/(dat[,i]*qden[[i]](dat[,i]))
      }
      else stop("family[[2]] must be in the range 1 to 6")
    }
  }

  if (log) wts <- log(wts)

  wts
}


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

  ## collect phi and any other parameters
  phi <- numeric(d)
  other_par <- vector(mode = "list", length=d)

  if (is.numeric(family)) {
    for (i in seq_len(d)) {
      whB <- match("beta", names(pars[[nms[i]]]), nomatch = 0L)
      if (family %in% c(1:3,6)) {
        phi[i] <- pars[[nms[i]]]$phi
        whP <- match("phi", names(pars[[nms[i]]]), nomatch = 0L)
      }

      if (length(c(whB, whP)) > 0) other_par <- pars[[nms[i]]][-c(whB, whP)]
      else other_par <- pars[[nms[i]]]

      # if (family[i] == 2) df[i] <- pars[[nms[i]]]$df
    }
  }
  else if (is(family[[1]], "causl_family")) {
    for (i in seq_len(d)) {
      prs <- pars[[nms[i]]]
      if ("phi" %in% family[[i]]$pars) phi[i] <- prs$phi
      other_par[[i]] <- prs[!(names(prs) %in% c("beta", "phi"))]
    }
  }
  else stop("'family' should be a valid numeric vector or list of causl_family objects")

  wts <- get_X_density(dat, eta=eta, phi=phi, qden=qden, family=family,
                       link=link, other_par=other_par)

  if (any(is.na(wts))) stop("Problem with weights")

  wts
}


##' @param careful should full, slower method be used?
## @param max_wt maximum weight obtained
##' @describeIn sim_inversion Rejection sampling code
sim_rejection <- function (out, proc_inputs, careful) {

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

    ## get range of bins for unbounded continuous variables
    rg <- matrix(nrow = sum(unb2_cts) + sum(unb_cts) + sum(b01), ncol=2)
    for (i in which(unb2_cts)) {
      rg[i,] <- range(Z0s[,i])
      rg[i,][1] <- floor(rg[[i]][1])
      rg[i,][2] <- ceiling(rg[[i]][2])
    }
    for (i in which(unb_cts)) {
      rg[i,] <- c(0, ceiling(max(Z0s[,i])))
    }
    for (i in which(b01)) {
      rg[i,] <- c(0,1)
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
  max_wt <- 1  # max value to divide by

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
    out[!OK,output] <- sim_copula(out[!OK,output,drop=FALSE], family=famCop,
                                  par = pars$cop, df=pars$cop$df, model_matrix=copMM)
    if (careful) OB <- rep(FALSE, nr)

    for (i in seq_along(LHS_Z)) {
      mms[[1]][[i]] <- model.matrix(formulas[[1]][[i]], data=out[!OK,,drop=FALSE])
      out[[LHS_Z[i]]][!OK] <- rescale_var(out[[LHS_Z[i]]][!OK], X=mms[[1]][[i]],
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
      out[[LHS_Y[i]]][!OK] <- rescale_var(out[[LHS_Y[i]]][!OK], X=mms[[3]][[i]],
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
      max_wt <- max(max(wts), max_wt)
      wts <- wts/max_wt
    }

    OK[!OK] <- runif(nr) < wts
    nr <- sum(!OK)
  }

  return(out)
}
