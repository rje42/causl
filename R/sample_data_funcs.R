gen_X_values <- function (n, famX, pars, LHS_X, dX) {
  ## list for densities used to simulate X's
  qden <- vector(mode="list", length=dX)
  out <- data.frame(matrix(0, ncol=dX, nrow=n))
  names(out) <- LHS_X

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

    ## obtain data for X's
    tmp <- sim_X(n, fam_x = famX[i], theta=theta)
    out[LHS_X[[i]]] <- tmp$x
    qden[[i]] <- tmp$qden
  }

  return(list(datX = out, qden = qden))
}

##' Simulate initial X values
##'
##' @param n number of observations
##' @param fam_x number for distribution family
##' @param theta parameters for model
##' @param offset optional mean correction
##'
##' @details Returns a list that includes a data frame containing a column
##' \code{x}, as well as the density that was used to generate it.  Possible
##' families are Gaussian (=1), t (=2), Exponential (=3), beta (=4)
##' Bernoulli/categorical (=5) and log-normal (=6).
##'
##' For the exponential distribution, \code{theta} is the mean.
##' Beta can take one or two parameters, and if there is only
##' one it is just repeated.
##'
##' The \code{offset} parameter alters the median for the normal and t-distributions,
##' or the median of the logarithm in the case of a log-normal.
##'
##' @return A list with two entries: \code{x} a vector of the simulated
##' values, and \code{qden}, which contains a function that evaluates to the
##' density of the distribution used to generate those values.
##'
##' @export
sim_X <- function(n, fam_x, theta, offset) {

  if (missing(offset)) offset <- 0

  if (fam_x == 1) {
    X <- rnorm(n, mean=offset, sd=sqrt(theta))
    qden <- function(x) dnorm(x, mean=offset, sd=sqrt(theta))
  }
  else if (fam_x == 6) {
    X <- exp(rnorm(n, mean=offset, sd=sqrt(theta)))
    qden <- function(x) dnorm(log(x), mean=offset, sd=sqrt(theta))/x
  }
  else if (fam_x == 2) {
    X <- theta[1]*(rt(n, df=theta[2]) + offset)
    qden <- function(x) dt((x-offset)/theta[1], df=theta[2])/theta[1]
  }
  else if (fam_x == 3) {
    X <- rgamma(n, shape=1, rate=1/theta)
    qden <- function(x) dgamma(x, shape=1, rate=1/theta)
  }
  else if (fam_x == 4) {
    if (length(theta) == 1) theta <- c(theta, theta)
    X <- rbeta(n, theta[1], theta[2])
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
    X <- sample(length(theta), size=n, replace=TRUE, prob=theta)-1
    qden <- function(x) theta[x+1]
  }
  else stop("X distribution must be normal (1), t (2), Gamma (3), Beta (4) or categorical (5)")

  return(list(x=X, qden=qden))
}

##' Rescale quantiles to arbitrary random variable.
##'
##' @param U vector of quantiles
##' @param X model matrix of covariates
##' @param pars list of parameters (see details)
##' @param family family of distributions to use
##' @param link link function
##'
##' @details \code{family} can be 1, 2, 3, 4 or 5 for Gaussian, t-distributed,
##' Gamma distributed, beta distributed or discrete respectively.
##' \code{pars} should be a list with entries \code{beta} and
##' \code{phi}, as well as possibly \code{par2} and \code{trunc} if the family
##' is set to 2 or 5.
##' \code{U} should have the same length as \code{X} has rows, and
##' \code{X} should have the same number of columns as the length of
##' \code{pars$beta}.
##'
##' @return vector of rescaled variables
##'
##' @export
rescaleVar <- function(U, X, pars, family=1, link) {

  ## get linear component
  eta <- X %*% pars$beta
  phi <- pars$phi

  ## make U normal, t or gamma
  if (family == 1) {
    if (link == "identity") Y <- qnorm(U, mean = eta, sd=sqrt(phi))
    else if (link == "log") Y <- qnorm(U, mean = exp(eta), sd=sqrt(phi))
    else if (link == "inverse") Y <- qnorm(U, mean = 1/eta, sd=sqrt(phi))
    else stop("invalid link function for Gaussian distribution")
  }
  else if (family == 2) {
    # Y <- sqrt(phi)*qt(U, df=pars$par2) + eta

    if (link == "identity") Y <- sqrt(phi)*qt(U, df=pars$par2) + eta
    else if (link == "log") Y <- sqrt(phi)*qt(U, df=pars$par2) + exp(eta)
    else if (link == "inverse") Y <- sqrt(phi)*qt(U, df=pars$par2) + 1/eta
    else stop("invalid link function for t-distribution")
  }
  else if (family == 3) {
    # Y <- qexp(U, rate = 1/(exp(eta)*sqrt(phi)))

    if (link == "log") mu <- exp(eta)
    else if (link == "identity") mu <- eta
    else if (link == "inverse") mu <- 1/eta
    else stop("invalid link function for Gamma distribution")

    Y <- qgamma(U, shape = 1/phi, rate = 1/(mu*phi))
  }
  else if (family == 4) {
    Y <- qbeta(U, shape1 = 1, shape2 = 1)
  }
  else if (family == 0 || family == 5) {
    if (link == "probit") Y <- 1*(eta + qnorm(U) > 0)
    else if (link == "logit") Y <- 1*(eta + qlogis(U) > 0)
    else stop("invalid link function for binomial distribution")

    # trunc <- pars$trunc
    # trnc <- 1
    #
    # stop("Not finished family==5 yet")
    # mat <- matrix(NA, length(U), length(trunc))
    # for (j in seq_along(trunc[[trnc]])) {
    #   mat[,j] <- 1*(U > trunc[[trnc]][j])
    # }
    # Y <- rowSums(mat)
  }
  else stop("family must be between 0 and 5")

  ### get Z values to correct families
  # nms <- names(dat)[grep("z", names(dat))]

  return(Y)
}


##' Simulate copula values
##'
##' @param dat data frame with empty columns
##' @param family numeric indicator of copula type
##' @param par mandatory parameters
##' @param par2 optional parameters
##' @param model_matrix design matrix for covariates
##'
##' @details Returns data frame containing columns \code{y}
##' and \code{z1, ..., zk}.
##'
##' The family variables are numeric and taken from \code{VineCopula}.
##' Use, for example, 1 for Gaussian, 2 for t, 3 for Clayton, 4 for Gumbel,
##' 5 for Frank, 6 for Joe and 11 for FGM copulas.
##'
##' @return A data frame of the same dimension as \code{dat} containing the
##' simulated values.
##'
##' @export
sim_CopVal <- function(dat, family, par, par2, model_matrix) {

  ## if more than one family given, must be a vine copula
  if (length(family) > 1) {
    if (choose(ncol(dat), 2) != length(family)) stop("family for copula has length > 1 but not equal to number of pairs of variables")

    dat <- sim_vinecop(dat, family=family, par=par$beta, par2=par2,
                       model_matrix=model_matrix)
    return(dat)
  }

  if (!(family %in% c(1:6,11)))  stop("family not supported")

  N <- nrow(dat)
  if (nrow(model_matrix) != N) stop("Incompatible dimension of dat and model_matrix")
  d <- ncol(dat)
  # col_nms <- colnames(dat)[-1] # paste("U", seq_len(d), sep="")

  if (d > 2 && family == 11) stop("Multidimensional fgmCopula not yet supported")
  dz <- choose(d,2)

  eta <- model_matrix %*% par$beta
  if (nrow(unique(eta)) == 1) only_one = TRUE
  else only_one = FALSE
  # if (missing(eta)) eta <- rep(par, nrow(dat))
  # if (length(eta) != nrow(dat)) stop("Incorrect number of parameters or rows of covariate matrix")

  ## transform to (-1,1) or other valid range
  if (family %in% 1:2 || family == 11) cors <- 2*expit(eta) - 1  # range (-1,1)
  #else if (family %in% 3:6) cors <- eta
  else if (family %in% 3) cors <- exp(eta) - 1  ## valid range (-1, infinity)
  else if (family %in% c(4,6)) cors <- exp(eta) + 1  ## valid range (1, infinity)
  else if (family %in% 5) cors <- eta   ## valid range all of real line

  # if (d > 2) stop("Not designed for multi-dimensional copulae yet")

  if (family == 1) {
    if (only_one) {
      Sigma <- diag(nrow=d)
      Sigma[upper.tri(Sigma)] <- cors[1,]
      Sigma <- t(Sigma)
      Sigma[upper.tri(Sigma)] <- cors[1,]

      # dat[, vnames] <- rnormCopula(N, Sigma=Sigma)
      dat[] <- rGaussCop(N, Sigma=Sigma)
    }
    else {
      Sigma <- array(diag(nrow = d), dim=c(d,d,N))

      Sigma[upper.tri(Sigma[,,1])] <- t(cors)
      Sigma <- aperm(Sigma, c(2,1,3))
      Sigma[upper.tri(Sigma[,,1])] <- t(cors)

      # if (d == 2) dat[, vnames] <- rnormCopula2(N, Sigma = Sigma)
      # else dat[, vnames] <- rnormCopula(N, Sigma = Sigma)
      dat[] <- rGaussCop(N, Sigma=Sigma)
    }
  }
  else if (family == 2) {
    if (only_one) {
      Sigma <- diag(nrow=d)
      Sigma[upper.tri(Sigma)] <- cors[1,]
      Sigma <- t(Sigma)
      Sigma[upper.tri(Sigma)] <- cors[1,]

      dat[] <- rtCop(N, Sigma=Sigma, df=par2)
    }
    else {
      Sigma <- array(diag(nrow = d), dim=c(d,d,N))

      Sigma[upper.tri(Sigma[,,1])] <- t(cors)
      Sigma <- aperm(Sigma, c(2,1,3))
      Sigma[upper.tri(Sigma[,,1])] <- t(cors)

      dat[] <- rtCop(N, Sigma = Sigma, df=par2)
    }
  }
  else if (family == 11) dat[] <- causl::rfgmCopula(N, d=2, cors)
  else if (family %in% 3:6) {

    ## have to do a loop for these guys
    if (family == 3) cop <- claytonCopula(dim=d)
    else if (family == 4) cop <- gumbelCopula(dim=d)
    else if (family == 5) cop <- frankCopula(dim=d)
    else if (family == 6) cop <- joeCopula(dim=d)
    else stop("Shouldn't get to here")

    if (length(eta) == dz*N) {
      ## here's the loop
      for (i in seq_len(N)) {
        cop <- setTheta(cop, cors[i,])
        dat[i, ] <- rCopula(1, cop)
      }
    }
    else if (length(eta) == dz) {
      cop <- setTheta(cop, cors)
      dat[] <- rCopula(N, cop)
    }
    else stop("Shouldn't get to here")
  }

  return(dat)
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

  wts <- rep(1, nrow(dat))

  for (i in seq_len(d)) {
    phi <- pars[[i]]$phi

    if (family[i] == 1) {
      if (link[i] == "identity") mu <- eta[[i]]
      else if (link[i] == "log") mu <- exp(eta[[i]])
      else if (link[i] == "inverse") mu <- 1/eta[[i]]
      else stop("not a valid link function for Gaussian distribution")

      if (is.numeric(qden[[i]])) wts <- wts*dnorm(dat[,i], mean=mu, sd=sqrt(phi))/qden[[i]]
      else wts <- wts*dnorm(dat[,i], mean=mu, sd=sqrt(phi))/qden[[i]](dat[,i])
    }
    else if (family[i] == 2) {
      if (link[i] == "identity") mu <- eta[[i]]
      else if (link[i] == "log") mu <- exp(eta[[i]])
      else if (link[i] == "inverse") mu <- 1/eta[[i]]
      else stop("not a valid link function for t-distribution")

      if (is.numeric(qden[[i]])) wts <- wts*dt((dat[,i] - mu)/sqrt(phi), df=pars[[i]]$par2)/(sqrt(phi)*qden[[i]])
      else wts <- wts*dt((dat[,i] - mu)/sqrt(phi), df=pars[[i]]$par2)/(sqrt(phi)*qden[[i]](dat[,i]))
    }
    else if (family[i] == 3) {
      if (link[i] == "identity") mu <- eta[[i]]
      else if (link[i] == "log") mu <- exp(eta[[i]])
      else if (link[i] == "inverse") mu <- 1/eta[[i]]
      else stop("not a valid link function for t-distribution")

      if (is.numeric(qden[[i]])) wts <- wts*dgamma(dat[,i], rate=1/(mu*phi), shape=1/phi)/qden[[i]]
      else wts <- wts*dgamma(dat[,i], rate=1/(mu*phi), shape=1/phi)/qden[[i]](dat[,i])
    }
    else if (family[i] == 4) {
      mu <- expit(eta[[i]])
      if (is.numeric(qden[[i]])) wts <- wts*dbeta(dat[,i], shape1=1+phi*mu, shape2=1+phi*(1-mu))/qden[[i]]
      else wts <- wts*dbeta(dat[,i], shape1=1+phi*mu, shape2=1+phi*(1-mu))/qden[[i]](dat[,i])
    }
    else if (family[i] == 5) {
      if (link[i] == "probit") mu <- qnorm(eta[[i]])
      else if (link[i] == "logit") mu <- expit(eta[[i]])
      else stop("invalid link function for binomial distribution")

      if (is.numeric(qden[[i]])) wts <- wts*dbinom(dat[,i], prob=mu, size=1)/qden[[i]]
      else wts <- wts*dbinom(dat[,i], prob=mu, size=1)/qden[[i]](dat[,i])
    }
    else if (family[i] == 6) {
      mu <- eta[[i]]
      if (is.numeric(qden[[i]])) wts <- wts*dnorm(log(dat[,i]), mean=mu, sd=sqrt(phi))/(dat[,i]*qden[[i]])
      else wts <- wts*dnorm(log(dat[,i]), mean=mu, sd=sqrt(phi))/(dat[,i]*qden[[i]](dat[,i]))
    }
    else stop("family[2] must be in the range 1 to 6")
  }

  if (any(is.na(wts))) stop("Problem with weights")

  wts
}

## Set up link functions
##
## @param link the input given to causalSamp()
## @param family the list of families for Z,X and Y variables
## @param vars a list of vectors of variable names with the same structure as \code{family}
##' @export
linkSetUp <- function(link, family, vars) {

  if (!missing(vars) && !all(lengths(vars) == lengths(family))) stop("length of variable names vector does not match number of families provided")

  ## set up output
  link_out <- lapply(family, function(x) {
    familyVals$family[match(x, familyVals$val)]
  })
  if (any(is.na(unlist(link_out)))) stop("Invalid family variable specified")
  link_out <- lapply(link_out, function(x) sapply(x, function(y) linksList[[y]][1]))

  if (!missing(vars)) {
    ## set names to match vars
    nms_list <- relist(unlist(vars), family)
    for (i in seq_along(link_out)) names(link_out[[i]]) <- nms_list[[i]]
  }

  ## if no link argument supplied, then just return this list
  if (missing(link) || is.null(link)) return(link_out)

  tmp <- unlist(linksList)

  ## now add in any modifications made in 'link'
  if (is.list(link)) {
    if (all(lengths(link) == lengths(family))) {
      ## if lengths the same, assume in same position as family variables
      for (i in seq_along(link_out)) {
        link_out[[i]][] <- tmp[pmatch(link[[i]], tmp)]
        if (any(is.na(link_out[[i]]))) stop("link not properly matched")
      }
      return(link_out)
    }
  }

  ## otherwise try to use names to deduce which is which
  if (!missing(vars)) {
    link <- unlist(link)
    if (is.null(names(link))) stop("names must be supplied for links in order to match with them")
    nms <- names(link)
    wh_set <- subsetmatch(as.list(nms), vars)
    if (any(is.na(wh_set))) stop("names must be supplied for links in order to match with them")

    ## now get particular entry in vector
    wh_val <- NA*wh_set
    for (i in seq_along(wh_set)) {
      wh_val[i] <- match(nms[i], vars[[wh_set[i]]])
      if (is.na(wh_val[i])) stop("some links not matched")
      link_out[[wh_set[i]]][wh_val[i]] <- tmp[pmatch(link[i], tmp)]
    }
    if (any(is.na(unlist(link_out)))) stop("some links not matched")
    #
    #
    # for (j in seq_along(vars)) {
    #   mask <- wh_set == j
    #
    #   for (i in seq_along(linksList)) {
    #     lki <- link[mask && fam_nm == nms[i]]
    #     linknm <- pmatch(lki, linksList[[i]])
    #     lki2 <- linksList[[i]][linknm]
    #     if (any(is.na(lki2))) stop("link ", paste(lki[which(is.na(lki2))], sep=", "), " not supported")
    #     lki[]
    #   }
    #   chk <- wh_set
    # }
  }
  else stop("variable names must be provided to match using them")

  # ## matching
  # fam_nm <- familyVals$family[familyVals$val==fams]
  # nms <- names(linksList)
  # for (i in seq_along(linksList)) {
  #   lki <- link[fam_nm == nms[i]]
  #   linknm <- pmatch(lki, linksList[[i]])
  #   lki2 <- linksList[[i]][linknm]
  #   if (any(is.na(lki2))) stop("link ", paste(lki[which(is.na(lki2))], sep=", "), " not supported")
  #

  return(link_out)
}

linksList <- list(
  gaussian = c("identity", "inverse", "log"),
  t = c("identity", "inverse", "log"),
  Gamma = c("log", "inverse", "identity"),
  beta = c("logit", "probit"),
  binomial = c("logit", "probit"),
  lognormal = c("exp", "identity")
)

familyVals <- data.frame(val=0:6,
                         family=c("binomial", "gaussian", "t", "Gamma", "beta", "binomial", "lognormal"))
