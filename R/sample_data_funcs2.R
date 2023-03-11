##' Function to obtain density from several linear models
##'
##' @param n number of samples
##' @param fam vector of families
##' @param eta matrix of covariates
##' @param link vector of link functions
##' @param pars list of parameters
##' @param LHS vector of names
##' @param sim logical: should data be generated, or just the density function returned?
##'
##' @export
get_X_dens <- function (n, fam, eta, link, pars, LHS, sim=TRUE) {
  d <- length(fam)
  if (missing(eta)) eta <- matrix(1, n, d)
  else if (ncol(eta) != d) stop("Mismatch in number of variables between 'fam' and 'eta'")
  if (length(LHS) != d) stop("'LHS' must contain name for each variable to be simulated")

  ## get list to output densities
  qden <- vector(mode="list", length=d)
  if (missing(n) && missing(eta)) stop("must supply either sample size or parameters")
  else if (missing(n)) n <- nrow(eta)

  out <- data.frame(matrix(0, ncol=length(LHS), nrow=n))
  names(out) <- LHS

  for (i in seq_along(LHS)) {
    ## get parameters for proposal
    if (fam[i] == 0 || fam[i] == 5) {
      if (!is.null(pars[[LHS[i]]][["p"]])) theta <- pars[[LHS[i]]]$p
      else theta <- 0.5
      fam[i] <- 5
    }
    else if (fam[i] == 1 || fam[i] == 2) {
      theta <- 2*pars[[LHS[i]]]$phi
    }
    else if (fam[i] == 6) {
      theta <- 1.5*pars[[LHS[i]]]$phi
    }
    else if (fam[i] == 3) {
      theta <- 2*pars[[LHS[i]]]$phi
    }
    else if (fam[i] == 4) {
      theta = c(1,1)
    }

    ## obtain data for X's
    tmp <- sim(n, fam = fam[i], theta=theta, sim=TRUE)
    out[LHS[[i]]] <- tmp$x
    qden[[i]] <- tmp$qden
  }

  return(list(datX = out, qden = qden))
}

##' Simulate initial X values
##'
##' @param n number of observations
##' @param fam number for distribution family
##' @param theta parameters for model
##' @param offset optional mean correction
##' @param sim logical: should values be simulated?
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
sim <- function(n, fam, theta, offset, sim=TRUE) {

  if (missing(offset)) offset <- 0

  if (fam == 1) {
    if (sim) X <- rnorm(n, mean=offset, sd=sqrt(theta))
    qden <- function(x) dnorm(x, mean=offset, sd=sqrt(theta))
  }
  else if (fam == 6) {
    if (sim) X <- exp(rnorm(n, mean=offset, sd=sqrt(theta)))
    qden <- function(x) dnorm(log(x), mean=offset, sd=sqrt(theta))/x
  }
  else if (fam == 2) {
    if (sim) X <- theta[1]*(rt(n, df=theta[2]) + offset)
    qden <- function(x) dt((x-offset)/theta[1], df=theta[2])/theta[1]
  }
  else if (fam == 3) {
    if (sim) X <- rgamma(n, shape=1, rate=1/theta)
    qden <- function(x) dgamma(x, shape=1, rate=1/theta)
  }
  else if (fam == 4) {
    if (length(theta) == 1) theta <- c(theta, theta)
    if (sim) X <- rbeta(n, theta[1], theta[2])
    qden <- function(x) dbeta(x, theta[1], theta[2])
  }
  else if (fam == 5) {
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
  else stop("X distribution must be normal (1), t (2), Gamma (3), Beta (4) or categorical (5)")

  if (sim) out <- list(x=X, qden=qden)
  else out <- list(qden)

  return(out)
}

##' Get means and linear forms for set of variables
##'
##' @param dat a data.frame containing all relevant variables
##' @param pars a list whose entries include every variable in \code{LHS}
##' @param full_form a merged formula whose left hand sides are in \code{LHS}
##' @param LHS a character vector of variables
##'
##' @description Both \code{eta_matrix} and \code{mu_matrix} produce a matrix
##' with \code{nrows(dat)} rows and a column for each entry in \code{LHS}.  In
##' the case of \code{eta_matrix} the entry is the corresponding linear form
##' \eqn{X\beta} and for \code{mu_matrix} it is the transformed \eqn{g^{-1}(X\beta)}
##' that is also \eqn{E [Y | X=x]}
##'
eta_matrix <- function (dat, pars, full_form, LHS) {
  mm <- model.matrix(full_form$formula, dat)
  wh <- full_form$wh
  d <- length(LHS)

  beta_m <- beta <- matrix(0, nrow=ncol(mm), ncol=d)
  for (i in seq_len(d)) {
    # beta_m[wh[i], i] <- 1
    beta[wh[[i]], i] <- pars[[LHS[[i]]]]$beta
  }

  eta <- mm %*% beta
  eta
}

##' @describeIn eta_matrix transformed to mean scale
##' @param link list of link functions, one for each entry of \code{LHS}
mu_matrix <- function(dat, pars, full_form, LHS, link) {
  mu <- eta_matrix(dat=dat, pars=pars, full_form=full_form, LHS=LHS)

  for (i in seq_along(LHS)) {
    mu[,i] <- stats::make.link(link[[LHS[[i]]]])$linkinv(mu[,i])
  }
  mu
}

## Convert linear responses into mean responses
##
## @param eta matrix of mean responses
## @param link list of vector of link functions
##
eta2mu <- function(eta, link) {
  if (length(link) != ncol(eta)) stop("Number of links does not match number of variables")

  for (i in seq_along(link)) {
    eta[,i] <- stats::make.link(link[[i]])$linkinv(eta[,i])
  }
  eta
}
