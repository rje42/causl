##' @importFrom methods is
NULL

##' Numbers for parametric families
##'
##' Data frames containing
##' \itemize{
##' \item `val`: an integer
##' \item `family`: a vector giving the associated parametric family for that integer.
##' }
##' The integer `val` may be used in place of the name of the parametric family
##' when specifying the `family` object.
##'
##' @format `family_vals` is a `data.frame` with 9 rows and 2 columns
##'
##' @export
family_vals <- data.frame(val=c(0:6,11,10),
                         family=c("binomial", "gaussian", "t", "Gamma", "beta", "binomial", "lognormal","ordinal","categorical"))

##' @describeIn family_vals Old name
##' @format `familyVals` is the same object as `family_vals`
##' @details
##' `familyVals` will be removed in version 1.0.0.
##'
##' @export
familyVals <- family_vals

##' @describeIn family_vals Values for copula families
##' @format `copula_vals` is a `data.frame` with 7 rows and 2 columns
##' @export
copula_vals <- data.frame(val=c(1:6,11),
                          family=c("gaussian", "t", "Clayton", "Gumbel", "Frank", "Joe", "FGM"))

##' Return causl_fam function from integer index
##'
##' @param val integer corresponding to distributional family
##'
##' @details
##' The functions `gaussian_causl_fam()` etc. represent the functions that are
##' returned by `get_family()`.
##'
##' A few function of this form can be defined by the user, and it should return
##' the following:
##'
##'   * `name`: the name of the relevant family;
##'   * `ddist`: a function returning the density of the distributions;
##'   * `qdist`: a function returning the quantiles from probabilities;
##'   * `rdist`: a function to sample values from the distribution;
##'   * `pdist`: a cumulative distribution function;
##'   * `pars`: a list of the names of the parameters used;
##'   * `default`: a function that returns a list of the default values for an
##'   observation and each of the parameters;
##'   * `link`: the specified link function.
##'
##' The function should also give the output the class `"causl_family"`, so that
##' it is interpreted appropriately.  Note that `ddist` should have a `log`
##' argument, to allow the log-likelihood to be evaluated.
##'
##'
##' @seealso [family_vals]
##'
get_family <- function (val) {
  if (is.numeric(val)) {
    fm <- match(val, family_vals$val)
    if (is.na(fm)) stop("Invalid family value")
  }
  else if (is.character(val)) {
    fm <- pmatch(val, family_vals$family)
    if (is.na(fm)) stop("Family not recognized")
  }

  fmly <- family_vals[fm,]$family
  if (is.numeric(val) && val > 5) stop("No function defined yet for this family")

  get(paste0(fmly, "_causl_fam"))
}

##' @describeIn get_family Gaussian distribution family
##' @param link link function
##' @export
gaussian_causl_fam <- function (link) {
  if (missing(link)) link = "identity"

  ## write functions
  dens <- function (x, mu, phi, log=FALSE) dnorm(x, mean=mu, sd=sqrt(phi), log=log)
  quan <- function (p, mu, phi) qnorm(p, mean=mu, sd=sqrt(phi))
  sim <- function (n, mu, phi) rnorm(n, mean=mu, sd=sqrt(phi))
  probs <- function (x, mu, phi) pnorm(x, mean=mu, sd=sqrt(phi))

  default <- function(theta) list(x=0, mu=0, phi=theta)

  ## define family
  out <- list(name="gaussian", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu", "phi"), default=default, link=link)
  class(out) <- "causl_family"

  return(out)
}

##' @describeIn get_family Student's t distribution family
##' @export
t_causl_fam <- function (link) {
  if (missing(link)) link = "identity"

  ## write functions
  dens <- function (x, mu, phi, par2, log=FALSE) {
    out <- dt((x-mu)/sqrt(phi), df=par2, log=log)
    if (log) out <- out - log(phi)/2
    else out <- out/sqrt(phi)

    return(out)
  }
  quan <- function (p, mu, phi, par2) qt(p, df=par2)*sqrt(phi) + mu
  sim <- function (n, mu, phi, par2) rt(n, df=par2)*sqrt(phi) + mu
  probs <- function (x, mu, phi, par2) pt((x-mu)/sqrt(phi), df=par2)

  default <- function (theta) list(x=0, mu=0, phi=theta[1], par2=theta[2])

  ## define family
  out <- list(name="t", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu", "phi", "par2"), default=default, link=link)
  class(out) <- "causl_family"

  return(out)
}

##' @describeIn get_family Gamma distribution family
##' @export
Gamma_causl_fam <- function (link) {
  if (missing(link)) link = "log"

  ## write functions
  dens <- function (x, mu, phi, log=FALSE) dgamma(x, rate=1/(mu*phi),
                                                  shape=1/phi, log=log)
  quan <- function (p, mu, phi) qgamma(p, rate=1/(mu*phi), shape=1/phi)
  sim <- function (n, mu, phi) rgamma(n, rate=1/(mu*phi), shape=1/phi)
  probs <- function (x, mu, phi) pgamma(x, rate=1/(mu*phi), shape=1/phi)

  default <- function(theta) list(x=1, mu=1, phi=theta, p=0.5)

  ## define family
  out <- list(name="Gamma", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu", "phi"), default=default, link=link)
  class(out) <- "causl_family"

  return(out)
}

##' @describeIn get_family binomial distribution family
##' @export
binomial_causl_fam <- function (link) {
  if (missing(link)) link = "logit"

  ## write functions
  dens <- function (x, mu, log=FALSE) dbinom(x, size=1, prob=mu, log=log)
  quan <- function (p, mu) qbinom(p, size=1, prob=mu)
  sim <- function (n, mu) rbinom(n, size=1, prob=mu)
  probs <- function (x, mu) pbinom(x, size=1, prob=mu)

  default <- function(theta) list(x=1, mu=0.5, p=0.5)

  ## define family
  out <- list(name="binomial", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu"), default=default, link=link)
  class(out) <- "causl_family"

  return(out)
}

##' @describeIn get_family beta distribution family
##' @export
beta_causl_fam <- function (link) {
  if (missing(link)) link = "logit"

  ## write functions
  dens <- function (x, mu, phi, log=FALSE) dbeta(x, shape1=1+phi*mu,
                                                 shape2=1+phi*(1-mu), log=log)
  quan <- function (p, mu, phi) qbeta(p, shape1=1+phi*mu, shape2=1+phi*(1-mu))
  sim <- function (n, mu, phi) rbeta(n, shape1=1+phi*mu, shape2=1+phi*(1-mu))
  probs <- function (x, mu, phi) pbeta(x, shape1=1+phi*mu, shape2=1+phi*(1-mu))

  default <- function(theta) list(x=1, mu=0.5, phi=2*(theta-1), p=0.5)

  ## define family
  out <- list(name="beta", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu", "phi"), default=default, link=link)
  class(out) <- "causl_family"

  return(out)
}

# lnormal_causl_fam <- function () {
#
# }

# ##' @describeIn get_family multinomial/categorical distribution family
# multinom_causl_fam <- function (link) {
#   if (missing(link)) link = "logit"
#
#   ## write functions
#   dens <- function (x, mu, log=FALSE) dmultinom(x, size=1, prob=mu, log=log)
#   quan <- function (p, mu) qmultinom(p, size=1, prob=mu)
#   sim <- function (n, mu) rmultinom(n, size=1, prob=mu)
#   probs <- function (x, mu) pmultinom(x, size=1, prob=mu)
#
#   default <- function(theta) list(x=1, mu=0.5, p=0.5)
#
#   ## define family
#   out <- list(name="multinomial", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
#               pars=c("mu"), default=default)
#   class(out) <- "causl_family"
#
#   return(out)
# }

##' @describeIn get_family multinomial/categorical distribution family
##' @export
categorical_causl_fam <- function (link) {
  if (missing(link)) link = "logit"

  set_mu_matrix <- function (mu, len) {
    if (is.matrix(mu)) {
      nd <- ncol(mu)
      if (!(nrow(mu) %in% c(1,len))) stop("'mu' must contain either 1 or length of other arg vectors of parameters")
      else if (nrow(mu) == 1) mu <- matrix(mu, len, length(mu), byrow = TRUE)
    }
    else {
      nd <- length(mu)
      mu <- matrix(mu, len, nd, byrow = TRUE)
    }

    return(mu)
  }

  ## write functions
  dens <- function (x, mu, log=FALSE) {
    x <- as.integer(x)

    ## get matrix form for mu
    mu <- set_mu_matrix(mu, length(x))

    nd <- ncol(mu)

    out <- mu[cbind(seq_along(x),x)]
    if (any(x > nd)) out[x > nd] <- NA
    if (log) out <- log(out)
    return(out)
  }
  quan <- function (p, mu) {
    ## get matrix form for mu
    mu <- set_mu_matrix(mu, length(p))

    qsum <- t(apply(mu, 1, cumsum))
    out <- rowSums(p > qsum) + 1
    return(out)
  }
  sim <- function (n, mu) {
    ## get matrix form for mu
    mu <- set_mu_matrix(mu, n)
    U <- runif(n)
    Y <- 1 + colSums(apply(mu, 1, cumsum) < rep(U, each=ncol(mu)))
    Y <- factor(Y, levels=seq_len(ncol(mu)))

    attr(Y, "quantile") <- U

    return(Y)
    # sample(length(mu), size=n, replace=TRUE, prob=mu)
  }
  probs <- function (x, mu) {
    x <- as.integer(x)

    ## get matrix form for mu
    mu <- set_mu_matrix(mu, length(x))

    qsum <- t(apply(mu, 1, cumsum))
    out <- qsum[cbind(seq_along(x),x)]
    #
    # qsum <- cumsum(mu)
    # out <- qsum[x]

    return(out)
  }

  default <- function(theta) list(x=1, mu=c(0.5,0.25,0.25), p=c(0.5,0.25,0.25), nlevel=3)

  ## define family
  out <- list(name="categorical", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu"), default=default, link=link)
  class(out) <- "causl_family"

  return(out)
}

##' @describeIn get_family ordinal categorical distribution family
##' @export
ordinal_causl_fam <- function (link) {
  if (missing(link)) link = "logit"

  set_mu_matrix <- function (mu, len) {
    if (is.matrix(mu)) {
      nd <- ncol(mu)
      if (!(nrow(mu) %in% c(1,len))) stop("'mu' must contain either 1 or length of other arg vectors of parameters")
      else if (nrow(mu) == 1) mu <- matrix(mu, len, length(mu), byrow = TRUE)
    }
    else {
      nd <- length(mu)
      mu <- matrix(mu, len, nd, byrow = TRUE)
    }

    return(mu)
  }

  ## write functions
  dens <- function (x, mu, log=FALSE) {
    x <- as.integer(x)

    ## get matrix form for mu
    mu <- set_mu_matrix(mu, length(x))

    nd <- ncol(mu)

    out <- mu[cbind(seq_along(x),x)]
    if (any(x > nd)) out[x > nd] <- NA
    if (log) out <- log(out)
    return(out)
  }
  quan <- function (p, mu) {
    ## get matrix form for mu
    mu <- set_mu_matrix(mu, length(p))

    qsum <- t(apply(mu, 1, cumsum))
    out <- rowSums(p > qsum) + 1
    return(out)
  }
  sim <- function (n, mu) {
    ## get matrix form for mu
    mu <- set_mu_matrix(mu, n)
    U <- runif(n)
    Y <- 1 + colSums(apply(mu, 1, cumsum) < rep(U, each=ncol(mu)))
    Y <- factor(Y, levels=seq_len(ncol(mu)))

    attr(Y, "quantile") <- U

    return(Y)
    # sample(length(mu), size=n, replace=TRUE, prob=mu)
  }
  probs <- function (x, mu) {
    x <- as.integer(x)

    ## get matrix form for mu
    mu <- set_mu_matrix(mu, length(x))

    qsum <- t(apply(mu, 1, cumsum))
    out <- qsum[cbind(seq_along(x),x)]
    #
    # qsum <- cumsum(mu)
    # out <- qsum[x]

    return(out)
  }


  default <- function(theta) {
    nl <- length(theta) + 1
    p <- rep(1/nl, nl)
    out <- list(x=factor(1, levels=seq(0,length(theta))+1), mu=p, p=p, nlevel=nl)
    return(out)
  }

  ## define family
  out <- list(name="ordinal", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu"), default=default, link=link)
  class(out) <- "causl_family"

  return(out)
}

##' Check if family is categorical
##'
##' @param x a family, either numerical, a name, or a `causl_family` object
##'
##' @details Returns a logical indicating if the object is the input object
##' represents a categorical or ordinal variable.  If it cannot represent a
##' family then `NA` is returned.
##'
##' @export
is_categorical <- function (x) {
  if (is.numeric(x)) return(x %in% 10:11)
  else if (is.character(x)) return(x %in% c("categorical", "ordinal"))
  else if (is(x[[1]], "causl_family")) return(sapply(x, function(y) y$name) %in% c("categorical", "ordinal"))
  else if (is(x, "causl_family")) return(x$name %in% c("categorical", "ordinal"))
  else return(NA)
}

##' @describeIn theta_to_p_cat for ordinal variables
theta_to_p_ord <- function (theta) {
  if (!is.matrix(theta)) theta <- matrix(theta, nrow=1)

  p <- expit(theta)
  for (i in rev(seq_len(ncol(p))[-1])) {
    p[,i] <- p[,i] - p[,i-1]
  }
  NArow <- apply(p, 1, function(x) any(x < 0))
  p[NArow, ] <- NA  ## negative probabilities not allowed

  return(p)
}

##' Transform categorical or ordinal parameters into probabilities
##'
##' @param theta provided log-linear parameters
##'
##' @details Returns the probabilities implied by given log-linear parameters.
##'
theta_to_p_cat <- function (theta) {
  if (!is.matrix(theta)) theta <- matrix(theta, nrow=1)

  p <- cbind(1, exp(theta))
  p <- p/rowSums(p)

  return(p)
}

##' Obtain list of family functions
##'
##' Obtain list of family functions from numeric or character representation
##'
##' @param family numeric or character vector of families
##' @param func_return function to apply to list of families
##'
##' @examples
##' family_list(c(1,3,5))
##' family_list(c("t","binomial"))
##'
##'
##' @export
family_list <- function (family, func_return=get_family) {
  if (is.numeric(family) || is.character(family)) {
    out <- lapply(family, func_return)
    names(out) <- names(family)
  }
  else if (all(sapply(family, is, "function"))) {
    for (i in seq_along(family)) {
      if (!is(family[[i]](), "causl_family")) stop("'family' should be a numeric or character vector, or a list of 'causl_family' objects")
    }
    out <- family
  }
  else stop("'family' should be a numeric or character vector, or a list of 'causl_family' objects")

  return(out)
}

