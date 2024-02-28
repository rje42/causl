##' @importFrom methods is
NULL

##' Numbers for parametric families
##'
##' Each function returns a data frame containing
##' \itemize{
##' \item `val`: an integer
##' \item `family`: a vector giving the associated parametric family for that integer.
##' }
##'
##' @export
familyVals <- data.frame(val=c(0:6,11,10),
                         family=c("binomial", "gaussian", "t", "Gamma", "beta", "binomial", "lognormal","ordinal","categorical"))

##' @describeIn familyVals Values for copula families
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
##'
##' @seealso [familyVals]
##'
get_family <- function (val) {
  if (is.numeric(val)) {
    fm <- match(val, familyVals$val)
    if (is.na(fm)) stop("Invalid family value")
  }
  else if (is.character(val)) {
    fm <- pmatch(val, familyVals$family)
    if (is.na(fm)) stop("Family not recognized")
  }

  fmly <- familyVals[fm,]$family
  if (is.numeric(val) && val > 5) stop("No function defined yet for this family")

  get(paste0(fmly, "_causl_fam"))
}

##' @describeIn get_family Gaussian distribution family
##' @param link link function
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
              pars=c("mu"), default=default)
  class(out) <- "causl_family"

  return(out)
}

##' @describeIn get_family beta distribution family
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
              pars=c("mu"), default=default)
  class(out) <- "causl_family"

  return(out)
}

##' @describeIn get_family ordinal categorical distribution family
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
              pars=c("mu"), default=default)
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

