## return the causl_family function given an integer index
get_family <- function (val) {
  fm <- match(val, familyVals$val)
  if (is.na(fm)) stop("Invalid family value")

  fmly <- familyVals[fm,]$family
  if (val > 2) stop("Only works for Gaussian and t variables")

  get(paste0(fmly, "_causl_fam"))
}

gaussian_causl_fam <- function () {
  ## write functions
  dens <- function (x, mu, phi, log=FALSE) dnorm(x, mean=mu, sd=sqrt(phi), log=log)
  quan <- function (p, mu, phi) qnorm(p, mean=mu, sd=sqrt(phi))
  sim <- function (n, mu, phi) rnorm(n, mean=mu, sd=sqrt(phi))
  probs <- function (x, mu, phi) pnorm(x, mean=mu, sd=sqrt(phi))

  default <- function(theta) list(x=0, mu=0, phi=theta)

  ## define family
  out <- list(name="gaussian", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu", "phi"), default=default)
  class(out) <- "causl_family"

  return(out)
}

t_causl_fam <- function () {
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
              pars=c("mu", "phi", "par2"), default=default)
  class(out) <- "causl_family"

  return(out)
}

gamma_causl_fam <- function () {
  ## write functions
  dens <- function (x, mu, phi, log=FALSE) dgamma(x, rate=1/(mu*phi),
                                                  shape=1/phi, log=log)
  quan <- function (p, mu, phi) qgamma(p, rate=1/(mu*phi), shape=1/phi)
  sim <- function (n, mu, phi) rgamma(n, rate=1/(mu*phi), shape=1/phi)
  probs <- function (x, mu, phi) pgamma(x, rate=1/(mu*phi), shape=1/phi)

  default <- function(theta) list(x=1, mu=1, phi=theta, p=0.5)

  ## define family
  out <- list(name="Gamma", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu", "phi"), default=default)
  class(out) <- "causl_family"

  return(out)
}

binomial_causl_fam <- function () {
  ## write functions
  dens <- function (x, mu, log=FALSE) dbinom(x, size=1, prob=mu, log=log)
  quan <- function (p, mu) qbinom(p, size=1, prob=mu)
  sim <- function (n, mu) rbinom(n, size=1, prob=mu)
  probs <- function (x, mu) pbinom(x, size=1, prob=mu)

  default <- function(theta) list(x=1, mu=1, p=0.5)

  ## define family
  out <- list(name="binomial", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu"), default=default)
  class(out) <- "causl_family"

  return(out)
}

beta_causl_fam <- function () {
  ## write functions
  dens <- function (x, mu, phi, log=FALSE) dbeta(x, shape1=1+phi*mu,
                                                 shape2=1+phi*(1-mu), log=log)
  quan <- function (p, mu, phi) qbeta(p, shape1=1+phi*mu, shape2=1+phi*(1-mu))
  sim <- function (n, mu, phi) rbeta(n, shape1=1+phi*mu, shape2=1+phi*(1-mu))
  probs <- function (x, mu, phi) pbeta(x, shape1=1+phi*mu, shape2=1+phi*(1-mu))

  default <- function(theta) list(x=1, mu=0.5, phi=2*(theta-1), p=0.5)

  ## define family
  out <- list(name="beta", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu", "phi"), default=default)
  class(out) <- "causl_family"

  return(out)
}

# lnormal_causl_fam <- function () {
#
# }
