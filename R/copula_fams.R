##' Copula family functions
##'
##' @param link link function
##'
##' @details
##' These return a list that contains, for each valid family:
##' \itemize{
##' - `name`: the name of the family;
##' - `ddist`: function  to evaluate the density;
##' - `rdist`: function to obtain samples;
##' - `pars`: character vector of the parameter names;
##' - `default`: list of default values;
##' - `link`: the chosen link function.
##' }
##'
##' @name causl_copula
NULL

##' @describeIn causl_copula Gaussian copula family
##' @export
gaussian_causl_cop <- function (link) {
  if (missing(link)) link = "tanh"

  ## write functions
  dens <- function (x, Sigma, log=FALSE) dGaussCop(x, Sigma=Sigma, log=log)
  sim <- function (n, mu, phi) rGaussCop(n, Sigma=Sigma)
  # probs <- function (x, mu, phi) pnorm(x, Sigma=Sigma)

  default <- function() list(x=c(0,0,0), Sigma=diag(3), d=3)

  ## define family
  out <- list(name="gaussian", ddist=dens, rdist=sim,
              pars=c("Sigma"), default=default, link=link)
  class(out) <- "causl_copula"

  return(out)
}

