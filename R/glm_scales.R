##' Obtain univariate densities
##'
##' Ultimately should also work for ordinal and categorical cases
##'
##' @inheritParams fit_causl
##' @inheritParams gen_dummy_dat
##'
process_discrete_dens <- function (dat, family, LHSs) {

  ## reorder variables so that discrete ones come last
  ## for discrete variables, plug in empirical probabilities
  disc <- (family == 5) | (family == 0)
  family[family == 0] <- 5
  trunc <- list()

  if (any(disc)) {
    wh_disc <- which(disc)

    ## tabulate discrete variables
    for (i in seq_along(wh_disc)) {
      trunc[[i]] <- tabulate(dat[[LHSs[wh_disc[i]]]] + 1)
      if (sum(trunc[[i]]) == 0) stop("tabulation of values failed")
      trunc[[i]] <- trunc[[i]]/sum(trunc[[i]])
    }

    ## then move discrete variables to the end
    wh_cnt <- which(!disc)
    new_ord <- c(wh_cnt, wh_disc)  # adjust for multiple copula formulae

    family <- family[new_ord]
  }
  else new_ord <- NULL

  return(list(family=family, trunc=trunc, order=new_ord))
}

##' Get univariate densities and uniform order statistics
##'
##' @param x vector of observations
##' @param eta,phi linear component and dispersion parameters
##' @param other_pars other parameters for certain families
##' @param family numeric indicator of family
##' @param link link function
##'
##' @details `fam` follows the usual numeric pattern: 1=normal,
##' 2=t-distribution and 3=Gamma with a log-link.
##'
##' @return A list with entries being the numeric vectors `u` (the
##' quantiles of the input values) and `ld` (the log density of each
##' observation).
##'
##'
glm_dens <- function (x, eta, phi, other_pars, family=1, link) {

  ## deal with case that family is a 'causl_family'
  if (is(family, "causl_family")) {
    if (!missing(link)) warning("Using link from 'causl_family' object")
    link <- family$link

    if (!(family$name %in% names(links_list))) {
      # if (!(link %in% family$links_list)) stop(paste0(link, " is not a valid link function for ", family$name, " family"))
      stop(paste0("Family ", family$name, " is not a valid and registered family"))
    }
    if (!(link %in% links_list[[family$name]])) stop(paste0(link, " is not a valid link function for ", family$name, " family"))

    ## get link function
    mu <- link_apply(eta, link, family$name)

    pars <- list(mu=mu)
    if ("phi" %in% family$pars) pars <- c(pars, list(phi=phi))
    if ("par2" %in% family$pars) pars <- c(pars, list(par2=other_pars$par2))

    lp <- do.call(family$ddist, c(list(x=x, log=TRUE), pars))
    u <- do.call(family$pdist, c(list(x=x), pars))
  }
  else if (is.numeric(family)) {
    if (missing(link)) link <- links_list[[family_vals[family_vals$val==family,"family"]]][1]

    ## get the densities for x
    if (family == 1) {
      if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else if (link=="log") mu <- exp(eta)
      else stop("Not a valid link function for Gaussian distribution")

      lp <- dnorm(x, mu, sd=sqrt(phi), log=TRUE)
      u <- pnorm(x, mu, sd=sqrt(phi))
    }
    else if (family == 2) {
      if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else if (link=="log") mu <- exp(eta)
      else stop("Not a valid link function for t-distribution")
      df <- other_pars$par2

      lp <- dt((x - mu)/sqrt(phi), df=df, log=TRUE) - log(sqrt(phi))
      u <- pt((x - mu)/sqrt(phi), df=df)
    }
    else if (family == 3) {
      if (link=="log") mu <- exp(eta)
      else if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else stop("Not a valid link function for gamma distribution")

      lp <- dgamma(x, shape=1/phi, scale=phi*mu, log=TRUE)
      u <- pgamma(x, shape=1/phi, scale=phi*mu)
    }
    else if (family == 5) {
      if (link=="logit") mu <- expit(eta)
      else if (link=="probit") mu <- pnorm(eta)
      else stop("Not a valid link function for Bernoulli distribution")

      lp <- x*log(mu) + (1-x)*log(1-mu)
      lp[is.nan(lp)] <- 0
      u <- x
    }
    else stop("Only Gaussian, t, gamma and Bernoulli distributions are allowed")
  }

  return(list(u=u, ld=lp))
}

##' @describeIn glm_dens old name
univarDens <- function (x, eta, phi, other_pars, family=1, link) {
  deprecate_soft("0.8.8", "univarDens()", "glm_dens()")
  glm_dens(x=x, eta=eta, phi=phi, other_pars=other_pars, family=family, link=link)
}

