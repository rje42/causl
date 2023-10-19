## version of make.link that includes 'exp' for lognormal
make.link <-
  function (link)
  {
    switch(link, logit = {
      linkfun <- function(mu) logit(mu) # .Call(C_logit_link, mu)
      linkinv <- function(eta) expit(eta) # .Call(C_logit_linkinv, eta)
      mu.eta <- function(eta) expit(eta)*expit(-eta) # .Call(C_logit_mu_eta, eta)
      valideta <- function(eta) TRUE
    }, probit = {
      linkfun <- function(mu) qnorm(mu)
      linkinv <- function(eta) {
        thresh <- -qnorm(.Machine$double.eps)
        eta <- pmin(pmax(eta, -thresh), thresh)
        pnorm(eta)
      }
      mu.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
      valideta <- function(eta) TRUE
    }, cauchit = {
      linkfun <- function(mu) qcauchy(mu)
      linkinv <- function(eta) {
        thresh <- -qcauchy(.Machine$double.eps)
        eta <- pmin(pmax(eta, -thresh), thresh)
        pcauchy(eta)
      }
      mu.eta <- function(eta) pmax(dcauchy(eta), .Machine$double.eps)
      valideta <- function(eta) TRUE
    }, cloglog = {
      linkfun <- function(mu) log(-log(1 - mu))
      linkinv <- function(eta) pmax(pmin(-expm1(-exp(eta)),
                                         1 - .Machine$double.eps), .Machine$double.eps)
      mu.eta <- function(eta) {
        eta <- pmin(eta, 700)
        pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
      }
      valideta <- function(eta) TRUE
    }, identity = {
      linkfun <- function(mu) mu
      linkinv <- function(eta) eta
      mu.eta <- function(eta) rep.int(1, length(eta))
      valideta <- function(eta) TRUE
    }, log = {
      linkfun <- function(mu) log(mu)
      linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
      mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
      valideta <- function(eta) TRUE
    }, exp = {
      linkfun <- function(mu) exp(mu)
      linkinv <- function(eta) log(eta)
      mu.eta <- function(eta) 1/eta
      valideta <- function(eta) all(is.finite(eta)) && all(eta >
                                                             0)
    },  sqrt = {
      linkfun <- function(mu) sqrt(mu)
      linkinv <- function(eta) eta^2
      mu.eta <- function(eta) 2 * eta
      valideta <- function(eta) all(is.finite(eta)) && all(eta >
                                                             0)
    }, `1/mu^2` = {
      linkfun <- function(mu) 1/mu^2
      linkinv <- function(eta) 1/sqrt(eta)
      mu.eta <- function(eta) -1/(2 * eta^1.5)
      valideta <- function(eta) all(is.finite(eta)) && all(eta >
                                                             0)
    }, inverse = {
      linkfun <- function(mu) 1/mu
      linkinv <- function(eta) 1/eta
      mu.eta <- function(eta) -1/(eta^2)
      valideta <- function(eta) all(is.finite(eta)) && all(eta !=
                                                             0)
    }, logit2 = {
      linkfun <- function(mu) logit((1+mu)/2)
      linkinv <- function(eta) 2*expit(eta) - 1
      mu.eta <- function(eta) 2*expit(eta)*expit(-eta)
      valideta <- function(eta) all(is.finite(eta))
    }, log1p = {
      linkfun <- function(mu) log1p(mu)
      linkinv <- function(eta) exp(eta) - 1
      mu.eta <- function(eta) exp(eta)
      valideta <- function(eta) all(is.finite(eta))
    }, log1m = {
      linkfun <- function(mu) log(mu - 1)
      linkinv <- function(eta) exp(eta) + 1
      mu.eta <- function(eta) exp(eta)
      valideta <- function(eta) all(is.finite(eta))
    }, stop(gettextf("%s link not recognised", sQuote(link)),
            domain = NA))
    environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(valideta) <- asNamespace("stats")
    structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                   valideta = valideta, name = link), class = "link-glm")
  }


normal <- function (link = "identity") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("identity", "inverse", "log")
  family <- "normal"
  if (linktemp %in% okLinks) {
    stats <- make.link(linktemp)
  } else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    stop(gettextf("link \"%s\" not available for %s family; available links are %s",
                  linktemp, family, paste(sQuote(okLinks), collapse = ", ")),
         domain = NA)
  }

  structure(list(family = family, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = function(mu) rep(1, length(mu)),
                 scale = function(u, pars) {
                   qnorm(u, mean=pars$mu, sd=sqrt(pars$phi))
                 },
                 simulate = function(n, pars) {
                   rnorm(n, mean=pars$mu, sd=sqrt(pars$phi))
                 },
                 mu.eta = stats$mu.eta,
                 validmu = function(mu) TRUE,
                 valideta = stats$valideta,
                 dispersion = TRUE,
                 par2 = FALSE))
}

student <- function (link = "identity") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("identity", "inverse", "log")
  family <- "normal"
  if (linktemp %in% okLinks) {
    stats <- make.link(linktemp)
  } else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    stop(gettextf("link \"%s\" not available for %s family; available links are %s",
                  linktemp, family, paste(sQuote(okLinks), collapse = ", ")),
         domain = NA)
  }

  structure(list(family = family, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = function(mu) rep(1, length(mu)),
                 scale = function(u, pars) {
                   sqrt(pars$phi)*qt(u, df=pars$par2) + pars$mu
                 },
                 simulate = function(n, pars) {
                   sqrt(pars$phi)*rt(n, df=pars$par2) + pars$mu
                 },
                 mu.eta = stats$mu.eta,
                 validmu = function(mu) TRUE,
                 valideta = stats$valideta,
                 dispersion = TRUE,
                 par2 = TRUE))
}


##' log-Gaussian distribution
##'
##' @param link link function to use
##'
lognormal <- function(link = "exp") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("exp", "identity")
  family <- "lognormal"
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for %s family; available links are %s",
                    linktemp, family, paste(sQuote(okLinks), collapse = ", ")),
           domain = NA)
    }
  }
  structure(list(family = family, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = function(mu) mu^2, dev.resids = function(y, mu, wt) wt *
                   ((y - mu)^2),
                 # aic = function(y, n, mu, wt, dev) {
                 #     nobs <- length(y)
                 #     nobs * (log(dev/nobs * 2 * pi) + 1) + 2 - sum(log(wt))
                 #   },
                 mu.eta = stats$mu.eta,
                 # initialize = expression({
                 #     n <- rep.int(1, nobs)
                 #     if (is.null(etastart) && is.null(start) && is.null(mustart) &&
                 #         ((family$link == "exp" && any(y == 0)) ||
                 #          (family$link == "identity" && any(y <= 0)))) stop("cannot find valid starting values: please specify some")
                 #     mustart <- log(y)
                 #   }),
                 validmu = function(mu) TRUE, valideta = stats$valideta),
            class = "family")
}
