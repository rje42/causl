##' Copula family functions
##'
##' @param link link function
##'
##' @details
##' `get_copula` returns the `causl_copula` that corresponds to the particular
##' integer given.  So far, `1` for Gaussian and `2` for t copulas are implemented.
##'
##' The `copula_fam` functions return a list that contains, for each valid family:
##' \itemize{
##' - `name`: the name of the family;
##' - `ddist`: function  to evaluate the density;
##' - `rdist`: function to obtain samples from copula;
##' - `pdist`: copula function
##' - `pars`: character vector of the parameter names;
##' - `default`: list of default values;
##' - `link`: the chosen link function.
##' }
##'
##' @name causl_copula
NULL


##' @describeIn causl_copula getter copula family
##' @param family_index integer representing copula family
##' @export
get_copula <- function(family_index, link = NULL){
  if(family_index == 1){
    if(is.null(link)) return(gaussian_causl_cop())
    else return(gaussian_causl_cop(link))
  }
  else if(family_index == 2){
    if(is.null(link)) return(t_causl_cop())
    else return(t_causl_cop(link))
  }
  else if (family_index > 2 & family_index < 10) stop("Not implemented yet")
  else stop("Not a valid index")
}


##' @describeIn causl_copula Gaussian copula family
##' @export
gaussian_causl_cop <- function (link) {
  if (missing(link)) link <- tanh

  ## write functions
  dens <- function (x, Sigma, log=FALSE) dGaussCop(x, Sigma=Sigma, log=log)
  sim <- function (n,Sigma, other_pars = NULL) rGaussCop(n, Sigma=Sigma)
  # probs <- function (x, mu, phi) pnorm(x, Sigma=Sigma)

  default <- function() list(x=c(0,0,0), Sigma=diag(3), d=3)

  ## define family
  out <- list(name="gaussian", ddist=dens, rdist=sim,
              pars=c("Sigma"), default=default, link=link)
  class(out) <- "causl_copula"

  return(out)
}

##' @describeIn causl_copula t copula family
##' @export
t_causl_cop <- function(link) {
  if (missing(link)) link <- tanh

  ## write functions
  dens <- function (x, Sigma, df, log=FALSE) dtCop(x, Sigma=Sigma, df=df, log=log)
  sim <- function (n, Sigma, df) rtCop(n, Sigma=Sigma, df=df)
  # probs <- function (x, mu, phi) pt(x, df=df, Sigma=Sigma)

  default <- function() list(x=c(0,0,0), Sigma=diag(3), d=3)

  ## define family
  out <- list(name="t", ddist=dens, rdist=sim,
              pars=c("Sigma"), default=default, link=link)
  class(out) <- "causl_copula"

  return(out)

}

##' @describeIn causl_copula empirical copula family
##' @export
emp_causl_cop <- function (link) {
  if (missing(link)) link <- identity

  ## write functions
  # dens <- function (x, Sigma, log=FALSE) dGaussCop(x, Sigma=Sigma, log=log)
  sim <- function (n, u) {
    if (!is.matrix(u)) u <- matrix(u, nrow=1)
    m <- nrow(u)
    return(u[sample(x=m, size=n, replace=TRUE),])
  }
  probs <- function (u) {
    if (!is.matrix(u)) u <- matrix(u, nrow=1)
    out <- numeric(nrow(u))

    for (j in seq(nrow(u))) {
      U2 <- rep(u[j,], each=nrow(pts))

      v <- U2 > pts
      emp <- apply(v, 1, all)
      out[j] <- mean(emp)
    }
    return(out)
  }

  default <- function() list(x=c(0,0,0), Sigma=diag(3), d=3)

  ## define family
  out <- list(name="empirical", rdist=sim, pdist=probs,
              pars=character(0), default=default, link=link)
  class(out) <- "causl_copula"

  return(out)
}


##' @describeIn causl_copula simulate from copula family
##' @param causl_copula family from which to simulate
##' @param beta_matrix matrix of regression coefficients
##' @param other_pars other parameters for some families
##' @param model_matrix matrix of regressors
##' @export
sim_cop <- function(causl_copula, beta_matrix, other_pars, model_matrix) {
  ## if causal_copula is numeric then obtain associated family
  if (is.numeric(causl_copula)) causl_copula <- get_copula(causl_copula)

  n <- nrow(model_matrix)
  p <- dim(beta_matrix)[2]


  us <- matrix(0, nrow = n, ncol = p)
  ## now, if there is a covariate, get levels of all the covariates
  if (nrow(unique(model_matrix)) > 1) {
    facts <- as.matrix(apply(model_matrix, 2, function(x) as.numeric(factor(x)))) - 1
    levs <- apply(facts, 2, max)+1
    b <- cumprod(c(1, levs[-length(levs)]))
    ovl_fact <- facts %*% b
    fact <- factor(ovl_fact)
  }
  else {
    fact <- factor(rep(1,n))
  }

  for (gp in levels(fact)) {
    gp_mem <- (fact == gp)
    n_gp <- sum(gp_mem)
    if (n_gp == 0) next

    eta <- c(model_matrix[match(TRUE, gp_mem), ,drop=FALSE] %*% beta_matrix)
    rhos <- causl_copula$link(eta)  # link is applied to eta, not beta
    Sigma <- matrix(1, nrow = p, ncol = p)
    Sigma[upper.tri(Sigma)] <- rhos
    Sigma[lower.tri(Sigma)] <- t(Sigma[upper.tri(Sigma)])
    eigen_values <- eigen(Sigma, only.values=TRUE)$values
    if(any(eigen_values < 0)) stop("specified non positive semi-definite correlation matrix")
    us[gp_mem,] <- causl_copula$rdist(n_gp, Sigma, other_pars)
  }

  return(us)
}

