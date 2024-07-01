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


##' @describeIn causl_copula getter copula family
##' @export
get_copula <- function(family_index, link = NULL){
  if(is.null(link)){
    if(family_index == 1){
      return(gaussian_causl_cop())
    }
    else if(family_index == 2){
      return(t_causl_cop())
    }
  }
  else{
    if(family_index == 1){
      return(gaussian_causl_cop(link))
    }
    else if(family_index == 2){
      return(t_causl_cop(link))
    }
  }


  
  
}

##' @describeIn causl_copula Gaussian copula family
##' @export
gaussian_causl_cop <- function (link) {
  if (missing(link)) link = tanh

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
t_causl_cop <- function(link){
  if (missing(link)) link = tanh
  
  ## write functions
  dens <- function (x, Sigma, log=FALSE) dGaussCop(x, Sigma=Sigma, log=log)
  sim <- function (n, Sigma, df) rtCop(n, Sigma=Sigma, df = df)
  # probs <- function (x, mu, phi) pnorm(x, Sigma=Sigma)
  
  default <- function() list(x=c(0,0,0), Sigma=diag(3), d=3)
  
  ## define family
  out <- list(name="t", ddist=dens, rdist=sim,
              pars=c("Sigma"), default=default, link=link)
  class(out) <- "causl_copula"
  
  return(out)
  
}



##' @describeIn causl_copula simulate from copula family
##' @export
sim_cop <- function(causl_copula, beta_matrix, other_pars, model_matrix){
  n <- nrow(model_matrix)

  p <- dim(beta_matrix)[2]
  link_beta <- causl_copula$link(beta_matrix)
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
    fact <- factor(rep(1,nrow(dat)))
  }
  
  for (gp in levels(fact)) {
    gp_mem <- (fact == gp)
    n_gp <- sum(gp_mem)
    if (n_gp == 0) next
    
    eta <- c(model_matrix[match(TRUE, gp_mem), ,drop=FALSE] %*% link_beta)
    SIGMA <- matrix(1, nrow = p, ncol = p)
    SIGMA[upper.tri(SIGMA)] <- eta
    SIGMA[lower.tri(SIGMA)] <- t(SIGMA[upper.tri(SIGMA)])
    eigen_values <- eigen(SIGMA)$values
    if(any(eigen_values < 0))stop("specified non positive semi-definite correlation matrix")
    us[gp_mem,] <- causl_copula$rdist(n_gp, SIGMA, other_pars)
  

  }
  
  return(us)
}


