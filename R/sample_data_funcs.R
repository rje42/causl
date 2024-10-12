##' Rescale quantiles to arbitrary random variable.
##'
##' @param U vector of quantiles
##' @param X model matrix of covariates
##' @param pars list of parameters (see details)
##' @param family family of distributions to use
##' @param link link function
##'
##' @details `family` can be 1, 2, 3, 4 or 5 for Gaussian, t-distributed,
##' Gamma distributed, beta distributed or discrete respectively, and 11 for
##' ordinal variables.
##' `pars` should be a list with entries `beta` and
##' `phi`, as well as possibly `par2`, `trunc` and `nlevel` if the family
##' is set to 2 or 5.
##' `U` should have the same length as `X` has rows, and
##' `X` should have the same number of columns as the length of
##' \code{pars$beta}.
##'
##' @return vector of rescaled variables
##'
##' @export
rescale_var <- function(U, X, pars, family=1, link) {

  # cat("fam=", family, ", link=", link, "\n", sep="")

  ## get linear component
  eta <- X %*% pars$beta
  phi <- pars$phi

  if (is.numeric(family)) {
    if (is.null(phi) && family %in% 1:3) stop("Variance parameter should be specified")

    ## make U normal, t or gamma
    if (family == 1) {
      if (link == "identity") Y <- qnorm(U, mean = eta, sd=sqrt(phi))
      else if (link == "log") Y <- qnorm(U, mean = exp(eta), sd=sqrt(phi))
      else if (link == "inverse") Y <- qnorm(U, mean = 1/eta, sd=sqrt(phi))
      else stop("invalid link function for Gaussian distribution")
    }
    else if (family == 6) {
      stop()
      if (link == "identity") Y <- qlnorm(U, meanlog = eta, sdlog=sqrt(phi))
      else if (link == "exp") Y <- qlnorm(U, meanlog = log(eta), sdlog=sqrt(phi))
      else stop("invalid link function for log-normal distribution")
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
      if (link == "logit") mu <- plogis(eta)
      else if (link == "probit") mu <- pnorm(eta)
      else stop("invalid link function for beta distribution")

      Y <- qbeta(U, shape1=1+phi*mu, shape2=1+phi*(1-mu))
    }
    else if (family == 0 || family == 5) {
      if (link == "probit") Y <- 1*(eta + qnorm(U) > 0)
      else if (link == "logit") Y <- 1*(eta + qlogis(U) > 0)
      else if (link == "log") Y <- 1*(eta - log(U) > 0)
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
    else if (family == 10) {
      if (link == "logit") mu <- theta_to_p_cat(eta)
      else stop("invalid link function for ordinal distribution")

      Y <- factor(1 + colSums(apply(mu, 1, cumsum) < rep(U, each=ncol(mu))), levels=seq_len(ncol(mu)))
    }
    else if (family == 11) {
      if (link == "logit") mu <- theta_to_p_ord(eta)
      else stop("invalid link function for ordinal distribution")

      Y <- factor(1 + colSums(apply(mu, 1, cumsum) < rep(U, each=ncol(mu))), levels=seq_len(ncol(mu)))
    }
    else stop("family must be between 0 and 5")
  }
  else if (is(family, "causl_family")) {
    mu <- link_apply(eta, link, family$name)

    if (family$name %in% c("binomial")) {
      if (link == "probit") mu <- 1*(eta + qnorm(U) > 0)
      else if (link == "logit") mu <- 1*(eta + qlogis(U) > 0)
      else stop("invalid link function for binomial distribution")
    }

    pars2 <- c(list(p=U, mu=mu), pars[names(pars) != "beta"])

    Y <- do.call(family$qdist, pars2)
  }
  else stop("'family' should be an integer or 'causl_family' object")

  ### get Z values to correct families
  # nms <- names(dat)[grep("z", names(dat))]

  return(Y)
}

##' @importFrom lifecycle deprecate_warn deprecate_soft
##' @describeIn rescale_var Old name, now deprecated
##' @export
rescaleVar <- function(U, X, pars, family=1, link) {
  deprecate_soft("0.8.0", "rescaleVar()", "rescale_var()")
  rescale_var(U, X, pars, family=family, link)
}

link_apply <- function(eta, link, family_nm) {
  if (is.function(link)) {
    mu <- link(eta)
  }
  else if (family_nm == "categorical") {
    if (link == "logit") mu <- theta_to_p_cat(eta)
    else stop("Not a valid link function for categorical variable")
  }
  else if (family_nm == "ordinal") {
    if (link == "logit") mu <- theta_to_p_ord(eta)
    else stop("Not a valid link function for ordinal variable")
  }
  else {
    if (link == "identity") mu <- eta
    else if (link == "log") mu <- exp(eta)
    else if (link == "inverse") mu <- 1/eta
    else if (link == "logit") mu <- plogis(eta)
    else if (link == "probit") mu <- pnorm(eta)
    else stop("Not a valid link function")
  }

  return(mu)
}


##' Rescale quantiles to conditional copula
##'
##' @param U matrix of quantiles
##' @param X model matrix of covariates
##' @param beta list of parameters (see details)
##' @param family variety of copula to use
##' @param par2 additional parameter for some copulas
## @param link link function
##'
##' @details The variable to be transformed must be in the final column of
##' `U`, with variables being conditioned upon in the earlier columns.
##'
##' `family` can be 1 for Gaussian, 2 for t, 3 for Clayton, 4 for
##' Gumbel, 5 for Frank, 6 for Joe and 11 for FGM copulas. Gamma distributed,
##' beta distributed or discrete respectively. `pars` should be a list
##' with entries `beta` and `phi`, as well as possibly `par2` if
##' \code{family=2}.
##' `U` should have the same length as `X` has rows, and `X`
##' should have the same number of columns as the length of \code{pars$beta}.
##'
##' @return vector of rescaled quantiles
##'
##' @importFrom copula normalCopula tCopula
## @inherit glm_sim
##'
##' @export
rescale_cop <- function(U, X, beta, family=1, par2) {

  if (!is.matrix(U)) stop("'U' should be a matrix")
  if (!is.matrix(X)) stop("'X' should be a matrix")
  if (nrow(U) != nrow(X)) stop("'U' and 'X' should have same number of rows")

  ## get linear component
  eta <- X %*% beta

  ## make U normal, t or gamma
  if (family == 1) {
    # if (link == "tanh")
    param <- 2*expit(eta) - 1
    # browser()
    # Y <- cVCopula(U, copula = normalCopula, param = param, inverse=TRUE)
    Y <- pnorm(qnorm(U[,2])*sqrt(1-param^2)+param*qnorm(U[,1]))
  }
  else if (family == 2) {
    # Y <- sqrt(phi)*qt(U, df=pars$par2) + eta
    param <- 2*expit(eta) - 1
    Y <- cVCopula(U, copula = tCopula, param = param, par2=par2, inverse=TRUE)
  }
  else if (family == 3) {
    param <- exp(eta) - 1
    Y <- cVCopula(U, copula = claytonCopula, param = param, inverse=TRUE)
  }
  else if (family == 4) {
    param <- exp(eta) + 1
    Y <- cVCopula(U, copula = gumbelCopula, param = param, inverse=TRUE)
  }
  else if (family == 5) {
    param <- eta
    Y <- cVCopula(U, copula = frankCopula, param = param, inverse=TRUE)
  }
  else if (family == 6) {
    param <- exp(eta) + 1
    Y <- cVCopula(U, copula = joeCopula, param = param, inverse=TRUE)
  }
  else stop("family must be between 0 and 5")

  ### get Z values to correct families
  # nms <- names(dat)[grep("z", names(dat))]

  return(Y[,ncol(Y)])
}

##' @describeIn rescale_cop Old name, now deprecated
##' @export
rescaleCop <- function(U, X, beta, family=1, par2) {
  deprecate_soft("0.8.0", "rescaleCop()", "rescale_cop()")
  rescale_cop(U, X, beta, family=family, par2)
}

##' Simulate copula values
##'
##' @param dat data frame with empty columns
##' @param family numeric indicator of copula type
##' @param par mandatory parameters
##' @param par2 optional parameters
##' @param model_matrix design matrix for covariates
##'
##' @details Returns data frame containing columns `y`
##' and \code{z1, ..., zk}.
##'
##' The family variables are numeric and taken from `VineCopula`.
##' Use, for example, 1 for Gaussian, 2 for t, 3 for Clayton, 4 for Gumbel,
##' 5 for Frank, 6 for Joe and 11 for FGM copulas.
##'
##' @return A data frame of the same dimension as `dat` containing the
##' simulated values.
##'
##' @export
sim_copula <- function(dat, family, par, par2, model_matrix) {

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

##' @describeIn sim_copula Old name, now deprecated
##' @export
sim_CopVal <- function(dat, family, par, par2, model_matrix) {
  deprecate_soft("0.8.0", "sim_CopVal()", "sim_copula()")
  sim_copula(dat, family, par, par2, model_matrix)
}


# #' @describeIn glm_sim Old name
# #' @inherit glm_sim
# sim_glm <- function (family, eta, phi, par2, link) {
#   glm_sim(family, eta, phi, other_pars=list(par2=par2), link=link)
# }

##' Simulate from a GLM
##'
##' Simulate values from some generalized linear models
##'
##' @inherit rejectionWeights
##' @inherit get_X_density
##' @param quantiles logical indicating whether to return quantiles
##' @param other_pars list of other parameters for specified family
##'
##' @export
glm_sim <- function (family, eta, phi, other_pars, link, quantiles=TRUE) {

  if (is.matrix(eta)) n <- nrow(eta)
  else n <- length(eta)
  if (missing(link)) {
    if (is(family, "causl_family")) link <- family$link
    else {
      ## get the default link for this family
      link <- family_vals[family_vals$val==family,2]
    }
  }

  if (is(family, "causl_family")) {
    if (!(family$name %in% names(links_list))) {
      ## this is not a registered family
      # if (!(link %in% family$links_list)) stop(paste0(link, " is not a valid link function for ", family$name, " family"))
      # stop(paste0("Family ", family$name, " is not a valid and registered family"))

      if (!(link %in% unlist(links_list)) &&
          !(link %in% names(family$custom_links))) stop("Link function not known")
    }
    else if (!(link %in% links_list[[family$name]])) stop(paste0(link, " is not a valid link function for ", family$name, " family"))

    if (family$name %in% c("categorical","ordinal")) {
      if (ncol(eta) != other_pars$nlevel - 1) stop("Invalid 'eta' input to glm_sim")
      if (link == "logit") {
        if (family$name == "categorical") {
          mu <- theta_to_p_cat(eta)
        }
        else if (family$name == "ordinal") {
          mu <- theta_to_p_ord(eta)
        }
      }
      else stop("We shouldn't get here")
    }
    else {
      if (link %in% unlist(links_list)) {
        if (link=="identity") mu <- eta
        else if (link=="inverse") mu <- 1/eta
        else if (link=="log") mu <- exp(eta)
        else if (link=="logit") mu <- expit(eta)
        else if (link=="probit") mu <- pnorm(eta)
        else stop("We shouldn't get here")
      }
      else {
        if (link %in% names(family$custom_links)) mu <- family$custom_links[[link]]$linkinv(eta)
        else stop("Link function not understood")
      }
    }

    pars <- list(mu=mu)
    if ("phi" %in% family$pars) pars <- c(pars, list(phi=phi))
    if ("par2" %in% family$pars) pars <- c(pars, list(par2=other_pars$par2))

    x <- do.call(family$rdist, c(list(n=n), pars))
    if (quantiles && !(family$name %in% c("categorical","ordinal"))) qx <- do.call(family$pdist, c(list(x=x), pars))
  }
  else if (is.numeric(family)) {

    ## old style approach
    ## get the densities for x
    if (family == 1 || family == 6) {
      if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else if (link=="log") mu <- exp(eta)
      else stop("Not a valid link function for the Gaussian distribution")

      x <- rnorm(n, mu, sd=sqrt(phi))
      qx <- pnorm(x, mu, sd=sqrt(phi))

      if (family == 6) {
        out <- exp(out)
        qx <- qx/out
      }
    }
    else if (family == 2) {
      par2 <- other_pars$par2

      if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else if (link=="log") mu <- exp(eta)
      else stop("Not a valid link function for the t-distribution")

      x <- rt(n, df=par2)*sqrt(phi) + mu
      qx <- pt((x - mu)/sqrt(phi), df=par2)
    }
    else if (family == 3) {
      if (link=="log") mu <- exp(eta)
      else if (link=="identity") mu <- eta
      else if (link=="inverse") mu <- 1/eta
      else stop("Not a valid link function for the gamma distribution")

      x <- rgamma(n, shape=1/phi, scale=phi*mu)
      qx <- pgamma(x, shape=1/phi, scale=phi*mu)
    }
    else if (family == 4) {
      if (link=="logit") mu <- expit(eta)
      else if (link=="probit") mu <- pnorm(eta)
      else stop("Not a valid link function for the beta distribution")

      th1 <- mu/phi
      th2 <- (1-mu)/phi

      x <- rbeta(n, th1, th2)
      qx <- pbeta(x, th1, th2)
    }
    else if (family == 5) {
      if (link=="logit") {
        mu <- expit(eta)
        z <- rlogis(n)
        qx <- plogis(z)
      }
      else if (link=="probit") {
        mu <- pnorm(eta)
        z <- rnorm(n)
        qx <- pnorm(z)
      }
      else if (link=="log") {
        mu <- exp(eta)
        z <- rexp(n)
        qx <- pexp(z)
      }
      else stop("Not a valid link function for the Bernoulli distribution")

      # x <- rbinom(n, size=1, prob=mu)
      x <- eta + z > 0
      # qx <- dbinom(x, size=1, prob=mu)
    }
    else if (family == 10) {
      if (link == "logit") mu <- theta_to_p_cat(eta)
      else stop("Not a valid link function for a categorical distribution")

      qx <- runif(n)
      qsum <- cumsum_mat(mu)

      x <- rowSums(qsum < qx) + 1
      x <- factor(x, levels=seq_len(ncol(mu)))
    }
    else if (family == 11) {
      if (link == "logit") mu <- theta_to_p_ord(eta)
      else stop("Not a valid link function for an ordinal distribution")

      qx <- runif(n)
      qsum <- cumsum_mat(mu)

      x <- rowSums(qsum < qx) + 1
      x <- factor(x, levels=seq_len(ncol(mu)))
    }
    else stop("Only Gaussian, t, beta, gamma, log-normal and categorical distributions are allowed")
  }
  else stop("family input should be an integer or 'causl_family' function")

  ## return quantile
  if (is.numeric(family) || !(family$name %in% c("categorical","ordinal"))) {
    if (quantiles) attr(x, "quantile") <- qx
  }

  return(x)
}

