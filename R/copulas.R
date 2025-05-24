##' Sample from multivariate copulas
##'
##' @param n sample size
##' @param Sigma in which each slice is a correlation matrix
##'
##' @details Quicker than rCopula.
##'
##' Note that `rfgmCopula` only works for \eqn{d = 2}.
##'
##' @return A vector of the simulated random variables.
##' @name sample_copulas
NULL

##' @describeIn sample_copulas Gaussian copula
##' @export
rGaussCop <- function(n, Sigma) {
  d <- nrow(Sigma)
  N <- length(Sigma)/d^2
  dim(Sigma) <- c(d,d,N)

  U <- apply(Sigma, 3, chol.default)
  dim(U) <- dim(Sigma)

  if (N < n) {
    U <- U[,,rep_len(seq_len(N), n),drop=FALSE]
  }

  x <- rnorm(d*n)
  x2 <- U*x[patternRepeat0(which = c(1,3), n = c(d,d,n))]
  x2 <- t(apply(x2, c(2,3), sum))

  ## return normalised values
  pnorm(x2)
}

##' @describeIn sample_copulas t-copula
##' @param df degrees of freedom
##' @export
rtCop <- function(n, Sigma, df) {
  if (length(dim(Sigma)) == 3) {
    d <- nrow(Sigma)
    out <- matrix(NA, n, d)
    for (i in seq_len(n)) out[i,] <- c(pt(rmvt(1, sigma=Sigma[,,i], df=df), df=df))

    out
  }
  else pt(rmvt(n, Sigma, df=df), df=df)
}


##' @describeIn sample_copulas FGM-copula
##'
##' @param d dimension of copula
##' @param alpha (vector of) parameter values
##'
##' @export
rfgmCopula <- function(n, d=2, alpha)
{
  stopifnot(d == 2, -1 <= alpha, alpha <= 1)

  if (n < 1) return(matrix(NA, 0, d))

  ## note: this is *only* for the case where S = {1,...,d}
  ##       => the FGM has only one parameter in this case
  ##       see Jaworski, Durante, Haerdle, Rychlik (2009, p. 19)
  U <- matrix(runif(n*d), nrow=n, ncol=d)
  B <- alpha * apply(1-2*U[,-d, drop=FALSE], 1, prod)
  C <- sqrt( (1 + B)^2 - 4 * B * U[,d])
  U[,d] <- 2 * U[,d] / (1 + B + C)

  U
}


##' Density of a multivariate copula
##'
##' @param x samples on (0,1)
##' @param Sigma collection of matrices
##' @param log logical: return log-density?
##' @param use_cpp logical: use the C routine?
##' @param N optional integer for number of covariance matrices
##'
##' @details Computes the density for data from a
##' Gaussian or t-copula.  Currently `use_cpp` only
##' works for `dGaussCop`.
##'
##' @name copula_density
NULL

##' @describeIn copula_density Gaussian copula
##' @export
dGaussCop <- function(x, Sigma, log=FALSE, use_cpp=TRUE, N) {
  d <- ncol(x)
  if (d <= 1) {
    if (!log) return(1)
    else return(0)
  }
  n <- nrow(x)
  if (missing(N)) N <- length(Sigma)/d^2
  dim(Sigma) <- c(d,d,N)

  if (any(x > 1) || any(x < 0)) {
    stop("x's outside valid range")
  }

  ## if matrix invalid, then return 0/-Inf
  if (any(is.na(Sigma))) {
    if (log) return(-Inf)
    else return(0)
  }

  ## use the C++ implementation
  if (use_cpp) {
    ## if all the same matrix, use single copula implementation
    if (N == 1) {
      dim(Sigma) <- c(d,d)
      qx <- qnorm(x)
      qx[qx > 1e100] <- 1e100
      qx[qx < -1e100] <- -1e100
      out <- c(dGcop(qx, Sigma, logd = TRUE))
      if (log) return(out)
      else return(exp(out))
    }
    else if (N != n) {
      Sigma <- Sigma[,,rep_len(seq_len(N), n)]
      if (N > n) warning("Some matrices not used")
      else if (N < n) warning("Recycling used")
    }
    # Sigma <- unlist(apply(Sigma, 3, list), recursive = FALSE)
    qx <- qnorm(x)
    qx[qx > 1e100] <- 1e100
    qx[qx < -1e100] <- -1e100
    out <- c(dGcop_sig(qx, Sigma, logd=TRUE))
    # out <- out - rowSums(dnorm(qnorm(x), log=TRUE))

    if (log) return(c(out))
    else return(c(exp(out)))
  }



  U <- tryCatch(apply(Sigma, 3, chol.default), error= function(e) e)
  if (inherits(U, "error")) {
    if (log) return(-Inf)
    else return(0)
  }
  dim(U) <- dim(Sigma)

  if (N < n) {
    U <- U[,,rep_len(seq_len(N), n),drop=FALSE]
  }
  U <- unlist(apply(U, 3, list), recursive=FALSE)
  x2 <- unlist(apply(qnorm(as.matrix(x)), 1, list), recursive = FALSE)
  rss <- colSums(mapply(backsolve, U, x2, transpose=TRUE)^2)
  dets <- sapply(U, function(v) sum(log(diag(v))))
  if (any(is.na(dets))) {
    if (log) return(-Inf)
    else return(0)
  }

  out <- - rss / 2 - dets - d*log(2*pi)/2
  out <- out - rowSums(dnorm(qnorm(x), log=TRUE))

  if (!log) out <- exp(out)

  out
}

##' @describeIn copula_density t-Copula density
##' @param df degrees of freedom
##' @export
dtCop <- function(x, Sigma, df, log=FALSE, use_cpp=TRUE) {
  if (missing(df)) stop("Degrees of freedom not specified")
  if (!is.null(ncol(x))) d <- ncol(x)
  else {
    d <- length(x)
    x <- matrix(x, ncol=d)
  }
  n <- length(x)/d
  N <- length(Sigma)/d^2
  dim(Sigma) <- c(d,d,N)

  U <- apply(Sigma, 3, chol.default)
  dim(U) <- dim(Sigma)

  if (N < n) {
    U <- U[,,rep_len(seq_len(N), n)]
  }
  U <- unlist(apply(U, 3, list), recursive=FALSE)
  x2 <- unlist(apply(qnorm(x), 1, list), recursive = FALSE)
  rss <- colSums(mapply(backsolve, U, x2, transpose=TRUE)^2)
  dets <- sapply(U, function(x) sum(log(diag(x))))

  out <- lgamma((d + df)/2) - (lgamma(df/2) + dets +
                                 d/2 * log(pi * df)) - 0.5 * (df + d) * log1p(rss/df)
  out <- out - rowSums(dt(qt(x, df=df), df=df, log=TRUE))

  if (!log) out <- exp(out)

  out
}


##' @param alpha parameter for copula
##'
##' @return numeric vector of densities
##'
##' @describeIn copula_density bivariate FGM copula
##' @export
dfgmCopula <- function(x, alpha, log=FALSE) {
  out <- 1 + alpha * (1 - 2 * x[,1]) * (1 - 2 * x[,2])
  if (log) return(log(out))
  else return(out)
}

# ##' Distribution Function of a Bivariate Copula
# ##'
# ##' Copy of `BiCopPDF` from `VineCopula` package,
# ##' but made faster by skipping checks.
# ##'
# ##' @seealso \link[package=VineCopula]{BiCopCDF}
# ##' @export BiCopPDF
# BiCopPDF <- function (u1, u2, family, par, par2 = 0, obj = NULL, check.pars = TRUE)
# {
#   n <- length(u1)
#   args <- list(par2 = rep(par2, length(par)), family=rep(family, length(par)))
#   # args <- preproc(c(as.list(environment()), call = match.call()),
#   #                 check_u, fix_nas, check_if_01, extract_from_BiCop, match_spec_lengths,
#   #                 check_fam_par)
#   list2env(args, environment())
#
#   if (length(par) == 1) {
#     coplik <- .C("PDF_seperate", as.integer(family), as.integer(n),
#                  as.double(u1), as.double(u2), as.double(par), as.double(par2),
#                  as.double(rep(0, n)), PACKAGE = "VineCopula")[[7]]
#   }
#   else {
#     coplik <- .C("PDF_seperate_vec", as.integer(family),
#                  as.integer(n), as.double(u1), as.double(u2), as.double(par),
#                  as.double(par2), as.double(rep(0, n)), PACKAGE = "VineCopula")[[7]]
#   }
#   # out <- reset_nas(coplik, args)
#   # out
#   coplik
# }


##' Density of a Mixed Copula
##'
##' @param x matrix of samples on (0,1)
##' @param m number of discrete variables
##' @param eta eta matrix
##' @param Sigma collection of matrices
##' @param log logical: return log=density?
##' @param use_cpp logical: use the C routine?
##'
##' @return numeric vector of densities
##' @importFrom mvtnorm pmvnorm
##'
##' @export
dGaussDiscCop <- function(x, m, Sigma, eta, log=FALSE, use_cpp=TRUE) {

  if(is.null(dim(x)) || length(dim(x)) != 2) stop("x must be a matrix-like object")
  if(is.null(dim(Sigma))) stop("Sigma must be a matrix-like object")

  ## if no truncation points given, then assume just a Gaussian copula
  # if(is.null(trunc)) return(dGaussCop(x, Sigma, log=log, use_cpp=use_cpp))
  ## check that trunc values are valid
  # if (any(is.na(unlist(trunc)))) stop("NA or NaN in trunc values")
  # if (any(unlist(trunc) > 1 | unlist(trunc) < 0)) stop("trunc values not in [0,1]")

  d <- ncol(x)
  # m <- length(trunc)
  # if (any(sapply(trunc, is.matrix))) {
  #   useC2 <- FALSE
  #   mv_trunc <- TRUE
  # }
  # else useC2 <- useC
  n <- nrow(x)
  N <- length(Sigma)/d^2
  dim(Sigma) <- c(d,d,N)

  if (any(x[,seq_len(d-m)] <= 0) ||
      any(x[,seq_len(d-m)] >= 1) ||
      any(x[,d-m+seq_len(m)] < 0)) {
    stop("x's outside valid range")
  }

  ## if matrix invalid, then return 0/-Inf
  if (any(is.na(Sigma))) {
    if (log) return(-Inf)
    else return(0)
  }


  ## check if all sigma matrices are the same
  Sigma2 <- rep(Sigma[,,1],dim(Sigma)[3])
  dim(Sigma2) <- dim(Sigma)
  same_mat <- prod(Sigma2 == Sigma)


  ## use the C++ implementation
  if (use_cpp) {
    ## transform to standard normal, but ensure that discrete variables are not transformed
    x[,seq_len(dim(Sigma)[2] - m)] <- qnorm(x[,seq_len(dim(Sigma)[2] - m)])

    ## if all the same matrix, use single copula implementation
    if (same_mat) {
      Sigma1 <- Sigma[,,1]
      out <- c(dGDcop2(x = x, sigma = Sigma1, eta = eta, q = m, logd = TRUE))
      if (log) return(out)
      else return(exp(out))
    }
    else {
      out <- c(dGDcop2_sig(x = x, sigma = Sigma, eta = eta, q = m, logd = TRUE))
      if (log) return(out)
      else return(exp(out))
    }

    # else if (N != n) {
    #   Sigma <- Sigma[,,rep_len(seq_len(N), n)]
    #   if (N > n) warning("Some matrices not used")
    #   else if (N < n) warning("Recycling used")
    # }
    # # Sigma <- unlist(apply(Sigma, 3, list), recursive = FALSE)
    # out <- c(dGcop_sig(qnorm(x), Sigma, logd=TRUE))
    # # out <- out - rowSums(dnorm(qnorm(x), log=TRUE))
    #
    # if (log) return(c(out))
    # else return(c(exp(out)))
  }
  # else {
  #   stop("Not implemented in R")
  # }

  SigmaC <- apply(Sigma, 3, schur, seq_len(m)+d-m, seq_len(m)+d-m, seq_len(d-m))
  dim(SigmaC) <- c(m,m,N)

  xd <- x[,d-m+seq_len(m),drop=FALSE]

  ## if different truncation values for each variable, then deal with this
  # if (mv_trunc) {
  #   truncM <- abind::abind(trunc, along=3)
  #   cumPM <- apply(truncM, c(1,3), function(x) c(0,cumsum(x))) # cumulative probabilities
  #   upper <- lower <- xd
  #   idx <- cbind(c(xd+1), seq_len(n), rep(seq_len(m), each=n))
  #   lower[] <- cumPM[idx]
  #   idx[,1] <- idx[,1] + 1
  #   upper[] <- cumPM[idx]
  #   lower <- qnorm(pmax(0,lower))  ## put on normal scale
  #   upper <- qnorm(pmin(1,upper))
  #   dim(lower) <- dim(upper) <- c(n,m)
  #   lower <- lapply(asplit(lower, 1), c)
  #   upper <- lapply(asplit(upper, 1), c)
  #
  #   if (N == n) {
  #     ## ALLOW THIS TO DEAL WITH A DISCRETE OUTCOME AS WELL
  #     SigmaC <- asplit(SigmaC, 3)
  #     SigmaM <- asplit(Sigma[d-m+seq_len(m),d-m+seq_len(m),,drop=FALSE], 3)
  #     out <- mapply(function(x,y,z) log(pmvnorm(x,y,sigma=z)), lower, upper, SigmaC) +
  #       - mapply(function(x,y,z) log(pmvnorm(x,y,sigma=z)), lower, upper, SigmaM)
  #   }
  #   else if (N == 1) {
  #     SigmaM <- Sigma[d-m+seq_len(m),d-m+seq_len(m),1]
  #     dim(SigmaC) <- dim(SigmaM) <- c(m,m)
  #     out <- mapply(function(x,y) log(pmvnorm(x,y,sigma=SigmaC)), lower, upper) +
  #       - mapply(function(x,y) log(pmvnorm(x,y,sigma=SigmaM)), lower, upper)
  #   }
  # }

  if (d > m) {
    rest <- dGaussCop(x=x[,seq_len(d-m),drop=FALSE], Sigma[seq_len(d-m),seq_len(d-m),,drop=FALSE], log = TRUE, use_cpp=use_cpp)
    out <- out + rest
  }

  if (!log) out <- exp(out)

  out
}

##' Vectorized conditional copula function
##'
##' @param U matrix of quantiles
##' @param copula family of copula to use
##' @param param vector of parameters
##' @param par2 Degrees of freedom for t-copula
##' @param inverse should inverse CDF be returned?
##'
##' @details Should have \code{nrow(U) = length(param)}.
##' @importFrom copula cCopula
##'
cVCopula <- function (U, copula, param, par2, inverse=FALSE) {
  ## check param has right length
  if (length(param) != nrow(U)) {
    if (length(param) == 1) param <- rep_len(param, nrow(U))
    else stop("'param' should have single entry or one for each row of 'U'")
  }
  ## get list of copulas
  if (missing(par2)) {
    cops <- lapply(param, copula)
  } else {
    cops <- lapply(param, function(x) copula(x, df=par2))
  }

  splU <- apply(U, 1, FUN = function(x) x, simplify = FALSE)
  out <- mapply(function (x,y) cCopula(x,y,inverse=inverse), splU, cops)

  if (is.matrix(out)) out <- t(out)

  return(out)
}
