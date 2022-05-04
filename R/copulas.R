##' Sample from Gaussian or t-copula
##'
##' @param n sample size
##' @param Sigma in which each slice is a correlation matrix
##'
##' @details Quicker than rCopula.
##'
##' Note that \code{rfgmCopula} only works for \eqn{d = 2}.
##'
##' @return A vector of the simulated random variables.
##'
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

##' @describeIn rGaussCop t-copula
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


##' @describeIn rGaussCop FGM-copula
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


##' Density of a Gaussian or t-Copula
##'
##' @param x samples on (0,1)
##' @param Sigma collection of matrices
##' @param log logical: return log=density?
##' @param useC logical: use the C routine?
##' @param N optional integer for number of covariance matrices
##'
##' @details Computes the density for data from a
##' Gaussian or t-copula.  Currently \code{useC} only
##' works for \code{dGaussCop}.
##'
##' @export
dGaussCop <- function(x, Sigma, log=FALSE, useC=TRUE, N) {
  d <- ncol(x)
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
  if (useC) {
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
    U <- U[,,rep_len(seq_len(N), n)]
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

##' @describeIn dGaussCop t-Copula Density
##' @param df degrees of freedom
##' @export
dtCop <- function(x, Sigma, df, log=FALSE, useC=TRUE) {
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


##' Density of a bivariate FGM copula
##'
##' @param u1,u2 vector of probabilities
##' @param alpha parameter for copula
##'
##' @return numeric vector of densities
##'
##' @export
dfgmCopula <- function(u1, u2, alpha) {
  1 + alpha * (1 - 2 * u1) * (1 - 2 * u2)
}

# ##' Distribution Function of a Bivariate Copula
# ##'
# ##' Copy of \code{BiCopPDF} from \code{VineCopula} package,
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
##' @param Sigma collection of matrices
##' @param trunc list of truncation points
##' @param log logical: return log=density?
##' @param useC logical: use the C routine?
##'
##' @return numeric vector of densities
##'
##' @export
dGaussDiscCop <- function(x, Sigma, trunc, log=FALSE, useC=TRUE) {

  if(is.null(dim(x))) stop("x must be a matrix-like object")
  if(is.null(dim(Sigma))) stop("Sigma must be a matrix-like object")

  ## if no truncation points given, then assume just a Gaussian copula
  if(is.null(trunc)) return(dGaussCop(x, Sigma, log=log, useC=useC))
  ## check that trunc values are valid
  if (any(is.na(unlist(trunc)))) stop("NA or NaN in trunc values")
  if (any(unlist(trunc) > 1 || unlist(trunc) < 0)) stop("trunc values not in [0,1]")

  d <- ncol(x)
  m <- length(trunc)
  n <- nrow(x)
  N <- length(Sigma)/d^2
  dim(Sigma) <- c(d,d,N)

  if (any(x[,seq_len(d-m)] >= 1) || any(x < 0)) {
    stop("x's outside valid range")
  }

  ## if matrix invalid, then return 0/-Inf
  if (any(is.na(Sigma))) {
    if (log) return(-Inf)
    else return(0)
  }

  ## use the C++ implementation
  if (useC) {
    ## transform to standard normal, but ensure that discrete variables are not transformed
    x[,dim(Sigma)[2]-length(trunc)] <- qnorm(x[,dim(Sigma)[2]-length(trunc)])

    ## if all the same matrix, use single copula implementation
    if (N == 1) {
      dim(Sigma) <- c(d,d)
      out <- c(dGDcop(x, Sigma, trunc=trunc, logd = TRUE))
      if (log) return(out)
      else return(exp(out))
    }
    else {
      out <- c(dGDcop_sig(x, Sigma, trunc=trunc, logd = TRUE))
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
  else {
    stop("Not implemented in R")
  }

  U <- tryCatch(apply(Sigma, 3, chol.default), error= function(e) e)
  if (inherits(U, "error")) {
    if (log) return(-Inf)
    else return(0)
  }
  dim(U) <- dim(Sigma)

  if (N < n) {
    U <- U[,,rep_len(seq_len(N), n)]
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

