##' Empirical CDF
##'
##' @param xo vector of observed values
##' @param inv logical: should inverse CDF be returned
##' @param zero value to use as inverse for zero
##'
##' @details `zero` can be `"min"` (for the minimum value of `xo`), `"m1"` (for
##' the minimum minus 1) or `"mInf"` (for `-Inf`).
##'
##' @export
emp_cdf <- function (xo, inv=FALSE, zero="min") {
  xo <- sort.int(xo)
  n <- length(xo)

  ## obtain function
  if (inv) {
    if (zero == "m1") z <- min(xo) - 1
    else if (zero == "min") z <- min(xo)
    else if (zero == "mInf") z <- -Inf
    else stop("invalid entry for 'zero'")
    xo <- c(z, xo)

    FUN <- function (x) {
      xo[locations(x, seq_len(n)/n)+1]
    }
  }
  else {
    FUN <- function (x) {
      locations(x, xo)/n
    }
  }

  return(FUN)
}


##' Empirical copula
##'
##' @param u matrix of integral probability transformed values
##' @param pts index of points at which to define copula
##' @param smoothing method of smoothing to use
##'
##' @export
emp_cop <- function (u, pts=u, smoothing=c("none", "checkerboard")) {

  smoothing <- smoothing[1]

  d <- ncol(u)
  if (d != ncol(pts)) stop("All points must have the same dimension")
  if (!(smoothing %in% c("none", "checkerboard"))) stop("Argument 'smoothing' must be 'none' or 'checkerboard'")

  # pts <- apply(pts, 2, sort.int)
  N <- nrow(pts)

  rks <- NA_real_*pts

  # ## recode this in a less stupid way
  # for (i in seq_len(d)) {
  #   rks[,i] <- colSums(outer(u[,i], pts[,i], `<=`))/N
  # }

  # pts <- rbind(0, pts)

  function (U, mode="copula", ...) {
    if (!is.matrix(U)) U <- matrix(U, nrow=1)
    if (mode == "copula") {
      out <- numeric(nrow(U))

      for (j in seq(nrow(U))) {
        U2 <- rep(U[j,], each=nrow(pts))

        v <- U2 > pts
        emp <- apply(v, 1, all)
        out[j] <- mean(emp)
      }
      return(out)
    }
    else if (mode == "h") {
      oargs <- list(...)


    }
  }
}


