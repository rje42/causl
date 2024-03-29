% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copulas.R
\name{dGaussCop}
\alias{dGaussCop}
\alias{dtCop}
\title{Density of a Gaussian or t-Copula}
\usage{
dGaussCop(x, Sigma, log = FALSE, useC = TRUE, N)

dtCop(x, Sigma, df, log = FALSE, useC = TRUE)
}
\arguments{
\item{x}{samples on (0,1)}

\item{Sigma}{collection of matrices}

\item{log}{logical: return log=density?}

\item{useC}{logical: use the C routine?}

\item{N}{optional integer for number of covariance matrices}

\item{df}{degrees of freedom}
}
\description{
Density of a Gaussian or t-Copula
}
\details{
Computes the density for data from a
Gaussian or t-copula.  Currently \code{useC} only
works for \code{dGaussCop}.
}
\section{Functions}{
\itemize{
\item \code{dtCop()}: t-Copula Density

}}
