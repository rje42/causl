% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_inversion.R, R/sim_rejection.R
\name{sim_inversion}
\alias{sim_inversion}
\alias{sim_rejection}
\title{Simulate for single time-step}
\usage{
sim_inversion(out, proc_inputs)

sim_rejection(out, proc_inputs, careful, control)
}
\arguments{
\item{out}{data frame for output}

\item{proc_inputs}{output of \code{process_inputs()}}

\item{careful}{should full, slower method be used?}

\item{control}{list of control parameters}
}
\description{
Simulate for single time-step
}
\details{
\code{sim_inversion} and \code{sim_rejection} correspond to
performing the sampling by inversion or using rejection sampling.
}
\section{Functions}{
\itemize{
\item \code{sim_rejection()}: Rejection sampling code

}}
