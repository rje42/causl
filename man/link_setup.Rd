% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/links.R
\name{link_setup}
\alias{link_setup}
\title{Set up link functions}
\usage{
link_setup(
  link,
  family,
  vars,
  sources = links_list,
  fam_list = list(family_vals)
)
}
\arguments{
\item{link}{input from \code{causl_model} object or similar}

\item{family}{the list of families for random variables}

\item{vars}{a list of vectors of variable names with the same structure as \code{family}}

\item{sources}{list of links for parametric families}

\item{fam_list}{list of data frames in the same format as \code{family_vals}}
}
\description{
\code{family} and \code{vars} should have the same structure if \code{vars} is
specified.
}
