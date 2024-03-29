% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_vinecop.R
\name{links_cop}
\alias{links_cop}
\title{List of link functions for copuulas}
\usage{
links_cop(family)
}
\arguments{
\item{family}{either the name of the family or its integer representative}
}
\description{
List of link functions for copuulas
}
\details{
This function returns the default link function for each possible copula.
These are given in the table:\tabular{rlll}{
   value \tab family \tab link \tab link name \cr
   1 \tab gaussian \tab logit((1+rho)/2) \tab logit2 \cr
   2 \tab t \tab logit((1+rho)/2) \tab logit2 \cr
   3 \tab Clayton \tab log(1+theta) \tab log1p \cr
   4 \tab Gumbel \tab log(theta - 1) \tab log1m \cr
   5 \tab Frank \tab logit(1 + theta) \tab log1p \cr
   6 \tab Joe \tab identity \tab identity \cr
   11 \tab FGM \tab logit((1+rho)/2) \tab logit2 \cr
}
}
