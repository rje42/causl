// [[Rcpp::depends("Rcpp")]]
#include <Rcpp.h>

// static double const log2pi = std::log(2.0 * M_PI);

// C++ function to compute density at points of Gaussian copula
// [[Rcpp::export]]
Rcpp::NumericMatrix cumsum_mat(Rcpp::NumericMatrix x) {

  for (int i=1; i < x.ncol(); i++) {
    for (int j=0; j < x.nrow(); j++) {
      x(j,i) += x(j,i-1);
    }
  }

  return x;
}


/*** R
N <- 100; p <- 3
x <- matrix(rnorm(N*p), N, p)
M <- rWishart(N, df=10, Sigma=diag(p))
M <- unlist(apply(M, 3, function(x) list(x/sqrt(diag(x)*rep(diag(x), each=nrow(x))))), recursive = FALSE)
dmvnrm_arma(x, M, TRUE)
*/
