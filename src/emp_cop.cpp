#include <Rcpp.h>

// [[Rcpp::export]]
void Rpr (Rcpp::NumericVector x, int max) {
  int n = x.length();
  if (max > n) max = n;

  for (int i = 0; i < max; i++) {
    Rprintf("%f ", x(i));
  }
  Rprintf("\n");

  return;
}


// C++ function to obtain position of one sorted vector within entries of another
// [[Rcpp::export]]
Rcpp::NumericVector locations (Rcpp::NumericVector u, Rcpp::NumericVector v) {
  int nu = u.length();
  int nv = v.length();

  Rcpp::NumericVector out(nu);
  u.push_back(v(nv-1)+1);  // add large entry to end of vector

  // Rpr(u,4);

  int j = 0;

  for (int i=0; i < nv; i++) {
    // obtain all the entries that are less than the ith entry of v
    while (u(j) <= v(i)) {
      out(j) = i;
      j++;
    }
    if (j >= nu) break;
  }
  while (j < nu) {
    // anything left is larger than all entries in v
    out(j) = nv;
    j++;
  }

  return out;
}
