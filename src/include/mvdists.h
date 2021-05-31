#include <R.h>

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);

arma::vec dGcop(arma::mat const &x,
                arma::mat const &sigma,
                bool const logd);

arma::vec dGcop_sig(arma::mat const &x,
                    arma::cube const &sigma,
                    bool const logd);

arma::vec dGDcop(arma::mat const &x,
                 arma::mat const &sigma,
                 Rcpp::List trunc,
                 bool const logd);


