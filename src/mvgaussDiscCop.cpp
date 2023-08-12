// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "include/mvtnorm.h"
#include "include/mvdists.h"

// #include <Rcpp.h>
// using namespace Rcpp;

void Rpr(double* x, int len) {
  // int len = x.length();
  for (int i = 0; i < len; i++) {
    Rprintf("%f ", x[i]);
  }
  Rprintf("\n");
  return;
}

void Rpri(int* x, int len) {
  // int len = x.length();
  for (int i = 0; i < len; i++) {
    Rprintf("%i ", x[i]);
  }
  Rprintf("\n");
  return;
}

// static double const log2pi = std::log(2.0 * M_PI);

using arma::mat;
using arma::uword;


// Calculates the efficient information \eqn{I_{bb}-I_{ba}I_{aa}^{-1}I_{ab}}.
//
// @param C complement to target
// @param A to condition upon
// @param B lower left part of matrix
// @return Numeric matrix.
arma::mat SchurC(const arma::mat C, const arma::mat A,
                 const arma::mat B){
  const arma::mat A_B = arma::solve(A,B.t(),arma::solve_opts::likely_sympd);
  const arma::mat C_A = join_horiz(C - B*A_B, A_B.t());
  return C_A;
}

// C++ function to compute density at points of Gaussian copula with discrete components
//  x :     matrix of values
//  sigma : correlation matrix of parameters
//  trunc : list of discretization points corresponding to final k dimensions of x/sigma
//  logd  : logical: return the log-density?
// [[Rcpp::export]]
arma::vec dGDcop2(arma::mat const &x,
                 arma::mat const &sigma,
                 Rcpp::List trunc,
                 bool const logd = false) {

  uword q = trunc.length();
  if (q == 0) return(dGcop(x, sigma, logd));

  uword const n = x.n_rows;
  //d = x.n_cols;
  arma::vec out(n);
  arma::rowvec z;
  uword const p = sigma.n_rows - trunc.length();

  // Rprintf("%i %i %i\n", p, sigma.n_rows, sigma.n_cols);

  arma::mat sigma0 = sigma.submat(0,0,p-1,p-1);
  arma::mat sigma1 = sigma.submat(p,p,sigma.n_rows-1, sigma.n_cols-1);
  arma::mat sigma10 = sigma.submat(p,0,sigma.n_rows-1, p-1);

  // Rprintf("%i %i %i\n", sigma0.n_elem, sigma1.n_elem, sigma10.n_elem);

  arma::mat sigma1_0 = SchurC(sigma1, sigma0, sigma10);
  arma::mat sigma1cov = sigma1_0.cols(0,q-1);
  arma::mat sigma1mn = sigma1_0.cols(q,sigma1_0.n_cols-1);

  // Rprintf("%i %i %i\n", sigma1_0.n_elem, sigma1cov.n_elem, sigma1mn.n_elem);


  // double const constants = -(double)d/2.0 * log2pi;

  // mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  arma::vec const eig = arma::eig_sym(sigma0);
  if (any(eig < 0)) {
    out = NA_REAL;
    // for (uword i = 0; i < n; i++) out(i) = nan;
    return out;
  }
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma0)));
  double const rootisum = arma::sum(log(rooti.diag()));

  // Rprintf("sum of log(inv_chol) = %3e\n", rootisum);

  // separate out continuous from discrete components
  for (uword i = 0; i < n; i++) {
    z = x.row(i);
    inplace_tri_mat_mult(z, rooti);
    // Rprintf("x(%i) = %3e %3e\n", i, z(0), z(1));
    out(i) = rootisum - 0.5 * (arma::dot(z, z) - arma::dot(x.row(i), x.row(i)));
  }
  // Rprintf("n = %i, p = %i, q = %i, x.n_rows = %i, x.n_cols = %i\n", n, p, q, x.n_rows, x.n_cols);
  arma::mat condmn = x.cols(0,p-1) * sigma1mn.t();
  arma::mat x2 = x.cols(p,x.n_cols-1);
  // x2 = x2 - condmn;  // this is the discrete part, with continuous mean subtracted

  // Rcpp::Rcout << "x2:\n" << x2 << std::endl;
  // Rcpp::Rcout << "condmn:\n" << condmn << std::endl;
  //
  // Rprintf("x2.n_rows = %i, x2.n_cols = %i\n", x2.n_rows, x2.n_cols);
  // Rprintf("x2(0,0) = %1.1e\n", x2(0,0));

  arma::mat wh = arma::zeros(n,q);

  //Rcpp::NumericVector upper(p), lower(p);
  arma::mat upper(n,q), lower(n,q);
  Rcpp::IntegerMatrix infin(n,q);

  // Rprintf("here...");
  //
  // Rprintf("q = %i, trunc.length() = %i\n", q, trunc.length());

  for (uword j = 0; j < q; j++) {
    // Rprintf("j = %i\n", j);
    Rcpp::NumericVector tmp = trunc[j];
    // Rprintf("tmp.length() = %i, tmp(0) = %2e, tmp(1) = %2e\n", tmp.length(), tmp(0), tmp(1));

    Rcpp::NumericVector tmp2 = Rcpp::cumsum(tmp);
    tmp2.push_front(0);
    tmp2 = qnorm(tmp2);
    // Rcpp::Rcout << "tmp2:\n" << tmp2 << std::endl;

    for (uword i = 0; i < n; i++) {
      // Rprintf("(i,j)=(%i,%i)...length(tmp) = %i...x2(i,j) = %2e, sqrt(sigma1cov(j,j))=%2e\n", i, j, tmp.length(), x2(i,j), sqrt(sigma1cov(j,j)));

      // lower and upper endpoints for intervals
      lower(i,j) = (tmp2(x(i,p+j)) + condmn(i,j))/sqrt(sigma1cov(j,j));
      upper(i,j) = (tmp2(x(i,p+j)+1) + condmn(i,j))/sqrt(sigma1cov(j,j));
      // indicates whether an end point is infinite
      infin(i,j) = 2 - 2*std::isinf(lower(i,j)) - std::isinf(upper(i,j));
      // if (infin(i,j) < 0) {
      //   Rcpp::stop("upper and lower end-points should not both be infinite");
      // }
    }
  }

  // lower.print();
  // upper.print();
  // infin.print();
  // Rcpp::Rcout << "infin:\n" << infin << std::endl;

  // Rprintf("1...");

  arma::mat correl1cov(q, q, arma::fill::eye);

  // Rprintf("2...");

  for (uword i = 0; i < q-1; i++) for (uword j = i+1; j < q; j++) {
    correl1cov(i,j) = correl1cov(j,i) = sigma1cov(i,j)/sqrt(sigma1cov(i,i)*sigma1cov(j,j));
  }

  // Rprintf("3...");

  // arma::mat cDec = arma::chol(sigma1cov);
  double tmp[q*(q-1)/2];

  // assign upper triangular elements to tmp
  int ct = 0;
  for (uword j=0; j < q; j++) {
    for (uword i=0; i < j; i++) {
      tmp[ct] = correl1cov(i,j);
      ct += 1;
    }
  }

  // initialize vector of zeros for FORTRAN function
  double zmn[q];
  for (uword j = 0; j < q; j++) {
    zmn[j] = 0.0;
  }
  int df = 0, rnd = 1, inform = 0;
  int maxpts = 25000;
  double abseps = 1E-3, releps = 0.0;

  // get limits and indicators of infinite limits
  for (uword i = 0; i < n; i++) {
    double lower2[q], upper2[q];
    int infin2[q];
    double error = 0.0, value = 0.0;

    for (uword j = 0; j < q; j++) {
      lower2[j] = lower(i,j);
      upper2[j] = upper(i,j);
      infin2[j] = infin(i,j);
    }

    // C_mvtdst2(&q);

    int q2 = q;
    // Rprintf("q = %i, df = %i\n", q2, df);
    // Rprintf("lower: %3e %3e\n", lower2[0]);
    // Rprintf("upper: %3e %3e\n", upper2[0]);
    // Rprintf("infin: %i %i\n", infin2[0]);

    // Rprintf("sigmaUT: ");
    // for (int nn=0; nn < q*(q+1)/2; nn++) Rprintf("%e ", tmp[nn]);
    // Rprintf("\n");
    // Rprintf("delta: %e %e\n", zmn[0], zmn[1]);
    // Rprintf("maxpts = %i, abseps = %e, releps = %e\n", maxpts, abseps, releps);
    // Rprintf("error = %e, value = %e, inform = %i, rnd = %i\n", error, value, inform, rnd);

    if (q == 1) {
      if (infin2[0] == 1) {
        value = 1.0 - R::pnorm(lower2[0], 0.0, 1.0, 1, 0);
      }
      else if (infin2[0] == 0) {
        value = R::pnorm(upper2[0], 0.0, 1.0, 1, 0);
      }
      else if (infin2[0] == 2) {
        value = R::pnorm(upper2[0], 0.0, 1.0, 1, 0) - R::pnorm(lower2[0], 0.0, 1.0, 1, 0);
      }
      else Rcpp::stop("infin2 should take value 0, 1 or 2");

    }
    else {
      C_mvtdst(&q2, // N
               &df,
               lower2,
               upper2,
               infin2,
               tmp,
               zmn,     // DELTA
               &maxpts,  // MAXPTS
               &abseps,  // ABSEPS
               &releps,  // RELEPS
               &error, // ERROR
               &value, // VALUE
               &inform,  // INFORM
               &rnd);  //RND
    }
    // Rprintf("value: %e\n", value);

    out(i) += log(value);
  }


  if (logd)
    return out;
  return exp(out);
}


// C++ function to compute density at points of Gaussian copula with discrete components
//  x :     matrix of values
//  sigma : correlation matrix of parameters
//  trunc : list of discretization points corresponding to final k dimensions of x/sigma
//  logd  : logical: return the log-density?
// [[Rcpp::export]]
arma::vec dGDcop(arma::mat const &x,
           arma::mat const &sigma,
           Rcpp::List trunc,
           bool const logd = false) {

  uword q = trunc.length();
  if (q == 0) return(dGcop(x, sigma, logd));

  uword const n = x.n_rows;
    //d = x.n_cols;
  arma::vec out(n);
  arma::rowvec z;
  uword const p = sigma.n_rows - trunc.length();

  Rprintf("%i %i %i\n", p, sigma.n_rows, sigma.n_cols);

  arma::mat sigma0 = sigma.submat(0,0,p-1,p-1);
  arma::mat sigma1 = sigma.submat(p,p,sigma.n_rows-1, sigma.n_cols-1);
  arma::mat sigma10 = sigma.submat(p,0,sigma.n_rows-1, p-1);

  Rprintf("%i %i %i\n", sigma0.n_elem, sigma1.n_elem, sigma10.n_elem);

  arma::mat sigma1_0 = SchurC(sigma1, sigma0, sigma10);
  arma::mat sigma1cov = sigma1_0.cols(0,q-1);
  arma::mat sigma1mn = sigma1_0.cols(q,sigma1_0.n_cols-1);

  Rprintf("%i %i %i\n", sigma1_0.n_elem, sigma1cov.n_elem, sigma1mn.n_elem);


  // double const constants = -(double)d/2.0 * log2pi;

  // mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  arma::vec const eig = arma::eig_sym(sigma0);
  if (any(eig < 0)) {
    out = NA_REAL;
    // for (uword i = 0; i < n; i++) out(i) = nan;
    return out;
  }
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma0)));
  double const rootisum = arma::sum(log(rooti.diag()));

 Rprintf("%3e\n", rootisum);

  for (uword i = 0; i < n; i++) {
    z = x.row(i);
    inplace_tri_mat_mult(z, rooti);
  Rprintf("%3e %3e\n", z(0), z(1));
    out(i) = rootisum - 0.5 * (arma::dot(z, z) - arma::dot(x.row(i), x.row(i)));
  }
  Rprintf("n = %i, p = %i, q = %i, x.n_rows = %i, x.n_cols = %i\n", n, p, q, x.n_rows, x.n_cols);
  arma::mat condmn = x.cols(0,p-1) * sigma1mn.t();
  arma::mat x2 = x.cols(p,x.n_cols-1);
  // x2 = x2 - condmn;

  Rprintf("x2.n_rows = %i, x2.n_cols = %i\n", x2.n_rows, x2.n_cols);
  Rprintf("x2(0,0) = %1.1e\n", x2(0,0));


  arma::mat wh = arma::zeros(n,q);

  //Rcpp::NumericVector upper(p), lower(p);
  arma::mat upper(n,q), lower(n,q);
  Rcpp::IntegerMatrix infin(n,q);

  Rprintf("here...");

  Rprintf("q = %i, trunc.length() = %i\n", q, trunc.length());

  for (uword j = 0; j < q; j++) {
    Rprintf("j = %i\n", j);
    Rcpp::NumericVector tmp = trunc[j];
    Rprintf("tmp.length() = %i, tmp(0) = %2e, tmp(1) = %2e\n", tmp.length(), tmp(0), tmp(1));

    Rcpp::NumericVector tmp2 = Rcpp::cumsum(tmp);
    tmp2.push_front(0);
    tmp2 = qnorm(tmp2);

    for (uword i = 0; i < n; i++) {
      Rprintf("(i,j)=(%i,%i)...length(tmp) = %i...x2(i,j) = %0d", i, j, tmp.length(), x2(i,j));
      lower(i,j) = (tmp2(x2(i,j)) + condmn(i,j))/sqrt(sigma1cov(j,j));
      upper(i,j) = (tmp2(x2(i,j)+1) + condmn(i,j))/sqrt(sigma1cov(j,j));
      infin(i,j) = 2 - 2*std::isinf(lower(i,j)) - std::isinf(upper(i,j));
    }
  }

  lower.print();
  upper.print();
  // infin.print();
  Rcpp::Rcout << "infin:\n" << infin << std::endl;

  Rprintf("1...");

  arma::mat correl1cov(q, q, arma::fill::eye);

  Rprintf("2...");

  for (uword i = 0; i < q-1; i++) for (uword j = i+1; j < q; j++) {
    correl1cov(i,j) = correl1cov(j,i) = sigma1cov(i,j)/sqrt(sigma1cov(i,i)*sigma1cov(j,j));
  }

  Rprintf("3...");

  // arma::mat cDec = arma::chol(sigma1cov);
  double tmp[q*(q-1)/2];

  // assign upper triangular elements to tmp
  int ct = 0;
  for (uword j=0; j < q; j++) {
    for (uword i=0; i < j; i++) {
      tmp[ct] = correl1cov(i,j);
      ct += 1;
    }
  }

  double zmn[q];
  for (uword j = 0; j < q; j++) {
    zmn[j] = 0.0;
  }
  int df = 0, rnd = 1, inform = 0;
  int maxpts = 25000;
  double abseps = 1E-3, releps = 0.0;

  for (uword i = 0; i < n; i++) {
    double lower2[q], upper2[q];
    int infin2[q];
    double error = 0.0, value = 0.0;

    for (uword j = 0; j < q; j++) {
      lower2[j] = lower(i,j);
      upper2[j] = upper(i,j);
      infin2[j] = infin(i,j);
    }

    // C_mvtdst2(&q);


    int q2 = q;
    // Rprintf("q = %i, df = %i\n", q2, df);
    // Rprintf("lower: %e %e\n", lower2[0], lower2[1]);
    // Rprintf("upper: %e %e\n", upper2[0], upper2[1]);
    // Rprintf("infin: %i %i\n", infin2[0], infin2[1]);
    //
    // Rprintf("sigmaUT: ");
    // for (int nn=0; nn < q*(q+1)/2; nn++) Rprintf("%e ", tmp[nn]);
    // Rprintf("\n");
    // Rprintf("delta: %e %e\n", zmn[0], zmn[1]);
    // Rprintf("maxpts = %i, abseps = %e, releps = %e\n", maxpts, abseps, releps);
    // Rprintf("error = %e, value = %e, inform = %i, rnd = %i\n", error, value, inform, rnd);

    C_mvtdst(&q2, // N
         &df,
         lower2,
         upper2,
         infin2,
         tmp,
         zmn,     // DELTA
         &maxpts,  // MAXPTS
         &abseps,  // ABSEPS
         &releps,  // RELEPS
         &error, // ERROR
         &value, // VALUE
         &inform,  // INFORM
         &rnd);  //RND

    Rprintf("value: %e\n", value);

    out(i) += log(value);
  }


  if (logd)
    return out;
  return exp(out);
}

// [[Rcpp::export]]
arma::vec dGDcop_sig(arma::mat const &x,
                 arma::cube const &sigma,
                 Rcpp::List trunc,
                 bool const logd = false) {

  // If no discrete component, then just call Gaussian only function
  uword q = trunc.length();
  if (q == 0) return(dGcop_sig(x, sigma, logd));

  // n = number of instances
  uword const n = x.n_rows;
  //d = x.n_cols;
  arma::vec out(n);
  arma::rowvec z;
  uword const p = sigma.n_rows - trunc.length();   // p = no of cts components
  // Rprintf("n = %i\n", n);

  // set up matrices for later use
  arma::mat wh = arma::zeros(n,q);

  //Rcpp::NumericVector upper(p), lower(p);
  arma::mat upper(n,q), lower(n,q);
  Rcpp::IntegerMatrix infin(n,q);
  Rcpp::NumericVector tmp2;
  Rcpp::List trunc2(trunc.length());

  for (uword j = 0; j < q; j++) {
    // Rprintf("j = %i\n", j);
    tmp2 = trunc[j];
    // Rprintf("tmp2(0) = %2e, tmp2(1) = %2e\n", tmp2(0), tmp2(1));
    Rcpp::NumericVector tmp3 = Rcpp::cumsum(tmp2);

    tmp3.push_front(0);
    tmp3 = qnorm(tmp3);
    // Rprintf("tmp2: %2e %2e\n", tmp3(0), tmp3(1));
    trunc2[j] = tmp3;
  }


  for (uword i = 0; i < n; i++) {
    // Rprintf("%i %i %i\n", p, sigma.n_rows, sigma.n_cols);
    // split the matrix into continuous/discrete pieces
    arma::mat sigma0 = sigma.slice(i).submat(0,0,p-1,p-1);
    arma::mat sigma1 = sigma.slice(i).submat(p,p,sigma.n_rows-1, sigma.n_cols-1);
    arma::mat sigma10 = sigma.slice(i).submat(p,0,sigma.n_rows-1, p-1);

    // Rprintf("%i %i %i\n", sigma0.n_elem, sigma1.n_elem, sigma10.n_elem);
    // get the conditional covariance given the continuous part
    arma::mat sigma1_0 = SchurC(sigma1, sigma0, sigma10);
    arma::mat sigma1cov = sigma1_0.cols(0,q-1);
    arma::mat sigma1mn = sigma1_0.cols(q,sigma1_0.n_cols-1);

    // Rprintf("%i %i %i\n", sigma1_0.n_elem, sigma1cov.n_elem, sigma1mn.n_elem);

    // double const constants = -(double)d/2.0 * log2pi;

    // mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    // check that continuous part has positive definite covariance
    arma::vec const eig = arma::eig_sym(sigma0);
    if (any(eig < 0)) {
      out(i) = NA_REAL;
      continue;
      // for (uword i = 0; i < n; i++) out(i) = nan;
    }
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma0)));
    double const rootisum = arma::sum(log(rooti.diag()));

    // Rprintf("%3e\n", rootisum);

    z = x.row(i);
    inplace_tri_mat_mult(z, rooti);
    // Rprintf("%3e %3e\n", z(0), z(1));
    out(i) = rootisum - 0.5 * (arma::dot(z, z) - arma::dot(x.row(i), x.row(i)));

    // Rprintf("n = %i, p = %i, q = %i, x.n_rows = %i, x.n_cols = %i\n", n, p, q, x.n_rows, x.n_cols);
    arma::mat condmn = x.cols(0,p-1) * sigma1mn.t();
    arma::mat x2 = x.cols(p,x.n_cols-1);
    // x2 = x2 - condmn;

    // Rprintf("x2.n_rows = %i, x2.n_cols = %i\n", x2.n_rows, x2.n_cols);
    // Rprintf("x2(0,0) = %1.1e\n", x2(0,0));

    // Rprintf("here...");

    // now go through getting upper and lower limits
    for (uword j = 0; j < q; j++) {
      // Rprintf("0...");
      tmp2 = trunc2[j];
      // Rprintf("%e %e\n", tmp2(0), tmp2(1));

      // Rprintf("(i,j)=(%i,%i), length(tmp2) = %i...x2(i,j) = %e, tmp2(x2(i,j)) = %e\n", i, j, tmp2.length(), x2(i,j), tmp2(x2(i,j)));
      lower(i,j) = (tmp2(x2(i,j)) + condmn(i,j))/sqrt(sigma1cov(j,j));
      upper(i,j) = (tmp2(x2(i,j)+1) + condmn(i,j))/sqrt(sigma1cov(j,j));
      infin(i,j) = 2 - 2*std::isinf(lower(i,j)) - std::isinf(upper(i,j));
    }

    // Rprintf("1...");

    arma::mat correl1cov(q, q, arma::fill::eye);

    // Rprintf("2...");

    // turn covariance into a correlation matrix
    for (uword j1 = 0; j1 < q-1; j1++) for (uword j2 = j1+1; j2 < q; j2++) {
      correl1cov(j1,j2) = correl1cov(j2,j1) = sigma1cov(j1,j2)/sqrt(sigma1cov(j1,j1)*sigma1cov(j2,j2));
    }

    // Rprintf("3...\n");

    // arma::mat cDec = arma::chol(sigma1cov);
    double tmp[q*(q-1)/2];

    // assign upper triangular elements of correlation to tmp
    int ct = 0;
    for (uword j=0; j < q; j++) {
      for (uword i=0; i < j; i++) {
        tmp[ct] = correl1cov(i,j);
        ct += 1;
      }
    }

    // set up elements for FORTRAN function
    double zmn[q];
    for (uword j = 0; j < q; j++) {
      zmn[j] = 0.0;
    }
    int df = 0, rnd = 1, inform = 0;
    int maxpts = 25000;
    double abseps = 1E-3, releps = 0.0;

    // set up limits and indicators of infinite limits
    double lower2[q], upper2[q];
    int infin2[q];
    double error = 0.0, value = 0.0;

    for (uword j = 0; j < q; j++) {
      lower2[j] = lower(i,j);
      upper2[j] = upper(i,j);
      infin2[j] = infin(i,j);
    }

    // C_mvtdst2(&q);

    int q2 = q;

    // Rprintf("q2 = %i\n", q2);
    // Rprintf("df = %i\n", df);
    // Rpr(lower2, q);
    // Rpr(upper2, q);
    // Rpri(infin2, q);
    //
    // Rprintf("tmp, zmn:\n");
    // Rpr(tmp, q*(q-1)/2);
    // Rpr(zmn, q);

    C_mvtdst(&q2, // N
             &df,
             lower2,
             upper2,
             infin2,
             tmp,
             zmn,     // DELTA
             &maxpts,  // MAXPTS
             &abseps,  // ABSEPS
             &releps,  // RELEPS
             &error, // ERROR
             &value, // VALUE
             &inform,  // INFORM
             &rnd);  //RND

    // Rprintf("value = %f\n", value);
    out(i) += log(value);
  }

  if (logd)
    return out;
  return exp(out);
}


// // C++ function to compute density at points of Gaussian copula
// // [[Rcpp::export]]
// arma::vec dGcop_sig(arma::mat const &x,
//                     arma::cube const &sigma,
//                     bool const logd = false) {
//   using arma::uword;
//   uword const n = x.n_rows;
//   arma::vec out(n);
//   // double const constants = -(double)d/2.0 * log2pi;
//
//   for (uword i = 0; i < n; i++) {
//     arma::vec eig = arma::eig_sym(sigma.slice(i));
//     if (any(eig < 0)) {
//       out(i) = NA_REAL;
//       continue;
//     }
//     arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma.slice(i))));
//     double const rootisum = arma::sum(log(rooti.diag()));
//     arma::rowvec z;
//
//     // Rprintf("%3e\n", rootisum);
//
//     z = x.row(i);
//     inplace_tri_mat_mult(z, rooti);
//     // Rprintf("%3e %3e\n", z(0), z(1));
//     out(i) = rootisum - 0.5 * (arma::dot(z, z) - arma::dot(x.row(i), x.row(i)));
//   }
//
//   if (logd)
//     return out;
//   return exp(out);
// }
//
//
// /*** R
// N <- 100; p <- 3
// x <- matrix(rnorm(N*p), N, p)
// M <- rWishart(N, df=10, Sigma=diag(p))
// M <- unlist(apply(M, 3, function(x) list(x/sqrt(diag(x)*rep(diag(x), each=nrow(x))))), recursive = FALSE)
// dmvnrm_arma(x, M, TRUE)
// */
