#include <hypergeo2.h>
#include <Rcpp.h>
#include "cbbinom_impl.h"
using namespace Rcpp;

// Macros

#ifndef CBBINOM_MACROS
#define CBBINOM_MACROS
#define GETV(x, i)      x(i % x.length())    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#endif // CBBINOM_MACROS

template <typename T1, typename T2>
Nullable<T2> nullable_getv(const Nullable<T1>& x, const int& idx) {
  if (x.isNull()) {
    return R_NilValue;
  }
  T1 x_vec = as<T1>(x);
  if (GETV(x_vec, idx) == R_NilValue) {
    return R_NilValue;
  }
  T2 out(1, GETV(x_vec, idx));
  return out;
}

// [[Rcpp::interfaces(r, cpp)]]

// Vectorized functions

// [[Rcpp::export]]
NumericVector cpp_pcbbinom(
    const NumericVector& q,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail,
    const bool& log_p,
    const Nullable<List>& prec
) {
  if (std::min({q.length(), size.length(),
               alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int n_max = std::max({
    q.length(),
    size.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(n_max);

  for (R_xlen_t idx = 0; idx < n_max; idx++) {
    p(idx) = pcbbinom_(
      GETV(q, idx),
      GETV(size, idx),
      GETV(alpha, idx),
      GETV(beta, idx),
      true,   // lower_tail
      false,  // log_p
      nullable_getv<List, IntegerVector>(prec, idx)
    );
  }

  if (lower_tail == false) {
    p = 1.0 - p;
  }
  if (log_p == true) {
    p = log(p);
  }
  return p;
}

// [[Rcpp::export]]
NumericVector cpp_qcbbinom(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail,
    const bool& log_p,
    const Nullable<List>& prec,
    const NumericVector& tol,
    const IntegerVector& max_iter
) {
  if (std::min({p.length(), size.length(),
               alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int n_max = std::max({
    p.length(),
    size.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector q(n_max);

  NumericVector p_ = clone(p);
  if (log_p == true) {
    p_ = exp(p_);
  }
  if (lower_tail == false) {
    p_ = 1 - p_;
  }

  NumericVector tol_ = clone(tol);
  IntegerVector max_iter_ = clone(max_iter);
  for (R_xlen_t idx = 0; idx < n_max; idx++) {
    q(idx) = qcbbinom_(
      GETV(p_, idx),
      GETV(size, idx),
      GETV(alpha, idx),
      GETV(beta, idx),
      true,   // lower_tail
      false,  // log_p
      nullable_getv<List, IntegerVector>(prec, idx),
      GETV(tol_, idx),
      GETV(max_iter_, idx)
    );
  }
  return q;
}

// [[Rcpp::export]]
NumericVector cpp_dcbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log,
    const Nullable<List>& prec
) {
  if (std::min({x.length(), size.length(),
               alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int n_max = std::max({
    x.length(),
    size.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector out(n_max);

  for (int i = 0; i < n_max; i++){
    out(i) = dcbbinom_(
      GETV(x, i),
      GETV(size, i),
      GETV(alpha, i),
      GETV(beta, i),
      false,
      nullable_getv<List, IntegerVector>(prec, i)
    );
  }
  if (log == true) {
    out = Rcpp::log(out);
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector cpp_rcbbinom(
  const int& n,
  const NumericVector& size,
  const NumericVector& alpha,
  const NumericVector& beta,
  const Nullable<List>& prec,
  const NumericVector& tol,
  const IntegerVector& max_iter
) {
  return cpp_qcbbinom(runif(n), size, alpha, beta, true, false,
                      prec, tol, max_iter);
}

/*** R
plot(seq(0, 11, 0.01),
     cpp_pcbbinom(q = seq(0, 11, 0.01), size = 10, alpha = 2, beta = 4,
                  lower_tail = TRUE, log_p = FALSE, prec = NULL))
# Density function
cpp_dcbbinom(x = 5, size = 10, alpha = 2, beta = 4, log = FALSE, prec = 20L)
# Distribution function
(test_val <- cpp_pcbbinom(q = 5, size = 10, alpha = 2, beta = 4,
                          lower_tail = TRUE, log_p = FALSE, prec = 20L))
# Quantile function
cpp_qcbbinom(p = test_val, size = 10, alpha = 2, beta = 4,
             lower_tail = TRUE, log_p = FALSE, prec = 20L,
             tol = 1e-6, max_iter = 10000L)
# Random generation
set.seed(1111L)
cpp_rcbbinom(n = 10L, size = 10, alpha = 2, beta = 4, prec = 20L,
             tol = 1e-6, max_iter = 10000L)
*/
