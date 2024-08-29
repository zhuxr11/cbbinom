#include "hypergeo.h"
#include "../inst/include/cbbinom/uniroot.h"
#include <Rcpp.h>
using namespace Rcpp;

// Macros

#ifndef CBBINOM_MACROS
#define CBBINOM_MACROS
#define GETV(x, i)      x(i % x.length())    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))
#endif // CBBINOM_MACROS


// Non-vectorized core functions

double pcbbinom_(
    const double& q,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& lower_tail,
    const bool& log_p,
    const double& tol,
    const int& max_iter
) {
  if (q < 0.0) {
    return 0.0;
  } else if (q > size + 1.0) {
    return 1.0;
  }
  // Coefficients
  double coef_log = R::lgammafn(size + 1.0) + R::lbeta(size + 1.0 - q + beta, alpha) -
    R::lgammafn(q) - R::lgammafn(size + 2.0 - q) - R::lbeta(alpha, beta);
  // Upper vector
  NumericVector U = {1.0 - q, size + 1.0 - q, size + 1.0 - q + beta};
  // Lower vector
  NumericVector L = {size + 2.0 - q, size + 1.0 - q + alpha + beta};
  // Final value
  double out = std::exp(coef_log + gen_hypergeo(U, L, 1.0, tol, max_iter, true, true));
  if (lower_tail == false) {
    out = 1 - out;
  }
  if (log_p == true) {
    out = std::log(out);
  }
  return out;
}

class PcbbinomEqn: public cbbinom::UnirootEqn
{
private:
  double size;
  double alpha;
  double beta;
  double tol;
  int max_iter;
  double p;
public:
  PcbbinomEqn(const double size_, const double alpha_, const double beta_,
              const double tol_, const int max_iter_, const double p_):
    size(size_), alpha(alpha_), beta(beta_),
    tol(tol_), max_iter(max_iter_), p(p_) {};
  double operator () (const double& x) const override {
    return pcbbinom_(x, this->size, this->alpha, this->beta,
                     true, false, this->tol, this->max_iter) - this->p;
  }
};

double qcbbinom_(
    double p,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& lower_tail,
    const bool& log_p,
    const double& p_tol,
    const int& p_max_iter,
    double root_tol,
    int root_max_iter
) {
  if (log_p == true) {
    p = std::exp(p);
  }
  if (VALID_PROB(p) == false) {
    warning("Wrong [p] as probability: %f, returning NA", p);
    return R_NaN;
  }
  if (lower_tail == false) {
    p = 1.0 - p;
  }
  PcbbinomEqn eqn_obj(size, alpha, beta, p_tol, p_max_iter, p);
  return cbbinom::cpp_uniroot(0.0, size + 1.0, -p, 1.0 - p,
                              &eqn_obj, &root_tol, &root_max_iter);
}

// Vectorized functions

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericVector cpp_pcbbinom(
    const NumericVector& q,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail,
    const bool& log_p,
    const NumericVector& tol,
    const IntegerVector& max_iter
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
      GETV(tol, idx),
      GETV(max_iter, idx)
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

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericVector cpp_qcbbinom(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail,
    const bool& log_p,
    const NumericVector& p_tol,
    const IntegerVector& p_max_iter,
    const NumericVector& root_tol,
    const IntegerVector& root_max_iter
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

  NumericVector root_tol_ = clone(root_tol);
  IntegerVector root_max_iter_ = clone(root_max_iter);
  for (R_xlen_t idx = 0; idx < n_max; idx++) {
    q(idx) = qcbbinom_(
      GETV(p_, idx),
      GETV(size, idx),
      GETV(alpha, idx),
      GETV(beta, idx),
      true,   // lower_tail
      false,  // log_p
      GETV(p_tol, idx),
      GETV(p_max_iter, idx),
      GETV(root_tol_, idx),
      GETV(root_max_iter_, idx)
    );
  }
  return q;
}

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericMatrix dcbblp(
    const NumericVector& x,
    const NumericVector& m,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& tol,
    const IntegerVector& max_iter
) {
  if (std::min({x.length(), m.length(),
               a.length(), b.length()}) < 1) {
    return NumericMatrix(0, 3);
  }

  int n_max = std::max({
    x.length(),
    m.length(),
    a.length(),
    b.length()
  });
  NumericMatrix f(n_max, 3);

  double h = 1e-6;
  for (int i = 0; i < n_max; i++){
    double xi = GETV(x, i);
    double mi = GETV(m, i);
    double ai = GETV(a, i);
    double bi = GETV(b, i);
    double ti = GETV(tol, i);
    double ri = GETV(max_iter, i);
    if (xi < 0) {
      f(i, 0) = R_NegInf;
      f(i, 1) = R_NegInf;
      f(i, 2) = h;
    } else if (xi <= EPSILON){
      xi = EPSILON * 1.0000001; // --> density at zero = density at eps
      f(i, 0) = pcbbinom_(xi + h, mi, ai, bi, true, true, ti, ri);
      f(i, 1) = pcbbinom_(xi, mi, ai, bi, true, true, ti, ri);
      f(i, 2) = h;
    } else if (xi > mi + 1) {
      f(i, 0) = 0;
      f(i, 1) = 0;
      f(i, 2) = h;
    } else if (xi >= h && xi <= mi + 1 - h){
      f(i, 0) = pcbbinom_(xi + h, mi, ai, bi, true, true, ti, ri);
      f(i, 1) = pcbbinom_(xi - h, mi, ai, bi, true, true, ti, ri);
      f(i, 2) = 2 * h;
    } else if (xi <= h) {
      f(i, 0) = pcbbinom_(xi + h, mi, ai, bi, true, true, ti, ri);
      f(i, 1) = pcbbinom_(xi, mi, ai, bi, true, true, ti, ri);
      f(i, 2) = h;
    } else {
      f(i, 0) = pcbbinom_(xi, mi, ai, bi, true, true, ti, ri);
      f(i, 1) = pcbbinom_(xi - h, mi, ai, bi, true, true, ti, ri);
      f(i, 2) = h;
    }
  }
  return(f);
}

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericVector cpp_rcbbinom(
  const int& n,
  const NumericVector& size,
  const NumericVector& alpha,
  const NumericVector& beta,
  const NumericVector& p_tol,
  const IntegerVector& p_max_iter,
  const NumericVector& root_tol,
  const IntegerVector& root_max_iter
) {
  return cpp_qcbbinom(runif(n), size, alpha, beta, true, false,
                      p_tol, p_max_iter, root_tol, root_max_iter);
}

/*** R
plot(seq(0, 11, 0.01),
     cpp_pcbbinom(q = seq(0, 11, 0.01), size = 10, alpha = 2, beta = 4,
                  lower_tail = TRUE, log_p = FALSE,
                  tol = 1e-6, max_iter = 10000L))
# Density function
cpp_dcbbinom(x = 5, size = 10, alpha = 2, beta = 4,
             log = FALSE, tol = 1e-6, max_iter = 10000L)
# Distribution function
(test_val <- cpp_pcbbinom(q = 5, size = 10, alpha = 2, beta = 4,
                          lower_tail = TRUE, log_p = FALSE,
                          tol = 1e-6, max_iter = 10000L))
# Quantile function
cpp_qcbbinom(p = test_val, size = 10, alpha = 2, beta = 4,
             lower_tail = TRUE, log_p = FALSE,
             p_tol = 1e-6, p_max_iter = 10000L,
             root_tol = 1e-6, root_max_iter = 10000L)
# Random generation
set.seed(1111L)
cpp_rcbbinom(n = 10L, size = 10, alpha = 2, beta = 4,
             p_tol = 1e-6, p_max_iter = 10000L,
             root_tol = 1e-6, root_max_iter = 10000L)
*/
