#include "hypergeo.h"
#include "../inst/include/cbbinom/uniroot.h"
#include "boost/math/differentiation/finite_difference.hpp"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH)]]

// Macros

#ifndef CBBINOM_MACROS
#define CBBINOM_MACROS
#define GETV(x, i)      x(i % x.length())    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))
#endif // CBBINOM_MACROS

template <typename T1, typename T2>
Nullable<T2> nullable_getv(const Nullable<T1>& x, const int& idx) {
  if (x.isNull()) {
    return R_NilValue;
  }
  T2 x_vec = as<NumericVector>(x);
  T2 out(1, GETV(x_vec, idx));
  return out;
}

// Non-vectorized core functions

double pcbbinom_(
    const double& q,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& lower_tail,
    const bool& log_p,
    const Nullable<NumericVector>& tol,
    const int& max_iter,
    const double& prec
) {
  if (q < 0.0) {
    if (log_p == true) {
      return R_NegInf;
    } else {
      return 0.0;
    }
  } else if (q > size + 1.0) {
    if (log_p == true) {
      return 0.0;
    } else {
      return 1.0;
    }
  }
  // Coefficients
  double coef_log = R::lgammafn(size + 1.0) + R::lbeta(size + 1.0 - q + beta, alpha) -
    R::lgammafn(q) - R::lgammafn(size + 2.0 - q) - R::lbeta(alpha, beta);
  // Upper vector
  NumericVector U = {1.0 - q, size + 1.0 - q, size + 1.0 - q + beta};
  // Lower vector
  NumericVector L = {size + 2.0 - q, size + 1.0 - q + alpha + beta};
  // Final value
  double out = std::exp(coef_log + gen_hypergeo(U, L, 1.0, tol, max_iter, prec, true, true));
  if (lower_tail == false) {
    out = 1 - out;
  }
  if (log_p == true) {
    out = std::log(out);
  }
  return out;
}

double dcbbinom_(
    const double& x,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& log,
    const Nullable<NumericVector>& tol,
    const int& max_iter,
    const double& prec
) {
  if ((x < 0.0) || (x > size + 1.0)) {
    if (log == true) {
      return R_NegInf;
    } else {
      return 0.0;
    }
  }
  // Compute derivative on log(pcbbinom_)
  auto f = [&size, &alpha, &beta, &tol, &max_iter, &prec](double x) {
    return pcbbinom_(x, size, alpha, beta, true, true, tol, max_iter, prec);
  };
  double out = boost::math::differentiation::finite_difference_derivative(f, x);
  if (out < 0.0) {
    warning("d[pcbbinom(x = %f, size = %f, alpha = %f, beta = %f, log = TRUE)]/dx \
              is negative: %f, which is set to 0; please consider using \
              higher [prec] level than %f\n",
            x, size, alpha, beta, out, prec);
    out = R_NegInf;
  } else {
    // Compute log result: log(d(log(pcbbinom_))/dx) + log(pcbbinom_)
    out = std::log(out) + f(x);
  }
  if (log == false) {
    out = std::exp(out);
  }
  return out;
}

class PcbbinomEqn: public cbbinom::UnirootEqn
{
private:
  double size;
  double alpha;
  double beta;
  Nullable<NumericVector> tol;
  int max_iter;
  double prec;
  double p;
public:
  PcbbinomEqn(const double size_, const double alpha_, const double beta_,
              const Nullable<NumericVector> tol_, const int max_iter_,
              const double prec_, const double p_):
    size(size_), alpha(alpha_), beta(beta_),
    tol(tol_), max_iter(max_iter_), prec(prec_), p(p_) {};
  double operator () (const double& x) const override {
    return pcbbinom_(x, this->size, this->alpha, this->beta,
                     true, false, this->tol, this->max_iter,
                     this->prec) - this->p;
  }
};

double qcbbinom_(
    double p,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& lower_tail,
    const bool& log_p,
    const Nullable<NumericVector>& p_tol,
    const int& p_max_iter,
    const double& p_prec,
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
  PcbbinomEqn eqn_obj(size, alpha, beta, p_tol, p_max_iter, p_prec, p);
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
    const Nullable<NumericVector>& tol,
    const IntegerVector& max_iter,
    const NumericVector& prec
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
      nullable_getv<NumericVector, NumericVector>(tol, idx),
      GETV(max_iter, idx),
      GETV(prec, idx)
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
    const Nullable<NumericVector>& p_tol,
    const IntegerVector& p_max_iter,
    const NumericVector& p_prec,
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
      nullable_getv<NumericVector, NumericVector>(p_tol, idx),
      GETV(p_max_iter, idx),
      GETV(p_prec, idx),
      GETV(root_tol_, idx),
      GETV(root_max_iter_, idx)
    );
  }
  return q;
}

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericVector cpp_dcbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log,
    const Nullable<NumericVector>& tol,
    const IntegerVector& max_iter,
    const NumericVector& prec
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
      true,
      nullable_getv<NumericVector, NumericVector>(tol, i),
      GETV(max_iter, i),
      GETV(prec, i)
    );
  }
  if (log == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericVector cpp_rcbbinom(
  const int& n,
  const NumericVector& size,
  const NumericVector& alpha,
  const NumericVector& beta,
  const Nullable<NumericVector>& p_tol,
  const IntegerVector& p_max_iter,
  const NumericVector& p_prec,
  const NumericVector& root_tol,
  const IntegerVector& root_max_iter
) {
  return cpp_qcbbinom(runif(n), size, alpha, beta, true, false,
                      p_tol, p_max_iter, p_prec, root_tol, root_max_iter);
}

/*** R
plot(seq(0, 11, 0.01),
     cpp_pcbbinom(q = seq(0, 11, 0.01), size = 10, alpha = 2, beta = 4,
                  lower_tail = TRUE, log_p = FALSE,
                  tol = NULL, max_iter = 10000L, prec = 0))
# Density function
cpp_dcbbinom(x = 5, size = 10, alpha = 2, beta = 4,
             log = FALSE, tol = NULL, max_iter = 10000L, prec = 20)
# Distribution function
(test_val <- cpp_pcbbinom(q = 5, size = 10, alpha = 2, beta = 4,
                          lower_tail = TRUE, log_p = FALSE,
                          tol = NULL, max_iter = 10000L, prec = 20))
# Quantile function
cpp_qcbbinom(p = test_val, size = 10, alpha = 2, beta = 4,
             lower_tail = TRUE, log_p = FALSE,
             p_tol = NULL, p_max_iter = 10000L, p_prec = 20,
             root_tol = 1e-6, root_max_iter = 10000L)
# Random generation
set.seed(1111L)
cpp_rcbbinom(n = 10L, size = 10, alpha = 2, beta = 4,
             p_tol = NULL, p_max_iter = 10000L, p_prec = 20,
             root_tol = 1e-6, root_max_iter = 10000L)
*/
