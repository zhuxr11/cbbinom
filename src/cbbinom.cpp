#include <hypergeo2.h>
#include <Rcpp.h>
#include "../inst/include/cbbinom/uniroot.h"
#include "boost/math/differentiation/finite_difference.hpp"
using namespace Rcpp;

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(hypergeo2)]]

// Macros

#ifndef CBBINOM_MACROS
#define CBBINOM_MACROS
#define GETV(x, i)      x(i % x.length())    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))
#define VALID_PARAM(x, size, alpha, beta) ((x >= 0.0) && (x <= size + 1.0) && (size >= 0.0) && (alpha > 0.0) && (beta > 0.0))
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

double gen_hypergeo(const double& q,
                    const double& size,
                    const double& alpha,
                    const double& beta,
                    const Nullable<IntegerVector>& prec) {
  NumericVector U = {1.0 - q, size + 1.0 - q, size + 1.0 - q + beta};
  NumericVector L = {size + 2.0 - q, size + 1.0 - q + alpha + beta};
  NumericVector z = {1.0};
  NumericVector out;
  if (prec.isNull()) {
    out = hypergeo2::genhypergeo_vec(
      List::create(U),
      List::create(L),
      z,
      R_NilValue,
      true,  // check_mode
      false,  // log
      "gmp"  // backend
    );
  } else {
    int prec_use = as<IntegerVector>(prec)(0);
    out = hypergeo2::genhypergeo_vec(
      List::create(U),
      List::create(L),
      z,
      List::create(prec_use),
      true,  // check_mode
      false,  // log
      "gmp"  // backend
    );
  }
  return out(0);
}

// Non-vectorized core functions

double pcbbinom_(
    const double& q,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& lower_tail,
    const bool& log_p,
    const Nullable<IntegerVector>& prec
) {
  if (VALID_PARAM(q, size, alpha, beta) == false) {
    warning("Invalid parameter set: q = %f, size = %f, alpha = %f, beta = %f; returing NaN",
            q, size, alpha, beta);
    return R_NaN;
  }
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
  double out = std::exp(coef_log) * gen_hypergeo(q, size, alpha, beta, prec);
  if (lower_tail == false) {
    out = 1 - out;
  }
  if (log_p == true) {
    out = std::log(out);
  }
  return out;
}

// Adapt boost::math::differentiation::finite_difference_derivative() to accept one-sided derivative
namespace boost {
namespace math {
namespace differentiation {

template<class F, class Real>
Real finite_difference_derivative_bound(const F f, Real x, Real low, Real high, Real* error = nullptr)
{
  using std::sqrt;
  using std::pow;
  using std::abs;
  using std::numeric_limits;

  const Real eps = (numeric_limits<Real>::epsilon)();
  // Error bound ~eps^6/7
  // Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
  Real h = pow(eps / 168, static_cast<Real>(1) / static_cast<Real>(7));
  h = detail::make_xph_representable(x, h);
  Real yh = static_cast<Real>(0);
  Real y2h = static_cast<Real>(0);
  Real y3h = static_cast<Real>(0);
  Real ymh = static_cast<Real>(0);
  Real ym2h = static_cast<Real>(0);
  Real ym3h = static_cast<Real>(0);

  Real y = f(x);
  // Rprintf("x = %f\n", x);
  // Rprintf("h = %f\n", h);
  // Rprintf("y = %f\n", y);
  if (x + 3 * h <= high) {
    yh = f(x + h);
    // Rprintf("yh = %f\n", yh);
    y2h = f(x + 2 * h);
    // Rprintf("y2h = %f\n", y2h);
    y3h = f(x + 3 * h);
    // Rprintf("y3h = %f\n", y3h);
  }
  if (x - 3 * h >= low) {
    ymh = f(x - h);
    // Rprintf("ymh = %f\n", ymh);
    ym2h = f(x - 2 * h);
    // Rprintf("ym2h = %f\n", ym2h);
    ym3h = f(x - 3 * h);
    // Rprintf("ym3h = %f\n", ym3h);
  }

  if ((x - 3 * h >= low) && (x + 3 * h <= high)) {
    Real y1 = yh - ymh;
    Real y2 = ym2h - y2h;
    Real y3 = y3h - ym3h;
    if (error)
    {
      // Mathematica code to generate fd scheme for 7th derivative:
      // Sum[(-1)^i*Binomial[7, i]*(f[x+(3-i)*h] + f[x+(4-i)*h])/2, {i, 0, 7}]
      // Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
      // Series[(f[x+4*h]-f[x-4*h] + 6*(f[x-3*h] - f[x+3*h]) + 14*(f[x-h] - f[x+h] + f[x+2*h] - f[x-2*h]))/2, {h, 0, 15}]
      Real y7 = (f(x + 4 * h) - f(x - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
      *error = abs(y7) / (140 * h) + 5 * (abs(yh) + abs(ymh))*eps / h;
    }
    return (y3 + 9 * y2 + 45 * y1) / (60 * h);
  } else if (x - 3 * h >= low) {
    // Backward differentiation
    if (error) {
      warning("[error] not implemented for backward differentiation");
    }
    return (11 * y - 18 * ymh + 9 * ym2h - 2 * ym3h) / (6 * h);
  } else if (x + 3 * h <= high) {
    // Forward differentiation
    if (error) {
      warning("[error] not implemented for forward differentiation");
    }
    return (-11 * y + 18 * yh - 9 * y2h + 2 * y3h) / (6 * h);
  } else {
    stop("Insufficient range: high - low < 6 * %f", h);
  }
}

}
}
}

double dcbbinom_(
    const double& x,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& log,
    const Nullable<IntegerVector>& prec
) {
  if (VALID_PARAM(x, size, alpha, beta) == false) {
    warning("Invalid parameter set: x = %f, size = %f, alpha = %f, beta = %f; returing NaN",
            x, size, alpha, beta);
    return R_NaN;
  }
  if ((x < 0.0) || (x > size + 1.0)) {
    if (log == true) {
      return R_NegInf;
    } else {
      return 0.0;
    }
  }
  // Compute derivative on pcbbinom_
  auto f = [&size, &alpha, &beta, &prec](double x) {
    return pcbbinom_(x, size, alpha, beta, true, false, prec);
  };
  double out = boost::math::differentiation::finite_difference_derivative_bound(f, x, 0.0, size + 1.0);
  if (out < 0.0) {
    String base_warning = "d[pcbbinom(q = %f, size = %f, alpha = %f, beta = %f)]/dq = %f < 0, which is set to 0, since probability density cannot be negative";
    if (prec.isNotNull()) {
      IntegerVector prec_ = as<IntegerVector>(prec);
      int prec_use = prec_(0);
      base_warning += "; you may use a higher [prec] level than %f";
      warning(base_warning.get_cstring(), x, size, alpha, beta, out, prec_use);
    } else {
      warning(base_warning.get_cstring(), x, size, alpha, beta, out);
    }
    out = 0.0;
  }
  return out;
  if (log == true) {
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
  Nullable<NumericVector> tol;
  int max_iter;
  Nullable<IntegerVector> prec;
  double p;
public:
  PcbbinomEqn(const double size_, const double alpha_, const double beta_,
              const Nullable<IntegerVector> prec_, const double p_):
    size(size_), alpha(alpha_), beta(beta_), prec(prec_), p(p_) {};
  double operator () (const double& x) const override {
    return pcbbinom_(x, this->size, this->alpha, this->beta,
                     true, false, this->prec) - this->p;
  }
};

double qcbbinom_(
    double p,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& lower_tail,
    const bool& log_p,
    const Nullable<IntegerVector>& prec,
    double tol,
    int max_iter
) {
  if (log_p == true) {
    p = std::exp(p);
  }
  if (VALID_PROB(p) == false) {
    warning("Wrong [p] as probability: %f; returning NaN", p);
    return R_NaN;
  }
  if (VALID_PARAM(0.0, size, alpha, beta) == false) {
    warning("Invalid parameter set: size = %f, alpha = %f, beta = %f; returing NaN",
            size, alpha, beta);
    return R_NaN;
  }
  if (lower_tail == false) {
    p = 1.0 - p;
  }
  PcbbinomEqn eqn_obj(size, alpha, beta, prec, p);
  return cbbinom::cpp_uniroot(0.0, size + 1.0, -p, 1.0 - p,
                              &eqn_obj, &tol, &max_iter);
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

// [[Rcpp::interfaces(r, cpp)]]
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

// [[Rcpp::interfaces(r, cpp)]]
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

// [[Rcpp::interfaces(r, cpp)]]
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
