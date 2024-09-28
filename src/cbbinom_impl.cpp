#include <hypergeo2.h>
#include <Rcpp.h>
#include "cbbinom_impl.h"
#include "../inst/include/cbbinom/uniroot.h"
#include "boost/math/differentiation/finite_difference.hpp"
using namespace Rcpp;

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(hypergeo2)]]

#ifndef CBBINOM_IMPL_MACROS
#define CBBINOM_IMPL_MACROS
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))
#define VALID_PARAM(x, size, alpha, beta) ((x >= 0.0) && (x <= size + 1.0) && (size >= 0.0) && (alpha > 0.0) && (beta > 0.0))
#endif // CBBINOM_IMPL_MACROS

double cbbinom_genhypergeo(const double& q,
                           const double& size,
                           const double& alpha,
                           const double& beta,
                           const Nullable<IntegerVector>& prec) {
  NumericVector U = {1.0 - q, size + 1.0 - q, size + 1.0 - q + beta};
  NumericVector L = {size + 2.0 - q, size + 1.0 - q + alpha + beta};
  double z = 1.0;
  double out = hypergeo2::genhypergeo_cpp(U, L, z, prec, true, false, "gmp");
  return out;
}

// [[Rcpp::interfaces(cpp)]]

// Non-vectorized core functions

// [[Rcpp::export]]
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
    warning("Invalid parameter set: q = %g, size = %g, alpha = %g, beta = %g; returing NaN",
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
  double out = std::exp(coef_log) * cbbinom_genhypergeo(q, size, alpha, beta, prec);
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

  // Rprintf("x = %g\n", x);
  const Real eps = (numeric_limits<Real>::epsilon)();
  // Error bound ~eps^6/7
  // Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
  Real h = pow(eps / 168, static_cast<Real>(1) / static_cast<Real>(7));
  h = detail::make_xph_representable(x, h);
  // Rprintf("h = %g\n", h);
  Real yh = static_cast<Real>(0);
  Real y2h = static_cast<Real>(0);
  Real y3h = static_cast<Real>(0);
  Real ymh = static_cast<Real>(0);
  Real ym2h = static_cast<Real>(0);
  Real ym3h = static_cast<Real>(0);

  Real y = f(x);
  // Rprintf("y = %g\n", y);
  if (x + 3 * h <= high) {
    yh = f(x + h);
    // Rprintf("yh = %g\n", yh);
    y2h = f(x + 2 * h);
    // Rprintf("y2h = %g\n", y2h);
    y3h = f(x + 3 * h);
    // Rprintf("y3h = %g\n", y3h);
  }
  if (x - 3 * h >= low) {
    ymh = f(x - h);
    // Rprintf("ymh = %g\n", ymh);
    ym2h = f(x - 2 * h);
    // Rprintf("ym2h = %g\n", ym2h);
    ym3h = f(x - 3 * h);
    // Rprintf("ym3h = %g\n", ym3h);
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
    stop("Insufficient range: high - low < 6 * %g", h);
  }
}

}
}
}

// [[Rcpp::export]]
double dcbbinom_(
    const double& x,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& log,
    const Nullable<IntegerVector>& prec
) {
  if (VALID_PARAM(x, size, alpha, beta) == false) {
    warning("Invalid parameter set: x = %g, size = %g, alpha = %g, beta = %g; returing NaN",
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
    String base_warning = "d[pcbbinom(q = %g, size = %g, alpha = %g, beta = %g)]/dq = %g < 0, which is set to 0, since probability density cannot be negative";
    if (prec.isNotNull()) {
      IntegerVector prec_ = as<IntegerVector>(prec);
      int prec_use = prec_(0);
      base_warning += "; you may use a higher [prec] level than %i";
      warning(base_warning.get_cstring(), x, size, alpha, beta, out, prec_use);
    } else {
      warning(base_warning.get_cstring(), x, size, alpha, beta, out);
    }
    out = 0.0;
  }
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

// [[Rcpp::export]]
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
    warning("Wrong [p] as probability: %g; returning NaN", p);
    return R_NaN;
  }
  if (VALID_PARAM(0.0, size, alpha, beta) == false) {
    warning("Invalid parameter set: size = %g, alpha = %g, beta = %g; returing NaN",
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
