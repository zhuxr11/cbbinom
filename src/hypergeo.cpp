#include "hypergeo.h"
#include <Rcpp.h>
using namespace Rcpp;

double prod(const NumericVector& x) {
  double out = 1.0;
  for (NumericVector::const_iterator i = x.cbegin(); i != x.cend(); i++) {
    out *= *i;
  }
  return out;
}

// [[Rcpp::interfaces(r, cpp)]]
//' Generalized hypergeometric function
//'
//' \code{gen_hypergeo} computes generalized hypergeometric function.
//'
//' @param U,L Numeric vectors for upper and lower values.
//' @param x Numeric (1L) as common ratio.
//' @param tol Numeric (1L) as convergence tolerance.
//' @param max_iter Integer (1L) as iteration limit.
//' @param check_mode Logical (1L) indicating whether the mode of \code{x}
//' should be checked for obvious convergence failures.
//' @param log Logical (1L) indicating whether result is given as log(result).
//'
//' @return Result of computation. Warnings are issued if failing to converge.
//'
//' @note Change log:
//' \itemize{
//'   \item{0.1.0 Xiurui Zhu - Initiate the function.}
//' }
//' @author Xiurui Zhu
//'
//' @export
//'
//' @examples
//' gen_hypergeo(U = c(1.1, 0.2, 0.3), L = c(10.1, 4 * pi), x = 1,
//'               max_iter = 10000L, tol = 1e-6, check_mode = TRUE, log = FALSE)
// [[Rcpp::export]]
double gen_hypergeo(const NumericVector& U,
                    const NumericVector& L,
                    const double& x,
                    const double& tol,
                    const R_xlen_t& max_iter,
                    const bool& check_mode,
                    const bool& log) {
  // Check inputs
  if (check_mode == true) {
    if (U.length() > L.length() + 1) {
      if (std::fabs(x) > 0) {
        warning("length(U) > length(L) + 1L: converge failure if abs(x) > 0; returning NaN");
        return R_NaN;
      }
    } else if (U.length() > L.length()) {
      if (std::fabs(x) > 1) {
        warning("length(U) > length(L): converge failure if abs(x) > 1; returning NaN");
        return R_NaN;
      }
    }
  }
  // Main code
  double out = x;
  double n = 1.0;
  double prev = 0.0;
  double term = x;
  bool conv = true;
  NumericVector U_ = clone(U);
  NumericVector L_ = clone(L);
  while (std::fabs(out - prev) > tol) {
    prev = out;
    term *= prod(U_) / prod(L_) * x / n;
    out += term;
    n++;
    if (n >= max_iter) {
      conv = false;
      break;
    }
    U_ = U_ + 1.0;
    L_ = L_ + 1.0;
  }
  if (conv == false) {
    warning("Generalized hypergeometric function fails to converge");
  }
  if (log == true) {
    return std::log(out);
  } else {
    return out;
  }
}

/*** R
gen_hypergeo(U = c(1.1, 0.2, 0.3), L = c(10.1, 4 * pi), x = 1,
             max_iter = 10000L, tol = 1e-6, check_mode = TRUE, log = FALSE)
*/
