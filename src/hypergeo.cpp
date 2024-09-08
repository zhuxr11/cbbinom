#include "hypergeo.h"
#include "boost/math/special_functions/hypergeometric_pFq.hpp"
#include "boost/multiprecision/mpfr.hpp"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH)]]

// Convert R vector to std::vector<cpp_bin_float<prec>>
template <typename T1, typename T2, typename T3>
std::vector<T3> conv_vec_prec(const T1& x) {
  std::vector<T3> out(x.length());
  for (R_xlen_t i = 0; i < x.length(); i++) {
    T2 x_i = x(i);
    out[i] = x_i;
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
//' @param prec \code{NULL} or (unsigned) integer (1L) as precision level during computation,
//' a.k.a the number of precise digits defined in \code{boost::multiprecision::mpfr_float}.
//' For \code{NULL}, double precision (default) is used.
//' @param check_mode Logical (1L) indicating whether the mode of \code{x}
//' should be checked for obvious convergence failures.
//' @param log Logical (1L) indicating whether result is given as log(result).
//'
//' @return Numeric (1L) as the result of computation (at \code{double} precision).
//' Warnings are issued if failing to converge.
//'
//' @note Change log:
//' \itemize{
//'   \item{0.1.0 Xiurui Zhu - Initiate the function.}
//'   \item{0.1.1 Xiurui Zhu - Use \code{boost} for higher accuracy.}
//' }
//' @author Xiurui Zhu
//'
//' @export
//'
//' @examples
//' gen_hypergeo(U = c(1.1, 0.2, 0.3), L = c(10.1, 4 * pi), x = 1,
//'              tol = NULL, max_iter = 10000L, prec = NULL, check_mode = TRUE, log = FALSE)
// [[Rcpp::export]]
double gen_hypergeo(const NumericVector& U,
                    const NumericVector& L,
                    const double& x,
                    const Nullable<IntegerVector>& prec,
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
  typedef typename boost::math::policies::policy<
    boost::math::policies::max_series_iterations<10000>
  > MyPol;
  double out = R_NaN;
  if (prec.isNull()) {
    double p_abs_error;
    out = boost::math::hypergeometric_pFq(U, L, x, &p_abs_error, MyPol());
  } else {
    IntegerVector prec_ = as<IntegerVector>(prec);
    unsigned int prec_use = (unsigned) prec_(0);
    typedef typename boost::multiprecision::number<boost::multiprecision::backends::mpfr_float_backend<0>> prec_float;
    boost::math::scoped_precision<prec_float> scoped(prec_use);
    prec_float x_ = x;
    std::vector<prec_float> U_ = conv_vec_prec<NumericVector, double, prec_float>(U);
    std::vector<prec_float> L_ = conv_vec_prec<NumericVector, double, prec_float>(L);
    prec_float p_abs_error;
    prec_float out_ = boost::math::hypergeometric_pFq(U_, L_, x_, &p_abs_error, MyPol());
    out = static_cast<double>(out_);
  }
  if (log == true) {
    return std::log(out);
  } else {
    return out;
  }
}

/*** R
gen_hypergeo(U = c(1.1, 0.2, 0.3), L = c(10.1, 4 * pi), x = 1,
             tol = NULL, max_iter = 10000L, prec = NULL, check_mode = TRUE, log = FALSE)
*/
