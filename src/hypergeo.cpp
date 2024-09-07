#include "hypergeo.h"
#include "boost/math/special_functions/hypergeometric_pFq.hpp"
#include "boost/multiprecision/mpfr.hpp"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH)]]

// Adapt some generalized hypergeometric functions from BH, with adjustable tolerance
namespace boost {
namespace math {
namespace detail {

template <class Seq, class Real, class Policy, class Terminal>
std::pair<Real, Real> hypergeometric_pFq_checked_series_impl(
    const Seq& aj,
    const Seq& bj,
    const Real& z,
    const Policy& pol,
    const Real& tol,
    const Terminal& termination,
    long long& log_scale
) {
  BOOST_MATH_STD_USING
  Real result = 1;
  Real abs_result = 1;
  Real term = 1;
  Real term0 = 0;
  // Real tol = boost::math::policies::get_epsilon<Real, Policy>();
  std::uintmax_t k = 0;
  Real upper_limit(sqrt(boost::math::tools::max_value<Real>())), diff;
  Real lower_limit(1 / upper_limit);
  long long log_scaling_factor = lltrunc(boost::math::tools::log_max_value<Real>()) - 2;
  Real scaling_factor = exp(Real(log_scaling_factor));
  Real term_m1;
  long long local_scaling = 0;
  // bool have_no_correct_bits = false;

  if ((aj.size() == 1) && (bj.size() == 0))
  {
    if (fabs(z) > 1)
    {
      if ((z > 0) && (floor(*aj.begin()) != *aj.begin()))
      {
        Real r = policies::raise_domain_error("boost::math::hypergeometric_pFq", "Got p == 1 and q == 0 and |z| > 1, result is imaginary", z, pol);
        return std::make_pair(r, r);
      }
      std::pair<Real, Real> r = hypergeometric_pFq_checked_series_impl(aj, bj, Real(1 / z), pol, Real(tol), termination, log_scale);
      Real mul = pow(-z, -*aj.begin());
      r.first *= mul;
      r.second *= mul;
      return r;
    }
  }

  if (aj.size() > bj.size())
  {
    if (aj.size() == bj.size() + 1)
    {
      if (fabs(z) > 1)
      {
        Real r = policies::raise_domain_error("boost::math::hypergeometric_pFq", "Got p == q+1 and |z| > 1, series does not converge", z, pol);
        return std::make_pair(r, r);
      }
      if (fabs(z) == 1)
      {
        Real s = 0;
        for (auto i = bj.begin(); i != bj.end(); ++i)
          s += *i;
        for (auto i = aj.begin(); i != aj.end(); ++i)
          s -= *i;
        if ((z == 1) && (s <= 0))
        {
          Real r = policies::raise_domain_error("boost::math::hypergeometric_pFq", "Got p == q+1 and |z| == 1, in a situation where the series does not converge", z, pol);
          return std::make_pair(r, r);
        }
        if ((z == -1) && (s <= -1))
        {
          Real r = policies::raise_domain_error("boost::math::hypergeometric_pFq", "Got p == q+1 and |z| == 1, in a situation where the series does not converge", z, pol);
          return std::make_pair(r, r);
        }
      }
    }
    else
    {
      Real r = policies::raise_domain_error("boost::math::hypergeometric_pFq", "Got p > q+1, series does not converge", z, pol);
      return std::make_pair(r, r);
    }
  }

  while (!termination(k))
  {
    for (auto ai = aj.begin(); ai != aj.end(); ++ai)
    {
      term *= *ai + k;
    }
    if (term == 0)
    {
      // There is a negative integer in the aj's:
      return std::make_pair(result, abs_result);
    }
    for (auto bi = bj.begin(); bi != bj.end(); ++bi)
    {
      if (*bi + k == 0)
      {
        // The series is undefined:
        result = boost::math::policies::raise_domain_error("boost::math::hypergeometric_pFq<%1%>", "One of the b values was the negative integer %1%", *bi, pol);
        return std::make_pair(result, result);
      }
      term /= *bi + k;
    }
    term *= z;
    ++k;
    term /= k;
    //std::cout << k << " " << *bj.begin() + k << " " << result << " " << term << /*" " << term_at_k(*aj.begin(), *bj.begin(), z, k, pol) <<*/ std::endl;
    result += term;
    abs_result += abs(term);
    //std::cout << "k = " << k << " term = " << term * exp(log_scale) << " result = " << result * exp(log_scale) << " abs_result = " << abs_result * exp(log_scale) << std::endl;

    //
    // Rescaling:
    //
    if (fabs(abs_result) >= upper_limit)
    {
      abs_result /= scaling_factor;
      result /= scaling_factor;
      term /= scaling_factor;
      log_scale += log_scaling_factor;
      local_scaling += log_scaling_factor;
    }
    if (fabs(abs_result) < lower_limit)
    {
      abs_result *= scaling_factor;
      result *= scaling_factor;
      term *= scaling_factor;
      log_scale -= log_scaling_factor;
      local_scaling -= log_scaling_factor;
    }

    if ((abs(result * tol) > abs(term)) && (abs(term0) > abs(term)))
      break;
    // if (abs_result * tol > abs(result))
    // {
    //   // Check if result is so small compared to abs_resuslt that there are no longer any
    //   // correct bits... we require two consecutive passes here before aborting to
    //   // avoid false positives when result transiently drops to near zero then rebounds.
    //   if (have_no_correct_bits)
    //   {
    //     // We have no correct bits in the result... just give up!
    //     result = boost::math::policies::raise_evaluation_error("boost::math::hypergeometric_pFq<%1%>", "Cancellation is so severe that no bits in the reuslt are correct, last result was %1%", Real(result * exp(Real(log_scale))), pol);
    //     return std::make_pair(result, result);
    //   }
    //   else
    //     have_no_correct_bits = true;
    // }
    // else
    //   have_no_correct_bits = false;
    term0 = term;
  }
  //std::cout << "result = " << result << std::endl;
  //std::cout << "local_scaling = " << local_scaling << std::endl;
  //std::cout << "Norm result = " << std::setprecision(35) << boost::multiprecision::mpfr_float_50(result) * exp(boost::multiprecision::mpfr_float_50(local_scaling)) << std::endl;
  //
  // We have to be careful when one of the b's crosses the origin:
  //
  if(bj.size() > BOOST_MATH_PFQ_MAX_B_TERMS)
    policies::raise_domain_error<Real>("boost::math::hypergeometric_pFq<%1%>(Seq, Seq, %1%)",
                                       "The number of b terms must be less than the value of BOOST_MATH_PFQ_MAX_B_TERMS (" BOOST_STRINGIZE(BOOST_MATH_PFQ_MAX_B_TERMS)  "), but got %1%.",
                                       Real(bj.size()), pol);

  unsigned crossover_locations[BOOST_MATH_PFQ_MAX_B_TERMS];

  unsigned N_crossovers = set_crossover_locations(aj, bj, z, crossover_locations);

  bool terminate = false;   // Set to true if one of the a's passes through the origin and terminates the series.

  for (unsigned n = 0; n < N_crossovers; ++n)
  {
    if (k < crossover_locations[n])
    {
      for (auto ai = aj.begin(); ai != aj.end(); ++ai)
      {
        if ((*ai < 0) && (floor(*ai) == *ai) && (*ai > static_cast<decltype(*ai)>(crossover_locations[n])))
          return std::make_pair(result, abs_result);  // b's will never cross the origin!
      }
      //
      // local results:
      //
      Real loop_result = 0;
      Real loop_abs_result = 0;
      long long loop_scale = 0;
      //
      // loop_error_scale will be used to increase the size of the error
      // estimate (absolute sum), based on the errors inherent in calculating
      // the pochhammer symbols.
      //
      Real loop_error_scale = 0;
      //boost::multiprecision::mpfi_float err_est = 0;
      //
      // b hasn't crossed the origin yet and the series may spring back into life at that point
      // so we need to jump forward to that term and then evaluate forwards and backwards from there:
      //
      unsigned s = crossover_locations[n];
      std::uintmax_t backstop = k;
      long long s1(1), s2(1);
      term = 0;
      for (auto ai = aj.begin(); ai != aj.end(); ++ai)
      {
        if ((floor(*ai) == *ai) && (*ai < 0) && (-*ai <= static_cast<decltype(*ai)>(s)))
        {
          // One of the a terms has passed through zero and terminated the series:
          terminate = true;
          break;
        }
        else
        {
          int ls = 1;
          Real p = log_pochhammer(*ai, s, pol, &ls);
          s1 *= ls;
          term += p;
          loop_error_scale = (std::max)(p, loop_error_scale);
          //err_est += boost::multiprecision::mpfi_float(p);
        }
      }
      //std::cout << "term = " << term << std::endl;
      if (terminate)
        break;
      for (auto bi = bj.begin(); bi != bj.end(); ++bi)
      {
        int ls = 1;
        Real p = log_pochhammer(*bi, s, pol, &ls);
        s2 *= ls;
        term -= p;
        loop_error_scale = (std::max)(p, loop_error_scale);
        //err_est -= boost::multiprecision::mpfi_float(p);
      }
      //std::cout << "term = " << term << std::endl;
      Real p = lgamma(Real(s + 1), pol);
      term -= p;
      loop_error_scale = (std::max)(p, loop_error_scale);
      //err_est -= boost::multiprecision::mpfi_float(p);
      p = s * log(fabs(z));
      term += p;
      loop_error_scale = (std::max)(p, loop_error_scale);
      //err_est += boost::multiprecision::mpfi_float(p);
      //err_est = exp(err_est);
      //std::cout << err_est << std::endl;
      //
      // Convert loop_error scale to the absolute error
      // in term after exp is applied:
      //
      loop_error_scale *= tools::epsilon<Real>();
      //
      // Convert to relative error after exp:
      //
      loop_error_scale = fabs(expm1(loop_error_scale, pol));
      //
      // Convert to multiplier for the error term:
      //
      loop_error_scale /= tools::epsilon<Real>();

      if (z < 0)
        s1 *= (s & 1 ? -1 : 1);

      if (term <= tools::log_min_value<Real>())
      {
        // rescale if we can:
        long long scale = lltrunc(floor(term - tools::log_min_value<Real>()) - 2);
        term -= scale;
        loop_scale += scale;
      }
      if (term > 10)
      {
        int scale = itrunc(floor(term));
        term -= scale;
        loop_scale += scale;
      }
      //std::cout << "term = " << term << std::endl;
      term = s1 * s2 * exp(term);
      //std::cout << "term = " << term << std::endl;
      //std::cout << "loop_scale = " << loop_scale << std::endl;
      k = s;
      term0 = term;
      long long saved_loop_scale = loop_scale;
      bool terms_are_growing = true;
      bool trivial_small_series_check = false;
      do
      {
        loop_result += term;
        loop_abs_result += fabs(term);
        //std::cout << "k = " << k << " term = " << term * exp(loop_scale) << " result = " << loop_result * exp(loop_scale) << " abs_result = " << loop_abs_result * exp(loop_scale) << std::endl;
        if (fabs(loop_result) >= upper_limit)
        {
          loop_result /= scaling_factor;
          loop_abs_result /= scaling_factor;
          term /= scaling_factor;
          loop_scale += log_scaling_factor;
        }
        if (fabs(loop_result) < lower_limit)
        {
          loop_result *= scaling_factor;
          loop_abs_result *= scaling_factor;
          term *= scaling_factor;
          loop_scale -= log_scaling_factor;
        }
        term_m1 = term;
        for (auto ai = aj.begin(); ai != aj.end(); ++ai)
        {
          term *= *ai + k;
        }
        if (term == 0)
        {
          // There is a negative integer in the aj's:
          return std::make_pair(result, abs_result);
        }
        for (auto bi = bj.begin(); bi != bj.end(); ++bi)
        {
          if (*bi + k == 0)
          {
            // The series is undefined:
            result = boost::math::policies::raise_domain_error("boost::math::hypergeometric_pFq<%1%>", "One of the b values was the negative integer %1%", *bi, pol);
            return std::make_pair(result, result);
          }
          term /= *bi + k;
        }
        term *= z / (k + 1);

        ++k;
        diff = fabs(term / loop_result);
        terms_are_growing = fabs(term) > fabs(term_m1);
        if (!trivial_small_series_check && !terms_are_growing)
        {
          //
          // Now that we have started to converge, check to see if the value of
          // this local sum is trivially small compared to the result.  If so
          // abort this part of the series.
          //
          trivial_small_series_check = true;
          Real d;
          if (loop_scale > local_scaling)
          {
            long long rescale = local_scaling - loop_scale;
            if (rescale < tools::log_min_value<Real>())
              d = 1;  // arbitrary value, we want to keep going
            else
              d = fabs(term / (result * exp(Real(rescale))));
          }
          else
          {
            long long rescale = loop_scale - local_scaling;
            if (rescale < tools::log_min_value<Real>())
              d = 0;  // terminate this loop
            else
              d = fabs(term * exp(Real(rescale)) / result);
          }
          if (d < boost::math::policies::get_epsilon<Real, Policy>())
            break;
        }
      } while (!termination(k - s) && ((diff > boost::math::policies::get_epsilon<Real, Policy>()) || terms_are_growing));

      //std::cout << "Norm loop result = " << std::setprecision(35) << boost::multiprecision::mpfr_float_50(loop_result)* exp(boost::multiprecision::mpfr_float_50(loop_scale)) << std::endl;
      //
      // We now need to combine the results of the first series summation with whatever
      // local results we have now.  First though, rescale abs_result by loop_error_scale
      // to factor in the error in the pochhammer terms at the start of this block:
      //
      std::uintmax_t next_backstop = k;
      loop_abs_result += loop_error_scale * fabs(loop_result);
      if (loop_scale > local_scaling)
      {
        //
        // Need to shrink previous result:
        //
        long long rescale = local_scaling - loop_scale;
        local_scaling = loop_scale;
        log_scale -= rescale;
        Real ex = exp(Real(rescale));
        result *= ex;
        abs_result *= ex;
        result += loop_result;
        abs_result += loop_abs_result;
      }
      else if (local_scaling > loop_scale)
      {
        //
        // Need to shrink local result:
        //
        long long rescale = loop_scale - local_scaling;
        Real ex = exp(Real(rescale));
        loop_result *= ex;
        loop_abs_result *= ex;
        result += loop_result;
        abs_result += loop_abs_result;
      }
      else
      {
        result += loop_result;
        abs_result += loop_abs_result;
      }
      //
      // Now go backwards as well:
      //
      k = s;
      term = term0;
      loop_result = 0;
      loop_abs_result = 0;
      loop_scale = saved_loop_scale;
      trivial_small_series_check = false;
      do
      {
        --k;
        if (k == backstop)
          break;
        term_m1 = term;
        for (auto ai = aj.begin(); ai != aj.end(); ++ai)
        {
          term /= *ai + k;
        }
        for (auto bi = bj.begin(); bi != bj.end(); ++bi)
        {
          if (*bi + k == 0)
          {
            // The series is undefined:
            result = boost::math::policies::raise_domain_error("boost::math::hypergeometric_pFq<%1%>", "One of the b values was the negative integer %1%", *bi, pol);
            return std::make_pair(result, result);
          }
          term *= *bi + k;
        }
        term *= (k + 1) / z;
        loop_result += term;
        loop_abs_result += fabs(term);

        if (!trivial_small_series_check && (fabs(term) < fabs(term_m1)))
        {
          //
          // Now that we have started to converge, check to see if the value of
          // this local sum is trivially small compared to the result.  If so
          // abort this part of the series.
          //
          trivial_small_series_check = true;
          Real d;
          if (loop_scale > local_scaling)
          {
            long long rescale = local_scaling - loop_scale;
            if (rescale < tools::log_min_value<Real>())
              d = 1;  // keep going
            else
              d = fabs(term / (result * exp(Real(rescale))));
          }
          else
          {
            long long rescale = loop_scale - local_scaling;
            if (rescale < tools::log_min_value<Real>())
              d = 0;  // stop, underflow
            else
              d = fabs(term * exp(Real(rescale)) / result);
          }
          if (d < boost::math::policies::get_epsilon<Real, Policy>())
            break;
        }

        //std::cout << "k = " << k << " result = " << result << " abs_result = " << abs_result << std::endl;
        if (fabs(loop_result) >= upper_limit)
        {
          loop_result /= scaling_factor;
          loop_abs_result /= scaling_factor;
          term /= scaling_factor;
          loop_scale += log_scaling_factor;
        }
        if (fabs(loop_result) < lower_limit)
        {
          loop_result *= scaling_factor;
          loop_abs_result *= scaling_factor;
          term *= scaling_factor;
          loop_scale -= log_scaling_factor;
        }
        diff = fabs(term / loop_result);
      } while (!termination(s - k) && ((diff > boost::math::policies::get_epsilon<Real, Policy>()) || (fabs(term) > fabs(term_m1))));

      //std::cout << "Norm loop result = " << std::setprecision(35) << boost::multiprecision::mpfr_float_50(loop_result)* exp(boost::multiprecision::mpfr_float_50(loop_scale)) << std::endl;
      //
      // We now need to combine the results of the first series summation with whatever
      // local results we have now.  First though, rescale abs_result by loop_error_scale
      // to factor in the error in the pochhammer terms at the start of this block:
      //
      loop_abs_result += loop_error_scale * fabs(loop_result);
      //
      if (loop_scale > local_scaling)
      {
        //
        // Need to shrink previous result:
        //
        long long rescale = local_scaling - loop_scale;
        local_scaling = loop_scale;
        log_scale -= rescale;
        Real ex = exp(Real(rescale));
        result *= ex;
        abs_result *= ex;
        result += loop_result;
        abs_result += loop_abs_result;
      }
      else if (local_scaling > loop_scale)
      {
        //
        // Need to shrink local result:
        //
        long long rescale = loop_scale - local_scaling;
        Real ex = exp(Real(rescale));
        loop_result *= ex;
        loop_abs_result *= ex;
        result += loop_result;
        abs_result += loop_abs_result;
      }
      else
      {
        result += loop_result;
        abs_result += loop_abs_result;
      }
      //
      // Reset k to the largest k we reached
      //
      k = next_backstop;
    }
  }

  return std::make_pair(result, abs_result);
}

} // namespace detail

template <class Seq, class Real, class Policy>
inline typename tools::promote_args<Real, typename Seq::value_type>::type hypergeometric_pFq(
    const Seq& aj,
    const Seq& bj,
    const Real& z,
    Real* p_abs_error,
    const Policy& pol,
    const Real& tol,
    const R_xlen_t& max_iter
) {
  typedef typename tools::promote_args<Real, typename Seq::value_type>::type result_type;
  typedef typename policies::evaluation<result_type, Policy>::type value_type;
  // typedef typename policies::normalise<
  //   Policy,
  //   policies::promote_float<false>,
  //   policies::promote_double<false>,
  //   policies::discrete_quantile<>,
  //   policies::assert_undefined<> >::type forwarding_policy;

  BOOST_MATH_STD_USING

  long long scale = 0;
  std::pair<value_type, value_type> r = boost::math::detail::hypergeometric_pFq_checked_series_impl(aj, bj, value_type(z), pol, value_type(tol), boost::math::detail::iteration_terminator(max_iter), scale);
  r.first *= exp(Real(scale));
  r.second *= exp(Real(scale));
  if (p_abs_error)
    *p_abs_error = static_cast<Real>(r.second) * boost::math::tools::epsilon<Real>();
  return policies::checked_narrowing_cast<result_type, Policy>(r.first, "boost::math::hypergeometric_pFq<%1%>(%1%,%1%,%1%)");
}

} // namespace math
} // namespace boost

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
//' @param tol \code{NULL} or Numeric (1L) as convergence tolerance.
//' If \code{NULL}, use the epsilon for the \code{value_type} determined by
//' \code{boost::math::tools::promote_args} (e.g. for \code{double} input,
//' usually use the epsilon for \code{long double}).
//' @param max_iter Integer (1L) as iteration limit.
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
                    const Nullable<NumericVector>& tol,
                    const R_xlen_t& max_iter,
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
  double tol_ = -1.0;
  if (tol.isNotNull()) {
    NumericVector tol_vec = as<NumericVector>(tol);
    tol_ = tol_vec(0);
  }
  // Main code
  typedef typename boost::math::policies::policy<> MyPol;
  double out = R_NaN;
  if (prec.isNull()) {
    double p_abs_error;
    out = boost::math::hypergeometric_pFq(U, L, x, &p_abs_error, MyPol(), tol_, max_iter);
  } else {
    IntegerVector prec_ = as<IntegerVector>(prec);
    unsigned int prec_use = (unsigned) prec_(0);
    typedef typename boost::multiprecision::number<boost::multiprecision::backends::mpfr_float_backend<0>> prec_float;
    boost::math::scoped_precision<prec_float> scoped(prec_use);
    prec_float x_ = x;
    std::vector<prec_float> U_ = conv_vec_prec<NumericVector, double, prec_float>(U);
    std::vector<prec_float> L_ = conv_vec_prec<NumericVector, double, prec_float>(L);
    prec_float p_abs_error;
    prec_float tol__ = tol_;
    prec_float out_ = boost::math::hypergeometric_pFq(
      U_, L_, x_, &p_abs_error, MyPol(), tol__, max_iter
    );
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
