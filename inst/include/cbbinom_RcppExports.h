// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_cbbinom_RCPPEXPORTS_H_GEN_
#define RCPP_cbbinom_RCPPEXPORTS_H_GEN_

#include <Rcpp.h>

namespace cbbinom {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("cbbinom", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("cbbinom", "_cbbinom_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in cbbinom");
            }
        }
    }

    inline NumericVector cpp_pcbbinom(const NumericVector& q, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& lower_tail, const bool& log_p, const NumericVector& tol, const IntegerVector& max_iter) {
        typedef SEXP(*Ptr_cpp_pcbbinom)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cpp_pcbbinom p_cpp_pcbbinom = NULL;
        if (p_cpp_pcbbinom == NULL) {
            validateSignature("NumericVector(*cpp_pcbbinom)(const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const bool&,const bool&,const NumericVector&,const IntegerVector&)");
            p_cpp_pcbbinom = (Ptr_cpp_pcbbinom)R_GetCCallable("cbbinom", "_cbbinom_cpp_pcbbinom");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_pcbbinom(Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(size)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(lower_tail)), Shield<SEXP>(Rcpp::wrap(log_p)), Shield<SEXP>(Rcpp::wrap(tol)), Shield<SEXP>(Rcpp::wrap(max_iter)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector cpp_qcbbinom(const NumericVector& p, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& lower_tail, const bool& log_p, const NumericVector& p_tol, const IntegerVector& p_max_iter, const NumericVector& root_tol, const IntegerVector& root_max_iter) {
        typedef SEXP(*Ptr_cpp_qcbbinom)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cpp_qcbbinom p_cpp_qcbbinom = NULL;
        if (p_cpp_qcbbinom == NULL) {
            validateSignature("NumericVector(*cpp_qcbbinom)(const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const bool&,const bool&,const NumericVector&,const IntegerVector&,const NumericVector&,const IntegerVector&)");
            p_cpp_qcbbinom = (Ptr_cpp_qcbbinom)R_GetCCallable("cbbinom", "_cbbinom_cpp_qcbbinom");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_qcbbinom(Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(size)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(lower_tail)), Shield<SEXP>(Rcpp::wrap(log_p)), Shield<SEXP>(Rcpp::wrap(p_tol)), Shield<SEXP>(Rcpp::wrap(p_max_iter)), Shield<SEXP>(Rcpp::wrap(root_tol)), Shield<SEXP>(Rcpp::wrap(root_max_iter)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericMatrix dcbblp(const NumericVector& x, const NumericVector& m, const NumericVector& a, const NumericVector& b, const NumericVector& tol, const IntegerVector& max_iter) {
        typedef SEXP(*Ptr_dcbblp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_dcbblp p_dcbblp = NULL;
        if (p_dcbblp == NULL) {
            validateSignature("NumericMatrix(*dcbblp)(const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const IntegerVector&)");
            p_dcbblp = (Ptr_dcbblp)R_GetCCallable("cbbinom", "_cbbinom_dcbblp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_dcbblp(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(m)), Shield<SEXP>(Rcpp::wrap(a)), Shield<SEXP>(Rcpp::wrap(b)), Shield<SEXP>(Rcpp::wrap(tol)), Shield<SEXP>(Rcpp::wrap(max_iter)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

    inline NumericVector cpp_rcbbinom(const int& n, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const NumericVector& p_tol, const IntegerVector& p_max_iter, const NumericVector& root_tol, const IntegerVector& root_max_iter) {
        typedef SEXP(*Ptr_cpp_rcbbinom)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cpp_rcbbinom p_cpp_rcbbinom = NULL;
        if (p_cpp_rcbbinom == NULL) {
            validateSignature("NumericVector(*cpp_rcbbinom)(const int&,const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const IntegerVector&,const NumericVector&,const IntegerVector&)");
            p_cpp_rcbbinom = (Ptr_cpp_rcbbinom)R_GetCCallable("cbbinom", "_cbbinom_cpp_rcbbinom");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_rcbbinom(Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(size)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(p_tol)), Shield<SEXP>(Rcpp::wrap(p_max_iter)), Shield<SEXP>(Rcpp::wrap(root_tol)), Shield<SEXP>(Rcpp::wrap(root_max_iter)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline double gen_hypergeo(const NumericVector& U, const NumericVector& L, const double& x, const double& tol, const R_xlen_t& max_iter, const bool& check_mode, const bool& log) {
        typedef SEXP(*Ptr_gen_hypergeo)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_gen_hypergeo p_gen_hypergeo = NULL;
        if (p_gen_hypergeo == NULL) {
            validateSignature("double(*gen_hypergeo)(const NumericVector&,const NumericVector&,const double&,const double&,const R_xlen_t&,const bool&,const bool&)");
            p_gen_hypergeo = (Ptr_gen_hypergeo)R_GetCCallable("cbbinom", "_cbbinom_gen_hypergeo");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_gen_hypergeo(Shield<SEXP>(Rcpp::wrap(U)), Shield<SEXP>(Rcpp::wrap(L)), Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(tol)), Shield<SEXP>(Rcpp::wrap(max_iter)), Shield<SEXP>(Rcpp::wrap(check_mode)), Shield<SEXP>(Rcpp::wrap(log)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

}

#endif // RCPP_cbbinom_RCPPEXPORTS_H_GEN_
