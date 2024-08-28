// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/cbbinom.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_pcbbinom
NumericVector cpp_pcbbinom(const NumericVector& q, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& lower_tail, const bool& log_p, const NumericVector& tol, const IntegerVector& max_iter);
static SEXP _cbbinom_cpp_pcbbinom_try(SEXP qSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP lower_tailSEXP, SEXP log_pSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type q(qSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_p(log_pSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pcbbinom(q, size, alpha, beta, lower_tail, log_p, tol, max_iter));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _cbbinom_cpp_pcbbinom(SEXP qSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP lower_tailSEXP, SEXP log_pSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_cbbinom_cpp_pcbbinom_try(qSEXP, sizeSEXP, alphaSEXP, betaSEXP, lower_tailSEXP, log_pSEXP, tolSEXP, max_iterSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cpp_qcbbinom
NumericVector cpp_qcbbinom(const NumericVector& p, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& lower_tail, const bool& log_p, const NumericVector& p_tol, const IntegerVector& p_max_iter, const NumericVector& root_tol, const IntegerVector& root_max_iter);
static SEXP _cbbinom_cpp_qcbbinom_try(SEXP pSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP lower_tailSEXP, SEXP log_pSEXP, SEXP p_tolSEXP, SEXP p_max_iterSEXP, SEXP root_tolSEXP, SEXP root_max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_p(log_pSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type p_tol(p_tolSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type p_max_iter(p_max_iterSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type root_tol(root_tolSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type root_max_iter(root_max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_qcbbinom(p, size, alpha, beta, lower_tail, log_p, p_tol, p_max_iter, root_tol, root_max_iter));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _cbbinom_cpp_qcbbinom(SEXP pSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP lower_tailSEXP, SEXP log_pSEXP, SEXP p_tolSEXP, SEXP p_max_iterSEXP, SEXP root_tolSEXP, SEXP root_max_iterSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_cbbinom_cpp_qcbbinom_try(pSEXP, sizeSEXP, alphaSEXP, betaSEXP, lower_tailSEXP, log_pSEXP, p_tolSEXP, p_max_iterSEXP, root_tolSEXP, root_max_iterSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// dcbblp
NumericMatrix dcbblp(const NumericVector& x, const NumericVector& m, const NumericVector& a, const NumericVector& b, const NumericVector& tol, const IntegerVector& max_iter);
static SEXP _cbbinom_dcbblp_try(SEXP xSEXP, SEXP mSEXP, SEXP aSEXP, SEXP bSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(dcbblp(x, m, a, b, tol, max_iter));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _cbbinom_dcbblp(SEXP xSEXP, SEXP mSEXP, SEXP aSEXP, SEXP bSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_cbbinom_dcbblp_try(xSEXP, mSEXP, aSEXP, bSEXP, tolSEXP, max_iterSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cpp_rcbbinom
NumericVector cpp_rcbbinom(const int& n, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const NumericVector& p_tol, const IntegerVector& p_max_iter, const NumericVector& root_tol, const IntegerVector& root_max_iter);
static SEXP _cbbinom_cpp_rcbbinom_try(SEXP nSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP p_tolSEXP, SEXP p_max_iterSEXP, SEXP root_tolSEXP, SEXP root_max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type p_tol(p_tolSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type p_max_iter(p_max_iterSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type root_tol(root_tolSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type root_max_iter(root_max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rcbbinom(n, size, alpha, beta, p_tol, p_max_iter, root_tol, root_max_iter));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _cbbinom_cpp_rcbbinom(SEXP nSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP p_tolSEXP, SEXP p_max_iterSEXP, SEXP root_tolSEXP, SEXP root_max_iterSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_cbbinom_cpp_rcbbinom_try(nSEXP, sizeSEXP, alphaSEXP, betaSEXP, p_tolSEXP, p_max_iterSEXP, root_tolSEXP, root_max_iterSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// gen_hypergeo
double gen_hypergeo(const NumericVector& U, const NumericVector& L, const double& x, const double& tol, const R_xlen_t& max_iter, const bool& check_mode, const bool& log);
static SEXP _cbbinom_gen_hypergeo_try(SEXP USEXP, SEXP LSEXP, SEXP xSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP check_modeSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type U(USEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const R_xlen_t& >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type check_mode(check_modeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_hypergeo(U, L, x, tol, max_iter, check_mode, log));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _cbbinom_gen_hypergeo(SEXP USEXP, SEXP LSEXP, SEXP xSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP check_modeSEXP, SEXP logSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_cbbinom_gen_hypergeo_try(USEXP, LSEXP, xSEXP, tolSEXP, max_iterSEXP, check_modeSEXP, logSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _cbbinom_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("NumericVector(*cpp_pcbbinom)(const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const bool&,const bool&,const NumericVector&,const IntegerVector&)");
        signatures.insert("NumericVector(*cpp_qcbbinom)(const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const bool&,const bool&,const NumericVector&,const IntegerVector&,const NumericVector&,const IntegerVector&)");
        signatures.insert("NumericMatrix(*dcbblp)(const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const IntegerVector&)");
        signatures.insert("NumericVector(*cpp_rcbbinom)(const int&,const NumericVector&,const NumericVector&,const NumericVector&,const NumericVector&,const IntegerVector&,const NumericVector&,const IntegerVector&)");
        signatures.insert("double(*gen_hypergeo)(const NumericVector&,const NumericVector&,const double&,const double&,const R_xlen_t&,const bool&,const bool&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _cbbinom_RcppExport_registerCCallable() { 
    R_RegisterCCallable("cbbinom", "_cbbinom_cpp_pcbbinom", (DL_FUNC)_cbbinom_cpp_pcbbinom_try);
    R_RegisterCCallable("cbbinom", "_cbbinom_cpp_qcbbinom", (DL_FUNC)_cbbinom_cpp_qcbbinom_try);
    R_RegisterCCallable("cbbinom", "_cbbinom_dcbblp", (DL_FUNC)_cbbinom_dcbblp_try);
    R_RegisterCCallable("cbbinom", "_cbbinom_cpp_rcbbinom", (DL_FUNC)_cbbinom_cpp_rcbbinom_try);
    R_RegisterCCallable("cbbinom", "_cbbinom_gen_hypergeo", (DL_FUNC)_cbbinom_gen_hypergeo_try);
    R_RegisterCCallable("cbbinom", "_cbbinom_RcppExport_validate", (DL_FUNC)_cbbinom_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_cbbinom_cpp_pcbbinom", (DL_FUNC) &_cbbinom_cpp_pcbbinom, 8},
    {"_cbbinom_cpp_qcbbinom", (DL_FUNC) &_cbbinom_cpp_qcbbinom, 10},
    {"_cbbinom_dcbblp", (DL_FUNC) &_cbbinom_dcbblp, 6},
    {"_cbbinom_cpp_rcbbinom", (DL_FUNC) &_cbbinom_cpp_rcbbinom, 8},
    {"_cbbinom_gen_hypergeo", (DL_FUNC) &_cbbinom_gen_hypergeo, 7},
    {"_cbbinom_RcppExport_registerCCallable", (DL_FUNC) &_cbbinom_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cbbinom(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
