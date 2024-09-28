#ifndef CBBINOM_IMPL_HEADER_GEN_H
#define CBBINOM_IMPL_HEADER_GEN_H

#include <Rcpp.h>
using namespace Rcpp;

double pcbbinom_(
    const double& q,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& lower_tail,
    const bool& log_p,
    const Nullable<IntegerVector>& prec
);

double dcbbinom_(
    const double& x,
    const double& size,
    const double& alpha,
    const double& beta,
    const bool& log,
    const Nullable<IntegerVector>& prec
);

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
);

#endif // CBBINOM_IMPL_HEADER_GEN_H
