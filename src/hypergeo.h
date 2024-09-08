#ifndef CBBINOM_HYPERGEO
#define CBBINOM_HYPERGEO

#include <Rcpp.h>
using namespace Rcpp;

double gen_hypergeo(const NumericVector& U,
                    const NumericVector& L,
                    const double& x,
                    const Nullable<IntegerVector>& prec,
                    const bool& check_mode,
                    const bool& log);

#endif
