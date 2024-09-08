
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Continuous beta-binomial distribution

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/cbbinom)](https://CRAN.R-project.org/package=cbbinom)
[![R-CMD-check](https://github.com/zhuxr11/cbbinom/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhuxr11/cbbinom/actions/workflows/R-CMD-check.yaml)
[![Download
stats](https://cranlogs.r-pkg.org/badges/grand-total/cbbinom)](https://CRAN.R-project.org/package=cbbinom)
<!-- badges: end -->

**Package**: [*cbbinom*](https://github.com/zhuxr11/cbbinom)
0.1.0.9000<br /> **Author**: Xiurui Zhu<br /> **Modified**: 2024-09-08
17:20:06<br /> **Compiled**: 2024-09-08 17:20:12

The goal of `cbbinom` is to implement continuous beta-binomial
distribution.

## Installation

### System requirements

If you are building from source (e.g. not installing binaries on
Windows), you need to prepare for two system requirements:
[GMP](https://gmplib.org/) and [MPFR](https://www.mpfr.org/), to
facilitate high-precision floating point types. You may follow the
installation instructions from their websites, or try the following
commands for quick (default) installation.

- For Windows with rtools40 (or newer) installed, please run:

``` bash
pacman -Syu mingw-w64-x86_64-make mingw-w64-x86_64-pkg-config gmp mpfr
```

- For macOS, please run:

``` bash
brew install gmp mpfr
```

- For linux, please build and install platform-specific libraries
  through their `configure` file, or install them using `conda`.

If the requirements are installed into their default paths (e.g. without
using the `--prefix` option), you are OK to go ahead installing the
package: `pkg-config` will take care finding them. However, if they are
not, you may first need to provide configuration arguments to
installation paths in one of the following ways (replacing
`<my_gmp_install>` and `<my_mpfr_install>` with string literals):

- Directly set `configure.args` in `install.packages()`:

``` r
my_gmp_install <- shQuote(<my_gmp_install>)
my_mpfr_install <- shQuote(<my_mpfr_install>)
install.packages(
  "cbbinom",
  configure.args = c(
    paste0("--with-gmp-include=", my_gmp_install, "/include"),
    paste0("--with-mpfr-include=", my_mpfr_install, "/include"),
    paste0("--with-gmp-lib=", my_gmp_install, "/lib"),
    paste0("--with-mpfr-lib=", my_mpfr_install, "/lib")
  )
)
```

- Set the following environment variables and run
  `install.packages("cbbinom")`:

``` bash
MY_GMP_INSTALL='<my_gmp_install>'
MY_MPFR_INSTALL='<my_mpfr_install>'
CBBINOM_GMP_INCLUDE="${MY_GMP_INSTALL}/include"
CBBINOM_MPFR_INCLUDE="${MY_MPFR_INSTALL}/include"
CBBINOM_GMP_LIB="${MY_GMP_INSTALL}/lib"
CBBINOM_MPFR_LIB="${MY_MPFR_INSTALL}/lib"
```

### The package

You can install the released version of `cbbinom` from
[CRAN](https://cran.r-project.org/) with:

``` r
install.packages("cbbinom")
```

Alternatively, you can install the developmental version of `cbbinom`
from [github](https://github.com/) with:

``` r
remotes::install_github("zhuxr11/cbbinom")
```

## Introduction to continuous beta-binomial distribution

The continuous beta-binomial distribution spreads the standard
probability mass of beta-binomial distribution at `x` to an interval
`[x, x + 1]` in a continuous manner. This can be validated via the
following plot, where we can see that the cumulative distribution
function (CDF) of the continuous beta-binomial distribution at `x + 1`
equals to that of the beta-binomial distribution at `x`.

``` r
library(cbbinom)
```

``` r
# The continuous beta-binomial CDF, shift by -1
cbbinom_plot_x <- seq(-1, 10, 0.01)
cbbinom_plot_y <- pcbbinom(
  q = cbbinom_plot_x,
  size = 10,
  alpha = 2,
  beta = 4,
  ncp = -1
)
# The beta-binomial CDF
bbinom_plot_x <- seq(0L, 10L, 1L)
bbinom_plot_y <- extraDistr::pbbinom(
  q = bbinom_plot_x,
  size = 10L,
  alpha = 2,
  beta = 4
)
ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_bar(
    data = data.frame(
      x = bbinom_plot_x,
      y = bbinom_plot_y
    ),
    stat = "identity"
  ) +
  ggplot2::geom_point(
    data = data.frame(
      x = cbbinom_plot_x,
      y = cbbinom_plot_y
    )
  ) +
  ggplot2::scale_x_continuous(
    n.breaks = diff(range(cbbinom_plot_x))
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "CDF(x)")
```

<img src="man/figures/README-cbbinom-vs-bbinom-p-1.png" width="100%" />

However, the central density at `x + 1/2` of the continuous
beta-binomial distribution may not equal to the corresponding
probability mass at `x`, especially around the summit and to the right
(since `alpha < beta`).

``` r
# The continuous beta-binomial CDF, shift by -1/2
cbbinom_plot_x_d <- seq(-1/2, 10 + 1/2, 0.01)
cbbinom_plot_y_d <- dcbbinom(
  x = cbbinom_plot_x_d,
  size = 10,
  alpha = 2,
  beta = 4,
  ncp = -1/2
)
# The beta-binomial CDF
bbinom_plot_x <- seq(0L, 10L, 1L)
bbinom_plot_y_d <- extraDistr::dbbinom(
  x = bbinom_plot_x,
  size = 10L,
  alpha = 2,
  beta = 4
)
ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_bar(
    data = data.frame(
      x = bbinom_plot_x,
      y = bbinom_plot_y_d
    ),
    stat = "identity"
  ) +
  ggplot2::geom_point(
    data = data.frame(
      x = cbbinom_plot_x_d,
      y = cbbinom_plot_y_d
    )
  ) +
  ggplot2::scale_x_continuous(
    n.breaks = diff(range(bbinom_plot_x))
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "CDF(x)")
```

<img src="man/figures/README-cbbinom-vs-bbinom-d-1.png" width="100%" />

For larger sizes, you may need higher precision than double for
accuracy, at the cost of computational speed.

``` r
cbbinom_plot_prec_x_p <- seq(0, 41, 0.1)
# Compute CDF at default (double) precision level
system.time(pcbbinom_plot_prec0_y <- pcbbinom(
  q = cbbinom_plot_prec_x_p,
  size = 40,
  alpha = 2,
  beta = 4,
  prec = NULL
))
#>    user  system elapsed 
#>    0.06    0.00    0.07
```

``` r
ggplot2::ggplot(data = data.frame(x = cbbinom_plot_prec_x_p,
                                  y = pcbbinom_plot_prec0_y),
                mapping = ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "CDF(x)")
```

<img src="man/figures/README-cbbinom-prec-p-1.png" width="100%" />

``` r
# Compute CDF at precision level 20
system.time(pcbbinom_plot_prec20_y <- pcbbinom(
  q = cbbinom_plot_prec_x_p,
  size = 40,
  alpha = 2,
  beta = 4,
  prec = 20L
))
#>    user  system elapsed 
#>    3.75    0.00    3.77
```

``` r
ggplot2::ggplot(data = data.frame(x = cbbinom_plot_prec_x_p,
                                  y = pcbbinom_plot_prec20_y),
                mapping = ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "CDF(x)")
```

<img src="man/figures/README-cbbinom-prec-p-2.png" width="100%" />

When computing probability density, it is advisable to use a higher
accuracy level than the one sufficient to compute the corresponding
cumulative probability (e.g. 20 -\> 25), since numerical derivative
taken to compute the probability density may be more sensitive to
computational errors.

``` r
cbbinom_plot_prec_x_d <- seq(0, 41, 0.1)
# Compute PDF at precision level 20
system.time(dcbbinom_plot_prec20_y <- dcbbinom(
  x = cbbinom_plot_prec_x_d,
  size = 40,
  alpha = 2,
  beta = 4,
  prec = 20L
))
#> Warning in cpp_dcbbinom(x = as.numeric(x - ncp), size = as.numeric(size), :
#> d[pcbbinom(q = 29.000000, size = 40.000000, alpha = 2.000000, beta =
#> 4.000000)]/dq = -0.000206 < 0, which is set to 0, since probability density
#> cannot be negative; you may use a higher [prec] level than 20
#>    user  system elapsed 
#>   27.22    0.02   27.26
```

``` r
ggplot2::ggplot(data = data.frame(x = cbbinom_plot_prec_x_d,
                                  y = dcbbinom_plot_prec20_y),
                mapping = ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "PDF(x)")
```

<img src="man/figures/README-cbbinom-prec-d-1.png" width="100%" />

``` r
# Compute PDF at precision level 25
system.time(dcbbinom_plot_prec25_y <- dcbbinom(
  x = cbbinom_plot_prec_x_d,
  size = 40,
  alpha = 2,
  beta = 4,
  prec = 25L
))
#>    user  system elapsed 
#>   38.67    0.03   38.86
```

``` r
ggplot2::ggplot(data = data.frame(x = cbbinom_plot_prec_x_d,
                                  y = dcbbinom_plot_prec25_y),
                mapping = ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "PDF(x)")
```

<img src="man/figures/README-cbbinom-prec-d-2.png" width="100%" />

## Examples of continuous beta-binomial distribution

As the probability distributions in `stats` package, `cbbinom` provides
a full set of density, distribution function, quantile function and
random generation for the continuous beta-binomial distribution.

``` r
# Density function
dcbbinom(x = 5, size = 10, alpha = 2, beta = 4)
#> [1] 0.12669
```

``` r
# Distribution function
(test_val <- pcbbinom(q = 5, size = 10, alpha = 2, beta = 4))
#> [1] 0.7062937
```

``` r
# Quantile function
qcbbinom(p = test_val, size = 10, alpha = 2, beta = 4)
#> [1] 5
```

``` r
# Random generation
set.seed(1111L)
rcbbinom(n = 10L, size = 10, alpha = 2, beta = 4)
#>  [1] 3.359039 3.038286 7.110936 1.311321 5.264688 8.709005 6.720415 1.164210
#>  [9] 3.868370 1.332590
```

For mathematical details, please check the details section of
`?cbbinom`.

## `Rcpp` implementation of `stats::uniroot()`

As a bonus, `cbbinom` also exports an `Rcpp` implementation of
`stats::uniroot()` function, which may come in handy to solve equations,
especially the monotonic ones used in quantile functions. Here is an
example to calculate `qnorm` from `pnorm` in `Rcpp`.

``` cpp
#include <iostream>
#include "cbbinom.h"
using namespace cbbinom;

// Define a functor as pnorm() - p
class PnormEqn: public UnirootEqn
{
private:
  double mu;
  double sd;
  double p;
public:
  PnormEqn(const double mu_, const double sd_, const double p_):
    mu(mu_), sd(sd_), p(p_) {}
  double operator () (const double& x) const override {
    return R::pnorm(x, this->mu, this->sd, true, false) - this->p;
  }
};

// Compute quantiles
int main() {
  double p = 0.975;  // Quantile
  PnormEqn eqn_obj(0.0, 1.0, 0.975);
  double tol = 1e-6;
  int max_iter = 10000;
  double q = cbbinom::cpp_uniroot(-1000.0, 1000.0, -p, 1.0 - p, &eqn_obj, &tol, &max_iter);
  std::cout << "Quantile at " << p << "is: " << q << std::endl;
  return 0;
}
```
