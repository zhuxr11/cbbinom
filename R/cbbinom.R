#' The Continuous Beta-Binomial Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' a continuous analog to the beta-binomial distribution with parameters \code{size},
#' \code{alpha} and \code{beta}. The usage and help pages are modeled
#' on the d-p-q-r families of functions for the commonly-used distributions
#' in the \code{stats} package.
#'
#' Derived from the continuous binomial distribution (Ilienko 2013), the continuous beta-binomial
#' distribution is defined as:
#' \deqn{P(x|n,\alpha,\beta)=\int_0^1\frac{B_{1-p}(n+1-x,x)}{B(n+1-x,x)}\frac{p^{\alpha-1}(1-p)^{\beta-1}}{B(\alpha,\beta)}dp,}
#' where \eqn{x} is the quantile, \eqn{n} is the size, \eqn{B_p(a,b)=\int_0^p{u^{a-1}(1-u)^{b-1}du}}
#' is the incomplete beta function.
#'
#' When simplified, the distribution becomes:
#' \deqn{P(x|n,\alpha,\beta)=\frac{\Gamma(n+1)B(n+1-x+\beta,\alpha)}{\Gamma(x)\Gamma(n+2-x)B(\alpha,\beta)}{}_3F_2(a;b;z),}
#' where \eqn{{}_3F_2(a;b;z)} is \link[hypergeo2:genhypergeo]{generalized hypergeometric function}, \eqn{a=\{1-x,n+1-x,n+1-x+\beta\}},
#' \eqn{b=\{n+2-x,n+1-x+\alpha+\beta\}}, \eqn{z=1}.
#'
#' Heuristically speaking, this distribution spreads the standard probability mass
#' at integer \code{x} to the interval \code{[x, x + 1]} in a continuous manner.
#' As a result, the distribution looks like a smoothed version of the standard,
#' discrete beta-binomial but shifted slightly to the right. The support of the continuous
#' beta-binomial is \code{[0, size + 1]}, and the mean is approximately
#' \code{size * alpha / (alpha + beta) + 1/2}.
#'
#' Supplying \code{ncp != 0} moves the support of beta-binomial to \code{[ncp, size + 1 + ncp]}. For example,
#' to build a continuous beta-binomial with approximately non-shifted mean, use \code{ncp = -0.5}.
#'
#' These functions are also available in \code{\link[Rcpp:Rcpp-package]{Rcpp}}
#' as \code{cbbinom::cpp_[d/p/q/r]cbbinom()}, and their non-vectorized versions
#' in \code{\link[Rcpp:Rcpp-package]{Rcpp}} as \code{cbbinom::[d/p/q/r]cbbinom_()}.
#' To use them, please use \code{[[Rcpp::depends(cbbinom)]]} and \code{#include <cbbinom.h>}.
#'
#' @param x,q vector of quantiles.
#' @param alpha,beta non-negative parameters of the Beta distribution.
#' @inheritParams stats::Binomial
#' @inheritParams stats::Beta
#' @param prec arguments passed on to \code{\link[hypergeo2]{genhypergeo}},
#' vectorized and recycled along with distribution parameters.
#' @param tol,max_iter arguments passed on to \code{\link[stats]{uniroot}},
#' vectorized and recycled along with distribution parameters.
#'
#' @return
#' \code{dcbbinom} gives the density, \code{pcbbinom} the distribution function,
#' \code{qcbbinom} the quantile function, and \code{rcbbinom} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN}, with a warning.
#'
#' The length of the result is determined by \code{n} for \code{rcbbinom},
#' and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result.
#' Only the first elements of the logical arguments are used.
#'
#' @note Change log:
#' \itemize{
#'   \item{0.1.0 Xiurui Zhu - Initiate the function.}
#'   \item{0.2.0 Xiurui Zhu - Re-implement distribution function with \code{\link[BH:BH-package]{BH}} package,
#'     add \code{NULL} default tolerance, and add precision parameters.}
#' }
#'
#' @references Ilienko, Andreii (2013). Continuous counterparts of Poisson and binomial
#' distributions and their properties. Annales Univ. Sci. Budapest., Sect. Comp.
#' 39: 137-147. \url{http://ac.inf.elte.hu/Vol_039_2013/137_39.pdf}
#'
#' @name cbbinom
#' @examples
#' # Density function
#' dcbbinom(x = 5, size = 10, alpha = 2, beta = 4)
#' # Distribution function
#' (test_val <- pcbbinom(q = 5, size = 10, alpha = 2, beta = 4))
#' # Quantile function
#' qcbbinom(p = test_val, size = 10, alpha = 2, beta = 4)
#' # Random generation
#' set.seed(1111L)
#' rcbbinom(n = 10L, size = 10, alpha = 2, beta = 4)
NULL

# Normalize precision
norm_prec <- function(prec) {
  if (is.null(prec) == TRUE) {
    NULL
  } else {
    res <- as.list(prec)
    lapply(res, function(x) {if (is.null(x) == TRUE) NULL else as.integer(x)})
  }
}

#' @aliases dcbbinom
#' @export
#' @rdname cbbinom
dcbbinom <- function(x, size, alpha = 1, beta = 1, ncp = 0, log = FALSE, prec = NULL) {
  p <- cpp_dcbbinom(x = as.numeric(x - ncp),
                    size = as.numeric(size),
                    alpha = as.numeric(alpha),
                    beta = as.numeric(beta),
                    log = TRUE,
                    prec = norm_prec(prec))
  p[is.na(p) == TRUE] <- -Inf
  if (log == FALSE) {
    p <- exp(p)
  }
  p
}

#' @aliases pcbbinom
#' @export
#' @rdname cbbinom
pcbbinom <- function(q, size, alpha = 1, beta = 1, ncp = 0,
                     lower.tail = TRUE, log.p = FALSE, prec = NULL) {
  cpp_pcbbinom(q = as.numeric(q - ncp),
               size = as.numeric(size),
               alpha = as.numeric(alpha),
               beta = as.numeric(beta),
               lower_tail = as.logical(lower.tail[[1L]]),
               log_p = as.logical(log.p[[1L]]),
               prec = norm_prec(prec))
}

#' @aliases qcbbinom
#' @export
#' @rdname cbbinom
qcbbinom <- function(p, size, alpha = 1, beta = 1, ncp = 0,
                     lower.tail = TRUE, log.p = FALSE, prec = NULL,
                     tol = 1e-6, max_iter = 10000L) {
  cpp_qcbbinom(p = as.numeric(p),
               size = as.numeric(size),
               alpha = as.numeric(alpha),
               beta = as.numeric(beta),
               lower_tail = as.logical(lower.tail[[1L]]),
               log_p = as.logical(log.p[[1L]]),
               prec = norm_prec(prec),
               tol = as.numeric(tol),
               max_iter = as.integer(max_iter)) + ncp
}

#' @aliases rcbbinom
#' @export
#' @rdname cbbinom
rcbbinom <- function(n, size, alpha = 1, beta = 1, ncp = 0, prec = NULL,
                     tol = 1e-6, max_iter = 10000L) {
  cpp_rcbbinom(n = as.integer(n[[1L]]),
               size = as.numeric(size),
               alpha = as.numeric(alpha),
               beta = as.numeric(beta),
               prec = norm_prec(prec),
               tol = as.numeric(tol),
               max_iter = as.integer(max_iter)) + ncp
}
