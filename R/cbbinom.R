#' The Continuous Beta-Binomial Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' a continuous analog to the beta-binomial distribution with parameters \code{size}
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
#' where \eqn{{}_3F_2(a;b;z)} is [generalized hypergeometric function][gen_hypergeo], \eqn{a=\{1-x,n+1-x,n+1-x+\beta\}},
#' \eqn{b=\{n+2-x,n+1-x+\alpha+\beta\}}, \eqn{z=1}.
#'
#' Heuristically speaking, this distribution spreads the standard probability mass
#' at integer \code{x} to the interval \code{[x, x + 1]} in a continuous manner.
#' As a result, the distribution looks like a smoothed version of the standard,
#' discrete beta-binomial but shifted slightly to the right. The support of the continuous
#' beta-binomial is \code{[0, size + 1]}, and the mean is approximately
#' \code{size * alpha / (alpha + beta) + 1/2}.
#'
#' Supplying \code{ncp} moves the support of beta-binomial to \code{[ncp, size + 1 + ncp]}, e.g.
#' for the continuous beta-binomial with non-shifted mean, use \code{ncp = -0.5}.
#'
#' @param x,q vector of quantiles.
#' @param alpha,beta non-negative parameters of the Beta distribution.
#' @inheritParams stats::Binomial
#' @inheritParams stats::Beta
#' @param tol,max_iter arguments passed on to \code{\link{gen_hypergeo}}.
#' @param p_tol,p_max_iter same as \code{tol}, \code{max_iter}.
#' @param root_tol,root_max_iter arguments passed on to \code{\link[stats]{uniroot}}.
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

#' @aliases dcbbinom
#' @section Numerical computation of the density function:
#' For simplicity, the density function is computed numerically through differentiation.
#' To achieve higher numerical resolution (given that \eqn{d\ln{u}/du>1,0<u<1}), it is computed as:
#' \deqn{p(x|n,\alpha,\beta)=\frac{\partial{P(x|n,\alpha,\beta)}}{\partial{x}}=\frac{\partial\exp[\ln{P(x|n,\alpha,\beta)}]}{\partial{x}}}
#' When simplified, it becomes:
#' \deqn{p(x|n,\alpha,\beta)=\frac{\partial\exp[\ln{P(x|n,\alpha,\beta)}]}{\partial\ln{P(x|n,\alpha,\beta)}}\frac{\partial\ln{P(x|n,\alpha,\beta)}}{\partial{x}}=\frac{\partial\ln{P(x|n,\alpha,\beta)}}{\partial{x}}P(x|n,\alpha,\beta),}
#' where the first term is computed numerically and the second term is the distribution function.
#' @export
#' @rdname cbbinom
dcbbinom <- function(x, size, alpha = 1, beta = 1, ncp = 0,
                     log = FALSE, tol = 1e-6, max_iter = 10000L) {
  lp <- dcbblp(x, size, alpha, beta, tol, max_iter)
  # To keep numerical resolution, compute log(p) originally
  p <- base::log((lp[, 1L, drop = TRUE] - lp[, 2L, drop = TRUE]) / lp[, 3L, drop = TRUE]) +
    pcbbinom(q = x, size = size, alpha = alpha, beta = beta,
             lower.tail = TRUE, log.p = TRUE,
             tol = tol, max_iter = max_iter)
  if (log == FALSE) {
    p <- exp(p)
  }
  p
}

#' @aliases pcbbinom
#' @export
#' @rdname cbbinom
pcbbinom <- function(q, size, alpha = 1, beta = 1, ncp = 0,
                     lower.tail = TRUE, log.p = FALSE,
                     tol = 1e-6, max_iter = 10000L) {
  cpp_pcbbinom(q = as.numeric(q - ncp),
               size = as.numeric(size),
               alpha = as.numeric(alpha),
               beta = as.numeric(beta),
               lower_tail = as.logical(lower.tail[[1L]]),
               log_p = as.logical(log.p[[1L]]),
               tol = as.numeric(tol),
               max_iter = as.integer(max_iter))
}

#' @aliases qcbbinom
#' @export
#' @rdname cbbinom
qcbbinom <- function(p, size, alpha = 1, beta = 1, ncp = 0,
                     lower.tail = TRUE, log.p = FALSE,
                     p_tol = 1e-6, p_max_iter = 10000L,
                     root_tol = 1e-6, root_max_iter = 10000L) {
  cpp_qcbbinom(p = as.numeric(p),
               size = as.numeric(size),
               alpha = as.numeric(alpha),
               beta = as.numeric(beta),
               lower_tail = as.logical(lower.tail[[1L]]),
               log_p = as.logical(log.p[[1L]]),
               p_tol = as.numeric(p_tol),
               p_max_iter = as.integer(p_max_iter),
               root_tol = as.numeric(root_tol),
               root_max_iter = as.integer(root_max_iter)) + ncp
}

#' @aliases rcbbinom
#' @export
#' @rdname cbbinom
rcbbinom <- function(n, size, alpha = 1, beta = 1, ncp = 0,
                     p_tol = 1e-6, p_max_iter = 10000L,
                     root_tol = 1e-6, root_max_iter = 10000L) {
  cpp_rcbbinom(n = as.integer(n[[1L]]),
               size = as.numeric(size),
               alpha = as.numeric(alpha),
               beta = as.numeric(beta),
               p_tol = as.numeric(p_tol),
               p_max_iter = as.integer(p_max_iter),
               root_tol = as.numeric(root_tol),
               root_max_iter = as.integer(root_max_iter)) + ncp
}
