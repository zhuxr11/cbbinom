library(testthat)
library(extraDistr)

# Parameters
test_size <- 10
test_alpha <- 2
test_beta <- 4
test_delta <- 1e-6

# Inputs
test_x <- 5
test_val <- pcbbinom(q = test_x, size = test_size, alpha = test_alpha, beta = test_beta)

testthat::test_that(
  "dcbbinom",
  {
    testthat::expect_equal(
      dcbbinom(x = test_x, size = test_size, alpha = test_alpha, beta = test_beta),
      (pcbbinom(q = test_x + test_delta, size = test_size,
                alpha = test_alpha, beta = test_beta, log.p = TRUE) -
         pcbbinom(q = test_x - test_delta, size = test_size,
                  alpha = test_alpha, beta = test_beta, log.p = TRUE)) /
        (2 * test_delta) * test_val
    )
  }
)

testthat::test_that(
  "pcbbinom",
  {
    testthat::expect_equal(
      test_val,
      extraDistr::pbbinom(q = test_x - 1, size = test_size,
                          alpha = test_alpha, beta = test_beta)
    )
  }
)

testthat::test_that(
  "qcbbinom",
  {
    testthat::expect_equal(
      qcbbinom(p = test_val, size = test_size, alpha = test_alpha, beta = test_beta),
      test_x
    )
  }
)

testthat::test_that(
  "rcbbinom",
  {
    testthat::expect_error(
      rcbbinom(n = 10L, size = test_size, alpha = test_alpha, beta = test_beta),
      NA
    )
  }
)

testthat::test_that(
  "gen_hypergeo",
  {
    testthat::expect_equal(
      gen_hypergeo(U = c(1 - test_x,
                         test_size + 1 - test_x,
                         test_size + 1 - test_x + test_beta),
                   L = c(test_size + test_alpha - test_x,
                         test_size + 1 - test_x + test_alpha + test_beta),
                   x = 1,
                   tol = 1e-6,
                   max_iter = 10000L,
                   check_mode = TRUE,
                   log = FALSE),
      101/5460
    )
  }
)
