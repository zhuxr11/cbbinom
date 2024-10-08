library(testthat)
library(extraDistr)

# Parameters
test_size <- 10
test_alpha <- 2
test_beta <- 4
test_delta <- 1e-3
test_prec <- 20L
test_tol <- 1e-6
test_max_iter <- 10000L

# Inputs
test_x <- 5
test_val <- pcbbinom(q = test_x, size = test_size, alpha = test_alpha, beta = test_beta)

testthat::test_that(
  "dcbbinom",
  {
    testthat::expect_equal(
      dcbbinom(x = test_x, size = test_size, alpha = test_alpha, beta = test_beta),
      (pcbbinom(q = test_x + test_delta, size = test_size,
                alpha = test_alpha, beta = test_beta, log.p = FALSE) -
         pcbbinom(q = test_x - test_delta, size = test_size,
                  alpha = test_alpha, beta = test_beta, log.p = FALSE)) /
        (2 * test_delta)
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
    # Verify that p, root_tol and root_max_iter are not changed
    test_val_log <- test_val_log_orig <- log(test_val)
    test_tol_orig <- test_tol
    test_max_iter_orig <- test_max_iter
    testthat::expect_equal(
      cpp_qcbbinom(p = test_val_log, size = test_size,
                   alpha = test_alpha, beta = test_beta,
                   lower_tail = TRUE, log_p = TRUE, prec = test_prec,
                   tol = test_tol, max_iter = test_max_iter),
      test_x
    )
    testthat::expect_identical(test_val_log, test_val_log_orig)
    testthat::expect_identical(test_tol, test_tol_orig)
    testthat::expect_identical(test_max_iter, test_max_iter_orig)
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
