library(testthat)
library(extraDistr)

test_x <- 5
test_val <- pcbbinom(q = test_x, size = 10, alpha = 2, beta = 4)

testthat::test_that(
  "dcbbinom",
  {
    testthat::expect_equal(
      dcbbinom(x = test_x, size = 10, alpha = 2, beta = 4),
      (pcbbinom(q = test_x + 1e-9, size = 10, alpha = 2, beta = 4) - test_val) / 1e-9,
      tolerance = 1e-6
    )
  }
)

testthat::test_that(
  "pcbbinom",
  {
    testthat::expect_equal(
      test_val,
      extraDistr::pbbinom(q = test_x - 1, size = 10, alpha = 2, beta = 4)
    )
  }
)

testthat::test_that(
  "qcbbinom",
  {
    testthat::expect_equal(
      qcbbinom(p = test_val, size = 10, alpha = 2, beta = 4),
      test_x
    )
  }
)

testthat::test_that(
  "rcbbinom",
  {
    testthat::expect_error(
      rcbbinom(n = 10L, size = 10, alpha = 2, beta = 4),
      NA
    )
  }
)

testthat::test_that(
  "gen_hypergeo",
  {
    testthat::expect_equal(
      gen_hypergeo(U = c(1 - test_x, 10 + 1 - test_x, 10 + 1 - test_x + 4),
                   L = c(10 + 2 - test_x, 10 + 1 - test_x + 2 + 4),
                   x = 1,
                   tol = 1e-6,
                   max_iter = 10000L,
                   check_mode = TRUE,
                   log = FALSE),
      101/5460
    )
  }
)
