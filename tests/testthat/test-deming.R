test_that("deming approaches OLS when x is measured exactly (large lambda)", {
  set.seed(11)
  n <- 400
  x <- rnorm(n)
  y <- 1 + 0.8 * x + rnorm(n, sd = 0.5)
  fit <- deming(x, y, lambda = 1e6, se = "none")
  b_ols <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
  expect_equal(unname(coef(fit)[2]), b_ols, tolerance = 1e-3)
})

test_that("deming approaches inverse regression when y is exact (small lambda)", {
  set.seed(12)
  n <- 400
  xstar <- rnorm(n)
  y <- 1 + 0.8 * xstar
  x <- xstar + rnorm(n, sd = 0.5)
  fit <- deming(x, y, lambda = 1e-6, se = "none")
  w <- sum((x - mean(x)) * (y - mean(y))) / sum((y - mean(y))^2)
  expect_equal(unname(coef(fit)[2]), 1 / w, tolerance = 1e-3)
})

test_that("deming with known error ratio recovers the true slope", {
  set.seed(21)
  n <- 2000
  xstar <- rnorm(n)
  sd_x_err <- 0.4
  sd_y_err <- 0.2
  x <- xstar + rnorm(n, sd = sd_x_err)
  y <- 1 + 0.9 * xstar + rnorm(n, sd = sd_y_err)
  fit <- deming(x, y, lambda = sd_y_err^2 / sd_x_err^2, se = "jackknife")
  expect_equal(unname(coef(fit)[2]), 0.9, tolerance = 0.05)
  expect_equal(unname(coef(fit)[1]), 1, tolerance = 0.05)
  expect_true(all(diag(vcov(fit)) > 0))
  ci <- confint(fit)
  expect_lt(ci[2, 1], 0.9)
  expect_gt(ci[2, 2], 0.9)
})

test_that("formula interface matches default interface", {
  set.seed(31)
  n <- 120
  xstar <- rnorm(n)
  d <- data.frame(y = 1 + 0.6 * xstar + rnorm(n, sd = 0.2),
                  x = xstar + rnorm(n, sd = 0.3))
  f1 <- deming(y ~ x, data = d, lambda = 0.5)
  f2 <- deming(d$x, d$y, lambda = 0.5)
  expect_equal(coef(f1), coef(f2))
  expect_equal(names(coef(f1)), c("(Intercept)", "x"))
})

test_that("invalid inputs are rejected", {
  expect_error(deming(1:10, 1:11, lambda = 1, se = "none"),
               "incompatible lengths")
  expect_error(deming(1:10, 1:10, lambda = -1, se = "none"), "lambda")
  expect_error(deming(rep(1, 10), 1:10, lambda = 1, se = "none"),
               "distinct values")
})
