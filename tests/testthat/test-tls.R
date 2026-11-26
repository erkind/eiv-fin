test_that("tls recovers an exact linear relationship", {
  x <- seq(1, 20, length.out = 25)
  y <- 2 + 0.5 * x
  fit <- tls(x, y, B = 0)
  expect_equal(unname(coef(fit)), c(2, 0.5), tolerance = 1e-8)
})

test_that("tls is less attenuated than OLS under noise in x", {
  set.seed(42)
  n <- 500
  xstar <- rnorm(n)
  x <- xstar + rnorm(n, sd = 0.5)
  y <- 1 + 1 * xstar + rnorm(n, sd = 0.5)
  fit <- tls(x, y, B = 0)
  b_ols <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
  expect_lt(abs(unname(coef(fit)[2]) - 1), abs(b_ols - 1))
  expect_equal(unname(coef(fit)[2]), 1, tolerance = 0.15)
})

test_that("formula interface matches default interface", {
  set.seed(7)
  n <- 100
  x <- rnorm(n)
  y <- 1 + 0.7 * x + rnorm(n)
  d <- data.frame(y = y, x = x)
  f1 <- tls(y ~ x, data = d, B = 0)
  f2 <- tls(d$x, d$y, B = 0)
  expect_equal(coef(f1), coef(f2))
  expect_equal(names(coef(f1)), c("(Intercept)", "x"))
})

test_that("bootstrap covariance is valid and supports confint", {
  set.seed(3)
  x <- rnorm(60)
  y <- 0.5 + 0.9 * x + rnorm(60)
  fit <- tls(x, y, B = 100, seed = 3)
  V <- vcov(fit)
  expect_true(all(diag(V) > 0))
  expect_true(isSymmetric(V))
  ci <- confint(fit)
  expect_equal(dim(ci), c(2L, 2L))
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("invalid inputs are rejected", {
  expect_error(tls(1:5, 1:6, B = 0), "incompatible lengths")
  expect_error(tls(rep(1, 10), 1:10, B = 0), "no variation")
  expect_error(vcov(tls(1:20, 2 + 0.5 * (1:20), B = 0)), "no variance estimate")
})
