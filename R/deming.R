# Deming regression with known error variance ratio

# Internal: Deming fit on raw vectors. lambda is the ratio
# sigma_e^2 / sigma_eta^2 of the error variance of y to the error variance
# of x. Moments are computed as sums of squares; the estimator is invariant
# to this choice.
deming_fit <- function(x, y, lambda) {
  mx <- mean(x)
  my <- mean(y)
  sxx <- sum((x - mx)^2)
  syy <- sum((y - my)^2)
  sxy <- sum((x - mx) * (y - my))
  if (sxy == 0) {
    stop("x and y are uncorrelated; the Deming slope is undefined")
  }
  disc <- (syy - lambda * sxx)^2 + 4 * lambda * sxy^2
  slope <- ((syy - lambda * sxx) + sqrt(disc)) / (2 * sxy)
  c("(Intercept)" = my - slope * mx, "x" = slope)
}

#' Deming Regression with Known Error Variance Ratio
#'
#' Fits the errors-in-variables model \eqn{y = \alpha + \beta x^* + e_y},
#' \eqn{x = x^* + e_x}, where \eqn{e_y} and \eqn{e_x} are independent and the
#' ratio \eqn{\lambda = \sigma_{e_y}^2 / \sigma_{e_x}^2} of the error variance
#' of \eqn{y} to the error variance of \eqn{x} is known (Deming 1943;
#' Cornbleet and Gochman 1979). The estimator is the maximum likelihood
#' estimator of the linear errors-in-variables model under normal errors.
#'
#' With \eqn{\lambda = 1} this is orthogonal regression. As
#' \eqn{\lambda \to \infty} (negligible error in \eqn{x}) the estimator
#' approaches ordinary least squares of \eqn{y} on \eqn{x}; as
#' \eqn{\lambda \to 0} (negligible error in \eqn{y}) it approaches the
#' inverse regression of \eqn{x} on \eqn{y}.
#'
#' Standard errors are computed by delete-one jackknife over the
#' observations, following the approach used for Deming regression in
#' method-comparison studies.
#'
#' @param x numeric vector of regressor observations (measured with error).
#' @param y numeric vector of response observations (measured with error).
#' @param lambda known ratio of the error variance of \eqn{y} to the error
#'   variance of \eqn{x}; positive and finite. Defaults to 1 (orthogonal
#'   regression).
#' @param se character; \code{"jackknife"} (default) computes delete-one
#'   jackknife standard errors, \code{"none"} skips variance estimation.
#' @param formula an object of class \code{\link{formula}} of the form
#'   \code{y ~ x} with exactly one regressor.
#' @param data an optional data frame containing the variables in
#'   \code{formula}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \code{\link{eivfit}}.
#'
#' @references
#' Deming, W. E. (1943). \emph{Statistical Adjustment of Data}. Wiley, New York.
#'
#' Cornbleet, P. J. and Gochman, N. (1979). Incorrect least-squares regression
#' coefficients in method-comparison analysis. \emph{Clinical Chemistry},
#' 25(3), 432-438.
#'
#' Fuller, W. A. (1987). \emph{Measurement Error Models}. Wiley, New York.
#'
#' @seealso \code{\link{tls}}, \code{\link{eivfit}}
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' xstar <- rnorm(n)
#' x <- xstar + rnorm(n, sd = 0.4)
#' y <- 1 + 0.8 * xstar + rnorm(n, sd = 0.2)
#' lambda <- 0.2^2 / 0.4^2
#' fit <- deming(x, y, lambda = lambda)
#' summary(fit)
#'
#' @export
deming <- function(x, y, lambda = 1, se = c("jackknife", "none")) {
  UseMethod("deming")
}

#' @rdname deming
#' @export
deming.default <- function(x, y, lambda = 1, se = c("jackknife", "none")) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) {
    stop("x and y have incompatible lengths")
  }
  if (length(x) < 3L) {
    stop("need at least 3 observations")
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("deming() requires fully observed finite x and y")
  }
  if (length(unique(x)) < 2L) {
    stop("x must have at least two distinct values")
  }
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda <= 0) {
    stop("lambda must be a single positive finite number")
  }
  se <- match.arg(se)
  fit <- deming_fit(x, y, lambda)
  vc <- NULL
  if (se == "jackknife") {
    n <- length(x)
    vals <- matrix(NA_real_, nrow = n, ncol = 2L)
    keep <- logical(n)
    for (i in seq_len(n)) {
      r <- tryCatch(deming_fit(x[-i], y[-i], lambda),
                    error = function(e) NULL)
      if (!is.null(r)) {
        vals[i, ] <- r
        keep[i] <- TRUE
      }
    }
    vals <- vals[keep, , drop = FALSE]
    m <- sum(keep)
    if (m >= 2L) {
      jse <- sqrt((m - 1) / m *
                    colSums(sweep(vals, 2L, colMeans(vals))^2))
      vc <- diag(c(jse[1L]^2, jse[2L]^2))
      dimnames(vc) <- list(names(fit), names(fit))
    }
  }
  out <- list(coefficients = fit,
              vcov = vc,
              sigma2 = NULL,
              lambda = lambda,
              n = length(x),
              B = NULL,
              method = "deming",
              call = match.call())
  class(out) <- "eivfit"
  out
}

#' @rdname deming
#' @export
deming.formula <- function(formula, data = NULL, lambda = 1,
                           se = c("jackknife", "none"), ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$lambda <- mf$se <- mf$... <- NULL
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  y <- stats::model.response(mf, "numeric")
  tt <- stats::delete.response(stats::terms(mf))
  attr(tt, "intercept") <- 0
  X <- stats::model.matrix(tt, mf)
  if (ncol(X) != 1L) {
    stop("deming() supports a single regressor")
  }
  xname <- colnames(X)
  res <- deming.default(X[, 1L], y, lambda = lambda, se = se)
  names(res$coefficients)[2L] <- xname
  if (!is.null(res$vcov)) {
    dimnames(res$vcov) <- list(names(res$coefficients),
                               names(res$coefficients))
  }
  res$call <- cl
  res
}
