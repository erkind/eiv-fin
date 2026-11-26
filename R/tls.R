# Total least squares (orthogonal) regression

# Internal: fit the TLS line/hyperplane on given data (no error handling).
# Returns coefficients (intercept first) and the estimated common error
# variance under the equal-variance assumption.
tls_fit <- function(X, y) {
  n <- nrow(X)
  k <- ncol(X)
  if (n <= k + 1L) {
    stop("need at least k + 2 observations for k regressors")
  }
  y <- as.numeric(y)
  xm <- colMeans(X)
  ym <- mean(y)
  D <- cbind(sweep(X, 2L, xm), y - ym)
  sv <- svd(D)
  v <- sv$v[, k + 1L]
  if (abs(v[k + 1L]) < sqrt(.Machine$double.eps)) {
    stop("degenerate TLS solution: response direction has no variation")
  }
  beta <- -v[seq_len(k)] / v[k + 1L]
  coef <- c("(Intercept)" = ym - sum(beta * xm),
            stats::setNames(as.numeric(beta), colnames(X)))
  list(coefficients = coef,
       sigma2 = sv$d[k + 1L]^2 / (n - k - 1L),
       n = n,
       k = k)
}

# Internal: single TLS fit on a matrix of regressors.
tls_matrix <- function(X, y) {
  X <- as.matrix(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }
  if (any(!is.finite(X)) || any(!is.finite(y))) {
    stop("tls() requires fully observed finite x and y")
  }
  tls_fit(X, y)
}

#' Total Least Squares (Orthogonal) Regression
#'
#' Fits the errors-in-variables model \eqn{y = \alpha + \beta x^* + e_y},
#' \eqn{x = x^* + e_x}, where both \eqn{x} and \eqn{y} are observed with
#' error. The estimator minimises the sum of squared perpendicular distances
#' to the fitted hyperplane and is computed from the singular value
#' decomposition of the centered data matrix (Golub and Van Loan 1980;
#' Fuller 1987, ch. 1). It assumes the error variances of all variables are
#' equal (the orthogonal regression case, error ratio 1).
#'
#' For unequal error variances use \code{\link{deming}} with the known ratio.
#'
#' @param x a numeric vector or numeric matrix of regressors observed with
#'   error (one column per regressor).
#' @param y numeric response vector observed with error.
#' @param B number of bootstrap replicates (resampling observations) used to
#'   estimate the covariance matrix of the coefficients. Use \code{B = 0} to
#'   skip the bootstrap.
#' @param seed optional integer seed for the bootstrap; passed to
#'   \code{\link{set.seed}} for reproducibility.
#' @param formula an object of class \code{\link{formula}} of the form
#'   \code{y ~ x1 + x2} (no intercept term is needed; one is always fitted).
#' @param data an optional data frame containing the variables in
#'   \code{formula}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \code{\link{eivfit}}.
#'
#' @references
#' Golub, G. H. and Van Loan, C. F. (1980). An analysis of the total least
#' squares problem. \emph{SIAM Journal on Numerical Analysis}, 17(6), 883-893.
#'
#' Fuller, W. A. (1987). \emph{Measurement Error Models}. Wiley, New York.
#'
#' @seealso \code{\link{deming}}, \code{\link{eivfit}}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' xstar <- rnorm(n)
#' x <- xstar + rnorm(n, sd = 0.3)
#' y <- 1 + 0.8 * xstar + rnorm(n, sd = 0.3)
#' fit <- tls(x, y, seed = 1)
#' fit
#' summary(fit)
#'
#' @export
tls <- function(x, y, B = 999, seed = NULL, ...) {
  UseMethod("tls")
}

#' @rdname tls
#' @export
tls.default <- function(x, y, B = 999, seed = NULL) {
  if (!is.numeric(B) || length(B) != 1L || B < 0) {
    stop("B must be a non-negative integer")
  }
  B <- as.integer(B)
  X <- as.matrix(x)
  if (ncol(X) == 0L) {
    stop("x must contain at least one column")
  }
  y <- as.numeric(y)
  if (length(y) != nrow(X)) {
    stop("x and y have incompatible lengths")
  }
  fit <- tls_matrix(X, y)
  vc <- NULL
  if (B > 0L) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    n <- nrow(X)
    draws <- matrix(NA_real_, nrow = B, ncol = length(fit$coefficients),
                    dimnames = list(NULL, names(fit$coefficients)))
    for (b in seq_len(B)) {
      idx <- sample.int(n, replace = TRUE)
      draw <- tryCatch(
        tls_matrix(X[idx, , drop = FALSE], y[idx])$coefficients,
        error = function(e) rep(NA_real_, length(fit$coefficients))
      )
      draws[b, ] <- draw
    }
    draws <- draws[stats::complete.cases(draws), , drop = FALSE]
    if (nrow(draws) >= 2L) {
      vc <- stats::cov(draws)
    }
  }
  out <- list(coefficients = fit$coefficients,
              vcov = vc,
              sigma2 = fit$sigma2,
              lambda = NULL,
              n = fit$n,
              B = B,
              method = "tls",
              call = match.call())
  class(out) <- "eivfit"
  out
}

#' @rdname tls
#' @export
tls.formula <- function(formula, data = NULL, B = 999, seed = NULL, ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$B <- mf$seed <- mf$... <- NULL
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  y <- stats::model.response(mf, "numeric")
  tt <- stats::delete.response(stats::terms(mf))
  attr(tt, "intercept") <- 0
  X <- stats::model.matrix(tt, mf)
  if (ncol(X) == 0L) {
    stop("formula must specify at least one regressor")
  }
  res <- tls.default(X, y, B = B, seed = seed)
  res$call <- cl
  res
}
