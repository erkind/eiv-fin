# Methods for fitted errors-in-variables models

#' Fitted Errors-in-Variables Model Objects
#'
#' Objects returned by \code{\link{tls}} and \code{\link{deming}}. The class
#' provides conventional extractor methods: \code{coef}, \code{vcov},
#' \code{confint}, \code{print} and \code{summary}.
#'
#' @param x an object of class \code{eivfit}.
#' @param object an object of class \code{eivfit}.
#' @param parm a specification of which parameters are to be given
#'   confidence intervals, either a vector of numbers or of names; defaults
#'   to all parameters.
#' @param level the confidence level required.
#' @param digits the number of significant digits to print.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{coef} returns the named vector of coefficient estimates.
#'   \code{vcov} returns the estimated covariance matrix of the
#'   coefficients. \code{confint} returns a matrix of confidence intervals.
#'   \code{summary} returns an object of class \code{summary.eivfit}.
#'
#' @seealso \code{\link{tls}}, \code{\link{deming}}
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' y <- 0.5 + 0.9 * x + rnorm(100, sd = 0.2)
#' fit <- tls(x, y, seed = 1)
#' coef(fit)
#' confint(fit)
#'
#' @name eivfit
NULL

#' @rdname eivfit
#' @export
coef.eivfit <- function(object, ...) {
  object$coefficients
}

#' @rdname eivfit
#' @export
vcov.eivfit <- function(object, ...) {
  if (is.null(object$vcov)) {
    stop("no variance estimate available; refit deming() with ",
         'se = "jackknife" or tls() with B > 0')
  }
  object$vcov
}

#' @rdname eivfit
#' @export
confint.eivfit <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  V <- vcov.eivfit(object)
  se <- sqrt(diag(V))
  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  ci <- cbind(cf - z * se, cf + z * se)
  dimnames(ci) <- list(names(cf),
                       paste0(format(100 * c(a, 1 - a), trim = TRUE), "%"))
  if (missing(parm)) {
    parm <- seq_along(cf)
  } else if (is.character(parm)) {
    parm <- match(parm, names(cf))
  }
  ci[parm, , drop = FALSE]
}

#' @rdname eivfit
#' @export
summary.eivfit <- function(object, ...) {
  cf <- coef(object)
  if (!is.null(object$vcov)) {
    se <- sqrt(diag(object$vcov))
    zval <- cf / se
    pval <- 2 * stats::pnorm(-abs(zval))
    tab <- cbind(Estimate = cf,
                 `Std. Error` = se,
                 `z value` = zval,
                 `Pr(>|z|)` = pval)
  } else {
    tab <- cbind(Estimate = cf)
  }
  out <- list(coefficients = tab,
              method = object$method,
              n = object$n,
              lambda = object$lambda,
              B = object$B,
              call = object$call)
  class(out) <- "summary.eivfit"
  out
}

#' @rdname eivfit
#' @export
print.summary.eivfit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  label <- switch(x$method,
                  tls = "Total least squares (orthogonal) regression",
                  deming = "Deming regression")
  cat(label, "\n")
  if (!is.null(x$lambda)) {
    cat(sprintf("Error variance ratio lambda = %g\n", x$lambda))
  }
  if (!is.null(x$B) && x$B > 0) {
    cat(sprintf("Bootstrap covariance based on B = %d replicates\n", x$B))
  }
  cat(sprintf("Observations: %d\n\n", x$n))
  stats::printCoefmat(x$coefficients, digits = digits)
  invisible(x)
}

#' @rdname eivfit
#' @export
print.eivfit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}
