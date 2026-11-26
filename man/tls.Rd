\name{tls}
\alias{tls}
\alias{tls.default}
\alias{tls.formula}
\title{Total Least Squares (Orthogonal) Regression}
\description{
Fits the errors-in-variables model \eqn{y = \alpha + \beta x^* + e_y},
\eqn{x = x^* + e_x}, where both \eqn{x} and \eqn{y} are observed with
error, by minimising the sum of squared perpendicular distances to the
fitted hyperplane (total least squares, computed from the singular value
decomposition of the centered data matrix). Assumes equal error variances
across variables; for unequal error variances use
\code{\link{deming}}.
}
\usage{
tls(x, y, B = 999, seed = NULL)

\method{tls}{formula}(formula, data = NULL, B = 999, seed = NULL, ...)
}
\arguments{
\item{x}{a numeric vector or numeric matrix of regressors observed with
error (one column per regressor).}
\item{y}{numeric response vector observed with error.}
\item{B}{number of bootstrap replicates (resampling observations) used to
estimate the covariance matrix of the coefficients. Use \code{B = 0} to
skip the bootstrap.}
\item{seed}{optional integer seed for the bootstrap; passed to
\code{\link{set.seed}} for reproducibility.}
\item{formula}{an object of class \code{\link{formula}} of the form
\code{y ~ x1 + x2} (no intercept term is needed; one is always fitted).}
\item{data}{an optional data frame containing the variables in
\code{formula}.}
\item{...}{further arguments passed to or from other methods.}
}
\details{
The estimator is computed via the singular value decomposition of the
centered matrix \eqn{[X, y]}; the solution direction is the right singular
vector associated with the smallest singular value (Golub and Van Loan
1980; Fuller 1987, ch. 1). The sampling covariance matrix of the
coefficients is estimated by the pairs bootstrap over observations.
}
\value{
An object of class \code{\link{eivfit}} with components
\code{coefficients}, \code{vcov} (when \code{B > 0}), \code{sigma2} (the
estimated common error variance), \code{n}, \code{B}, \code{method} and
\code{call}.
}
\references{
Golub, G. H. and Van Loan, C. F. (1980). An analysis of the total least
squares problem. \emph{SIAM Journal on Numerical Analysis}, 17(6), 883-893.

Fuller, W. A. (1987). \emph{Measurement Error Models}. Wiley, New York.
}
\seealso{
\code{\link{deming}}, \code{\link{eivfit}}
}
\examples{
set.seed(1)
n <- 200
xstar <- rnorm(n)
x <- xstar + rnorm(n, sd = 0.3)
y <- 1 + 0.8 * xstar + rnorm(n, sd = 0.3)
fit <- tls(x, y, seed = 1)
fit
summary(fit)
}
\keyword{regression}
