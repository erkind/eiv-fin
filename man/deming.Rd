\name{deming}
\alias{deming}
\alias{deming.default}
\alias{deming.formula}
\title{Deming Regression with Known Error Variance Ratio}
\description{
Fits the errors-in-variables model \eqn{y = \alpha + \beta x^* + e_y},
\eqn{x = x^* + e_x}, where \eqn{e_y} and \eqn{e_x} are independent and the
ratio \eqn{\lambda = \sigma_{e_y}^2 / \sigma_{e_x}^2} of the error variance
of \eqn{y} to the error variance of \eqn{x} is known. The estimator is the
maximum likelihood estimator of the linear errors-in-variables model under
normal errors (Deming 1943; Cornbleet and Gochman 1979).
}
\usage{
deming(x, y, lambda = 1, se = c("jackknife", "none"))

\method{deming}{formula}(formula, data = NULL, lambda = 1,
        se = c("jackknife", "none"), ...)
}
\arguments{
\item{x}{numeric vector of regressor observations (measured with error).}
\item{y}{numeric vector of response observations (measured with error).}
\item{lambda}{known ratio of the error variance of \eqn{y} to the error
variance of \eqn{x}; positive and finite. Defaults to 1 (orthogonal
regression).}
\item{se}{character; \code{"jackknife"} (default) computes delete-one
jackknife standard errors, \code{"none"} skips variance estimation.}
\item{formula}{an object of class \code{\link{formula}} of the form
\code{y ~ x} with exactly one regressor.}
\item{data}{an optional data frame containing the variables in
\code{formula}.}
\item{...}{further arguments passed to or from other methods.}
}
\details{
The slope estimator is
\eqn{\hat\beta = (s_{yy} - \lambda s_{xx} + \sqrt{(s_{yy} - \lambda s_{xx})^2 + 4\lambda s_{xy}^2}) / (2 s_{xy})}.
With \eqn{\lambda = 1} this is orthogonal regression. As
\eqn{\lambda \to \infty} (negligible error in \eqn{x}) the estimator
approaches ordinary least squares of \eqn{y} on \eqn{x}; as
\eqn{\lambda \to 0} (negligible error in \eqn{y}) it approaches the
inverse regression of \eqn{x} on \eqn{y}.

Standard errors are computed by delete-one jackknife over the
observations, following the approach used for Deming regression in
method-comparison studies.
}
\value{
An object of class \code{\link{eivfit}} with components
\code{coefficients}, \code{vcov} (when \code{se = "jackknife"}),
\code{lambda}, \code{n}, \code{method} and \code{call}.
}
\references{
Deming, W. E. (1943). \emph{Statistical Adjustment of Data}. Wiley, New York.

Cornbleet, P. J. and Gochman, N. (1979). Incorrect least-squares regression
coefficients in method-comparison analysis. \emph{Clinical Chemistry},
25(3), 432-438.

Fuller, W. A. (1987). \emph{Measurement Error Models}. Wiley, New York.
}
\seealso{
\code{\link{tls}}, \code{\link{eivfit}}
}
\examples{
set.seed(1)
n <- 300
xstar <- rnorm(n)
x <- xstar + rnorm(n, sd = 0.4)
y <- 1 + 0.8 * xstar + rnorm(n, sd = 0.2)
lambda <- 0.2^2 / 0.4^2
fit <- deming(x, y, lambda = lambda)
summary(fit)
}
\keyword{regression}
