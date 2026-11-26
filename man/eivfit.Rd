\name{eivfit}
\alias{eivfit}
\alias{coef.eivfit}
\alias{vcov.eivfit}
\alias{confint.eivfit}
\alias{summary.eivfit}
\alias{print.eivfit}
\alias{print.summary.eivfit}
\title{Fitted Errors-in-Variables Model Objects}
\description{
Objects returned by \code{\link{tls}} and \code{\link{deming}}. The class
provides conventional extractor methods: \code{coef}, \code{vcov},
\code{confint}, \code{print} and \code{summary}.
}
\usage{
\method{coef}{eivfit}(object, ...)

\method{vcov}{eivfit}(object, ...)

\method{confint}{eivfit}(object, parm, level = 0.95, ...)

\method{summary}{eivfit}(object, ...)

\method{print}{eivfit}(x, digits = max(3L, getOption("digits") - 3L), ...)

\method{print}{summary.eivfit}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{object}{an object of class \code{eivfit}.}
\item{parm}{a specification of which parameters are to be given confidence
intervals, either a vector of numbers or of names; defaults to all
parameters.}
\item{level}{the confidence level required.}
\item{x}{an object of class \code{eivfit} or \code{summary.eivfit}.}
\item{digits}{the number of significant digits to print.}
\item{...}{further arguments passed to or from other methods.}
}
\value{
\code{coef} returns the named vector of coefficient estimates. \code{vcov}
returns the estimated covariance matrix of the coefficients (an error is
signalled if no variance estimate was requested at fit time).
\code{confint} returns a matrix with columns giving lower and upper
confidence limits for each parameter, based on a normal approximation.
\code{summary} returns an object of class \code{summary.eivfit} containing
a coefficient table.
}
\seealso{
\code{\link{tls}}, \code{\link{deming}}
}
\examples{
set.seed(1)
x <- rnorm(100)
y <- 0.5 + 0.9 * x + rnorm(100, sd = 0.2)
fit <- tls(x, y, seed = 1)
coef(fit)
confint(fit)
}
\keyword{regression}
