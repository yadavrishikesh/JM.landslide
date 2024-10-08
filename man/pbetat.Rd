% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other-supporting-functions.R
\name{pbetat}
\alias{pbetat}
\title{Distribution Function of Scaled Beta Distribution}
\usage{
pbetat(x, u.thr, shape1, shape2, log.p = TRUE)
}
\arguments{
\item{x}{A numeric vector of quantiles for which the CDF is to be evaluated.}

\item{u.thr}{A positive numeric value specifying the upper bound of the support of the Beta distribution.}

\item{shape1}{A positive numeric value representing the first shape parameter of the Beta distribution.}

\item{shape2}{A positive numeric value representing the second shape parameter of the Beta distribution.}

\item{log.p}{Logical; if \code{TRUE}, the log of the CDF is returned. Defaults to \code{TRUE}.}
}
\value{
A numeric vector of CDF values corresponding to the quantiles in \code{x}.
}
\description{
This function computes the cumulative distribution function (CDF) of a Beta distribution scaled to have its support in the interval \eqn{(0, u.thr)}.
}
\details{
The standard Beta distribution has support on the interval \eqn{(0, 1)}. This function rescales the distribution so that
its support lies in \eqn{(0, u.thr)}. The CDF is computed by applying a change of variables to the standard Beta
distribution function.

If \code{log.p} is \code{TRUE}, the log of the CDF is returned.
}
\examples{
# Evaluate the CDF of the scaled Beta distribution
x <- seq(0.01, 5, length.out = 100)
u.thr <- 5
shape1 <- 2
shape2 <- 3
cdf_values <- pbetat(x, u.thr, shape1, shape2, log.p = FALSE)
plot(x, cdf_values, type = "l", main = "Scaled Beta CDF", ylab = "CDF", xlab = "x")

}
\seealso{
\code{\link{pbeta}} for the standard Beta distribution CDF.
}
