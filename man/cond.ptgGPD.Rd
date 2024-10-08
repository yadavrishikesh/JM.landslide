% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other-supporting-functions.R
\name{cond.ptgGPD}
\alias{cond.ptgGPD}
\title{Conditional Cumulative Distribution Function of Mixture of Truncated Gamma and GPD}
\usage{
cond.ptgGPD(x, mu, u.thr, k, xi, sigma, sigma.GP, thr.prob)
}
\arguments{
\item{x}{A numeric value specifying the quantile at which the CDF is evaluated.}

\item{mu}{A numeric vector representing the location parameter values for the conditional CDF.}

\item{u.thr}{A numeric value specifying the threshold. Values of \code{x} greater than \code{u.thr} are treated as GPD, while those less than \code{u.thr} follow the truncated Gamma distribution.}

\item{k}{A positive numeric value representing the shape parameter of the truncated Gamma distribution.}

\item{xi}{A numeric value representing the shape parameter of the GPD.}

\item{sigma}{A numeric vector specifying the scale parameter of the truncated Gamma distribution.}

\item{sigma.GP}{A positive numeric value representing the scale parameter of the GPD.}

\item{thr.prob}{A numeric value in (0, 1) specifying the probability mass of the GPD component in the mixture.}
}
\value{
A numeric vector of the same length as \code{mu}, containing the conditional CDF values of the mixture distribution at the specified quantile \code{x}.
}
\description{
This function calculates the conditional cumulative distribution function (CDF) for a mixture of truncated Gamma and Generalized Pareto Distributions (GPD), given the parameter vector \code{mu}.
The CDF is computed based on whether the input value \code{x} is greater than or less than the threshold \code{u.thr}.
}
\details{
The mixture distribution is modeled conditionally based on the threshold \code{u.thr}. Specifically:
\itemize{
\item If \code{x > u.thr}, the conditional CDF is computed using the GPD with a certain probability mass, while the truncated Gamma distribution contributes when \code{exp(mu)} is less than \code{u.thr}.
\item If \code{x <= u.thr}, the conditional CDF is calculated using only the truncated Gamma distribution.
}
}
\examples{
# Example usage of cond.ptgGPD function
set.seed(123)
x <- 5
mu <- rnorm(10)
u.thr <- 4
k <- 2
xi <- 0.5
sigma <- rep(1, length(mu))
sigma.GP <- 0.8
thr.prob <- 0.3
cdf_values <- cond.ptgGPD(x, mu, u.thr, k, xi, sigma, sigma.GP, thr.prob)
print(cdf_values)
}
\seealso{
\code{\link[evd]{pgpd}} for the CDF of the standard GPD,
\code{\link{pgammat}} for the CDF of the truncated Gamma distribution.
}
