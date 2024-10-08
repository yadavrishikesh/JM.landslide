% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other-supporting-functions.R
\name{rtgGPD}
\alias{rtgGPD}
\title{Simulate Random Numbers from a Mixture of Truncated Gamma and GPD Distribution}
\usage{
rtgGPD(n, u.thr, k, xi, sigma, sigma.GP, ind)
}
\arguments{
\item{n}{An integer specifying the total number of random samples to generate.}

\item{u.thr}{A numeric value specifying the threshold that separates the Gamma and GPD components.}

\item{k}{A positive numeric value representing the shape parameter of the truncated Gamma distribution.}

\item{xi}{A numeric value representing the shape parameter of the GPD.}

\item{sigma}{A numeric vector specifying the scale parameter for the truncated Gamma distribution.}

\item{sigma.GP}{A positive numeric value representing the scale parameter of the GPD.}

\item{ind}{A logical vector of length \code{n} indicating which samples should be generated from the GPD (\code{TRUE}) and which should be generated from the truncated Gamma distribution (\code{FALSE}).}
}
\value{
A numeric vector of length \code{n}, containing the random samples drawn from the mixture of truncated Gamma and GPD distributions.
}
\description{
This function generates random samples from a mixture of a truncated Gamma and a Generalized Pareto Distribution (GPD). The samples are drawn based on an indicator vector, where elements specify which distribution (Gamma or GPD) to sample from.
}
\details{
The function generates \code{n} random samples from a mixture distribution where:
\itemize{
\item Samples corresponding to \code{ind = TRUE} are drawn from the GPD with threshold \code{u.thr}, scale parameter \code{sigma.GP}, and shape parameter \code{xi}.
\item Samples corresponding to \code{ind = FALSE} are drawn from a truncated Gamma distribution with shape \code{k}, scale \code{sigma}, and truncated at \code{u.thr}.
}
}
\examples{
# Example usage of rtgGPD function
set.seed(123)
n <- 100
u.thr <- 5
k <- 2
xi <- 0.3
sigma <- 1
sigma.GP <- 0.8
ind <- runif(n) > 0.7  # 70\% Gamma, 30\% GPD
samples <- rtgGPD(n, u.thr, k, xi, sigma, sigma.GP, ind)
hist(samples, breaks = 30, main = "Samples from Mixture of Truncated Gamma and GPD")
}
\seealso{
\code{\link[evd]{rgpd}} for sampling from the Generalized Pareto Distribution,
\code{\link{rgammat}} for sampling from the truncated Gamma distribution.
}
