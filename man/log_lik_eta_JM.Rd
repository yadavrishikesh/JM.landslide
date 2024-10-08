% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_lik_latent.R
\name{log_lik_eta_JM}
\alias{log_lik_eta_JM}
\title{Log-Likelihood for at the count data level in Joint Model}
\usage{
log_lik_eta_JM(Y, eta)
}
\arguments{
\item{Y}{A numeric vector representing the response variable in the model.}

\item{eta}{A numeric vector representing the linear predictor or latent variable associated with the model.}
}
\value{
A numeric vector representing the log-likelihood for each value of \code{eta}.
}
\description{
This function calculates the log-likelihood for the count data.
}
\examples{
\dontshow{if (FALSE) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Example usage (not meant to be run directly):
Y <- rpois(10, lambda = 3)
eta <- rnorm(10)
loglik_values <- log_lik_eta_JM(Y, eta)
\dontshow{\}) # examplesIf}
}
