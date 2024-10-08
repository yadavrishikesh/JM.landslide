% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_cov_coeff.R
\name{beta2_sim_thr_ind}
\alias{beta2_sim_thr_ind}
\title{Simulate the Coefficients \code{beta2} for Threshold Indicator Models}
\usage{
beta2_sim_thr_ind(
  mu,
  intercept2,
  beta,
  kappa_mu,
  W2,
  Z2,
  Z2.crossprd,
  hyper_fixed
)
}
\arguments{
\item{mu}{A numeric vector representing the latent mean parameter of the model.}

\item{intercept2}{A numeric value representing the intercept term in the model.}

\item{beta}{A numeric vector of initial values for \code{beta2}.}

\item{kappa_mu}{A numeric value representing the precision parameter for the latent process \code{mu}.}

\item{W2}{A numeric vector indicating the latent (iCAR) random effects.}

\item{Z2}{A numeric matrix of covariates corresponding to the latent process \code{mu}.}

\item{Z2.crossprd}{A pre-computed cross-product matrix of \code{Z2} \(\code{t(Z2) Z2}\) for computational efficiency.}

\item{hyper_fixed}{A list containing hyperparameters, specifically the precision of the normal prior for \code{beta2}.}
}
\value{
A numeric vector of proposed values for \code{beta2}.
}
\description{
This function simulates the regression coefficients \code{beta2} for threshold indicator models using a normal full conditional distribution. It performs a precision-based update by considering both prior and data contributions.
}
\examples{
# Sample latent variables and covariates
mu <- rnorm(100)
intercept2 <- 0.5
beta <- rnorm(10)
kappa_mu <- 1.5
W2 <- rnorm(100)
Z2 <- matrix(rnorm(1000), ncol = 10)
Z2.crossprd <- crossprod(Z2) # pre-computed cross-product
hyper_fixed <- list(beta2 = 0.1)

# Simulate beta2 for threshold indicator model
beta2_sim <- beta2_sim_thr_ind(mu, intercept2, beta, kappa_mu, W2, Z2, Z2.crossprd, hyper_fixed)
}
