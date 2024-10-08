% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_hyperparameters.R
\name{lb_ub_GP}
\alias{lb_ub_GP}
\title{Lower and Upper Bounds for Generalized Pareto Distribution (GPD) Parameters}
\usage{
lb_ub_GP(mark_dist)
}
\arguments{
\item{mark_dist}{A character string indicating the type of GPD distribution. Supported values are:
\itemize{
\item \code{"bGPD"}: Beta-Generalized Pareto Distribution.
\item \code{"tgGPD"}: Truncated Generalized Pareto Distribution.
}}
}
\value{
A list containing two elements:
\itemize{
\item \code{lb}: A numeric vector of the lower bounds for the parameters of the specified GPD distribution.
\item \code{ub}: A numeric vector of the upper bounds for the parameters of the specified GPD distribution.
}
}
\description{
Provides the lower and upper bounds for the parameters of certain Generalized Pareto Distributions based on the specified distribution type.
}
\details{
The function returns parameter bounds depending on the type of GPD distribution:
\itemize{
\item For both \code{"bGPD"} and \code{"tgGPD"}, the lower bounds are \code{c(0, 0)}, and the upper bounds are \code{c(Inf, 1)}.
}
}
\examples{
# Get parameter bounds for a beta-GPD distribution
lb_ub_GP("bGPD")

# Get parameter bounds for a truncated GPD distribution
lb_ub_GP("tgGPD")
}
