% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{sim_data_cs}
\alias{sim_data_cs}
\title{generate simulated cross-section data for testing purposes}
\usage{
sim_data_cs(N = 500, x_coeff = c(10, 6, 3), z_coeff = c(3, 1.5, 0.5),
  sigma_u = 2, sigma_v = 3, ineff = -1, aslist = F)
}
\arguments{
\item{N}{number of observations}

\item{x_coeff}{coefficients for explanatory variables}

\item{z_coeff}{coefficients for exogeneous varaibles}

\item{sigma_u}{variance of the inefficiency term}

\item{sigma_v}{variance of the random noise}

\item{ineff}{production (-1) or cost inefficiency form (1),}

\item{aslist}{TRUE/FALSE to return list or array}
}
\description{
generate simulated cross-section data for testing purposes
}
