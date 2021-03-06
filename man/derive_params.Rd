% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_params.R
\name{mean.list}
\alias{mean.list}
\alias{median.list}
\alias{variance.list}
\title{mean.list}
\usage{
\method{mean}{list}(params, distr = "tnorm")

\method{median}{list}(params, distr = "tnorm")

variance.list(params, distr = "tnorm")
}
\arguments{
\item{params}{named list of distribution parameters}

\item{distr}{name of distribution, Default: 'tnorm'}
}
\value{
numeric
}
\description{
Derive mean, median a variance of various distributions.
}
\examples{
\dontrun{
if(interactive()){
 mean(list(mu = 2, sigma = 3))
 mean(list(shape = 2, scale = 3), distr = "gamma")
 }
}
}
