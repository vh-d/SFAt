# derive mean, median and variance of various distributions from their parameters ----------------------------

#' @title mean.list
#' @description Derive mean, median a variance of various distributions.
#' @param params named list of distribution parameters
#' @param distr name of distribution, Default: 'tnorm'
#' @return numeric
# @details
#' @examples
#' \dontrun{
#' if(interactive()){
#'  mean(list(mu = 2, sigma = 3))
#'  mean(list(shape = 2, scale = 3), distr = "gamma")
#'  }
#' }
#' @export
#' @rdname derive_params
mean.list <- function(params, distr = "tnorm") {
  switch(distr,
         tnorm = params$mu + params$sigma*dnorm(-params$mu/params$sigma)/(1 - pnorm(-params$mu/params$sigma)), # lower bound 0 is assumed
         gamma = params$shape * params$scale)
}

#' @title median.list
#' @export
#' @rdname derive_params
median.list <- function(params, distr = "tnorm") {
  switch(distr,
         tnorm = params$mu + params$sigma*qnorm(0.5*pnorm(-params$mu/params$sigma) + 0.5)) # lower bound 0 is assumed
}

#' @title variance.list
#' @export
#' @rdname derive_params
variance.list <- function(params, distr = "tnorm") {
  switch(distr,
         gamma = params$shape * params$scale^2)
}

# # tests
# require(SFAt)
# mean(list(mu = -2, sigma = 2))
# mean(rtnorm(10000, mean = -2, sd = 2, lower = 0))
#
# mean(list(shape = 2, scale = 10), distr = "gamma")
# mean(rgamma(10000, shape = 2, rate = 1/10))
