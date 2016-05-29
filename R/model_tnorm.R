# NORMAL / T-NORMAL - HOMOSCEDASTIC - CROSS SECTION DATA ----------------------------

par_cs_tnorm <- c(mu = 0.0, lnsigma2_u = 0.5, lnsigma2_v = 0.5)

# likelihood function normal/truncated-normal distributional assumption
# params: beta, mu, log(sigma_u^2), log(sigma_v^2)

ll_cs_tnorm <- function(params, y, X, ineff, deb) {

  nbetas <- ncol(X)

  if (length(params) != nbetas + 3) {
    stop("Incorrect nuber of parameters. ", nbetas + 3, " needed, but ", length(params), " supplied.")
  }

  beta_coef <- params[1:nbetas]

  mu <- params[nbetas + 1]
  lnsigma2_u <- params[nbetas + 2]
  lnsigma2_v <- params[nbetas + 3]

  sigma2_u <- exp(lnsigma2_u)
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(beta_coef), "\n",
               "Mu: ", mu,
               "Sigma_u: ", sigma_u,
               "Sigma_v: ", sigma_v, "\n")

  epsilon <- -ineff * as.vector(y - X %*% beta_coef)
  lambda <- sigma_u / sigma_v
  sigma <- sqrt(sigma2_u + sigma2_v)

  if (deb) cat("lambda: ", lambda,
               "sigma: ", sigma, "\n")

  N <- length(y)

  ll <-
    - N * (log(sigma) + 0.5*log(2*pi) + log(pnorm(mu/sigma_u))) +
    sum(
      log(pnorm((mu / (sigma*lambda)) - ((epsilon*lambda) / sigma))) +
        -0.5*((epsilon + mu) / sigma)^2
    )

  # todo: Kumbhakar 2014 has this function that is actually much more stable
  # check derivation of the likelihood functions
  # ll <-
  #   - N * (log(sigma) + log(pnorm(mu/sigma_u))) +
  #   sum(
  #     log(pnorm((mu / (sigma*lambda)) - ((epsilon*lambda) / sigma))) +
  #       -((epsilon + mu) / sigma)^2
  #   )

  if (deb) cat("Likelihood: ", ll, "\n")

  return(-ll)
}

# INEFFICIENCY TERM ESTIMATION FUNCTION ------
u_cs_tnorm <- function(object, estimator) {
  # extract sigmas from model object
  mu <- object$parameters[length(object$coeff) + 1]
  lnsigma2_u <- object$parameters[length(object$coeff) + 2]
  lnsigma2_v <- object$parameters[length(object$coeff) + 3]

  # derive the rest of the parameters
  sigma2_u <- exp(lnsigma2_u)
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  sigma2 <- sigma2_u + sigma2_v

  mu_ast <- (object$ineff * object$residuals * sigma2_u + sigma2_v * mu) / sigma2
  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  u <- switch(estimator,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_ast*(dnorm(mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast)))

  if (is.null(u)) stop(paste0("Unknown type ", estimator," of conditional mean estimator."))

  return(u)
}