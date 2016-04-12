# NORMAL / T-NORMAL - HOMOSCEDASTIC - CROSS SECTION DATA ----------------------------

par_cs_tnorm <- c(sigma_u = 1, sigma_v = 1)

ll_cs_tnorm <- function(params, y, X, ineff, deb) {

  nbetas <- ncol(X)

  if (length(params) != nbetas + 3) {
    stop("Incorrect nuber of parameters. ", nbetas + 3, " needed, but ", length(params), " supplied.")
  }

  beta_coef <- params[1:nbetas]

  mu <- params[nbetas + 1]
  sigma_u <- params[nbetas + 2]
  sigma_v <- params[nbetas + 3]

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(beta_coef), "\n",
               "Mu: ", mu,
               "Sigma_u: ", sigma_u,
               "Sigma_v: ", sigma_v, "\n")

  if (sigma_u <= 0 | sigma_v <= 0) return(1e12)

  epsilon <- -ineff * as.vector(y - X %*% beta_coef)
  lambda <- sigma_u / sigma_v
  sigma <- sqrt(sigma_u^2 + sigma_v^2)

  if (deb) cat("lambda: ", lambda,
               "sigma: ", sigma, "\n")

  N <- length(y)

  ll <-
    - N * (log(sigma) + 0.5*log(2*pi) + log(pnorm(mu/sigma_u))) +
    sum(
      log(pnorm((mu / (sigma*lambda)) - ((epsilon*lambda) / sigma))) +
        -0.5*((epsilon + mu) / sigma)^2
    )

  # Kumbhakar 2014 has this function that is actually much more stable
  # ll <-
  #   - N * (log(sigma) + log(pnorm(mu/sigma_u))) +
  #   sum(
  #     log(pnorm((mu / (sigma*lambda)) - ((epsilon*lambda) / sigma))) +
  #       -((epsilon + mu) / sigma)^2
  #   )

  if (deb) cat("Likelihood: ", ll, "\n")

  return(-ll)
}
u_cs_tnorm <- function(object, type) {
  # extract sigmas from model object
  mu <- object$parameters[length(object$coefficients) + 1]
  sigma <- object$parameters[length(object$coefficients) + 2]
  lambda <- object$parameters[length(object$coefficients) + 3]

  sigma_u <- lambda * sigma / sqrt(1 + lambda^2)
  sigma2_u <- sigma_u^2

  sigma2 <- sigma^2
  sigma2_v <- sigma2 - sigma2_u
  sigma_v <- sqrt(sigma2_v)

  mu_ast <- (object$ineff * object$residuals * sigma2_u + sigma2_v * mu) / sigma2
  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  u <- switch(type,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_ast*(dnorm(mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast)))
  if (is.null(u)) stop(paste0("Unknown type ", type ," of conditional mean estimator."))

  return(u)
}