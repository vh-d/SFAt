# NORMAL / EXPONENTIAL - HOMOSCEDASTIC - CROSS SECTION DATA ----------------------------

par_cs_exp <- c(lnsigma2_u = 0, lnsigma2_v = 0)

# likelihood function normal/exponential distributional assumption
# params: beta, log(sigma_u^2), log(sigma_v^2)

ll_cs_exp <- function(params, y, X, ineff, deb) {

  nbetas <- ncol(X)

  if (length(params) != nbetas + 2) {
    stop("Incorrect nuber of parameters. ", nbetas + 2, " needed, but ", length(params), " supplied.")
  }

  beta_coef <- params[1:nbetas]

  lnsigma2_u <- params[nbetas + 1]
  lnsigma2_v <- params[nbetas + 2]

  sigma2_u <- exp(lnsigma2_u)
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(beta_coef), "\n",
               "Sigma_u: ", sigma_u,
               "Sigma_v: ", sigma_v, "\n")

  epsilon <- -ineff * (y - X %*% beta_coef)

  N <- length(y)

  ll <-
    - N * log(sigma_u) +
    0.5 * N * (sigma2_v / sigma2_u) +
    sum(log(pnorm(-(epsilon + (sigma2_v / sigma_u)) / sigma_v))) +
    sum(epsilon) / sigma_u

  return(-ll)
}

u_cs_exp <- function(object, estimator) {

  # extract sigmas from model object
  lnsigma2_u <- object$parameters[length(object$coeff) + 1]
  lnsigma2_v <- object$parameters[length(object$coeff) + 2]

  # derive the rest of the parameters
  sigma2_u <- exp(lnsigma2_u)
  sigma2_v <- exp(lnsigma2_v)

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  sigma2 <- sigma2_u + sigma2_v

  mu_ast <- object$ineff * object$residuals - sigma2_v / sigma_u
  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  u <- switch(estimator,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_v*dnorm(mu_ast/sigma_v)/pnorm(mu_ast/sigma_v))
  if (is.null(u)) stop(paste0("Unknown type ", estimator," of conditional mean estimator."))

  return(u)
}