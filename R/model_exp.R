# NORMAL / EXPONENTIAL - HOMOSCEDASTIC - CROSS SECTION DATA ----------------------------

par_cs_exp <- c(sigma2_u = 1, sigma2_u = 1)

# likelihood function normal/exponential distributional assumption
# params: beta, sigma_u^2, sigma_v^2
ll_cs_exp <- function(params, y, X, ineff, deb) {
  nbetas <- ncol(X)

  if (length(params) != nbetas + 2) {
    stop("Incorrect nuber of parameters. ", nbetas + 2, " needed, but ", length(params), " supplied.")
  }

  beta_coef <- params[1:nbetas]

  sigma2_u <- params[nbetas + 1]
  sigma2_v <- params[nbetas + 2]

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(beta_coef), "\n",
               "Sigma_u^2: ", sigma2_u,
               "Sigma_v^2: ", sigma2_v, "\n")

  if (sigma2_u <= 0 | sigma2_v <= 0) return(1e12)

  epsilon <- -ineff * (y - X %*% beta_coef)
  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)

  N <- length(y)

  ll <-
    - N * log(sigma_u) +
    0.5 * N * (sigma2_v / sigma2_u) +
    sum(log(pnorm(-(epsilon + (sigma2_v / sigma_u)) / sigma_v))) +
    sum(epsilon) / sigma_u

  return(-ll)
}

u_cs_exp <- function(object, type) {

  # extract sigmas from model object
  sigma2_u <- object$parameters[length(object$coefficients) + 1]
  sigma2_v <- object$parameters[length(object$coefficients) + 2]
  sigma_v <- sqrt(sigma2_v)
  sigma_u <- sqrt(sigma2_u)

  mu_ast <- object$ineff * object$residuals - sigma2_v / sigma_u

  u <- switch(type,
              # ME = pmax(mu_ast, 0),
              # BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_v*dnorm(mu_ast/sigma_v)/pnorm(mu_ast/sigma_v))
  if (is.null(u)) stop(paste0("Unknown type ", type ," of conditional mean estimator."))

  return(u)
}