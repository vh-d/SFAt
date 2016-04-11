par_cs_hnorm <- c(sigma2_u = 1, sigma2_u = 1)

# likelihood function normal/half-normal distributional assumption
# params: beta, sigma_u^2, sigma_v^2
ll_cs_hnorm <- function(params, y, X, ineff, deb) {
  nbetas <- ncol(X)

  if (length(params) != nbetas + 2) {
    stop("Incorrect nuber of parameters. ", nbetas + 2, " needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]

  sigma2_u <- params[nbetas + 1]
  sigma2_v <- params[nbetas + 2]

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(p_beta), "\n",
               "Sigma_u^2: ", sigma2_u,
               "Sigma_v^2: ", sigma2_v, "\n")

  if (sigma2_u <= 0 | sigma2_v <= 0) return(abs(min(sigma2_v, sigma2_u)-0.001)*1e12)

  epsilon <- -ineff*as.vector(y - X %*% p_beta)
  sigma <- sqrt(sigma2_u + sigma2_v)

  N <- length(y)

  ll <-
    -N*log(sigma) +
    -(N*log(sqrt(2 / pi))) + # the const term
    +sum(log(pnorm(-(epsilon * sqrt(sigma2_u) / sqrt(sigma2_v)) / sigma))) - 0.5 / (sigma2_u + sigma2_v) * sum(epsilon^2)

  if (deb) {
    if (is.nan(ll)) print(paste(sigma))
    cat("Log-ll: ", ll, "\n")
  }

  return(-ll)
}

# Jondrow et al. (1982) as in Parmeter-Kumbhakar
u_cs_hnorm <- function(object, type) {
  # extract sigmas from model object
  sigma2_u <- object$parameters[length(object$coefficients) + 1]
  sigma2_v <- object$parameters[length(object$coefficients) + 2]

  sigma2 <- sigma2_u + sigma2_v
  mu_ast <- object$ineff * object$residuals * sigma2_u / sigma2
  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  u <- switch(type,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_ast*(dnorm(mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast)))
  if (is.null(u)) stop(paste0("Unknown type ", type ," of conditional mean estimator."))

  return(u)
}