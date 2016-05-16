# NORMAL / T-NORMAL MODEL- CROSS SECTION DATA ----------------------------

par_cs_tnorm <- c(sigma_u = 1, sigma_v = 1)

# Advanced version of (Battese and Coelli, 1995) and (Huang and Liu, 1994) models
# heterogeneity in efficiency term: endogeneous location parameter mu
# implemented as Hadri et al. 2003
# parameters: beta_coef, delta_coef, sigma2_v, sigma2_u
ll_cs_tnorm_bc95 <- function(params, y, X, Z, ineff, deb) {

  # extract parameters from parameter vector
  nbetas <- ncol(X) # number of beta coeffs
  ndeltas <- ncol(Z) # number of delta coeffs

  if (length(params) != nbetas + ndeltas + 2) {
    stop("Incorrect nuber of parameters. ",
         nbetas, "+", ndeltas, "+ 2 needed, but ", length(params), " supplied.")
  }

  beta_coef <- params[1:nbetas]
  delta_coef <- params[(nbetas + 1):(nbetas + ndeltas)]

  sigma2_u <- params[(nbetas + ndeltas + 1)] # variance of inefficiency term
  sigma2_v <- params[(nbetas + ndeltas + 2)] # variance of symmetric error term

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(beta_coef), "\n",
               "Deltas: ", paste(delta_coef), "\n",
               "Sigma2_u: ", sigma2_u,
               "Sigma2_v: ", sigma2_v, "\n")

  # penalize negative variance parameters
  if (sigma2_v <= 0 | sigma2_u <= 0 ) {
    if (deb) {
      cat("One of the sigmas is not positive...", sigma2_v, " or ", sigma2_u)
    }
    return(1e12)
  }

  N <- length(y)

  sigma2 <- sigma2_v + sigma2_u # variance of composite error term

  sigma_u <- sqrt(sigma2_u)
  sigma_v <- sqrt(sigma2_v)
  sigma <- sqrt(sigma2)

  sigma_ast <- sigma_u * sigma_v / sigma

  Zdelta <- as.vector(Z %*% delta_coef) # fitted means of inefficiency term

  eps <- as.vector(y - (X %*% beta_coef)) # composite error terms

  mu_ast <- (sigma2_v*Zdelta - sigma2_u*eps) / sigma2

  temp_expr1 <- -0.5 * N * (log(2*pi) + log(sigma2))
  temp_expr2 <- -0.5 * sum((-ineff*eps + Zdelta)^2)/sigma2
  temp_expr3 <- sum(log(pnorm(mu_ast / sigma_ast)) - log(pnorm(Zdelta / sigma_u)))

  if (deb) cat("Log-likelihood: ",
               temp_expr1, " + ", temp_expr2, " + ",
               temp_expr3, " + ",
               # temp_expr4,
               "\n")

  # result <- temp_expr1 + temp_expr2 + temp_expr3 + temp_expr4
  result <- temp_expr1 + temp_expr2 + temp_expr3

  if (deb) cat("Loglikelihood: ", result,  "\n")

  return(-result)
}

u_cs_tnorm_bc95 <- function(object, type) {
  # extract sigmas from model object
  sigma2_u <- object$parameters[length(object$coefficients) + length(object$coefficients_Z) + 1]
  sigma2_v <- object$parameters[length(object$coefficients) + length(object$coefficients_Z) + 2]

  mu <- as.vector(object$data$Z %*% object$coefficients_Z)
  sigma2 <- sigma2_u + sigma2_v
  mu_ast <- (object$ineff * object$residuals * sigma2_u + sigma2_v * mu) / sigma2
  sigma_ast <- sqrt(sigma2_u * sigma2_v / sigma2)

  u <- switch(type,
              ME = pmax(mu_ast, 0),
              BC = -log(exp(-mu_ast + 0.5 * sigma_ast^2)*(pnorm(-sigma_ast + mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast))),
              JLMS = mu_ast + sigma_ast*(dnorm(mu_ast/sigma_ast)/pnorm(mu_ast/sigma_ast)))
  if (is.null(u)) stop(paste0("Unknown type ", type ," of conditional mean estimator."))

  return(u)
}
