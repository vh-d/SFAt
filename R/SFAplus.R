# LIKELIHOOD FUNCTIONS -----------------------------------------------------

# likelihood function normal/half-normal distributional assumption
ll_hnorm <- function(params, y, X, deb) {
  nbetas <- ncol(X)

  # +/- scale to control for cost/production frontiers
  sc <- 1

  if (length(params) != nbetas + 1) {
    stop("Incorrect nuber of parameters. ", nbetas + 1, " needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]

  p_sigma2_u <- params[nbetas + 1]
  p_sigma2_v <- params[nbetas + 2]

  epsilon <- y - X %*% p_beta
  sigma <- sqrt(p_sigma2_u + p_sigma2_v)

  N <- length(y)

  ll <-
    -N*log(sigma) +
    -(N*log(sqrt(2) / sqrt(pi))) + # the const term
    +sum(log(pnorm(-sc*(epsilon * (sqrt(p_sigma2_u) / sqrt(p_sigma2_u))) / sigma))) +
    -0.5 * (p_sigma2_u + p_sigma2_v) * sum(epsilon^2)

  return(-ll)
}

# likelihood function normal/exponential distributional assumption
ll_exp <- function(params, y, X, deb) {
  nbetas <- ncol(X)

  # +/- scale to control for cost/production frontiers
  sc <- 1

  if (length(params) != nbetas + 1) {
    stop("Incorrect nuber of parameters. ", nbetas + 1, " needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]

  p_sigma2_u <- params[nbetas + 1]
  p_sigma2_v <- params[nbetas + 2]

  epsilon <- y - X %*% p_beta
  sigma_u <- sqrt(p_sigma2_u)
  sigma_v <- sqrt(p_sigma2_v)

  N <- length(y)

  ll <-
    - N * log(sigma_u) +
    0.5 * N * (p_sigma2_v / p_sigma2_u) +
    sum(log(pnorm(-(sc*epsilon + (p_sigma2_v / sigma_u)) / sigma_v))) +
    sc / sigma_u * sum(epsilon)

  return(-ll)
}

# likelihood function normal/truncated-normal distributional assumption
ll_tnorm <- function(params, y, X, deb) {

  nbetas <- ncol(X)

  # +/- scale to control for cost/production frontiers
  sc <- 1

  if (length(params) != nbetas + 1) {
    stop("Incorrect nuber of parameters. ", nbetas + 1, " needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]

  p_sigma2_u <- params[nbetas + 1]
  p_sigma2_v <- params[nbetas + 2]

  p_mu <- params[nbetas + 3]

  epsilon <- y - X %*% p_beta
  sigma_u <- sqrt(p_sigma2_u)
  sigma_v <- sqrt(p_sigma2_v)

  N <- length(y)

  ll <-
    -(N / 2* log(1 / (2 * pi)) + N * log(1 / (sqrt(sigmau2 + sigmav2)))
      - log(pnorm((p_mu/sqrt(sigmau2))))
      + sum(log(pnorm(((1 - (sigmau2/(sigmau2 + sigmav2)))
                       - sc * (sigmau2 / (sigmau2 + sigmav2)) * epsilon)
                      / (sqrt(sigmau2 * (1 - sigmau2 / (sigmau2 + sigmav2)))))))
      - 1 / (2 * (sigmau2 + sigmav2)) * sum((epsilon + sc * p_mu)^2))

}

# Battise-Coelli model
# parameters: p_beta, p_delta, p_sigma_w, p_sigma_v
ll_bc95_tnorm <- function(params, y, X, Z, deb) {

  # extract parameters from parameter vector
  nbetas <- ncol(X) # number of beta coeffs
  ndeltas <- ncol(Z) # number of delta coeffs

  if (length(params) != nbetas + ndeltas + 2) {
    stop("Incorrect nuber of parameters. ", nbetas, "+", ndeltas, "+ 2 needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]
  p_delta <- params[(nbetas + 1):(nbetas + ndeltas)]

  sigma2_v <- params[(nbetas + ndeltas + 1)] # variance of inefficiency term
  sigma2_w <- params[(nbetas + ndeltas + 2)] # variance of symmetric error term

  if (deb) cat("Total of ", length(params), " parameters: \n",
               "Betas: ", paste(p_beta), "\n",
               "Deltas: ", paste(p_delta), "\n",
               "Sigma2_v: ", sigma2_v,
               "Sigma2_w: ", sigma2_w, "\n")

  # penalize negative variance parameters
  if (sigma2_w <= 0 | sigma2_v <= 0 ) {
    if (deb) {
      cat("One of the sigmas is not positive...", sigma2_w, " or ", sigma2_v)
    }
    return(1e12)
  }

  sigma2 <- sigma2_w + sigma2_v # variance of composite error term

  sigma_v <- sqrt(sigma2_v)
  sigma_w <- sqrt(sigma2_w)
  sigma <- sqrt(sigma2)

  sigma_ast <- sigma_v * sigma_w / sigma

  Zdelta <- as.vector(Z %*% p_delta) # fitted means of inefficiency term

  eps <- as.vector(y - (X %*% p_beta)) # composite error terms

  mu_ast <- (sigma2_w*Zdelta - sigma2_v*eps) / sigma2

  temp_expr1 <- -0.5 * length(y) * (log(2*pi) + log(sigma2))
  temp_expr2 <- -0.5 * sum((eps + Zdelta)^2)/sigma2
  temp_expr3 <- - sum(log(pnorm(Zdelta / sigma_v))
                      - log(pnorm(mu_ast / sigma_ast)))

  if (deb) cat("Log-likelihood: ", temp_expr1, " + ", temp_expr2, " + ", temp_expr3, "\n")

  if (is.na(temp_expr3) | is.infinite(temp_expr3)) {
    if (deb) {
      cat("Infinite term3...", temp_expr3, "\n")
    }
    temp_expr3 <-
      pmax(log(pnorm(Zdelta / sqrt(sigma2_v))), -6e12) -
      pmax(log(pnorm(mu_ast / sigma_ast)), -6e12)

    temp_expr3 <- - sum(temp_expr3)

    if (deb) cat(temp_expr3)

    # return(1e9)
  }

  result <- temp_expr1 + temp_expr2 + temp_expr3
  return(-result)
}


# MAIN FUNCTION -----------------------------------------------------------

#' stochastic frontier analysis
#' @export
sfa.fit <- function(y, X,
                    Z = NULL,
                    intercept = TRUE,
                    intercept_Z = TRUE,
                    model = "panel",
                    dist = "tnorm",
                    start_val = NULL,
                    par_mu = NULL, # set if the mu parameter is known
                    form = "cost",
                    opt_method = "BFGS",
                    deb = T, # TRUE for debug reports
                    control_opt = NULL,
                    ...) {

  n_betas <- ncol(X)
  n_deltas <- ncol(Z)

  if (intercept) {
    n_betas <- n_betas + 1
    X <- cbind(1, X)
  }

  if (intercept_Z) {
    n_deltas <- n_deltas + 1
    Z <- cbind(1, Z)
  }

  # fit by OLS for LR test and starting values
  lmfit <- lm.fit(y = y, x = X)
  if (deb) print(summary(lmfit))

  if (is.null(start_val)) {
    start_val = lmfit
  }

  ll_fn_call <- parse(text = paste0("ll", "_", model, "_", dist))

  initpar = c(lmfit$coefficients[1:n_betas], rep(0, n_deltas), 1, 1)

  # lowerb <- c(rep(-Inf, n_betas + n_deltas), 0.000001, 0.000001)

  est <- optim(initpar,
               fn = eval(ll_fn_call),
               # lower = lowerb,
               method = opt_method,
               control = control_opt,
               hessian = T,
               y = y,
               X = X,
               Z = Z,
               deb = deb)

  return(est)
}


# FORMULA FUNCTION --------------------------------------------------------

#' stochastic frontier analysis
#' @export
sfa <- function(formula,
                data = NULL,
                intercept = TRUE,
                intercept_Z = TRUE,
                dist = "hnormal",
                start_val = NULL,
                par_mu = NULL,
                form = "cost",
                opt_method = "BFGS", ...){

  formula_ext <- Formula(formula)
  formula_length <- length(formula_ext)


  y <- as.vector(
    model.frame(
      formula(formula_ext,
              lhs = 1, rhs = 0)))

  X <- as.matrix(
    model.frame(
      formula(formula_ext,
              lhs = 0, rhs = 1)))

  # exdogenous variables
  if (formula_length[2] > 1) {
    Z <- as.matrix(
      model.frame(
        formula(formula_ext,
                lhs = 0, rhs = 2)))
  } else {
    Z = NULL
  }

  sfa.fit(y = y,
          X = X,
          Z = Z,
          intercept = intercept,
          intercept_Z = intercept_Z,
          dist = dist,
          start_val = start_val,
          par_mu = par_mu,
          form = form,
          method = method,
          ...)
}
