# LIKELIHOOD FUNCTIONS -----------------------------------------------------

# likelihood function normal/half-normal distributional assumption
# params: beta, sigma_u^2, sigma_v^2
ll_cs_hnorm <- function(params, y, X, ineff, deb) {
  nbetas <- ncol(X)

  if (length(params) != nbetas + 2) {
    stop("Incorrect nuber of parameters. ", nbetas + 2, " needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]

  p_sigma2_u <- params[nbetas + 1]
  p_sigma2_v <- params[nbetas + 2]

  epsilon <- -ineff*(y - X %*% p_beta)
  sigma <- sqrt(p_sigma2_u + p_sigma2_v)

  N <- length(y)

  ll <-
    -N*log(sigma) +
    -(N*log(sqrt(2) / sqrt(pi))) + # the const term
    +sum(log(pnorm(-(epsilon * (sqrt(p_sigma2_u) / sqrt(p_sigma2_u))) / sigma))) +
    -0.5 * (p_sigma2_u + p_sigma2_v) * sum(epsilon^2)

  return(-ll)
}

# likelihood function normal/exponential distributional assumption
# params: beta, sigma_u^2, sigma_v^2
ll_cs_exp <- function(params, y, X, ineff, deb) {
  nbetas <- ncol(X)

  if (length(params) != nbetas + 2) {
    stop("Incorrect nuber of parameters. ", nbetas + 2, " needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]

  p_sigma2_u <- params[nbetas + 1]
  p_sigma2_v <- params[nbetas + 2]

  epsilon <- -ineff * (y - X %*% p_beta)
  sigma_u <- sqrt(p_sigma2_u)
  sigma_v <- sqrt(p_sigma2_v)

  N <- length(y)

  ll <-
    - N * log(sigma_u) +
    0.5 * N * (p_sigma2_v / p_sigma2_u) +
    sum(log(pnorm(-(epsilon + (p_sigma2_v / sigma_u)) / sigma_v))) +
    sum(epsilon) / sigma_u

  return(-ll)
}

# likelihood function normal/truncated-normal distributional assumption
# params: beta, sigma^2, lambda, mu
ll_cs_tnorm <- function(params, y, X, ineff, deb) {

  nbetas <- ncol(X)

  if (length(params) != nbetas + 3) {
    stop("Incorrect nuber of parameters. ", nbetas + 3, " needed, but ", length(params), " supplied.")
  }

  p_beta <- params[1:nbetas]

  p_mu <- params[nbetas + 1]

  p_sigma <- params[nbetas + 2]
  p_lambda <- params[nbetas + 3]

  epsilon <- -ineff * (y - X %*% p_beta)
  sigma_u <- p_lambda * p_sigma / sqrt(1 + p_lambda^2)

  N <- length(y)

  ll <-
    - N * (0.5 * log(2*pi) + log(p_sigma) + log(pnorm(p_mu/sigma_u))) +
    sum(
      log(pnorm((p_mu / (p_sigma * p_lambda)) - ((epsilon * p_lambda) / p_sigma))) +
        - 0.5 * ((epsilon + p_mu) / p_sigma)^2
    )

  return(-ll)
}


# likelihood function normal/gamma distributional assumption
# ll_gamma <- function() {
#
# }


# Advanced version of (Battese and Coelli, 1995) and (Huang and Liu, 1994) models
# heterogeneity in efficiency term: endogeneous location parameter mu
# implemented as Hadri et al. 2003
# parameters: p_beta, p_delta, p_sigma_w, p_sigma_v
ll_cs_tnorm_bc95 <- function(params, y, X, Z, ineff, deb) {

  # extract parameters from parameter vector
  nbetas <- ncol(X) # number of beta coeffs
  ndeltas <- ncol(Z) # number of delta coeffs

  if (length(params) != nbetas + ndeltas + 2) {
    stop("Incorrect nuber of parameters. ",
         nbetas, "+", ndeltas, "+ 2 needed, but ", length(params), " supplied.")
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

  N <- length(y)

  sigma2 <- sigma2_w + sigma2_v # variance of composite error term

  sigma_v <- sqrt(sigma2_v)
  sigma_w <- sqrt(sigma2_w)
  sigma <- sqrt(sigma2)

  sigma_ast <- sigma_v * sigma_w / sigma

  Zdelta <- as.vector(Z %*% p_delta) # fitted means of inefficiency term

  eps <- as.vector(y - (X %*% p_beta)) # composite error terms

  mu_ast <- (sigma2_w*Zdelta - sigma2_v*eps) / sigma2

  temp_expr1 <- -0.5 * N * (log(2*pi) + log(sigma2))
  temp_expr2 <- -0.5 * sum((-ineff*eps + Zdelta)^2)/sigma2
  temp_expr3 <- - sum(log(pnorm(Zdelta / sigma_v)))
  temp_expr4 <- + sum(log(pnorm(mu_ast / sigma_ast)))


  if (deb) cat("Log-likelihood: ",
               temp_expr1, " + ", temp_expr2, " + ",
               temp_expr3, " + ", temp_expr4, "\n")

  if (is.na(temp_expr3) | is.infinite(temp_expr3)) {
    if (deb) {
      cat("Infinite term3...", temp_expr3, "\n")
    }
    temp_expr3 <- pmin(log(pnorm(Zdelta / sqrt(sigma2_v))), 6e12)
    temp_expr3 <- - sum(temp_expr3)

    if (deb) cat(temp_expr3, "\n")
  }

  if (is.na(temp_expr4) | is.infinite(temp_expr4)) {
    if (deb) {
      cat("Infinite term3...", temp_expr4, "\n")
    }
    temp_expr4 <- pmax(log(pnorm(mu_ast / sigma_ast)), -6e12)
    temp_expr4 <- sum(temp_expr4)

    if (deb) cat(temp_expr4, "\n")
  }

  result <- temp_expr1 + temp_expr2 + temp_expr3 + temp_expr4
  return(-result)
}


# MAIN FUNCTION -----------------------------------------------------------

#' sfa.fit
#'
#' Fits stochastic frontier analysis (SFA) model
#'
#' @param y dependent (production/cost) variable.
#' @param X variables of the production/cost function.
#' @param Z exogneous determinants of the mean inefficiency location.
#' @param intercept TRUE/FALSE if the intercept should be included in the main formula.
#' @param intercept_Z TRUE/FALSE if the intercept should be included in the inefficiency location formula.
#' @param structure "cs" for cross-section or "panel" for panel data model.
#' @param dist distribution of inefficiency term ("hnorm", "exp", "tnorm").
#' @param spec specifies what model of endogeneous inefficiency term should be used (currently only bc95 for cross-section implemented).
#' @param start_val starting value of model parameters to be passed to optimization routine.
#' @param ineff -1 (or 1) for production (or cost) function.
#' @param opt_method optimization method.
#' @param deb debug mode (TRUE/FALSE).
#' @param control_opt list of options for optimization routine.
#'
#' @return list
#' @export
sfa.fit <- function(y, X,
                    Z = NULL,
                    intercept = TRUE,
                    intercept_Z = TRUE,
                    structure = "cs",
                    dist = "tnorm",
                    spec = "bc95",
                    start_val = NULL,
                    ineff = -1,
                    opt_method = "BFGS",
                    deb = F, # TRUE for debug reports
                    control_opt = NULL) {

  n_betas <- ncol(X)
  n_deltas <- if (is.null(Z)) 0 else ncol(Z)
  x_names <- colnames(X)

  if (intercept) {
    n_betas <- n_betas + 1
    X <- cbind(1, X)
    x_names <- c("intercept", x_names)
  }

  if (deb) {
    cat(ifelse(intercept == T, "X", "no X"),  "intercept, ",
        n_betas, " X coefficient parameters, ",
        n_deltas, " Z coefficient parameters,", "\n")
  }


  # STARTING VALUES

  # lmfit <- lm.fit(y = y, x = X)

  if (is.null(start_val)) {
    if (is.null(spec)) {

      # fit OLS for starting values of betas and sigma
      lmfit <- lm(y ~ X - 1)
      if (deb) print(summary(lmfit))

      if (dist == "tnorm") {
        start_val <- c(lmfit$coefficients,
                       mu = 0,
                       sigma_u = 1,
                       sigma_v = 1)
      }
      else {
        start_val <- c(lmfit$coefficients,
                       sigma_u = 1,
                       sigma_v = 1)
      }
    } else {

      # fit OLS for starting values of betas and sigma
      lmfit <- lm(y ~ X + Z - 1)
      if (deb) print(summary(lmfit))

      start_val <-
        c(lmfit$coefficients[1:n_betas],
          if (intercept_Z) 0 else NULL,
          lmfit$coefficients[(n_betas + 1) : (n_betas + n_deltas)],
          1, 1)

      if (intercept_Z) {
        # n_deltas <- n_deltas + 1
        Z <- cbind(1, Z)
      }

    }
  }

  if (deb) print(start_val)

  ll_fn_call <- parse(text = paste0("ll", "_",
                                    structure, "_",
                                    dist,
                                    if (is.null(spec)) NULL else paste0("_", spec)))


  # lowerb <- c(rep(-Inf, n_betas + n_deltas), 0.000001, 0.000001)

  if (is.null(Z)) {
   est <- optim(start_val,
               fn = eval(ll_fn_call),
               method = opt_method,
               control = control_opt,
               hessian = T,
               y = y,
               X = X,
               ineff = ineff,
               deb = deb)
  } else {
     est <- optim(start_val,
               fn = eval(ll_fn_call),
               method = opt_method,
               control = control_opt,
               hessian = T,
               y = y,
               X = X,
               Z = Z,
               ineff = ineff,
               deb = deb)
  }

  fit <- list(coeffs = est$par,
              N = length(y),
              ineff_name = ifelse(ineff == -1, "production", "cost"),
              call = list(ineff = ineff,
                          intercept = intercept,
                          intercept_Z = intercept_Z,
                          dist = dist,
                          spec = spec,
                          structure = structure),
              loglik = est$val,
              hessian = est$hessian,
              lmfit = lmfit)

  class(fit) <- c("sfa")

  return(fit)
}


# FORMULA FUNCTION --------------------------------------------------------

#' stochastic frontier analysis
# this formula interface is not ready yet
sfa <- function(formula,
                data = NULL,
                intercept = TRUE,
                intercept_Z = TRUE,
                dist = "hnormal",
                start_val = NULL,
                par_mu = NULL,
                ineff = -1,
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
          model = model,
          start_val = start_val,
          form = form,
          method = method,
          ...)
}


# SUMMARY FUNCTION --------------------------------------------------------
#' @export
summary.sfa <- function(sfa_model) {
var_beta <- diag(solve(est$hessian))

tstats <- est$coeff/sqrt(var_beta)

pvalues <- 2 * pt(q = abs(tstats),
                  df = est$N - length(est$par)+3,
                  lower.tail = FALSE)

coef.table <- round(cbind(est$coeff,
                          sqrt(var_beta),
                          tstats,
                          pvalues),
                    3)
colnames(coef.table) <- c("Estimate", "Std. Error","t-stat", "Pr(>|t|)")
row.names(coef.table) <- names(est$coeff)

print(coef.table)

ll_sfa <- -est$loglik
ll_ols <- logLik(lmfit)
LR_test_stat <- 2*(ll_sfa - ll_ols)[1]
chisq_df <- (4+3+2) - attributes(logLik(lmfit))$df
p_value <- pchisq(LR_test_stat, chisq_df, lower.tail = FALSE)

print(LR_test_stat)
print(p_value)

}
