#' @export
sim_data <- function(N = 500,
                     x_coeff = c(10, 6, 3),
                     z_coeff = c(5, 1.5, -0.9),
                     sigma_u = 2,
                     sigma_v = 3,
                     ineff = -1,
                     aslist = F) {

  require(msm)

  x_ncols <- length(x_coeff) - 1
  z_ncols <- length(z_coeff) - 1

  X <- matrix(rtnorm(n = N*x_ncols,
                     mean = 5,
                     sd = 3,
                     lower = 0),
              nrow = N,
              ncol = x_ncols)

  colnames(X) <- paste0("x", 1:x_ncols)

  if (z_ncols > 0 ) {
    Z <- matrix(rnorm(n = N*z_ncols,
                      mean = 3,
                      sd = 2),
                nrow = N,
                ncol = z_ncols)

    colnames(Z) <- paste0("z", 1:z_ncols)

    u <- rtnorm(N,
                mean = as.vector(cbind(1,  Z) %*% z_coeff),
                lower = 0,
                sd = sigma_u)
  } else {

    u <- rtnorm(N,
                mean = z_coeff,
                lower = 0,
                sd = sigma_u)
    Z <- NULL
  }

  v <- rnorm(N,
             mean = 0,
             sd = sigma_v)

  eps <- ineff*u + v

  # y <- 10 + 6*x1 + 3*x2 + eps
  y <- as.vector( cbind(1, X) %*% x_coeff) + eps

  if (aslist) {
    return(list(y = y, X = X, Z = Z, u = u, v = v, eps = eps))
  } else
    return(cbind(y, X, Z, u, v, eps))
}

#' @export
sim_panel_data <- function(obs = 10,
                           k = 20,
                           x_coeff = c(10, 6, 3),
                           z_coeff = c(5, 1.9, -0.9),
                           sigma_u = 2,
                           sigma_v = 3,
                           ineff = -1,
                           aslist = F) {

  N <- obs * k
  require(msm)

  x_ncols <- length(x_coeff) - 1
  z_ncols <- length(z_coeff) - 1

  X <- matrix(rtnorm(n = N*x_ncols,
                     mean = 5,
                     sd = 3,
                     lower = 0),
              nrow = N,
              ncol = x_ncols)

  colnames(X) <- paste0("x", 1:x_ncols)

  mu <- rep(rnorm(k, z_coeff[1], 6), each = k)

  if (z_ncols > 0 ) {
    Z <- matrix(rnorm(n = N*z_ncols,
                      mean = 3,
                      sd = 2),
                nrow = N,
                ncol = z_ncols)

    colnames(Z) <- paste0("z", 1:z_ncols)

    u <- rtnorm(N,
                mean = mu + as.vector(Z %*% z_coeff[-1]),
                lower = 0,
                sd = sigma_u)
  } else {

    u <- rtnorm(N,
                mean = mu,
                lower = 0,
                sd = sigma_u)
    Z <- NULL
  }

  v <- rnorm(N,
             mean = 0,
             sd = sigma_v)

  eps <- ineff*u + v

  # y <- 10 + 6*x1 + 3*x2 + eps
  y <- as.vector( cbind(1, X) %*% x_coeff) + eps

  if (aslist) {
    return(list(y = y, X = X, Z = Z, u = u, v = v, eps = eps))
  } else
    return(cbind(y, X, Z, u, v, eps))
}

sim_data_cs_homo <- function(N, ineff = -1) {

    require(msm)

  x1 <- rnorm(N, 8)
  x2 <- rnorm(N, 15)

  u <- rtnorm(N, mean = 2, lower = 1, sd = 3)
  v <- rnorm(N, 0, 1)
  eps <- ineff*u + v

  y <- 10 + 6*x1 + 3*x2 + eps

  return(cbind(y, x1, x2, z1, z2, u, v, eps))
}
