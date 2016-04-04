require(data.table)
require(msm)

sim_data <- function(N = 500,
                     x_coeff = c(10, 6, 3),
                     z_coeff = c(4, 1.5, -0.9),
                     sigma_u = 2,
                     sigma_v = 3,
                     ineff = -1) {

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
                      mean = 4,
                      sd = 3),
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

  return(cbind(y, X, Z, u, v, eps))
}

sim_data_cs_homo <- function(N, ineff = -1) {
  # require(data.table)

  # N = 100

  x1 <- rnorm(N, 8)
  x2 <- rnorm(N, 15)

  u <- rtnorm(N, mean = 2, lower = 1, sd = 3)
  v <- rnorm(N, 0, 1)
  eps <- ineff*u + v

  y <- 10 + 6*x1 + 3*x2 + eps

  return(cbind(y, x1, x2, z1, z2, u, v, eps))
}
