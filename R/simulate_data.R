require(data.table)
require(msm)
# @import data.table

#' @export
sim_data <- function(N, ineff = -1) {
  # require(data.table)

  # N = 100

  x1 <- rnorm(N, 8)
  x2 <- rnorm(N, 15)
  z1 <- rnorm(N, 2, 1)
  z2 <- rnorm(N, 5, 1)

  # x1x2 <- x1*x2

  u <- rtnorm(N, mean = 5 + z1 - 0.3*z2, lower = 1, sd = 3)

  v <- rnorm(N, 0, 1)

  y <- 10 + 6*x1 + 3*x2 + ineff*u + v

  return(cbind(y, x1, x2, z1, z2, u, v))
}
