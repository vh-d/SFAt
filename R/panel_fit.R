#' fit SFA models on panel data
#' @param y - dependent variable (typically log of total product or cost)
#' @param X - explanatory variables (typically log of inputs or translog form)
#' @param K - matrix of panel data indeces
#' @param ineff - -1L/1L for production/cost inefficiency
#' @param spec - panel data model specification (currently only CSS, 1990)
#' @param deb - debug option
#' @param debll - debug option for log-likelohood functions
#' @export
panel.sfa.fit <- function(y,
                          X,
                          # CM = NULL,
                          # CV_u = NULL,
                          # CV_v = NULL,
                          K = NULL,
                          ineff = -1L,
                          spec = c("css90"),
                          # intercept = list(f = TRUE,
                          #                  cm = TRUE,
                          #                  cv_u = TRUE,
                          #                  cv_v = TRUE),
                          # sv = list(f = NULL,
                          #           cm = NULL,
                          #           cv_u = NULL,
                          #           cv_v = NULL),
                          # ll = NULL,
                          # opt_method = "SANN",
                          # opt_control = NULL,
                          deb = F, # TRUE for debug reports
                          debll = F) {

  result <- do.call(
    paste0("fit_", spec, "_model"),
    args = list(y = y,
                X = X,
                K = K,
                ineff = ineff)
  )

  result$data <- list(y = y,
                      X = X,
                      K = K)

  result$ineff = ineff

  class(result) <- "psfa"

  return(result)
}

