# LR test function --------------------------------------------------------
#' LR test for SFA class
#' @export
lrtest.SFA <- function(object) {

  LR_test_stat <- 2*(object$loglik - logLik(object$lmfit))
  LR_chisq_df <- length(object$parameters) - attributes(logLik(object$lmfit))$df
  if (LR_chisq_df > 1) {
    LR_pvalue <-
      0.25*pchisq(LR_test_stat, LR_chisq_df-2, lower.tail = FALSE) +
      0.5*pchisq(LR_test_stat, LR_chisq_df-1, lower.tail = FALSE) +
      0.25*pchisq(LR_test_stat, LR_chisq_df, lower.tail = FALSE)
  } else {
    LR_pvalue <-
      0.5*pchisq(LR_test_stat, LR_chisq_df-1, lower.tail = FALSE) +
      0.5*pchisq(LR_test_stat, LR_chisq_df, lower.tail = FALSE)
  }

  cat("Likelihood ratio test:",
      "======================",
      "Null-model: OLS",
      sep = "\n")

  cat(paste0("LR Chisq: ", round(LR_test_stat, 3)),
      paste0("Chisq Df: ", round(LR_chisq_df, 3)),
      paste0("Pr(>Chisq): ", round(LR_pvalue, 3)),
      sep = "\n")
}

#' @export
lrtest <- function(object, ...){
  UseMethod("lrtest")
}