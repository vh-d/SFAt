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

  ans <- list(LR_test_stat = LR_test_stat,
              LR_chisq_df = LR_chisq_df,
              LR_pvalue = LR_pvalue)

  class(ans) <- "lrtest.SFA"

  return(ans)
}


#' Print LR test results for SFA class
#' @export
print.lrtest.SFA <- function(object) {
  cat("Likelihood ratio test:",
      "======================",
      "Null-model: OLS",
      sep = "\n")

  cat(paste0("LR Chisq: ", format(object$LR_test_stat)),
      paste0("Chisq Df: ", format(object$LR_chisq_df)),
      paste0("Pr(>Chisq): ", format(object$LR_pvalue)),
      sep = "\n")

}