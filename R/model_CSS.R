#' Cornwell, Schmidt, and Sickles (1990) model
#' as is Parmeter & Kumbhakar 2014, pp.78:
#' 1. use within estimator to fit a fixed effect panel data model
#' 1. estimate residuals
#' 1. regress residuals on time polynomial
#' 1. predict \alpha_{it}
#' 1. compute inefficiency as a difference between max(\alpha_{i}) - \alpha_{it}
fit_css_model <- function(y,
                          X,
                          K,
                          ineff,
                          deb = T) {
  require(plm)

  formula_frontier <- as.formula(
    paste0(c("y ~ ",
             paste0(colnames(pdata$X), collapse = " + "),
             " - 1"),
           collapse = ""))

  dataframe <- pdata.frame(data.frame(K, y, X), index = c("k", "t"))

  if (deb) print(head(dataframe))

  step1fit <- plm(formula_frontier, data = dataframe, model = "within")
  alpha_0 <- fixef(step1fit)

  if (deb) {
    print(summary(step1fit))
    print(alpha_0)
  }

  step2fit <- lm(formula = res ~ as.factor(k):t + as.factor(k):I(t^2),
                 data    = data.frame(res = step1fit$residuals,
                                         pdata$K))

  if (deb) {
    print(summary(step2fit))
  }

  result <- list(
    fe = alpha_0,
    tvar = step2fit$coefficients,
    step1fit = step1fit,
    step2fit = step2fit
  )

  # lmfit <- summary(lmfit)

  return(result)
}

uit_css <- function(x) {
  predict(x$step2fit) + rep(x$fe, 10)
}

# require(SFAt)
# require(data.table)
# require(ggplot2)
# pdata <- sim_data_panel(z_mean = 10:30, sigma_u = 1, sigma_v =1, z_coeff = 0, aslist = T)
# temp <- fit_css_model(pdata$y, X = pdata$X, K = pdata$K, ineff = -1)
# plot(temp$fe)
#
# uit <- data.table(pdata$K, uit = uit_css(temp))
# uit[, ui := max(uit), by = t]
# uit[, ineff := ui - uit]
# # ui <- aggregate(x = uit$uit, by = list(uit$t), FUN = max)
# # ineff <- rep(ui$x, times = 20) - uit$uit
# # hist(ineff)
#
# ggplot(uit,
#       aes(x = t,
#           y = ineff,
#           color = factor(k))) +
#   geom_line()
#
# ggplot(uit,
#        aes(x = ineff)) +
#   geom_histogram()
#
