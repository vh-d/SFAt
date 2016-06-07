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

  formula_frontier <- as.formula(paste0(c("y ~ ", paste0(colnames(pdata$X), collapse = " + ")), collapse = ""))
  dataframe <- pdata.frame(data.frame(K, y, X), index = c("k", "t"))
  head(dataframe)

  plmfit <- plm(formula_frontier, data = dataframe, model = "within")

  if (deb) {
    print(summary(plmfit))
    print(fixef(plmfit))
  }

  rfit <- plmfit$residuals

  lmfit <- lm(rfit ~ as.factor(k) + as.factor(k):t  + as.factor(k):I(t^2) - 1, data = as.data.frame(pdata$K))

  # lmfit <- summary(lmfit)

  return(lmfit)
}

# pdata <- sim_data_panel(aslist = T)
#
# temp <- fit_css_model(pdata$y, X = pdata$X, K = pdata$K, ineff = -1)
#
# print(summary(temp))
