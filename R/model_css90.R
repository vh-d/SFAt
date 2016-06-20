#' Cornwell, Schmidt, and Sickles (1990) model
# as is Parmeter & Kumbhakar 2014, pp.78:
# 1. use within estimator to fit a fixed effect panel data model
# 2. estimate residuals
# 3. regress residuals on time polynomial
# 4. predict \alpha_{it}
# 5. compute inefficiency as a difference between max(\alpha_{i}) - \alpha_{it}
fit_css90_model <- function(y,
                            X,
                            K,
                            ineff,
                            deb = T) {

  # ---- step 1: ----
  require(plm) #this dependence may be removed in the future

  # construct formula for step 1 regression
  formula_frontier <- as.formula(
    paste0(c("y ~ ",
             paste0(colnames(X), collapse = " + "),
             " - 1"),
           collapse = ""))

  if (deb) print(formula_frontier)

  dataframe <- pdata.frame(data.frame(K, y, X), index = colnames(K))

  if (deb) print(head(dataframe))

  #  fit the panel data using within estimator
  step1fit <- plm(formula_frontier,
                  data = dataframe,
                  model = "within")

  alpha_0 <- fixef(step1fit)

  if (deb) {
    cat("======== Step 1 regression: =======", "\n")
    print(summary(step1fit))
    cat("\n")

    cat("======== Fixed effects: =======", "\n")
    print(alpha_0)
    cat("\n")
  }

  # ---- step 2: ----
  # regress residuals againts time polynomial
  step2fit <- lm(formula = res ~ as.factor(k):t + as.factor(k):I(t^2),
                 data    = data.frame(res = step1fit$residuals,
                                      K))

  if (deb) {
    cat("======== Step 2 regression: =======", "\n")
    print(summary(step2fit))
    cat("\n")
  }

  # ---- return ----
  result <- list(
    model = "css90",
    fe = alpha_0,
    tvar = step2fit$coefficients,
    step1fit = step1fit,
    step2fit = step2fit
  )

  return(result)
}

# fit inefficiency (u term) for CSS 1990 model
predict_uit_css90 <- function(x) {
  predict(x$step2fit) + x$fe[attr(x$step1fit$model, "index")$k]
}

# summary statistics for CSS 1990
summary_css90 <- function(x) {
  require(data.table)

  u_it <- data.table(x$data$K,
                     y = x$data$y,
                     x$data$X,
                     uit = predict_uit_css90(x))

  u_it[, ui := max(uit), by = t]
  u_it[, ineff := ui - uit]

  mean_u <- mean(u_it$ineff)
  sigma_u <- sd(u_it$ineff)

  sigma_v <- sd(x$step1fit$residuals)

  ans <- list(
    data = u_it,
    sigma_u = sigma_u,
    sigma_v = sigma_v,
    u = u_it$ineff)

  return(ans)
}


