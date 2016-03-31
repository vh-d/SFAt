# source("./R/simulate_data.R")

require(msm)
require(SFAplus)


require(data.table)

dtDT <- as.data.table(SFAplus::sim_data(1000))

# require(GGally)
# ggpairs(dtDT)
plot(dtDT)
hist(dtDT$y)
hist(dtDT$u)
# FIT MODEL ---------------------------------------------------------------

est <- SFAplus::sfa.fit(y = dtDT$y,
                        X = as.matrix(dtDT[, .(1, x1, x2)]),
                        Z = as.matrix(dtDT[, .(1, z1, z2)]),
                        intercept = F,
                        intercept_Z = F,
                        model  = "bc95",
                        deb = F,
                        opt_method = "BFGS",
                        # control_opt = list(maxit = 300),
                        trace = T)

print(res)

# OLS FIT ----------------------------------------------------------
# for LR test and starting values of parameters


lmfit <- lm(y ~ x1 + x2 + z1 + z2, data = dtDT)
summary(lmfit)



# BOUNDS ------------------------------------------------------------------
# set lower bounds, important for sigmas
lowerb <- c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0.0001, 0.0001)

# initpar = c(lmfit$coefficients[1:4], 0, lmfit$coefficients[5:6], 1, 1)
# OR
initpar = c(rep(0, 4+3), 1, 1)


# estimate by BFGS method with bounds
est <- optim(initpar, fn = ll, method = "BFGS", control = list(trace = F, REPORT = 1, maxit = 300), hessian = T,
             nbeta = 4, ngamma = 3,
             y = dtDT$y,
             X = as.matrix(dtDT[, .(1, x1, x2, x1x2)]),
             Z = as.matrix(dtDT[, .(1, z1, z2)]))
print(est)


# SUMMARY -----------------------------------------------------------------

# compute statistics
var_beta <- abs(diag(solve(-est$hess)))
tvalue <- est$par/sqrt(var_beta)
coef.table <- cbind(est$par, sqrt(var_beta), tvalue)
colnames(coef.table) <- c("Estimate", "Std. Error","t value")
row.names(coef.table) <- names(est$par)
ll_sfa <- -est$value
ll_ols <- logLik(lmfit)
LRtest <- 2*(ll_sfa - ll_ols)[1]
chisq_df <- (4+3+2) - attributes(logLik(lmfit))$df
p_value <- pchisq(LRtest, chisq_df, lower.tail = FALSE)



# benchmark with sfa package ----------------------------------------------
require(sfa)
summary(sfa::sfa(y ~ x1 + x2 + x1x2 + z1 + z2, data= dtDT))


# BENCHMARK WITH FRONTIER PACKAGE methods -----------------------------------------
require(frontier)

frontest <- frontier::sfa(y ~ x1 + x2 + x1x2 | z1 + z2,
                          data = dtDT,
                          ineffDecrease = T,
                          muBound = 0,
                          maxit = 100)

summary(frontest)
lrtest(frontest)
frontier:::lrtest.frontier()
