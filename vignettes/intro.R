## ----Init----------------------------------------------------------------
require(SFAt)
# require(msm) # package for drawing from truncated-normal distribution
# require(data.table) # for more convenient data manipulation

## ------------------------------------------------------------------------
data_ex1 <- as.data.frame(sim_data_cs(1000, 
                                      x_coeff = c(50, 4, 2), 
                                      z_coeff = c(20), 
                                      sigma_u = 20, sigma_v = 5))


## ------------------------------------------------------------------------
# require(GGally)
# ggpairs(dtDT)

plot(data_ex1)
hist(data_ex1$y)
hist(data_ex1$u)
hist(data_ex1$eps)


## ------------------------------------------------------------------------

model_ex1 <- sfa.fit(y = data_ex1$y,
                     X = cbind(a = data_ex1$x1, b = data_ex1$x2),
                     CM = NULL,    
                     dist = "tnorm",
                     spec  = NULL,
                     deb = F,
                     opt_method = "L-BF"
                     # control_opt = list(maxit = 300)
)

print(summary(model_ex1))

## ------------------------------------------------------------------------
ineff_fit <- inefficiencyTerm(model_ex1, estimator = "JLMS")
hist(ineff_fit)

## ------------------------------------------------------------------------
pr <- predict(model_ex1, type = "frontier")
hist(pr, xlim = c(0, max(data_ex1$y)))
hist(data_ex1$y, xlim = c(0, max(data_ex1$y)))

## ------------------------------------------------------------------------
eff_scaled <- (pr - ineff_fit)/pr
hist(eff_scaled)

## ------------------------------------------------------------------------
ineff_scaled <- 1-eff_scaled
hist(ineff_scaled)

## ------------------------------------------------------------------------
hist(efficiency(model_ex1))

## ------------------------------------------------------------------------
data_ex2_cm <- as.data.frame(sim_data_cs(1000, 
                                         x_coeff = c(50, 4, 2), 
                                         z_coeff = c(20, -2, 2), 
                                         sigma_u = 10, sigma_v = 5))


## ------------------------------------------------------------------------

model_ex2 <- sfa.fit(y = data_ex2_cm$y,
                     X = cbind(a = data_ex2_cm$x1, b = data_ex2_cm$x2),
                     CM = cbind(c = data_ex2_cm$z1, d = data_ex2_cm$z2),
                     dist = "tnorm",
                     spec  = "bc95",
                     deb = F,
                     opt_method = "SANN",
                     control_opt = list(maxit = 6000))

print(model_ex2$coefficients)
print(model_ex2$coefficients_Z)
print(summary(model_ex2))

## ------------------------------------------------------------------------
ex2_ineff_term <- predict(model_ex2, 
                          type = "inefficiency", 
                          estimator = "JLMS")

hist(ex2_ineff_term)

## ------------------------------------------------------------------------
ex2_efficiencies <- predict(model_ex2, 
                            type = "efficiency", 
                            estimator = "JLMS")

hist(ex2_efficiencies)

