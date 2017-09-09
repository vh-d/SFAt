## ----Init----------------------------------------------------------------
require(SFAt)

## ----Simulate data-------------------------------------------------------

EX1 <- as.data.frame(sim_data_cs(500))


## ----Scatter plot matrix, fig.width = 10, fig.height = 10----------------
# require(GGally)
# ggpairs(dtDT)

plot(EX1)

## ----Plot individual variables-------------------------------------------
hist(EX1$y)
hist(EX1$u)
hist(EX1$eps)


## ----OSL model-----------------------------------------------------------

lmfit <- lm(y ~ x1 + x2 + z1 + z2, data = EX1)
summary(lmfit)


## ----SFAt model----------------------------------------------------------
sfat <- SFAt::sfa.fit(y = EX1$y,
                     X = as.matrix(EX1[, c("x1", "x2")]),
                     CM = as.matrix(EX1[, c("z1", "z2")]),
                     dist = "tnorm",
                     deb = F,
                     optim_method = "BFGS",
                     optim_control = list(maxit = 3e4, trace = T)
)
summary(sfat)

## ----frontier package----------------------------------------------------
require(frontier)

frontest <- frontier::sfa(y ~ x1 + x2 | z1 + z2,
                          data = EX1,
                          ineffDecrease = T)

print(summary(frontest))
print(lrtest(frontest))


## ----SFAt BC95-----------------------------------------------------------
data("front41Data")

sfa2 <- SFAt::sfa.fit(y = log(front41Data$output), 
                      X = cbind(capital = log(front41Data$capital), 
                                labour = log(front41Data$labour)), 
                      dist = "tnorm",
                      optim_method = "CG",
                      optim_control = list(trace = T),
                      deb = F)

print(summary(sfa2))

## ----frontier package BC95-----------------------------------------------
sfa2_front <- frontier::sfa( log( output ) ~ log( capital ) + log( labour ),
                             truncNorm = T,
                             data = front41Data )

print(summary(sfa2_front))
print(lrtest(sfa2_front))

