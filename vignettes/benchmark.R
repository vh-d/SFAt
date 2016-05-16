## ----Init----------------------------------------------------------------
require(SFAt)
# require(msm) # package for drawing from truncated-normal distribution
require(data.table) # for more convenient data manipulation

## ------------------------------------------------------------------------

dtDT <- as.data.table(sim_data_cs(500))


## ---- fig.width = 10, fig.height = 10------------------------------------
# require(GGally)
# ggpairs(dtDT)

plot(dtDT)

## ------------------------------------------------------------------------
hist(dtDT$y)
hist(dtDT$u)
hist(dtDT$eps)


## ------------------------------------------------------------------------

lmfit <- lm(y ~ x1 + x2 + z1 + z2, data = dtDT)
summary(lmfit)


## ------------------------------------------------------------------------
sfat <- SFAt::sfa.fit(y = dtDT$y,
                     X = as.matrix(dtDT[, .(x1, x2)]),
                     CM = as.matrix(dtDT[, .(z1, z2)]),
                     spec = "bc95",
                     deb = F,
                     opt_method = "BFGS"
                     # control_opt = list(maxit = 300, trace = T)
)
summary(sfat)

## ------------------------------------------------------------------------
require(frontier)

frontest <- frontier::sfa(y ~ x1 + x2 | z1 + z2,
                          data = dtDT,
                          ineffDecrease = T,
                          truncNorm = T,
                          muBound = 0,
                          maxit = 100)

print(summary(frontest))
print(lrtest(frontest))


## ------------------------------------------------------------------------
data("front41Data")

sfa2 <- SFAt::sfa.fit(y = log(front41Data$output), 
                         X = cbind(capital = log(front41Data$capital), 
                                   labour = log(front41Data$labour)), 
                         dist = "tnorm",
                         spec = NULL,
                         deb = F)

print(summary(sfa2))

## ------------------------------------------------------------------------
sfa2_front <- frontier::sfa( log( output ) ~ log( capital ) + log( labour ),
                             truncNorm = T,
                             data = front41Data )

print(summary(sfa2_front))
print(lrtest(sfa2_front))

