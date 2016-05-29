SFAt: Stochastic frontier analysis on cross-section and panel data
==================================================================

This package provides a collection of methods for doing SFA (Stochastic frontier analysis) based on traditional statistical techniques.

It was inspired by some existing software for SFA (R packages *frontier*, *sfa* and *Benchmarking*, *Stata* commands `frontier`, `xtfrontier` and `sfcross`/`sfpanel`) but it aims to be more comprehansive, modular/extendable in the future. 

**This package is in beta version. It is functional but the interface is under development. Therefor it is not ready for regular use.**

Installation
------------

The package is not on CRAN yet. Install it from its GitHub repo.

```{r}
library(devtools)
install_github("vh-d/SFAt")
```

To-do:
------

- more verbose summary output
- fix efficiency prediction for u computed via E(u_i)
- panel data models
    - time-invariant model
    - time-varying
      - decay model
      - fixed effects in conditional mean equations in (Battese-Coelli, 1995) model
      - fixed effects in conditional variance equations in (Battese-Coelli, 1995) model
- robust standard errors
- SFA() funtion with Formula interface 
- add tests
- rewrite likelihood functions in C++ (using Rcpp)
- explore `nloptr`, `ucminf` and `maxLik` packages for MLE optimizations  

Contributions and suggestions are welcome!