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
Features
--------

**cross-section models**
- distributions: 
  - normal/half-normal
  - normal/truncated normal
  - normal/exponential
- homoskedasticity and heteroskedasticity for both symmetric and inefficiency terms
- conditional mean (BC, 1995) model

**panel data models**
- panel data conditional mean (BC, 1995) model can be estimated via cross-section model specification
- CSS, 1990 model

Known issues
-----------
- truncated-normal model have to be specified as bc95 model (without passing data in CM) 
- fix efficiency prediction for u computed via E(u_i)

To-do:
------

- more verbose summary output
- panel data models
    - time-invariant model
    - time-varying
      - decay model
      - fixed effects in conditional mean equations in (Battese-Coelli, 1995) model (can be implemented via cross-section model already)
      - fixed effects in conditional variance equations in (Battese-Coelli, 1995) model (can be implemented via cross-section model already)
- robust standard errors
- SFA() funtion with Formula interface 
- add tests
- rewrite likelihood functions in C++ (using Rcpp), altough the speed increase is small based on first experiments
- explore `nloptr`, `ucminf` and `maxLik` packages for MLE optimizations  

Contributions and suggestions are welcome!