SFAt: Stochastic frontier analysis on cross-section and panel data
==================================================================

This package provides a collection of methods for doing SFA (Stochastic frontier analysis) based on traditional statistical techniques.

Its development was inspired by some existing software for SFA (R packages *frontier*, *sfa* and *Benchmarking*, *Stata* commands `frontier`, `xtfrontier` and `sfcross`/`sfpanel`).

**This package is in beta version. It is mostly functional but the interface and some features are still under development. Therefore it is not ready for regular use yet.**

Contributions and suggestions are welcome!

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

To-do:
------

- [ ] more verbose summary output
- panel data models
    - [ ] time-invariant model
    - [ ] time decay model
    - [x] fixed effects in conditional mean equations in (Battese-Coelli, 1995) model (can be implemented via cross-section model already)
    - [x] fixed effects in conditional variance equations in (Battese-Coelli, 1995) model (can be implemented via cross-section model already)
- [ ] explore `optimx` package for MLE optimization  
- [x] add analytic gradient functions
- [ ] analytic hessians
- [ ] LR test for two SFA objects
- [ ] variance of truncated normal distribution