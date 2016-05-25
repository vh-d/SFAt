SFAt: Stochastic frontier analysis on cross-section and panel data
==================================================================

This package provides a collection of methods for doing SFA (Stochastic frontier analysis) based on traditional statistical techniques.

It was inspired by some existing software for SFA (R packages frontier, sfa and Benchmarking, STATA commands frontier, xtfrontier and sfcross/sfpanel) but it aims to be more comprehansive, modular/extendable in the future. 

**This package is in alpha version and therefore is not ready for regular use.**

Installation
------------

The package is not on CRAN yet. Install it from its GitHub repo.

```{r}
library(devtools)
install_github("vh-d/SFAt")
```

To-do:
------

- panel data models
    - firm specific intercepts in conditional mean equations in (Battese-Coelli, 1995) model
    - firm specific intercepts in conditional variance equations in (Battese-Coelli, 1995) model
- SFA() funtion with Formula interface 
- add tests
- rewrite likelihood functions in C++ (using Rcpp)
- explore `likelihood`, `ucmin` and `mle` packages for MLE optimizations  


Contributions and suggestions are welcome!