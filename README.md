SFAt: Stochastic frontier analysis on cross-section and panel data
==================================================================

*This package is in alpha version and therefore is not ready for regular use.*

This package provides a collection of methods for doing SFA (Stochastic frontier analysis) based on traditional statistical techniques.

It was inspired by some existing software for SFA (R packages frontier, sfa and Benchmarking, STATA commands frontier, xtfrontier and sfcross/sfpanel) but it aims to be
more comprehansive, modular/extendable and open (open source).

Contributors are welcome!

to-do:
------

- panel data models
    - firm specific intercepts in conditional mean equations in (Battese-Coelli, 1995) model
    - firm specific intercepts in conditional variance equations in (Battese-Coelli, 1995) model
- sfa() funtion with Formula interface 
- rewrite likelihood functions in C++ (using Rcpp)
- explore *likelihood*, *ucminf* and *mle2* packages for optimizations in MLE 