# Bayesian Model Averaging with INLA

This repository contains the code for the following paper:

* V. Gómez-Rubio, R. S. Bivand and H. Rue (2019). *Bayesian model averaging with the integrated nested Laplace approximation*.

This paper describes how to do Bayesian model averaging with INLA, which
is illustated using the the SAC spatial econometrics model.

These are the files provided:

* `loaddata.RData`, turnout in the 2011 ellection in Italy dataset plus adjacency matrix.

* `utils.R`, some configuration options and functiions required by the main file.

* `analysis.R`, main `R` file (i.e., the one you want to run). `n.xy` controls 
  the number of grid points in each dimension and `model.to.run` the model
  to run (`0` for no covariates and `1` for including covariates).

The number of models fit with INLA is `n.xy * n.xy` so set this parameter
depending on the computer you are using. A value of `20` seems to work
on recent laptops, but set it to `5` for testing. `mclapply` is
used for parallel computing, which may not work on some computers. 
You can set the number of cores and the number of threads used by each 
INLA proces in `utils.R`.


**NOTE:** The original source of the data is Ward and Gleditsch (2008)
and it is available from [http://ksgleditsch.com/srm_book.html](http://ksgleditsch.com/srm_book.html).

**References**

Michael Ward and Kristian Skrede Gleditsch. 2008. *Spatial Regression Models*. Thousand Oaks, CA: Sage.
