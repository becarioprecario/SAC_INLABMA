# Analysis of turnout data in Italy

library(rgdal)
library(sp)
library(RColorBrewer)
library(spdep)

# LOad data
load("load_data.RData")

# Models with INLABMA
library(INLA)
library(INLABMA)

# Fit models
library(parallel)
options(mc.cores = 4)

# Set hthis for parallel computing
inla.setOption(num.threads = 2)

# Get mliks
get.mliks <- function(objlist) {
  as.vector(unlist(mclapply(objlist, function(X) {
   res <- try(X$mlik[1, 1]) #+ 2 * log(0.5))# Add 2 * log-prior
    if(!is.numeric(res))
      res <- NA
    res
  })))
}

# Re-weight
reweight <- function(xx) {
  ww <- exp(xx - max(xx, na.rm = TRUE))
  ww <- ww / sum(ww, na.rm = TRUE)
}

# Statistics
posterior <- function(xx, ww) {
  p.mean <- sum(xx * ww, na.rm = TRUE)

  p.var <- sum(xx^2 * ww, na.rm = TRUE) - p.mean^2

  return(list(mean = p.mean, p.var = p.var, p.sqrt = sqrt(p.var)))
}


