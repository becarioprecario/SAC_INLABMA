# Analysis of turnout data in Italy

library(rgdal)
library(sp)
library(RColorBrewer)
library(spdep)
# Models with INLABMA
library(INLA)
library(INLABMA)

# Load functions and stuff
source("utils.R")

# LOad data
load("load_data.RData")


#Index for spatial random effects
turnout$idx <- 1:nrow(turnout)

# Do not use Gaussian error
zero.variance <- list(prec = list(initial = 15, fixed = TRUE))


# Number of grid points in each dimension
n.xy <-5 
model.to.run <- 0 #0: No covs; 1: with covs

# Formulas
if(!model.to.run) {
  form <- TURNOUT01 ~ 1
  sac <- sac0
  pdffile <- "ww0.pdf"
  savefile <- "results0.RData"
} else {
  form <- TURNOUT01 ~ 1 + GDPCAP
  sac <- sac1
  pdffile <- "ww1.pdf"
  savefile <- "results1.RData"
}

# Grid of points
#xx.rho <- seq(0.85, 0.95, length.out = 5) #by =  0.005)
#xx.lambda <- seq(0.00, 0.3, length.out = 10) #by = 0.005)
# Using ML estimates
xx.rho <- seq(sac$rho - 3 * sac$rho.se, sac$rho + 3 * sac$rho.se,
  length.out = n.xy)
xx.lambda <- seq(sac$lambda - 3 * sac$lambda.se,
  sac$lambda + 3 * sac$lambda.se, length.out = n.xy)

#rl.step <- c(0.001, 0.001)
#xx.rho <- seq(0.88, 0.90, by =  rl.step[1])
#xx.lambda <- seq(0.29, 0.31, by = rl.step[2])
rholambda <- expand.grid(rho = xx.rho, lambda = xx.lambda)

#data.frame(rho = rep(xx.rho, length(xx.lambda)),
#  lambda = rep(xx.lambda, each = length(xx.rho))
#)

# VAlues
rholambda


# No covariates
sacinla <- mclapply(1:nrow(rholambda), function(idx){
 
  rl <- rholambda[idx, ]
  print(idx)

  res <- try(sac.inla(form, d = as.data.frame(turnout), W.rho = W, W.lambda = W,
    fhyper = list(prec = list(param = c(0.01, 0.01))),
    rho = rl$rho,
    lambda = rl$lambda,
    family = "gaussian", impacts = FALSE,
    control.fixed = list(prec.intercept = 0.001),
    control.family = list(hyper = zero.variance),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, cpo = TRUE),
    control.inla = list(print.joint.hyper = TRUE, #), 
      strategy = "laplace", tolerance = 1e-10, h = 0.001),
    # Start close to ML estimate
    control.mode = list(theta = log(0.2), restart = TRUE),
    improve = FALSE,
    verbose = FALSE
  ))

  logdet <- res$logdet
  res <- try(inla.rerun(res))

  res$logdet <- logdet
  res$mlik <- res$mlik + 0.5 * logdet

  return(res)

})

rholambda$mliks <- get.mliks(sacinla)

rholambda$w <- reweight(rholambda$mliks)


# Figures of weights and mliks
pdf(file = pdffile, width = 8, height = 5)
par(mfrow = c(1, 2))
hist(rholambda$w, xlab = "weight", main = NULL)
hist(rholambda$mliks, xlab = "marginal log-likelihood", main = NULL)
dev.off()

posterior(rholambda$rho, rholambda$w)
posterior(rholambda$lambda, rholambda$w)


# Merge models
idx <- 1:nrow(rholambda) # which(!is.na(rholambda$w) & rholambda$w > 0.0000001)

mergemodel <- inla.merge(sacinla[idx], rholambda$w[idx])

summary(mergemodel)

# Transform marginals of precision
inla.zmarginal(inla.tmarginal(function (x) exp(-x),
  mergemodel$internal.marginals.hyperpar[[1]]))

sacinlaml <- sac.inla(form, d = as.data.frame(turnout), W.rho = W, W.lambda = W,
    fhyper = list(prec = list(param = c(0.01, 0.01))),
    rho = sac$rho,
    lambda = sac$lambda,
    family = "gaussian", impacts = FALSE,
    control.fixed = list(prec.intercept = 0.001),
    control.family = list(hyper = zero.variance),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, cpo = TRUE),
    control.inla = list(print.joint.hyper = TRUE,
      strategy = "laplace", tolerance = 1e-10, h = 0.001),
    improve = FALSE,
    verbose = FALSE
  )

logdet <- sacinlaml$logdet
sacinlaml <- inla.rerun(sacinlaml)
sacinlaml$logdet <- logdet
sacinlaml$mlik <- sacinlaml$mlik + 0.5 * logdet
summary(sacinlaml)
inla.zmarginal(
  inla.tmarginal(function(x) exp(-x),
    sacinlaml$internal.marginals.hyperpar[[1]])
)

rholambda

par(mfrow = c(1, 2))
plot(xx.rho, by(rholambda$mliks, rholambda$rho, mean),
  xlab = expression(rho), ylab  = "marginal likelihood")
plot(xx.lambda, by(rholambda$mliks, rholambda$lambda, mean),
  xlab = expression(lambda), ylab = "marginal likelihood")

save(file = savefile, list = c("rholambda", "mergemodel", "sacinlaml"))

