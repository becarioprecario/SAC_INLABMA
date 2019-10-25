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
n.xy <-  20
model.to.run <- 1 #0: No covs; 1: with covs
create.fig <- FALSE # Create plot of mliks and weights?
save.results <- TRUE

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

# Using 'internal scale': log( (1 + x) / (1 - x))
#  Variance in internal scale is: serr(x)^2 * (2 / [(1 + x) (1 - x)])^2
#    using Delta method.
f.trans <- function(x) {
  log( (1 + x) / (1 - x))
}

f.trans.inv <- function(y) {
  2 / (1 + exp (-y)) - 1
}

err.trans <- function(x, se.x) {
  sqrt( (se.x^2) * (2 / ( (1 + x) * (1 - x)))^2 )
}

xx.rho <- f.trans.inv(
  seq(f.trans(sac$rho) - 3 * err.trans(sac$rho, sac$rho.se),
    f.trans(sac$rho) + 3 * err.trans(sac$rho, sac$rho.se),
    length.out = n.xy)
)

xx.lambda <- f.trans.inv(
  seq(f.trans(sac$lambda) - 3 * err.trans(sac$lambda, sac$lambda.se),
    f.trans(sac$lambda) + 3 * err.trans(sac$lambda, sac$lambda.se),
    length.out = n.xy)
)


rholambda <- expand.grid(rho = xx.rho, lambda = xx.lambda)

# Log-prior in internal scale: lg-unif + change of scale
#  log(1 / 2) + log( 2 * exp(theta) / (1 + exp(theta))^2)

# Internal scale
theta1 <- f.trans(rholambda$rho)
theta2 <- f.trans(rholambda$lambda)

#rholambda$logprior <- -log(2) + log(2) + theta1 - 2 * log(1 + exp(theta1)) +
#  -log(2) + log(2) + theta2 - 2 * log(1 + exp(theta2))
# Simplfied expression
rholambda$logprior <- theta1 - 2 * log(1 + exp(theta1)) + 
  theta2 - 2 * log(1 + exp(theta2))

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

# Marginal log-likelihood
rholambda$mliks <- get.mliks(sacinla)

# Weights by re-scaling mliks + log-prior
rholambda$w <- reweight(rholambda$mliks + rholambda$logprior)


# Figures of weights and mliks
if(create.fig) {
  pdf(file = pdffile, width = 8, height = 5)
  par(mfrow = c(1, 2))
  hist(rholambda$w, xlab = "weight", main = NULL)
  hist(rholambda$mliks, xlab = "marginal log-likelihood", main = NULL)
  dev.off()
}
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


# Save marginals to compute impacts
if(model.to.run) {
  margs.impacts <- lapply(sacinla, function (X) {
    X$marginals.fixed
  })
} else {
  margs.impacts <- NA
}


if(save.results) 
  save(file = savefile, list = c("rholambda", "mergemodel", "sacinlaml",
    "margs.impacts"))

