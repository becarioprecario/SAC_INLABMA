# MCMC and plots


library(INLA)
library(spatialreg)
library(coda)

load("load_data.RData")

if(FALSE) { #DO not re-run
set.seed(1)
res0 <- spBreg_sac(TURNOUT01 ~ 1, data = as.data.frame(turnout),
  listw = it.lw, control = list(ndraw = 100000L, nomit = 10000L, thin = 10L,
  prior = list(nu = 0.01, d0 = 0.01, a1 = 1, a2 = 1, rho = sac0$rho,
    lambda = sac0$lambda,
    Tbeta = diag(1) * 1000)))
summary(res0)
raftery.diag(res0)

set.seed(1)
res1 <- spBreg_sac(TURNOUT01 ~ GDPCAP, data = as.data.frame(turnout),
  listw = it.lw, control = list(ndraw = 100000L, nomit = 10000L, thin = 10L,
  prior = list(nu = 0.01, d0 = 0.01, a1 = 1, a2 = 1, 
    rho = sac1$rho, lambda = sac1$lambda, Tbeta = diag(2) * 1000)))
raftery.diag(res1)
summary(res1)

save(file = "MCMC.RData", list = c("res0", "res1"))
} #DO NOT re-run 

# Plots
# Load results
load("MCMC.RData")



# NO covariates
load("results0.RData")




# Density estimates
pdf(file = "../marginals0.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
# Rho
plot(density(res0[, "rho"], bw = 0.02), xlab = expression(rho),
  ylab = "density", main = "")
lines(density(rholambda$rho, weights = rholambda$w, bw = 0.02), lty = 2)
abline(v = sac0$rho, lty = 3)
legend("topleft", legend = c("MCMC", "BMA", "ML"), lty = 1:3, bty = "n")
# Lambda
plot(density(res0[, "lambda"], bw = 0.03), xlab = expression(lambda),
  ylab = "density", main = "")
lines(density(rholambda$lambda, weights = rholambda$w, bw = 0.03), lty = 2)
abline(v = sac0$lambda, lty = 3)
legend("topleft", legend = c("MCMC", "BMA", "ML"), lty = 1:3, bty = "n")

dev.off()

# Density estimates of other parameters
pdf(file = "../params0.pdf",  width = 10, height = 5)
par(mfrow = c(1, 2))
#beta_0
plot(density(res0[, "(Intercept)"], bw = 0.5), xlab = expression(beta[0]),
  ylab = "density", main = "", ylim = c(0, 4.25))
lines(mergemodel$marginals.fixed[[1]], lty = 2)
abline(v = sac0$coefficients[1], lty = 3)
lines(sacinlaml$marginals.fixed[[1]], lty = 4)
legend("topleft", legend = c("MCMC", "BMA", "ML", "INLA-ML"),
  lty = 1:4, bty = "n")

# Variance
plot(density(res0[, "sige"], bw = 0.08), xlab = expression(tau^{-1}),
  ylab = "density", main = "", ylim = c(0, 1.8))
lines(inla.tmarginal( function (x) 1 / x, 
  mergemodel$marginals.hyperpar[[1]]), lty = 2)
abline(v = sac0$s2, lty = 3)
lines(inla.tmarginal( function (x) 1 / x,
  sacinlaml$marginals.hyperpar[[1]]), lty = 4)
legend("topleft", legend = c("MCMC", "BMA", "ML", "INLA-ML"),
  lty = 1:4, bty = "n")

dev.off()


# Joint posterior of rho and lambda
library(MASS)
library(ggtern) # For kde2d.weighted
z.inla <- kde2d.weighted(rholambda$rho, rholambda$lambda, h = 0.07,
  w = rholambda$w)
z.mcmc <- kde2d(res0[, "rho"], res0[, "lambda"], h = 0.07)

pdf(file = "../joint0.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
contour(z.inla, xlim = c(0.85, 1), ylim = c(-0.2, 0.4), main = "BMA with INLA",
  xlab = expression(rho), ylab = expression(lambda))
points(sac0$rho, sac0$lambda, pch = 19)
contour(z.mcmc, xlim = c(0.85, 1), ylim = c(-0.2, 0.4), main = "MCMC",
  xlab = expression(rho), ylab = expression(lambda))
points(sac0$rho, sac0$lambda, pch = 19)
dev.off()

# WITH covariates
load("results1.RData")

# Density estimates
pdf(file = "../marginals1.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
# Rho
plot(density(res1[, "rho"], bw = 0.02), xlab = expression(rho),
  ylab = "density", main = "", ylim = c(0, 9.5))
lines(density(rholambda$rho, weights = rholambda$w, bw = 0.02), lty = 2)
abline(v = sac1$rho, lty = 3)
legend("topleft", legend = c("MCMC", "BMA", "ML"), lty = 1:3, bty = "n")
# Lambda
plot(density(res1[, "lambda"], bw = 0.035), xlab = expression(lambda),
  ylab = "density", main = "", ylim = c(0, 3.5))
lines(density(rholambda$lambda, weights = rholambda$w, bw = 0.035), lty = 2)
abline(v = sac0$lambda, lty = 3)
legend("topleft", legend = c("MCMC", "BMA", "ML"), lty = 1:3, bty = "n")

dev.off()


# Joint posterior of rho and lambda
library(MASS)
library(ggtern) # For kde2d.weighted
z.inla <- kde2d.weighted(rholambda$rho, rholambda$lambda, h = 0.1,
  w = rholambda$w)
z.mcmc <- kde2d(res1[, "rho"], res1[, "lambda"], h = 0.1)

pdf(file = "../joint1.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
contour(z.inla, xlim = c(0.75, 0.95), ylim = c(-0.1, 0.45), main = "BMA with INLA",
  xlab = expression(rho), ylab = expression(lambda))
points(sac1$rho, sac1$lambda, pch = 19)
contour(z.mcmc, xlim = c(0.75, 0.95), ylim = c(-0.1, 0.45), main = "MCMC",
  xlab = expression(rho), ylab = expression(lambda))
points(sac1$rho, sac1$lambda, pch = 19)
dev.off()


# Density estimates of other parameters
pdf(file = "../params1.pdf",  width = 10, height = 5)
par(mfrow = c(1, 3))
#beta_0
plot(density(res1[, "(Intercept)"], bw = 0.8), xlab = expression(beta[0]),
  ylab = "density", main = "", ylim = c(0, 1.15), xlim = c(0, 20))
lines(mergemodel$marginals.fixed[[1]], lty = 2)
abline(v = sac1$coefficients[1], lty = 3)
lines(sacinlaml$marginals.fixed[[1]], lty = 4)
legend("topright", legend = c("MCMC", "BMA", "ML", "INLA-ML"),
  lty = 1:4, bty = "n")

#beta_1
plot(density(res1[, "GDPCAP"], bw = 0.005), xlab = expression(beta[1]),
  ylab = "density", main = "", ylim = c(0, 40))
lines(mergemodel$marginals.fixed[[2]], lty = 2)
abline(v = sac1$coefficients[2], lty = 3)
lines(sacinlaml$marginals.fixed[[2]], lty = 4)
legend("topright", legend = c("MCMC", "BMA", "ML", "INLA-ML"),
  lty = 1:4, bty = "n")


# Variance
plot(density(res1[, "sige"], bw = 0.08), xlab = expression(tau^{-1}),
  ylab = "density", main = "", ylim = c(0, 1.8))
lines(inla.tmarginal( function (x) 1 / x, 
  mergemodel$marginals.hyperpar[[1]]), lty = 2)
abline(v = sac1$s2, lty = 3)
lines(inla.tmarginal( function (x) 1 / x,
  sacinlaml$marginals.hyperpar[[1]]), lty = 4)
legend("topright", legend = c("MCMC", "BMA", "ML", "INLA-ML"),
  lty = 1:4, bty = "n")

dev.off()



# Impacts
# IMPORTANT!!!! No impacts for res0 

# Required to compute theimpacts
trMatc <- trW(W, type="mult")

# ML
impacts.ML <- impacts(sac1, tr = trMatc)
impacts.ML


# INLA - ML
totimp.marg <- inla.tmarginal(function (x) 1 / (1 - sac1$rho) * x,
  sacinlaml$marginals.fixed[[2]])
inla.zmarginal(totimp.marg)

# Average of trace
IrhoW.inv <- solve(Diagonal(nrow(turnout), 1) - sac1$rho * W)
avgtraceW <- mean(diag(IrhoW.inv))
dirimp.marg <- inla.tmarginal(function(x)  {avgtraceW * x},
  sacinlaml$marginals.fixed[[2]])
inla.zmarginal(dirimp.marg)

# Sum of off diagonal elements / n
avgoffW <- sum(IrhoW.inv) / nrow(turnout) - avgtraceW 
idirimp.marg <- inla.tmarginal(function(x)  {avgoffW * x},
  sacinlaml$marginals.fixed[[2]])
inla.zmarginal(idirimp.marg)


# MCMC
impacts.MCMC <- impacts(res1, tr = trMatc)
summary(impacts.MCMC, short = TRUE, zstats = TRUE)


# BMA with INLA
library(INLABMA)
# Direct impacts
margs <-  mclapply(1:nrow(rholambda),  function(i) {
    IrhoW.inv <- solve(Diagonal(nrow(turnout), 1) - rholambda$rho[i] * W)
    avgtraceW <- mean(diag(IrhoW.inv))
    aux <- inla.tmarginal(function(x)  {avgtraceW * x}, margs.impacts[[i]][[2]])
    return(aux)
})
dirimp.margBMA  <- INLABMA:::fitmargBMA(margs, rholambda$w)
inla.zmarginal(dirimp.margBMA)

# INdirect impacts
margs <-  mclapply(1:nrow(rholambda),  function(i) {
    IrhoW.inv <- solve(Diagonal(nrow(turnout), 1) - rholambda$rho[i] * W)
    avgoffW <- (sum(IrhoW.inv) - sum(diag(IrhoW.inv))) / nrow(turnout)
    aux <- inla.tmarginal(function(x)  {avgoffW * x}, margs.impacts[[i]][[2]])
    return(aux)
})
idirimp.margBMA  <- INLABMA:::fitmargBMA(margs, rholambda$w)
inla.zmarginal(idirimp.margBMA)


# TOTAL impacts
margs <-  mclapply(1:nrow(rholambda),  function(i) {
    aux <- inla.tmarginal(function(x)  {x / (1 - rholambda$rho[i])},
      margs.impacts[[i]][[2]])
    return(aux)
})
totimp.margBMA  <- INLABMA:::fitmargBMA(margs, rholambda$w)
inla.zmarginal(totimp.margBMA)


# PLot of impacts post. marginals
# Density estimates of other parameters
pdf(file = "../impacts.pdf",  width = 10, height = 5)
par(mfrow = c(1, 3))

#Direct impacts
plot(density(impacts.MCMC$sres$direct[, 1], bw = 0.0055), xlab = "GDPCAP",
  ylab = "density", main = "", ylim = c(0, 25), xlim = c(0, 0.20))
lines(dirimp.margBMA, lty = 2)
abline(v = impacts.ML$direct, lty = 3)
lines(dirimp.marg, lty = 4)
legend("topright", legend = c("MCMC", "BMA", "ML", "INLA-ML"),
  lty = 1:4, bty = "n")

# Indirect impacts
plot(density(impacts.MCMC$sres$indirect[, 1], bw = 0.02), 
  xlab = "GDPCAP",
  ylab = "density", main = "", ylim = c(0, 5.5))
lines(idirimp.margBMA, lty = 2)
abline(v = impacts.ML$indirect, lty = 3)
lines(idirimp.marg, lty = 4)
legend("topright", legend = c("MCMC", "BMA", "ML", "INLA-ML"),
  lty = 1:4, bty = "n")

# Total impacts
plot(density(impacts.MCMC$sres$total[, 1], bw = 0.02), 
  xlab = "GDPCAP",
  ylab = "density", main = "", ylim = c(0, 5.5))
lines(totimp.margBMA, lty = 2)
abline(v = impacts.ML$total, lty = 3)
lines(totimp.marg, lty = 4)
legend("topright", legend = c("MCMC", "BMA", "ML", "INLA-ML"),
  lty = 1:4, bty = "n")


dev.off()



if(FALSE) {
# MCMC
library(SEMCMC)
# Trick to estimate model
turnout$zero <- 0
mcmc0 <- SEMCMC(TURNOUT01 ~ 1 + zero, as.data.frame(turnout), W = as.matrix(W),
  model = "sac",
  n.burnin = 2000, n.iter = 1000, n.thin = 10)

mcmc1 <- SEMCMC(form1, as.data.frame(turnout), W = as.matrix(W),
  model = "sac",
  n.burnin = 2000, n.iter = 1000, n.thin = 10)
} # FALSE

