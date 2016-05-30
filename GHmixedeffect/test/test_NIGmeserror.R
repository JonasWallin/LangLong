##
# simple test that verifies that the model can correctly idenitfy the parameters, 
# from simulated data for NIG measurement error
# the model is
# Y_{ij} = \simga\sqrt{V}Z, where
# V      \sim NIG(\nu,\nu)  
#
##
rm(list=ls())
graphics.off()
library(GHmixedeffect)
library(rGIG)
library(MASS)
nu <- 3
sigma <- 2.1
n.pers <- 2
n.obs  <- 500

V_list     <- list()
Y_list     <- list()
Vin        <- list()
for(i in 1:n.pers)
{
  V <- rGIG(rep(-0.5,n.obs), rep(nu, n.obs), rep(nu, n.obs))
  V_list[[i]] <- V
  Y_list[[i]] <- sqrt(V) * rnorm(n = n.obs, 0, sd = sigma)
  Vin[[i]] <- rep(1, n.obs)
}

input <- list(Y  = Y_list, 
              Vs = Vin,  
              sigma_eps   = as.matrix(1.),
              nu          = as.matrix(1.),
              Niter       = 5000)
res <- estimateNIGnoise(input)

x11()
par(mfrow=c(2,1))
plot(res$sigmaVec, type='l', col='red')
n_ <- length(res$sigmaVec)
lines(c(1, n_), c(sigma[1],sigma[1]))
lines(c(1, n_), c(sigma[2],sigma[2]))
plot(res$nuVec, type='l', col='red')
lines(c(1, n_), c(nu[1],nu[1]))