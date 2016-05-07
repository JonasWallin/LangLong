##
# simple test that verifies that the model can correctly idenitfy the parameters, 
# from simulated data
#
##
rm(list=ls())
graphics.off()
library(GHmixedeffect)
library(rGIG)
library(MASS)
n.pers <- 1000 #number of patients
n.obs  <- 10 #number of obs per patient

COV_beta <- matrix(c(0.2,0.1,0.1,0.2), ncol = 2, nrow = 2)
sd_Y    <- 0.1 # error of the noise

Br_list <- list()
betar <- c(0.9,0.4)
betaf <- c(1.)
mu   <- c(0.2, -0.2)
nu <- 0.5
betar_list <- list()
Bf_list    <- list()
V_list     <- list()
Y_list     <- list()
for(i in 1:n.pers)
{
  Bf_list[[i]]    <- as.matrix(runif(n = n.obs))
  Br_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  V <- rGIG(-0.5, nu, nu)
  V_list[[i]] <- V
  betar_list[[i]] <- betar - mu * 1  + V * mu +
                     sqrt(V) * mvrnorm(n = 1, mu  =c(0,0), Sigma = COV_beta)
  Y_list[[i]]        <- rnorm(n = n.obs,
                              Br_list[[i]]%*%betar_list[[i]]
                            + Bf_list[[i]]%*%betaf, sd = sd_Y)

}

input <- list(Y = Y_list, 
              B_random    = Br_list, 
              B_fixed     = Bf_list, 
              Sigma       = COV_beta, 
              beta_random = c(0,0.), 
              beta_fixed  = c(0.), 
              sigma_eps   = sd_Y,
              mu          = 0*as.matrix(mu),
              nu          = as.matrix(1.),
              Niter       = 5000,
              noise       = "NIG")

beta_mat <- t(matrix(unlist(betar_list), nrow= 2, ncol = n.pers))
x11()
par(mfrow=c(2,1))
hist(beta_mat[,1],100)
hist(beta_mat[,2],100)
res <- estimateME(input)
print(res$mixedeffect$Sigma)
print(res$mixedeffect$beta_random)
print(res$mixedeffect$mu)
print(res$mixedeffect$nu)
print(res$mixedeffect$beta_fixed)
x11()
par(mfrow=c(3,1))
plot(res$mu[,1], type='l', col='red', ylim=c(min(res$mu), max(res$mu)))
lines(res$mu[,2], col='red')
n_ <- length(res$mu[,1])
lines(c(1, n_), c(mu[1],mu[1]))
lines(c(1, n_), c(mu[2],mu[2]))
plot(res$betaVec[,1], type='l', col='red', ylim=c(min(res$betaVec), max(res$betaVec)))
lines(res$betaVec[,2], col='red')
n_ <- length(res$mu[,1])
lines(c(1, n_), c(betar[1], betar[1]))
lines(c(1, n_), c(betar[2], betar[2]))

plot(res$nuVec, type='l', col='red' )
n_ <- length(res$mu[,1])
lines(c(1, n_), c(nu, nu))
