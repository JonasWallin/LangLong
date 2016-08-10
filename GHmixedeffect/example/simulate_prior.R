##
# simple script that simulates from the prior models
#
##


graphics.off()
library(GHmixedeffect)


n.pers <- 2000
n.obs  <- 10

sd_beta <- 0.1
sd_Y    <- 0.1

B_list <- list()
beta <- c(0.9,0.9)
beta_list <- list()
Y_list <- list()
Vin        <- list()
for(i in 1:n.pers)
{
  B_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  Vin[[i]] <- rep(1, n.obs)
}
mError_list <- list(Vs = Vin, noise = "Normal")
mixed_list  <- list(B_random = B_list, Sigma = diag(2), noise = "Normal")
sim_res_Nor <- simulateMixed(mixed_list)

x11()
par(mfrow=c(2, 2))
for(i in 1:2){
  hist(sim_res_Nor$U[i,],50, prob=T)
  x_ <- seq(min(sim_res_Nor$U[i,]), 
          max(sim_res_Nor$U[i,]),
          length=100)
  lines(x_,
      dnorm(x_, sd = mixed_list$Sigma[i,i]),
      col='red')
}
##
# for NIG
##
mixed_list$noise = "NIG"
mixed_list$mu = as.matrix(c(-1, 2.))
mixed_list$nu = 1. 
sim_res_NIG <- simulateMixed(mixed_list)
for(i in 1:2){
  hist(sim_res_NIG$U[i,],50, prob=T)
  x_ <- seq(min(sim_res_NIG$U[i,]), 
            max(sim_res_NIG$U[i,]),
            length=100)
  f <- dnig(x_, -mixed_list$mu[i], mixed_list$mu[i], mixed_list$nu, mixed_list$Sigma[i,i])
  lines(x_, f, col = 'red')
}


###
# measurement error simulate
###
mError_list$sigma <- 1.4
mError_list$Y <- mError_list$Vs
sim_res_Normal <- simulateNoise(mError_list)
Y <- unlist(sim_res_Normal$Y)
x11()
par(mfrow=c(2, 1))
hist(Y,50, prob=T)
x_ <- seq(min(Y), 
          max(Y),
          length=100)
lines(x_,
      dnorm(x_, sd = mError_list$sigma),
      col='red')
##
# nig meas error
##
mError_list$noise = "NIG"
mError_list$nu <- 0.78
sim_res_NIG <- simulateNoise(mError_list)
Y <- unlist(sim_res_NIG$Y)
hist(Y,50, prob=T)
x_ <- seq(min(Y), 
          max(Y),
          length=100)
f <- dnig(x_, 0,0, mError_list$nu, mError_list$sigma)
lines(x_, f, col = 'red')