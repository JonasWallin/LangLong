#####
# test function for estimating parameters for stochastic proccesses mixed effect
# #TODO: add a small regularization on the random effect otherwise it can be come singular
#
####
graphics.off()
library(LANGlong)
library(methods)

nIter <- 5000
n.pers <- 100
nSim  <- 1
n.obs  <- 1000
n <- 100

use.matern <- TRUE
kappa.true <- 2
tau.true <- 15
###
# setup of objects
#
###

Y <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()
Vin <- list()
for(i in 1:n.pers)
{
  B_random[[i]] <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  B_fixed[[i]]  <- as.matrix(rnorm(n = n.obs, 0,sd=2))
  Y[[i]] <- rep(1,n.obs)
  locs[[i]] <- seq(from=0,to=10,length.out=n.obs)
  Vin[[i]] <- rep(1, n.obs)
}
mError_list <- list(Vs = Vin, noise = "NIG", sigma = 0.1, nu = 1)
mixedEffect_list  <- list(B_random = B_random,
                          B_fixed  = B_fixed,
                          beta_random = as.matrix(c(2,-1)),
                          beta_fixed  = as.matrix(c(.1)),
                          Sigma = diag(c(0.1, 0.2)),
                          noise = "Normal",
                          Sigma_epsilon=1)

if(use.matern){
  operator_list <- create_operator(locs, n, name = "Matern")
  operator_list$type  <- "Matern"
  operator_list$kappa <- kappa.true
  operator_list$tau   <- tau.true
} else {
  operator_list <- create_operator(locs, n, name = "fd2")
  operator_list$type  <- "fd2"
  operator_list$tau   <- tau.true
}

processes_list = list(noise = "Normal",
                      nu  = 0.,
                      mu  = 0.)
processes_list$V <- list()
for(i in 1:length(locs))
{
  processes_list$V[[i]] <- operator_list$h
}
###
# simulation
#
###
mixedEffect_list_in = mixedEffect_list

sim_res <- simulateLongPrior( Y                 = Y,
                              locs              = locs,
                              mixedEffect_list  = mixedEffect_list_in,
                              measurment_list   = mError_list,
                              processes_list    = processes_list,
                              operator_list     = operator_list)

if(use.matern){
  operator_list$kappa <- kappa.true*2
  operator_list$tau   <- tau.true/3
} else {
  operator_list$tau   <- tau.true/2
}

#mixedEffect_list$beta_random = as.matrix(c(0.,0.))
#mixedEffect_list$Sigma = diag(c(1,1.))
#mError_list$sigma = as.matrix(c(1.))
processes_list$X <- sim_res$X
if(1){
res <- estimateLong(Y                = sim_res$Y,
                    nIter            = nIter,
                    nSim             = nSim,
                    locs             = locs,
                    mixedEffect_list = mixedEffect_list,
                    measurment_list  = mError_list,
                    processes_list   = processes_list,
                    operator_list    = operator_list)

if(use.matern){
  par(mfrow=c(1,2))
  plot(res$operator_list$tauVec,col=2,type="l")
  lines(rep(tau.true,nIter))
  plot(res$operator_list$kappaVec,col=2,type="l")
  lines(rep(kappa.true,nIter))

} else {
  plot(res$operator_list$tauVec,col=2,type="l")
  lines(rep(tau.true,nIter))
  }
}
x11()

x_t <-  res$obs_list[[1]]$A%*%sim_res$X[[1]] + res$mixedEffect_list$Br[[1]]%*%(sim_res$U[,1] + mixedEffect_list$beta_random)
x_t <- x_t + res$mixedEffect_list$Bf[[1]]%*%mixedEffect_list$beta_fixed
x_t <- as.matrix(x_t)
x_ <- res$obs_list[[1]]$A%*%res$X[[1]] + res$mixedEffect_list$Br[[1]]%*%(res$mixedEffect_list$U[,1] + res$mixedEffect_list$beta_random)
x_ <- x_ + res$mixedEffect_list$Bf[[1]]%*%res$mixedEffect_list$beta_fixed
x_ <- as.matrix(x_)
plot(sim_res$Y[[1]], ylim = c(min(unlist(c(x_,x_t,sim_res$Y[[1]]))),
                              max(x_,x_t,sim_res$Y[[1]])))
lines(x_, col = 'green')
lines(x_t, col = 'red')
x11()
#plot(sim_res$Y[[1]])
plot(sim_res$X[[1]], col = 'red' , ylim = c(min(res$X[[1]],sim_res$X[[1]]), max(res$X[[1]],sim_res$X[[1]])))
lines(res$X[[1]], col = 'green')
lines(res$mixedEffect_list$Br[[1]]%*%(sim_res$U[,1] + res$mixedEffect_list$beta_random),col='red')
lines(res$mixedEffect_list$Br[[1]]%*%(res$mixedEffect_list$U[,1] + res$mixedEffect_list$beta_random),col='green')
#plot tracj of mixed effect
x11()
par(mfrow = c(2,1))
n_ = length(res$mixedEffect_list$betarVec[,2])
plot(res$mixedEffect_list$betarVec[,1], type='l', col='red')
lines(c(1, n_), c(mixedEffect_list$beta_random[1],mixedEffect_list$beta_random[1]))
plot(res$mixedEffect_list$betarVec[,2], type='l', col='red')
lines(c(1, n_), c(mixedEffect_list$beta_random[2],mixedEffect_list$beta_random[2]))
print(res$mixedEffect_list$Sigma)

x11()
plot(res$tauVec, type='l', col='red')
