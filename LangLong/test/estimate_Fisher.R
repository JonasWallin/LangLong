#####
# test function for estimating parameters for stochastic proccesses mixed effect
# #TODO: add a small regularization on the random effect otherwise it can be come singular
#
####
graphics.off()
library(LANGlong)
library(methods)

nIter <- 50
n.pers <- 100
nBurnin <- 5
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
  operator_list$kappa <- kappa.true
  operator_list$tau   <- tau.true
} else {
  operator_list$tau   <- tau.true
}

#mixedEffect_list$beta_random = as.matrix(c(0.,0.))
#mixedEffect_list$Sigma = diag(c(1,1.))
#mError_list$sigma = as.matrix(c(1.))
processes_list$X <- sim_res$X
obs_list <- list()
for(i in 1:length(locs))
  obs_list[[i]] <- list(A = spde.A(locs[[i]],
                                   operator_list$loc),
                        Y=sim_res$Y[[i]],
                        locs = locs[[i]])
in_list <- list(Y                = sim_res$Y,
                nIter            = nIter,
                nSim             = nSim,
                nBurnin          = nBurnin,
                pSubsample       = 0.1,
                locs             = locs,
                obs_list         = obs_list,
                mixedEffect_list = mixedEffect_list,
                measurementError_list  = mError_list,
                processes_list   = processes_list,
                operator_list    = operator_list,
                silent           = 0)
res <- estimateFisher(in_list)

