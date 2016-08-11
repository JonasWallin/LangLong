#simulate V
#simulate X
#simulate Y
#use second order RW
rm(list=ls())
library(LANGlong)
library(INLA)
library(rGIG)
library(methods)
graphics.off()
npers <- 2
nobs  <- 100
nIter <- 15000
n     <- 100 #n grid points

nu_true <- 5
mu_true <- 10
nu_guess <- 2
mu_guess <- 0
theta <- list()
theta$sigma <- 0.1 # meas error
theta$tau   <- 0.5
theta$nu    <- nu_true
theta$mu    <- mu_true
locs   <- list()
for(i in 1:nobs)
{ 
  locs[[i]]   <- seq(0, 1, length = nobs)
}

output_sim <- simulateLong.R(locs, 
               theta,
               noise = "NIG",
               operatorType = "fd2",
               n = n)
operator_list <- create_operator(locs, n, name = "fd2")

obs_list <- list()
X <- list()
V <- list()
for(i in 1:length(locs)){
  obs_list[[i]] <- list(A =  spde.A(x = operator_list$loc, loc = locs[[i]]), 
                        Y=output_sim$Y[[i]], 
                        locs = locs[[i]])
  X[[i]] <- rep(0, n) 
  V[[i]] <- operator_list$h
}


mError_list <- list(noise = "Normal",
                    sigma = theta$sigma)

operator_list$tau <-theta$tau
processes_list <- list(nu = nu_guess, mu = mu_guess, X = output_sim$X, V = output_sim$V, noise = "NIG")
input <- list( obs_list         = obs_list,
               operator_list    = operator_list,
               processes_list   = processes_list,
               nIter            = nIter,     # iterations to run the stochastic gradient
               nSim             = 2,
               nBurnin          = 100,   # steps before starting gradient estimation
               silent           = 0, # print iteration info)
               step0            = 1,
               alpha            = 0.01,
               measurment_list   = mError_list
              )
output <- estimateProcess_cpp(input)
x11()
par(mfrow=c(3,1))
plot(locs[[1]],output_sim$Y[[5]])
lines(output_sim$xloc, output_sim$X[[5]])
lines(operator_list$mesh1d$loc, output$Xs[[5]],col='red',lty='dashed')
plot(output$mu_vec)
plot(output$nu_vec)
K = sqrt(theta$tau)*operator_list$Q
Z = as.vector(K%*%output_sim$X[[1]])
sum(Z- ( -operator_list$h + output_sim$V[[1]])*-97.7582)