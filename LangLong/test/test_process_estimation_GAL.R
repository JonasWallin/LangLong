# testing that the parameters converge for GAL distribution
rm(list=ls())
library(LANGlong)
library(INLA)
library(rGIG)
library(methods)
graphics.off()
npers <- 100
nobs  <- 100
nIter <- 20000
n     <- 100 #n grid points
noise <- "GAL"
nu_true <- 50
mu_true <- 5
nu_guess <- 1.1
mu_guess <- 1
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
               noise = noise,
               operatorType = "fd2",
               n = n)
operator_list <- create_operator(locs, n, name = "fd2")

obs_list <- list()
X <- list()
V <- list()
for(i in 1:length(locs)){
  obs_list[[i]] <- list(A = inla.mesh.1d.A(operator_list$mesh1d, locs[[i]]), 
                        Y=output_sim$Y[[i]], 
                        locs = locs[[i]])
  X[[i]] <- rep(0, n) 
  V[[i]] <- operator_list$h
}
operator_list$tau <-theta$tau
processes_list <- list(nu = nu_guess, mu = mu_guess, X = output_sim$X, V = output_sim$V, noise = noise)
input <- list( obs_list         = obs_list,
               operator_list    = operator_list,
               processes_list   = processes_list,
               nIter            = nIter,     # iterations to run the stochastic gradient
               nSim             = 2,
               nBurnin          = 10,   # steps before starting gradient estimation
               silent           = 0, # print iteration info)
               step0            = 1,
               alpha            = 0.01,
               sigma            = theta$sigma
              )
output <- estimateProcess_cpp(input)
x11()
par(mfrow=c(3,1))
plot(locs[[1]],output_sim$Y[[5]])
lines(output_sim$xloc, output_sim$X[[5]])
lines(operator_list$mesh1d$loc, output$X[[5]],col='red',lty='dashed')
plot(output$mu_vec)
plot(output$nu_vec)
K = sqrt(theta$tau)*operator_list$Q
Z = as.vector(K%*%output_sim$X[[1]])
sum(Z- ( -operator_list$h + output_sim$V[[1]])*-97.7582)

print(paste("nu - nu_true =",output$nu-nu_true))
print(paste("mu - mu_true =",output$mu-mu_true))