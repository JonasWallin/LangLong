library(INLA)
library(LANGlong)
library(methods)


n.pers <- 10
n.obs  <- 100
n <- 100
Y   <- list()
locs <- list()
B_random <- list()
B_fixed  <- list()
Vin <- list()
for(i in 1:n.pers)
{
  B_random[[i]] <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  B_fixed[[i]]  <- as.matrix(rnorm(n = n.obs, 0,sd=2))
  Y[[i]] <- rep(1,n.obs)
  locs[[i]] <- 1:n.obs
  Vin[[i]] <- rep(1, n.obs)
}
mError_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1)
mixedEffect_list  <- list(B_random = B_random, 
                    B_fixed  = B_fixed,
                    beta_random = as.matrix(c(2,-1)), 
                    beta_fixed  = as.matrix(c(1.)), 
                    Sigma = diag(0.1, 0.2), 
                    noise = "Normal")

operator_list <- create_operator(locs, n, name = "Matern")
operator_list$type  <- "matern"
operator_list$kappa <- -4
operator_list$tau   <- 10


processes_list = list(noise = "Normal", 
                      nu  = 2., 
                      mu = 1.)
processes_list$V <- list() 
for(i in 1:length(locs))
{
  processes_list$V[[i]] <- operator_list$h
}

sim_res <- simulateLongPrior( Y               = Y, 
                  locs             = locs,
                  mixedEffect_list      = mixedEffect_list,
                  measurment_list = mError_list,
                  processes_list  = processes_list,
                  operator_list   = operator_list)

plot(locs[[1]], sim_res$Y[[1]])
lines(operator_list$mesh1d$loc, sim_res$X[[1]], col = 'red')