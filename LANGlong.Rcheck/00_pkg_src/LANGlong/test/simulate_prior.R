library(LANGlong)


n.pers <- 10
n.obs  <- 11

Y   <- list()
loc <- list()
B_random <- list()
B_fixed  <- list()
Vin <- list()
for(i in 1:n.pers)
{
  B_random[[i]] <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  B_fixed[[i]]  <- rnorm(n = n.obs, 0,sd=2)
  Y[[i]] <- rep(1,n.obs)
  loc[[i]] <- 1:n.obs
  Vin[[i]] <- rep(1, n.obs)
}
mError_list <- list(Vs = Vin, noise = "Normal", sigma = 0.1)
mixed_list  <- list(B_random = B_random, 
                    B_fixed  = B_fixed,
                    beta_random = c(2,-1), 
                    beta_fixed  = c(1.), 
                    Sigma = diag(0.1), 
                    noise = "Normal")
processes_list = list(noise = "Normal")
operator_list  = list(type = "Matern" )
simulateLongPrior( Y               = Y, 
                               loc             = loc,
                               mixed_list      = mixed_list,
                               measurment_list = mError_list,
                               processes_list  = processes_list,
                               operator_list   = operator_list)