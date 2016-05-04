###
# exmaple for testing estimation of mixed effect.
###
library(GHmixedeffect)
n.pers <- 10
n.obs  <- 10

sd_beta <- 0.2
sd_Y    <- 0.1

B_list <- list()
beta <- c(0.9,0.4)
beta_list <- list()
Y_list <- list()
for(i in 1:n.pers)
{
  B_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  beta_list[[i]] <- rnorm(n = length( beta), beta, sd = sd_beta)
  Y_list[[i]]        <- rnorm(n = n.obs, B_list[[i]]%*%beta_list[[i]], sd = sd_Y)
}

input <- list(Y = Y_list, 
              B = B_list, 
              Sigma = sd_beta*diag(2), 
              beta = c(0.,0.), 
              sigma_eps = 0.1,
              Niter = 1000)
res <- estimateME(input)
print(res$beta)
print(res$sigma_eps)
print(res$mixedeffect$Sigma)