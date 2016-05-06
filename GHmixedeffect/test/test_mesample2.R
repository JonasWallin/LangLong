###
# exmaple for testing estimation of mixed effect.
# using five covaraites due to previous errors
###
library(GHmixedeffect)
n.pers <- 20
n.obs  <- 20

sd_beta <- 0.2
sd_Y    <- 0.1

B_list <- list()
beta <- c(0.9,0.4,0.2,0.2,0.1)
beta_list <- list()
Y_list <- list()
for(i in 1:n.pers)
{
  B_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs, runif(n.obs), rnorm(n.obs), rnorm(n.obs) )
  beta_list[[i]] <- rnorm(n = length( beta), beta, sd = sd_beta)
  Y_list[[i]]        <- rnorm(n = n.obs, B_list[[i]]%*%beta_list[[i]], sd = sd_Y)
}

input <- list(Y = Y_list, 
              B_random = B_list, 
              Sigma = sd_beta*diag(5), 
              beta_ranom = rep(0,5), 
              sigma_eps = 0.1,
              Niter = 2000)
res <- estimateME(input)
print(res$beta)
print(res$sigma_eps)
print(res$mixedeffect$Sigma)