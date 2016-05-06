###
# exmaple for testing estimation of random effect.
# with both fixed and random effect
###
library(GHmixedeffect)
n.pers <- 100
n.obs  <- 10

sd_beta <- 0.2
sd_Y    <- 0.1

Br_list <- list()
beta_r <- c(0.9,0.4)
beta_f <- 1
betar_list <- list()
Bf_list <- list()
Y_list <- list()
for(i in 1:n.pers)
{
  Bf_list[[i]]    <- as.matrix(runif(n = n.obs))
  Br_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  betar_list[[i]] <- beta_r + rnorm(n = length( beta_r), sd = sd_beta)
  Y_list[[i]]        <- rnorm(n = n.obs, 
                              Br_list[[i]]%*%betar_list[[i]] + Bf_list[[i]]%*%beta_f,
                               sd = sd_Y)
}

input <- list(Y = Y_list, 
              B_random = Br_list, 
              B_fixed  = Bf_list,
              Sigma = sd_beta*diag(2), 
              beta_random = c(0.,0.), 
              sigma_eps = 0.1,
              Niter = 100)
res <- estimateME(input)
print(res$mixedeffect$Sigma)
print(res$mixedeffect$beta_random)
print(res$mixedeffect$beta_fixed)