###
# exmaple for testing estimation of mixed effect.
###
library(GHmixedeffect)
library(MASS)
n.pers <- 40
n.obs  <- 10
Sigma_beta <- matrix(c(0.4,0.12,0.12,0.2), nrow = 2, ncol =  2)
sd_Y    <- 0.1

B_list <- list()
beta <- c(0.9,0.4)
beta_list <- list()
Y_list <- list()
beta_vec <- c()
for(i in 1:n.pers)
{
  B_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  beta_list[[i]] <- mvrnorm(n =  1, mu = beta, Sigma  = Sigma_beta)
  Y_list[[i]]        <- rnorm(n = n.obs, B_list[[i]]%*%beta_list[[i]], sd = sd_Y)
  beta_vec <- rbind(beta_vec, beta_list[[i]] )
}

input <- list(Y = Y_list, 
              B = B_list, 
              Sigma = diag(2), 
              beta = c(0.9,0.4), 
              sigma_eps = 0.5,
              Niter = 20)
res <- estimateME(input)
print(res$beta)
print(res$sigma_eps)
print(res$mixedeffect$Sigma)
print(cov(t(res$mixedeffect$U)))