##
# simple test that verifies that the model can correctly idenitfy the parameters, 
# from simulated data
#
##
rm(list=ls())
graphics.off()
library(GHmixedeffect)
library(rGIG)
library(MASS)
n.pers <- 1000 #number of patients
n.obs  <- 10 #number of obs per patient

COV_beta <- matrix(c(0.2,0.1,0.1,0.2), ncol = 2, nrow = 2)
sd_Y    <- 0.1 # error of the noise

Br_list <- list()
betar <- c(0.9,0.4)
mu   <- -c(0.2, -0.2)
nu <- .7
betar_list <- list()
V_list    <- list()
Y_list    <- list()
for(i in 1:n.pers)
{
  Br_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  V <- rGIG(-0.5, nu, nu)
  V_list[[i]] <- V
  betar_list[[i]] <- betar - mu * 1  + V * mu +
                     sqrt(V) * mvrnorm(n = 1, mu  =c(0,0), Sigma = COV_beta)
  Y_list[[i]]        <- rnorm(n = n.obs,
                              Br_list[[i]]%*%betar_list[[i]], sd = sd_Y)

}

input <- list(Y = Y_list, 
              B_random    = Br_list, 
              Sigma       = COV_beta, 
              beta_random = betar, 
              sigma_eps = sd_Y,
              mu = as.matrix(mu),
              nu = as.matrix(nu),
              Niter = 10,
              noise = "NIG")

beta_mat <- t(matrix(unlist(betar_list), nrow= 2, ncol = n.pers))
#x11()
#hist(beta_mat[,1],100)
#hist(beta_mat[,2])
res <- estimateME(input)
print(res$mixedeffect$Sigma)
print(cov(beta_mat))
print(cov(t(res$mixedeffect$U)))