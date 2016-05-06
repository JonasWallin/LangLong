##
# simple test that verifies that the model can correctly idenitfy the parameters, 
# from simulated data
#
##
library(rGIG)
library(MASS)
n.pers <- 2000 #number of patients
n.obs  <- 10 #number of obs per patient

COV_beta <- matrix(c(0.2,0.1,0.1,0.2), ncol = 2, nrow = 2)
sd_Y    <- 0.1 # error of the noise

B_list <- list()
beta <- c(0.9,0.4)
mu   <- c(0.2, 0)
nu <- 1

beta_list <- list()
V_list    <- list()
Y_list    <- list()
for(i in 1:n.pers)
{
  B_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  V <- rGIG(0.5, nu, nu)
  V_list[[i]] <- V
  beta_list[[i]] <- beta - mu * 1  + V*mu + sqrt(V) * mvrnorm(n = 1, mu  =c(0,0), Sigma = COV_beta)
  Y_list[[i]]        <- rnorm(n = n.obs, B_list[[i]]%*%beta_list[[i]], sd = sd_Y)
}
beta_mat <- t(matrix(unlist(beta_list), nrow= 2, ncol = n.pers))
x11()
hist(beta_mat[,1],100)
#hist(beta_mat[,2])
