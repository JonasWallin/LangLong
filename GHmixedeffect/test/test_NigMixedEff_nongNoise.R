###
# exmaple for testing estimation of random effect.
###
graphics.off()
library(GHmixedeffect)
library(rGIG)
nu <- 3
nu_mix <- 1
n.pers <- 1000
n.obs  <- 10

sd_beta <- 0.1
sd_Y    <- 0.1

B_list <- list()
beta <- c(0.9,0.9)
beta_list <- list()
Y_list <- list()
Vin        <- list()
for(i in 1:n.pers)
{
  B_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  V_mix <- rGIG( -0.5, nu_mix, nu_mix)
  beta_list[[i]] <- beta +  sqrt(V_mix)*rnorm(n = length( beta), 0, sd = sd_beta)
  V <- rGIG(rep(-0.5,n.obs), rep(nu, n.obs), rep(nu, n.obs))
  Y_list[[i]]        <-  B_list[[i]]%*%beta_list[[i]] + sqrt(V)*rnorm(n = n.obs, 0, sd = sd_Y)
  Vin[[i]] <- rep(1, n.obs)
}
meas_list <- list(Vs = Vin, sigma.eps = 1., nu = 1.)
input <- list(Y = Y_list, 
              B_random = B_list, 
              Sigma = sd_beta*diag(2),
              sigma_eps = 0,
              beta_random = c(0.,0.), 
              Niter = 3000,
              meas_list = meas_list,
              meas_noise = 'NIG',
              noise      = 'NIG')
res <- estimateME(input)
print(res$beta)
print(res$mixedeffect$Sigma)
x11()
par(mfrow=c(2,2))
plot(res$sigma_eps, type='l', col='red')
n_ <- length(res$sigma_eps)
lines(c(1, n_), c(sd_Y[1],sd_Y[1]))
lines(c(1, n_), c(sd_Y[2],sd_Y[2]))
plot(res$nu_measerror, type='l', col='red')
lines(c(1, n_), c(nu[1],nu[1]))

plot(res$betaVec[,2], type='l', col='red')
lines(c(1, n_), c(beta[2],beta[2]))
plot(res$nuVec, type='l', col='red')
#lines(c(1, n_), c(beta[2],beta[2]))
