library(methods)
library(INLA)
library(LANGlong)

library(Matrix)
graphics.off()
n.obs  = 80
n.rep  = 40
n.proc = 100
n.sim  = 2
N      = 10
stepsize = 1/(n.rep*n.obs)
use.R = FALSE
B = list()
locs <- list()
P <- seq(0,500, length.out = n.obs)
for(i in 1:n.rep){
  locs[[i]] <- P
  B[[i]] = cBind(rep(1, length(P)),locs[[i]]/2000)
}
theta <- list(kappa = 0, sigma = 0.1, tau =1e4, beta = c(2.0,-0.001), lambda = 10, mu = 0, sigma_r = 1)
output <- simulateLongGH(locs,
                         theta,
                         B,
                         n = n.proc,
                         noise="Normal", 
                         operatorType = 'FD2',
                         randomIntercept = TRUE) 


plot(locs[[1]],output$Y[[1]])
#theta0 = list(kappa = log(theta$kappa), sigma = 2*log(theta$sigma), tau = log(theta$tau),
#              beta = theta$beta, lambda = log(theta$lambda), mu  =theta$mu)

theta0 = list(kappa = log(theta$kappa), 
              sigma = log(1), 
              tau = log(1e3), 
              beta = c(0,0), 
              lambda = 1, 
              mu  =0.,
              sigma_r = 1)


#result <- estimateLongGH(output$Y, B, locs, theta = theta0, 
#                         N=N, stepsize = stepsize,
#                         n=n.proc, nsim = n.sim, GAL = 0, operatorType = 'fd2',
#                         mixedEffect = TRUE)  

#plotLongGh(result)
#res <- sample_posteriror(result, sim = 2)