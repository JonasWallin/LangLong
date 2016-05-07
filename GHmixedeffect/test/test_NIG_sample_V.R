###
# exmaple for testing estimation of random effect.
# with both fixed and random effect
###
graphics.off()
rm(list=ls())
library(GHmixedeffect)

f <- function(v, sd_Y, sd, y, nu, mu  = 0){ 
  
    dens <- exp(-(y -(- mu + v*mu))^2/(2 * (sd_Y^2 + sd^2 * v ))) / sqrt(sd_Y^2 + sd^2 * v)
    dens <- dens * v^(-1.5) * exp(- nu/2 * (v + v^-1))
    dens[v == 0] <- 0
    return(dens)
}

n.pers <- 1
n.obs  <- 1

sd_beta <- 0.2
sd_Y    <- 0.1
nu <- 0.56
mu = 1.
Br_list <- list()
beta_r <- c(0.)
betar_list <- list()
Y_list <- list()
for(i in 1:n.pers)
{
  Br_list[[i]]    <- as.matrix(1)
  betar_list[[i]] <- 1
  Y_list[[i]]        <- 2#rnorm(n = n.obs, 
                          #    Br_list[[i]]%*%betar_list[[i]],
                           #    sd = sd_Y)
}

input <- list(Y = Y_list, 
              B_random = Br_list,
              Sigma = sd_beta^2*diag(1), 
              beta_random = beta_r, 
              sigma_eps = 0.1,
              Niter = 20000,
              nu = nu)
res <- samplePostV(input)

x <- seq( 0, 60, length = 1000)
x11()
par(mfrow=c(2,1))
g <- function(x){ f(x, sd_Y, sd_beta, Y_list[[1]], nu)} 
dens <- g(x)
hist(res,100, prob=T)
lines(x, dens/integrate(g,0,Inf)$value, col='red' )
#plot(x, dens/integrate(g,0,Inf)$value, col='red')

input <- list(Y = Y_list, 
              B_random = Br_list,
              Sigma = sd_beta^2*diag(1), 
              beta_random = beta_r, 
              sigma_eps = 0.1,
              Niter = 20000,
              nu = nu,
              mu = as.matrix(mu))
res <- samplePostV(input)

g <- function(x){ f(x, sd_Y, sd_beta, Y_list[[1]], nu, mu = mu)} 
dens <- g(x)
hist(res,100, prob=T)
lines(x, dens/integrate(g,0,Inf)$value, col='red' )
#plot(x, dens/integrate(g,0,Inf)$value, col='red')