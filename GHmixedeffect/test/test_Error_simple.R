###
# simple test to just see if estimation method does not crash
###
graphics.off()
library(testthat)
library(GHmixedeffect)
library(rGIG)
nu <- 3
n.pers <- 10
n.obs  <- 200

sd_beta <- 0.01
sd_Y    <- 0.1

B_list <- list()
beta <- c(0.9,0.4)
beta_list <- list()
Y_list <- list()
Vin        <- list()
for(i in 1:n.pers)
{
  B_list[[i]]    <- cbind(rep(1, n.obs), (1:n.obs) / n.obs )
  beta_list[[i]] <-beta+  rnorm(n = length( beta), 0, sd = sd_beta)
  V <- rGIG(rep(-0.5,n.obs), rep(nu, n.obs), rep(nu, n.obs))
  Y_list[[i]]        <-  B_list[[i]]%*%beta_list[[i]] + sqrt(V)*rnorm(n = n.obs, 0, sd = sd_Y)
  Vin[[i]] <- rep(1, n.obs)
}
meas_list <- list(Vs = Vin, sigma.eps = 1., nu = 1., noise = "NIG")

mixedEffect_list <- list(B_random = B_list, 
                         Sigma = sd_beta*diag(2), 
                         beta_random = c(0.,0.),
                         noise = "Normal")
input <- list(Y = Y_list, 
              mixedEffect_list = mixedEffect_list,
              measurementError_list = meas_list,
              nSim = 2,
              alpha = 0.3,
              step0 = 1,
              Niter = 20,
              silent = 1)



test_that("simple Normal error does not crash",
{
  meas_list$noise = "normal"
  res <- estimateME(input)
})

test_that("simple NIG error does not crash",
{
  meas_list$noise = "NIG"
  res <- estimateME(input)
})
test_that("simple IG error does not crash",
{
  meas_list$noise = "IG"
  res <- estimateME(input)
})