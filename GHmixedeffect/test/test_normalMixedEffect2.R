###
# exmaple for testing estimation of mixed effect.
# using five covaraites, and no fixed effect
###
library(testthat)
library(GHmixedeffect)
n.pers <- 100
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
  beta_list[[i]] <- beta + sd_beta*rnorm(n = length( beta), 0, sd = 1)
  Y_list[[i]]        <- rnorm(n = n.obs, B_list[[i]]%*%beta_list[[i]], sd = sd_Y)
}

meas_list <- list(sigma_eps <- 0.1, noise = "Normal")
mixedEffect_list <- list(B_random = B_list, 
                         Sigma = sd_beta*diag(5), 
                         beta_random =rep(0,5),
                         noise = "Normal")
input <- list(Y = Y_list, 
              mixedEffect_list = mixedEffect_list,
              measurementError_list = meas_list,
              nSim = 2,
              alpha = 0.3,
              step0 = 1,
              Niter = 1000)

res <- estimateME(input)


test_that("simple Gaussian-Gaussian random effect (5 cov)",
{
  expect_equal( res$mixedEffect_list$beta_random, beta, tolerance  = 0.1)
})
test_that("simple Gaussian-Gaussian measuerment sigma (5 cov)",
{
expect_equal( res$measurementError_list$sigma, sd_Y, tolerance  = 0.1)
})
test_that("simple Gaussian-Gaussian measuerment Sigma (5 cov)",
{
  expect_equal( mean(diag(res$mixedEffect_list$Sigma)-sd_beta^2), 0, tolerance  = 0.05)
})