###
# exmaple for testing estimation of random effect.
# with both fixed and random effect,
# and the distribution is normal
###
library(testthat)
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
meas_list <- list(sigma_eps <- 0.1, noise = "Normal")
mixedEffect_list <- list(B_random = Br_list, 
                         B_fixed  = Bf_list,
                         Sigma = sd_beta*diag(2), 
                         beta_random = c(0.,0.),
                         noise = "Normal")
input <- list(Y = Y_list, 
              mixedEffect_list = mixedEffect_list,
              measurementError_list = meas_list,
              nSim = 2,
              alpha = 0.3,
              step0 = 1,
              Niter = 100)
res <- estimateME(input)
print(res$mixedEffect_list$Sigma)
test_that("simple Gaussian-Gaussian random",
          {
            expect_equal(res$mixedEffect_list$beta_random, beta_r, tolerance  = 0.1)
          })
test_that("simple Gaussian-Gaussian fixed",
{
  expect_equal(res$mixedEffect_list$beta_fixed, beta_f, tolerance  = 0.1)
})
#cat('result :')
#cat("beta_random true = ", beta_r,'\n')
#cat("diff = ",res$mixedEffect_list$beta_random - beta_r,"\n")
#cat("beta_random true = ", beta_f,'\n')
#cat("diff = ",res$mixedEffect_list$beta_fixed - beta_f,'\n')

res$Niter     <- 10
res$nSim      <- 100
res$nBurnin   <- 10
res$Y         <- Y_list
res$sigma_eps <- res$sigma_eps[length(res$sigma_eps)]
res$meas_list <- res$measerror
output <- EstimateFisherInformation(res)
#output$mixedeffect$Cov_theta
#output$measerror$Cov_theta