##
# test parameter estimation with fixed V (known variance components)
##
library(testthat)
library(GHmixedeffect)
library(rGIG)

noises = c("IG", "NIG")
for( k in 1:length(noises)){
  noise <- noises[k]
  nu <- 3
  n.pers <- 10
  n.obs  <- 400
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
    if(noise == "NIG"){
      V <- rGIG(rep(-0.5,n.obs), rep(nu, n.obs), rep(nu, n.obs))
    }else if(noise == "IG"){
      V <- 1/rgamma(n.obs, nu, nu)
    }
    Y_list[[i]]        <-  B_list[[i]]%*%beta_list[[i]] + sqrt(V)*rnorm(n = n.obs, 0, sd = sd_Y)
    Vin[[i]] <- V
  }
  meas_list <- list(Vs = Vin, sigma.eps = 1., nu = 1., noise = noise)
  
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
                Niter = 300,
                silent = 1)
  name <- paste("known V, noise = ",noise,", testing latent parameters",sep="")
  test_that(name,
  {
    res <- estimateME_Vfixed(input)
    expect_equal(res$measurementError_list$nu, nu, tolerance  = 0.05)
  })
}