library(INLA)

#' Estimating longitudal model
#'
#' @param locs list of location of observations
#' @param theta list with covariates mu, kappa, sigma_eps, tau
#' @param n discretization size
#' @param Niter - max number of iterations
#' @param burnin
#' @param long.percentage percentage of longitudinal samples to use in each iteration.
#' @param nsim
#' @param B list of matrix with covariates
#' @param noise - either Normal, Generalized assymetric Laplace, or Normal inverse Gaussian
#' @param silent - print redsults
#' @param V - list of inital guess of Variance components
#' @param U - initaul guess of random intercept
#' @param mixedEffectList - list contaning $B (mixed effects) 
#' @param Bmixed - the covariates for mixed effects (if mixeEffectList is not null, otherwise uses only mixed effects)
#' @param operatorType - type of operator: 1) Matern 2) Finite Difference 
#' 
#' 
estimateLongGH <- function(Y, 
                           B, 
                           loc, 
                           theta = NULL, 
                           stepsize = 1e-4, 
                           n = 100, 
                           Niter = 100, 
                           burnin = 10,
                           long.percentage = 100,
                           commonsigma = TRUE,
                           nsim = 1, 
                           noise  = c("Normal","NIG","GAL"),
                           mNoise = c("Normal", "NIG"),
                           silent=FALSE,
                           V = NULL,
                           mixedEffect = FALSE,
                           mixedType   = c("Normal", "NIG"),
                           Bmixed      = NULL,
                           mixedEffectList = NULL,
                           operatorType = "Matern",
                           U = NULL
                           )
{
  noise <- match.arg(noise)
  mixedType <- match.arg(mixedType)
  mixedType <- match.arg(mixedType)
  mNoise    <- match.arg(mNoise)
  if(missing(Y))
    stop("Must supply data")

  if(missing(B))
    stop("Must supply list with covariates")
  
  if(missing(loc))
    stop("Must supply list with locations")
  
  operator_List <- create_operator(locs, n, name = operatorType)
  mesh1d <- operator_List$mesh1d
  
  obs_ <- list()
  for(i in 1:length(locs))
  {
    obs_[[i]] <- list(A = inla.mesh.1d.A(mesh1d, locs[[i]]), B = B[[i]],Y=Y[[i]], locs = locs[[i]])
  }
  Nlong = length(locs)
  Nlong = max(min(round(long.percentage*Nlong/100),Nlong),1)
  
  if(is.null(V))
  {
    if(noise == "Normal"){
      for(i in 1:length(locs))
      {
        V[[i]] <- operator_List$h
      }
    }else{
      for(i in 1:length(locs))
      {
        V[[i]] <- c(rep(1, length(operator_List$h)))
      }
    }
        
  }
  if(is.null(U))
    U <- as.vector(rep(0,length(locs)))
  if(is.null(theta))
    theta <- list(tau = as.matrix(0) , 
                  sigma = 0, 
                  asigma = 0,
                  bsigma = 0,
                  beta = as.matrix(rep(0, dim(B[[1]])[2] )),
                  kappa = 0,
                  lambda = 10,
                  mu     = 0,
                  sigma_r = 0)
  
  if(is.null(mixedEffectList) ==FALSE){
    mixedEffect_list <- mixedEffectList
    mixedEffect_list$on =TRUE
  }else{
    mixedEffect_list <- list(on = mixedEffect)
    if(mixedEffect == TRUE & is.null(Bmixed) == TRUE)
    {
      mixedEffect_list$B = B 
    }else{
      mixedEffect_list$B = Bmixed
    }
    if(mixedEffect == TRUE)
    {
      if(is.null(theta$Sigma) == TRUE){
        mixedEffect_list$Sigma = diag(dim(mixedEffect_list$B[[1]])[2])
      }else{mixedEffect_list$Sigma = theta$Sigma}
    }
  }
  if(noise == "Normal"){
    noise.i = -1
  } else if(noise=="GAL") {
    noise.i = 1
  } else {
    noise.i = 0
  }
  
  theta <- estimateLongGH_cpp(obs_, 
                              operator_List, 
                              theta, 
                              mixedEffect_list, 
                              stepsize, 
                              Niter, 
                              nsim, 
                              burnin, 
                              noise.i,
                              commonsigma, 
                              silent, 
                              V, 
                              U,
                              Nlong)

  output              <- list(theta = theta, 
                              noise = noise, 
                              commonsigma = commonsigma,
                              mixedEffectList = theta$mixedeffect)
  output$mixedEffectList$on = mixedEffect
  output$obs_          <- obs_
  output$operator_List <- operator_List
  return(output)
}

#' posterior samples generates mean and variance of X, and V
#'
#' @param output estimateLongGH
#' 
sample_posteriror<- function(output, sim = 100)
{
  noise.i = 0
  if(output$noise == "Normal"){
    noise.i = -1
  }else if(output$noise == "GAL"){
    noise.i = 1
  } 
  sim_res <- samplePosteriorGH(output$obs_, output$operator_List, output$theta,output$mixedEffectList, output$theta$V, sim, noise.i,output$commonsigma)
  return(sim_res)
}


#' plot parameter tracjetories 
#'
#' @param output estimateLongGH
#' 
plotLongGh <-function(output, beta = c(1))
{
  theta <- output$theta
  par(mfrow=c(3,2))
  plot(theta$beta_vec[,beta], main=expression(beta))
  if(output$operator_List$type == 'matern')
    plot(theta$kappa_vec, main = expression(kappa))
  if(output$commonsigma == 1){
    plot(theta$sigma_vec, main = expression(sigma))  
  }else{
    plot(theta$asigma_vec, main = expression(asigma))
    plot(theta$bsigma_vec, main = expression(bsigma))
  }
  
  plot(theta$tau_vec,   main = expression(tau))
  
  if(output$noise != "Normal")
    plot(theta$lambda_vec,   main = expression(lambda))
  
}

plottraj <- function(nr, output, samples)
{
  X_mean <- samples$X_mean[[nr]]
  plot(output$operator_List$mesh1d$loc, X_mean, xlab = 't', ylab = 'X(t)', type = 'l')
  
}
#' simualte from the prior model using output from parameter estimate
#' @param output estimateLongGH
sample_prior_from_result_LongGH <- function(output)
{
  theta_in <- output$theta
  theta_in$GAL <- result$GAL 
  return(simulateLongGH_cpp(output$obs_, output$operator_List, theta_in))
}