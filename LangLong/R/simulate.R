

#' Simulating Gaussian longitudal model
#'
#' @param locs list of location of observations
#' @param theta lost with covariates beta, kappa, sigma_eps, sigma
#' @param B list of matrix with covariates
#' @param n number of gridpoints to approximate the random series
#' @param operaterType
simulateLong <- function(locs, 
                         theta, 
                         B=NULL, 
                         n = 100, 
                         operatorType = "Matern"
                         )
{

  if(is.null(B))
  {
    B=list()
    for(i in 1:length(locs))
      B[[i]] <- as.matrix(rep(1, length(locs[[i]])))
    theta$beta <- c(0)
  }

  if(length(locs) != length(B))
  {
    cat('locs and B should be equal length lists\n')
    return(-1)
  }
  
  operator_List <- create_operator(locs, n, name = operatorType)
  mesh1d <- operator_List$mesh1d

  obs_ <- list()
  for(i in 1:length(locs))
  {
    obs_[[i]] <- list(A = inla.mesh.1d.A(mesh1d, locs[[i]]), B = B[[i]])
  }
  return(simulateLong_cpp(obs_, operator_List, theta))
}

#' Simulating NIG or GAL longitudal model
#'
#' @param locs list of location of observations
#' @param theta lost with covariates mu, beta, kappa, sigma_eps, tau, lambda 
#' @param B list of matrix with covariates
#' @param n number of base points for operator
#' @param GAL - 1: if Generalized assymetric Laplace else Normal inverse Gaussian
#' @param randomIntercept - does the data have random intercept 
#' @param operatorType
simulateLongGH <- function(locs, 
                           theta, 
                           B=NULL, 
                           n = 100, 
                           noise = c("Normal","NIG","GAL"),
                           randomIntercept = FALSE,
                           operatorType = "Matern")
{
  noise <- match.arg(noise)
  if(is.null(B))
  {
    B=list()
    for(i in 1:length(locs))
      B[[i]] <- as.matrix(rep(1, length(locs[[i]])))
    theta$beta <- c(0)
  }
  theta$randomIntercept = randomIntercept
  if(length(locs) != length(B))
  {
    cat('locs and B should be equal length lists\n')
    return(-1)
  }
  
  operator_List <- create_operator(locs, n, name = operatorType)
  mesh1d <- operator_List$mesh1d
  
  obs_ <- list()
  for(i in 1:length(locs))
  {
    obs_[[i]] <- list(A = inla.mesh.1d.A(mesh1d, locs[[i]]), B = B[[i]])
  }
  theta[['noise']] <- noise
  return(simulateLongGH_cpp(obs_, operator_List, theta))
}

