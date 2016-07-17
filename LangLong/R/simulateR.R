library(rGIG)

create_matrices_ <- function(locs, n)
{
  min_l <- min(locs[[1]])
  max_l <- max(locs[[1]])
  if(length(locs) > 1){
    for(i in 2:length(locs))
    {
      min_l <- min(min_l, min(locs[[i]]))
      max_l <- max(max_l, max(locs[[i]]))
    }
  }
  
  P <- seq(min_l, max_l, length.out = n)
  mesh1d <- inla.mesh.1d(P)
  MatrixBlock <- inla.mesh.1d.fem(mesh1d)
  output <- list(MatrixBlock = MatrixBlock, mesh1d = mesh1d)
  return(output)
}

#' Simulating longitudal model
#'
#' @param locs list of location of observations
#' @param theta lost with covariates mu, kappa, sigma_eps, sigma
#' @param B list of matrix with covariates
simulateLong.R <- function(loc, 
                           theta, 
                           n = 100, 
                           B = NULL,
                           noise = c("Normal","NIG","GAL"),
                           mixedEffect = FALSE,
                           Bmixed      = NULL,
                           mixedEffectList = NULL,
                           operatorType = "Matern")
{
  
  nrep = length(loc)
  if(is.null(B))
  {
    B=list()
    for(i in 1:nrep)
      B[[i]] <- as.matrix(rep(1, length(loc[[i]])))
    theta$beta <- c(0)
  }
  
  if(nrep != length(B))
  {
    cat('locs and B should be equal length lists\n')
    return(-1)
  }
  
  
  operator_List <- create_operator(loc, n, name = operatorType)


  tau = as.double(theta$tau)
  sigma = theta$sigma
  beta = as.double(theta$beta)
  if(length(sigma) == 1){
    sigma <- rep(sigma,nrep)
  }
  
  if(operatorType=="Matern"){
    kappa = theta$kappa
    K = sqrt(tau) %*% (operator_List$G + kappa*operator_List$C)
    Q = (K%*%operator_List$Ci%*%K)
    R = chol(Q)  
  }else{
    K = sqrt(tau)*operator_List$Q
    
    Ci = as(sparseMatrix(i=1:n,j=1:n,x=1/operator_List$h,dims=c(n, n)), "CsparseMatrix")
    Q = Matrix::t(K) %*% Ci %*% K
    R = chol(Q)
  }
  
  A <- list()
  for(i in 1:length(loc))
  {
    A[[i]] <- inla.mesh.1d.A(operator_List$mesh1d, loc[[i]])
    
  }
  
  Y = list()
  X = list()
  V = list()
  for(i in 1:nrep){
    if(noise == "Normal"){
      X[[i]] = solve(R,rnorm(dim(R)[1]))
    }else if (noise == "NIG"){
      V[[i]] =  rGIG(rep(-0.5, n), 
                     rep( theta$nu^2, n), 
                     (operator_List$h * theta$nu)^2)
      Z <- (- operator_List$h  + V[[i]]) * theta$mu + sqrt(V[[i]]) * rnorm(n)
      X[[i]] <- solve(K, Z)
    }else{
      print("GAL  NOT IMPLIMENTED YET")
    }
      
    Y[[i]] = (B[[i]]%*%beta + A[[i]]%*%X[[i]] + sigma[i]*rnorm(dim(A[[i]])[1]))@x
  }
  return(list(Y=Y, X=X, xloc = operator_List$mesh1d$loc, A=A, V= V)) 
}