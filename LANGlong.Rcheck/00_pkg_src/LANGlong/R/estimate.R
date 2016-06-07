

#' Estimating longitudal model
#'
#' @param locs list of location of observations
#' @param theta lost with covariates mu, kappa, sigma_eps, tau
#' @param B list of matrix with covariates
estimateLong <- function(Y, B, loc, theta = NULL, stepsize = 1e-4, n = 100,N=100,nsim = 1, silent=FALSE)
{
  if(missing(Y))
    stop("Must supply data")

  if(missing(B))
    stop("Must supply list with covariates")
  
  if(missing(loc))
    stop("Must supply list with locations")
  
  
  tmp          <- create_matrices_(locs, n)
  MatrixBlock  <- tmp$MatrixBlock
  mesh1d       <- tmp$mesh1d
  operator_List <- list(C = as(as(MatrixBlock$c0,"CsparseMatrix"), "dgCMatrix"),
                        Ci = as(sparseMatrix(i=1:n,j=1:n,x=1/(MatrixBlock$c0@x),
                                             dims=c(n, n)), "CsparseMatrix"),
                        G = MatrixBlock$g1,
                        B.kappa = as.matrix(rep(1, dim(MatrixBlock$c0)[1])),
                        kappa = 0
  )
  
  obs_ <- list()
  for(i in 1:length(locs))
  {
    obs_[[i]] <- list(A = inla.mesh.1d.A(mesh1d, locs[[i]]), B = B[[i]],Y=Y[[i]])
  }
  
  if(is.null(theta))
    theta <- list(tau = as.matrix(0) , 
                  sigma = 0, 
                  beta = as.matrix(rep(0, dim(B[[1]])[2] )),
                  kappa = 0
                  )
  return(estimateLong_cpp(obs_, operator_List, theta, stepsize,N,nsim, silent))
}