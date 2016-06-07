estimateLong.R <- function(Y, B, loc, theta = NULL, stepsize = 1e-4, n = 100,N=100,nsim = 1)
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
  C = as(as(MatrixBlock$c1,"CsparseMatrix"), "dgCMatrix")
  Ci = as(sparseMatrix(i=1:n,j=1:n,x=1/(MatrixBlock$c0@x),
                       dims=c(n, n)), "CsparseMatrix")
  G = MatrixBlock$g1

  
  A <- list()
  for(i in 1:length(locs)){
    A[[i]] <- inla.mesh.1d.A(mesh1d, locs[[i]])
  }
  
  if(is.null(theta))
    theta <- list(tau = as.matrix(0) , 
                  sigma = 0, 
                  beta = as.matrix(rep(0, dim(B[[1]])[2] )),
                  kappa = 0)
  
  tau = as.double(theta$tau)
  sigma = theta$sigma
  beta = as.double(theta$beta)
  kappa = theta$kappa
  
  d = dim(G)[1]
  nrep = length(loc)
  x11()
  for(iter in 1:N){
    cat("sigma = ",exp(0.5*sigma), ", beta = ", beta, ", tau = ", exp(tau), ", kappa = ", exp(kappa),"\n")
    K = G + exp(kappa)*C
    Q = exp(tau)*t(K)%*%Ci%*%K
    dtau = nsim*nrep*d/2    
    dsigma = 0
    dbeta = 0
    for(i in 1:nrep){
      Qi = Q + t(A[[i]])%*%A[[i]]/exp(sigma)
      b = t(A[[i]])%*%(Y[[i]] - B[[i]]%*%beta)/exp(sigma)
      R = chol(Qi)
      for(ii in 1:nsim){
        z = rnorm(d)
        X = solve(R,solve(t(R),b)+z)  
        Ysim = Y[[i]] - B[[i]]%*%beta - A[[i]]%*%X
        dsigma = dsigma - length(Ysim)/2 + exp(-sigma)*sum(Ysim^2)/2
        dbeta = as.double(dbeta + exp(-sigma)*t(B[[i]])%*%Ysim)
        dtau = dtau - t(X)%*%Q%*%X/2
      }
      if(i==1 && iter%%10==0){
        plot(locs[[i]],Y[[i]])
        lines(locs[[i]],B[[i]]%*%beta + A[[i]]%*%X,col=2)
        Sys.sleep(0)
      }
    }
    sigma = sigma + stepsize*dsigma/nsim
    beta = beta + stepsize * dbeta/nsim
    tau = as.double(tau + stepsize * dtau/nsim)
  } 
  return(list(sigma=exp(sigma),beta=beta,tau=exp(tau)))
}
