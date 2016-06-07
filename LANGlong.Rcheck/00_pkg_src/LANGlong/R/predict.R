#' Compute posterior mean and credible bands 
#'
#' @param i number of the patient to predict.
#' @param results file from estimateLongGH 
#' @param n number of basis functions in FEM approx.
#' @param type type of prediction, 'filtering' or 'smoothing' supported.
library(spam)
predictLong <- function(i, resullt, n = -1, type="filter")
{
  if(n == -1)
    n <- result$operator_List$mesh1d$n
  tmp          <- create_matrices_(result$obs_[[i]]$loc, n)
  MatrixBlock  <- tmp$MatrixBlock
  mesh1d       <- tmp$mesh1d
  C = as(as(MatrixBlock$c0,"CsparseMatrix"), "dgCMatrix")
  Ci = as(sparseMatrix(i=1:n,j=1:n,x=1/(MatrixBlock$c0@x),
                                             dims=c(n, n)), "CsparseMatrix")
  G = MatrixBlock$g1

  K = G + result$theta$kappa*C
  Q = result$theta$tau* spam::t(K)%*%Ci%*%K
  A = inla.mesh.1d.A(mesh1d, locs[[i]])
  if(type=="smooth"){
    Q.post = Q + spam::t(A)%*%A/(result$theta$sigma^2)
    mu.post = solve(Q.post,spam::t(A)%*%(Y[[i]]-B[[i]]%*%result$theta$beta)/(result$theta$sigma^2))
    vars <- diag(inla.qinv(Q.post))
  } else if(type == "filter"){
    n = length(mesh1d$loc)
    mu.post = rep(0,n)
    vars = rep(0,n)
    for(j in 1:length(mesh1d$loc)){
      obs.ind = locs[[i]]<= max(mesh1d$loc[1:j])
      Ai = A[obs.ind,,drop=FALSE]
      Qp = Q + spam::t(Ai)%*%Ai/(result$theta$sigma^2)
      mp = solve(Qp,
                 as.vector(spam::t(Ai)%*%(Y[[i]][obs.ind,,drop=F]-B[[i]][obs.ind, ,drop=F]%*%result$theta$beta)/(result$theta$sigma^2)
                ))
      mu.post[j]= mp[j]
      vars[j] = diag(inla.qinv(Qp))[j]
    }
  } else {
    stop("Only smoothing and filtering supported.")
  }
  
  out <- list(mu = as.double( mu.post),
              upper = as.double(mu.post + 1.96*sqrt(vars)),
              lower = as.double( mu.post - 1.96*sqrt(vars)),
              vars = as.double(vars),
              loc = mesh1d$loc,
              Y = Y[[i]],
              Yloc = locs[[i]])
  

  
  return(out)
}