library(INLA)

create_operator <- function(locs, n, name = "matern")
{
  if(tolower(name) == "matern"){
    return(create_matrices_Matern(locs, n))
  }else{
    return(create_matrices_FD2(locs, n))
  }
  
}
#' creates matrices for Matern 1D operator
#'
#' @return operator_List list to to use in simulation and estimation object
create_matrices_Matern <- function(locs, n)
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
  operator_List <- list(type = 'matern', 
                        mesh1d = mesh1d,
                        C = as(as(MatrixBlock$c0,"CsparseMatrix"), "dgCMatrix"),
                        Ci = as(sparseMatrix(i=1:n,j=1:n,x=1/(MatrixBlock$c0@x),
                                             dims=c(n, n)), "CsparseMatrix"),
                        G = MatrixBlock$g1,
                        h = MatrixBlock$c0@x,
                        B.kappa = as.matrix(rep(1, dim(MatrixBlock$c0)[1])),
                        kappa = 0,
                        loc   = mesh1d$loc
  )
  return(operator_List)
}
#' creates matrices for Finite difference operator, one sided 
#'
create_matrices_FD2 <- function(locs, n)
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
  vec_toeplitz <- rep(0, length=n)
  h <- (P[2] - P[1])
  vec_toeplitz[1] <- -1 / h # -1/h
  vec_toeplitz[2] <- 1  / h  # 1/h
  Operator_1D <- Matrix(toeplitz(vec_toeplitz), sparse=T)
  Operator_1D[upper.tri(Operator_1D)] = 0
  Operator_2D <- (h*Operator_1D) %*% Operator_1D #mult with h for the first operator is scaled
  # due to W(t+h) - W(t) = N(0,h), not (W(t+h) - W(t))/h = N(0,h)
  operator_List <- list(type   = 'fd2',
                        mesh1d = mesh1d,
                        Q      = as(Operator_2D, "dgCMatrix"), 
                        h      = rep(h, n),
                        loc   = mesh1d$loc
                        )
  return(operator_List)
}