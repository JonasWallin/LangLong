##
# testing gradient estimation
##
library(MASS) 
library(solvercpp)
n <- 10
Niter <- 10
stepsize <- 1/n
Sigma_beta <- matrix(c(0.4,0.12,0.12,0.2), nrow = 2, ncol =  2)

beta_vec <- mvrnorm(n =  n, mu = c(0,0), Sigma  = Sigma_beta)
#cov(beta_vec)
Sigma <- diag(2)
Sigma_vech <- vech_cpp(Sigma)
D <- duplicatematrix_cpp(Sigma)
for(i in 1:Niter){
  stepsize <- 1
  iSigma <- solve(Sigma)
  dSigma_vech <- 0.5 * t(D) %*% kronecker(iSigma, iSigma) %*% vec_cpp(t(beta_vec)%*%beta_vec - n*Sigma)
  e_min <- -1
  Sigma_vech_old <- Sigma_vech
  while(e_min <= 0){
    Sigma_vech <- Sigma_vech_old + stepsize * solve(n*0.5*t(D) %*% kronecker(iSigma, iSigma) %*%D, dSigma_vech)
    print(solve(n*0.5*t(D) %*% kronecker(iSigma, iSigma) %*%D, dSigma_vech))
    Sigma <- veci_cpp( D %*% Sigma_vech, dim(Sigma)[1] )
    e_min <- min(eigen(Sigma)$values)
    stepsize <- 0.5 * stepsize
  }
  print(stepsize * 2)
  print(Sigma- t(beta_vec)%*%beta_vec/n)
}
if(0){
Sigma <-t(beta_vec)%*%beta_vec/n
Sigma_vech <- vech_cpp(Sigma)
for(i in 1:3){
h <- 10^-7
e <- c(0,0,0)
e[i] <- h
loglik <-  0.5 * (-n*log(abs(det(Sigma)) )  -  sum( diag(solve(Sigma) %*% t(beta_vec)%*%beta_vec) ))
Sigma_eps <- veci_cpp( D %*% (Sigma_vech + e), dim(Sigma)[1] )
loglik_eps <-   0.5 * (-n*log(abs(det(Sigma_eps)) )  -  sum(diag( solve(Sigma_eps) %*% t(beta_vec)%*%beta_vec) ))
print((loglik - loglik_eps)/h)
iSigma <- solve(Sigma)
}
dSigma_vech <- 0.5 * t(D) %*% kronecker(iSigma, iSigma) %*% vec_cpp(t(beta_vec)%*%beta_vec - n*Sigma)
print(dSigma_vech)      
}
