#
# test sampling N_C(mu, Q)
#
Q  <- diag(4)
mu <- rep(4, 4)
n_sample <- 10000
Xs  <- matrix(ncol=4, nrow=n_sample)
for(i in 1:n_sample)
  Xs[i, ] <- sample_Nccpp(mu, Q)

print(paste('max diff Q = ',max(abs(cov(Xs)-solve(Q)))))
print(paste('max diff mu = ',max(abs(colMeans(Xs)- mu))))

Sigma <- toeplitz(c(4,1,1,0,0))
mu <- c(1, 2, 3, 4, 5)
Q <- solve(Sigma)
b <- Q%*%mu
Xs  <- matrix(ncol=5, nrow=n_sample)
for(i in 1:n_sample)
  Xs[i, ] <- sample_Nccpp(b, Q)

print(paste('max diff Q = ',max(abs(cov(Xs)-Sigma))))
print(paste('max diff mu = ',max(abs(colMeans(Xs)- mu))))