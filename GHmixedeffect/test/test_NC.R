#
# testing sampling N_C(mu, Q)
#
library(testthat)
library(GHmixedeffect)
Q  <- diag(4)
mu <- rep(4, 4)
n_sample <- 10000
Xs  <- matrix(ncol=4, nrow=n_sample)
for(i in 1:n_sample)
  Xs[i, ] <- sample_Nccpp(mu, Q)


test_that(" testing sample_N_C",
{
  expect_equal(max(abs(cov(Xs)-solve(Q))), 0, tolerance  = 0.1)
  expect_equal(max(abs(colMeans(Xs)- mu)), 0, tolerance  = 0.1)
})

Sigma <- toeplitz(c(4,1,1,0,0))
mu <- c(1, 2, 3, 4, 5)
Q <- solve(Sigma)
b <- Q%*%mu
Xs  <- matrix(ncol=5, nrow=n_sample)
for(i in 1:n_sample)
  Xs[i, ] <- sample_Nccpp(b, Q)

test_that(" testing sample_N_C v2",
{
  expect_equal(max(abs(cov(Xs)-Sigma)), 0, tolerance  = 0.1)
  expect_equal(max(abs(colMeans(Xs)- mu)), 0, tolerance  = 0.1)
})