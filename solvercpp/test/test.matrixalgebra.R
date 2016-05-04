library(testthat)
library(solvercpp)
library(Matrix)
context("matrixAlgebra")

n <- sample(3:20, 1)

test_that("matrixAlgebra, duplicatematrix", {
  A <- matrix(rnorm(n*n), n, n)
  AtA <- t(A)%*%A
  vech_AtA <- vech_cpp(AtA)
  vec_AtA <- vec_cpp(AtA)
  D <- duplicatematrix_cpp(AtA)

  expect_equal(c(D%*%vech_AtA),vec_AtA)
})

test_that("matrixAlgebra, veci", {
  A <- matrix(rnorm(n*n), n, n)
  AtA <- t(A)%*%A
  AtA2 <- veci_cpp(vec_cpp(AtA), n)
  expect_equal(AtA2, AtA)
})
