
library(testthat)
library(solvercpp)
library(Matrix)
context("Solver")

############################################
###   #TEST FOR LU
############################################
test_that("Solver, solve LU", {
  A <- matrix(rnorm(9),3,3)
  expect_equal(logdetWithLU(A),log(abs(det(A))),tolerance=1e-7)
})
test_that("Solver, solve LU", {
  A <- matrix(rnorm(9),3,3)
  b <- rnorm(3)
  expect_equal(solveWithLU(A,b),solve(A,b),tolerance=1e-7)
})
test_that("trace, solve LU", {
  A <- matrix(rnorm(9),3,3)
  B <- matrix(rnorm(9),3,3)
  trace  <- (sum(diag(B%*%solve(A))))
  trace_LU <- (traceWithLU(A, B))
  expect_equal(trace, trace_LU,tolerance=1e-7)
})
test_that("trace sparse, solve LU", {
  n <- 11
  A <- matrix(rnorm(n*n),n,n)
  B = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  trace  <- (sum(diag(B%*%solve(A))))
  trace_LU <- (traceWithLUSparse(A, as(B, "dgCMatrix")))
  expect_equal(trace, trace_LU,tolerance=1e-7)
})
#####################################
#####################################

test_that("Solver, solve Cholesky", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  v = seq(-5, 5, length=n)
  datalist = list(Q=as(Q, "dgCMatrix"), v=v, type=0, trace.iter=10,
                  tol=0.001,solver.max.iter=1000, operation=0)
  res = solver_tests(datalist)
  X = as.vector(solve(Q,v))
  expect_equal(res$X,X,tolerance=1e-7)
})




test_that("Solver, solve iterative", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  v = seq(-5, 5, length=n)
  datalist = list(Q=as(Q, "dgCMatrix"), v=v, type=1, trace.iter=1000,
                  tol=0.00001,solver.max.iter=1000, operation=0)
  res = solver_tests(datalist)
  X = as.vector(solve(Q,v))
  expect_equal(res$X,X,tolerance=1e-7)
})


test_that("Solver, trace Cholesky", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  M = Matrix(toeplitz(c(1, 0.5, rep(0, n-2))))
  datalist = list(Q=as(Q, "dgCMatrix"), M=as(M, "dgCMatrix"), type=0,
                  trace.iter=1000, tol=0.00001,
                  solver.max.iter=1000, operation=1)
  res = solver_tests(datalist)
  tr = sum(diag(solve(Q,M)))
  expect_equal(res$trace,tr,tolerance=1e-7)
})

test_that("Solver, trace iter", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  M = Matrix(toeplitz(c(1, 0.5, rep(0, n-2))))
  datalist = list(Q=as(Q, "dgCMatrix"), M=as(M, "dgCMatrix"), type=1,
                  trace.iter=10000, tol=1e-10,
                  solver.max.iter=10000, operation=1)
  res = solver_tests(datalist)
  tr = sum(diag(solve(Q,M)))

  #N = 100000
  #U = matrix(sample(c(-1,1),n*N,replace=1),n,N)
  #QU = solve(Q,U)
  #MQU = M%*%QU
  #tr2 = sum(U*MQU)/N
  expect_equal(res$trace,tr,tolerance=0.1)
})

test_that("Solver, trace2 Cholesky", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  M1 = Matrix(toeplitz(c(1, 0.5, rep(0, n-2))))
  M2 = Matrix(toeplitz(c(1, 0.2, rep(0, n-2))))
  datalist = list(Q=as(Q, "dgCMatrix"), M1=as(M1, "dgCMatrix"),
                  M2=as(M2, "dgCMatrix"),type=0,
                  trace.iter=1000, tol=0.00001,
                  solver.max.iter=1000, operation=2)
  res = solver_tests(datalist)
  tr = sum(diag(solve(Q,M1)%*%solve(Q,M2)))
  expect_equal(res$trace2,tr,tolerance=1e-7)
})


test_that("Solver, trace2 iter", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  M1 = Matrix(toeplitz(c(1, 0.5, rep(0, n-2))))
  M2 = Matrix(toeplitz(c(1, 0.2, rep(0, n-2))))
  datalist = list(Q=as(Q, "dgCMatrix"), M1=as(M1, "dgCMatrix"),
                  M2=as(M2, "dgCMatrix"),type=1,
                  trace.iter=10000, tol=0.00001,
                  solver.max.iter=10000, operation=2)
  res = solver_tests(datalist)
  tr = sum(diag(solve(Q,M1)%*%solve(Q,M2)))
  expect_equal(res$trace2,tr,tolerance=0.1)
})


test_that("Solver, Qinv_diag", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  datalist = list(Q=as(Q, "dgCMatrix"),type=0,
                  trace.iter=1000, tol=0.00001,
                  solver.max.iter=1000, operation=3)
  res = solver_tests(datalist)
  vars = diag(solve(Q))
  expect_equal(res$vars,vars,tolerance=1e-7)
})

test_that("Solver, logdet", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  datalist = list(Q=as(Q, "dgCMatrix"),type=0,
                  trace.iter=1000, tol=0.00001,
                  solver.max.iter=1000, operation=4)
  res = solver_tests(datalist)
  ld = log(det(Q))
  expect_equal(res$logdet,ld,tolerance=1e-7)
})


test_that("Solver, Qinv", {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  v = seq(-5, 5, length=n)
  datalist = list(Q=as(Q, "dgCMatrix"), v=v, type=0, trace.iter=10,
                  tol=0.001,solver.max.iter=1000, operation=5)
  res = solver_tests(datalist)
  Q_solve <- solve(Q)
  index <-res$Qinv>0
  expect_equal(Q_solve[index],res$Qinv[index],tolerance=1e-7)
})
