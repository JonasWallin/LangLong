#include <Rcpp.h>
#include "MatrixAlgebra.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
Eigen::VectorXd vech_cpp(Eigen::MatrixXd A)
{
  return(vech(A));
}

// [[Rcpp::export]]
Eigen::MatrixXd veci_cpp(Eigen::VectorXd vecA, int n)
{
  return(veci(vecA, n, n));
}

// [[Rcpp::export]]
Eigen::VectorXd vec_cpp(Eigen::MatrixXd A)
{
  return(vec(A));
}
// [[Rcpp::export]]
Eigen::MatrixXi duplicatematrix_cpp(Eigen::MatrixXd A)
{
  return(duplicatematrix(A.rows()));
}

