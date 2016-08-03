#include "operatorMatrix.h"
#include "error_check.h"

using namespace std;

double fd2Operator::trace_variance( const Eigen::SparseMatrix<double,0,int>& A){
  Eigen::VectorXd  obs_loc = A * loc;
  obs_loc.array() -= m_loc;
  // covariance / h
  // covariance = t^3/3 (integrated brownian moition
  return(obs_loc.array().pow(3).sum() / (3 * h_average * tau));

}