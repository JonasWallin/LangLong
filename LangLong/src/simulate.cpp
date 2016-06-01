
#include <Rcpp.h>
#include <random>
#include <chrono>
#include "Qmatrix.h"
#include "solver.h"
using namespace Rcpp;



// [[Rcpp::export]]
List simulateLong_cpp(List obs_, List operator_, List theta_)
{
  
  List obs_list      = clone(obs_);
  List operator_list = clone(operator_);
  List theta_list    = clone(theta_);
  
  /*
    CREATING THE OPERATOR AND SOLVER
  */
  Qmatrix* Kobj;   
  Kobj     = new MaternMatrixOperator;
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  cholesky_solver Solver;
  Solver.init(Kobj->d, 0, 0, 0);
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd>(theta_list["kappa"]);
  Kobj->vec_to_theta(kappa.log()); 
  /*
    Simulation object
  */
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  Eigen::VectorXd  z;
  z.resize( Kobj->d);
  
  /*
    Building Q matrix
  */
  Eigen::SparseMatrix<double,0,int> iC  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(operator_["Ci"]);
  Eigen::SparseMatrix<double,0,int> Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
  Q  = Q  * iC;
  double tau = theta_list["tau"];
  Q  = tau*Q  * Kobj->Q;
  Solver.analyze(Q);
  Solver.compute(Q);
  
  
  Eigen::VectorXd  X, mu;
  mu.resize( Kobj->d);
  mu.setZero( Kobj->d);
  
  //
  List out(obs_list.length());
  int counter = 0;
  double sigma_eps = Rcpp::as<double>(theta_list["sigma"]);
  Eigen::VectorXd beta      = Rcpp::as<Eigen::VectorXd>(theta_list["beta"]);
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
     List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    
    
    for(int i =0; i < Kobj->d; i++)
      z[i] =  normal(random_engine);
    
    X = Solver.rMVN(mu, z);
      
    Eigen::SparseMatrix<double,0,int> A  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    Eigen::MatrixXd B  = Rcpp::as<Eigen::MatrixXd >(obs_tmp["B"]);
    Eigen::VectorXd Y = A*X;
    Y += B * beta;
    for(int i =0; i < Y.size(); i++)
      Y[i] += normal(random_engine) * sigma_eps;
    
    obs_list[counter++] = Y;
  }
  
  return(obs_list);
}



// [[Rcpp::export]]
List testSimulateX_cpp( List operator_, List theta_)
{
  
  List operator_list = clone(operator_);
  List theta_list    = clone(theta_);
  
  /*
    CREATING THE OPERATOR AND SOLVER
  */
  Qmatrix* Kobj;   
  Kobj     = new MaternMatrixOperator;
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  cholesky_solver Solver;
  Solver.init(Kobj->d, 0, 0, 0);
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd>(theta_list["kappa"]);
  Kobj->vec_to_theta(kappa.log()); 
  
  /*
    Simulation object
  */
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  Eigen::VectorXd  z;
  z.resize( Kobj->d);
  
  /*
    Building Q matrix
  */
  Eigen::SparseMatrix<double,0,int> iC  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(operator_["Ci"]);
  Eigen::SparseMatrix<double,0,int> Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
  Q  = Q  * iC;
  double tau = theta_list["tau"];
  Q  = tau* (Q * Kobj->Q);
  Solver.analyze(Q);
  Solver.compute(Q);
  
  
  Eigen::VectorXd  X, mu;
  mu.resize( Kobj->d);
  mu.setZero(Kobj->d);
  
  //
    
  for(int i =0; i < Kobj->d; i++)
    z[i] =  normal(random_engine);
  
  
  X = Solver.rMVN(mu, z);
  List out_list;
  out_list["X"] = X;
  Q /= tau;
  out_list["Q"] = Q;
  return(out_list);
}
// [[Rcpp::export]]
List testSimulateX2_cpp( List operator_, List theta_)
{
  
  List operator_list = clone(operator_);
  List theta_list    = clone(theta_);
  
  /*
    CREATING THE OPERATOR AND SOLVER
  */
  Qmatrix* Kobj;   
  Kobj     = new MaternMatrixOperator;
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  cholesky_solver Solver;
  Solver.init(Kobj->d, 0, 0, 0);
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd>(theta_list["kappa"]);
  Kobj->vec_to_theta(kappa.log()); 
  
  /*
    Simulation object
  */
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  Eigen::VectorXd  z;
  z.resize( Kobj->d);
  
  /*
    Building Q matrix
  */
  Eigen::SparseMatrix<double,0,int> iC  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(operator_["Ci"]);
  Eigen::SparseMatrix<double,0,int> Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
  Eigen::SparseMatrix<double,0,int> A = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(operator_["A"]);
  Eigen::VectorXd Y = Rcpp::as<Eigen::VectorXd >(operator_["Y"]);

  Q  = Q  * iC;
  double tau = theta_list["tau"];
  Q  = tau* (Q * Kobj->Q);
  Q = Q + (A.transpose()*A)/Rcpp::as<double>(theta_list["sigma"]);
  Solver.analyze(Q);
  Solver.compute(Q);
  
  
  Eigen::VectorXd  X, mu;
  mu.resize( Kobj->d);
  mu.setZero(Kobj->d);
  
    
  for(int i =0; i < Kobj->d; i++)
    z[i] =  normal(random_engine);
  
  mu = A.transpose()*Y/Rcpp::as<double>(theta_list["sigma"]);
  X = Solver.rMVN(mu, z);
  List out_list;
  out_list["X"] = X;
  out_list["Q"] = Q;
  out_list["XtQX"] = 0.5 * X.dot(X.transpose() *  Q);
  return(out_list);
}


