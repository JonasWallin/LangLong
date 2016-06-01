
#include <Rcpp.h>
#include <random>
#include <chrono>
#include <string>
#include "Qmatrix.h"
#include "solver.h"
#include "operator_helper.h"
#include "GIG.h"
using namespace Rcpp;


List simulateLongV_cpp(List obs_list, List operator_list, List theta_list)
{
  /*
    CREATING THE OPERATOR AND SOLVER
  */
  std::string type_operator = Rcpp::as<std::string>(operator_list["type"]); 
  Qmatrix* Kobj;   
  //Kobj = new MaternMatrixOperator;
  operator_select(type_operator, &Kobj);
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  cholesky_solver Solver;
  Solver.init(Kobj->d, 0, 0, 0);
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd>(theta_list["kappa"]);
  Kobj->vec_to_theta(kappa.log()); 
  gig rgig;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  
  
  Eigen::VectorXd h  = Rcpp::as<Eigen::VectorXd>(operator_list["h"]);
  double lambda      = Rcpp::as<double>(theta_list["lambda"]);
  int    GAL         = Rcpp::as<double>(theta_list["GAL"]); 
  //
  List out;
  List out_V(obs_list.length());
  int counter = 0;
  
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
     List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    
    Eigen::VectorXd V = sampleV_pre(rgig, h, lambda ,GAL );

    out_V[counter]   = V;
  }
  return(out_V);
}

// [[Rcpp::export]]
List simulateLongGH_cpp(List obs_list, List operator_list, List theta_list)
{
  
  /*
    CREATING THE OPERATOR AND SOLVER
  */
  std::string type_operator = Rcpp::as<std::string>(operator_list["type"]); 
  Qmatrix* Kobj;   
  //Kobj = new MaternMatrixOperator;
  operator_select(type_operator, &Kobj);
  
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  cholesky_solver Solver;
  Solver.init(Kobj->d, 0, 0, 0);
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd>(theta_list["kappa"]);
  Kobj->vec_to_theta(kappa.log()); 
  Eigen::VectorXd h      = Rcpp::as<Eigen::VectorXd>(operator_list["h"]);
  /*
    Simulation object
  */
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  gig rgig;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  
  Eigen::VectorXd  z;
  z.resize( Kobj->d);
 
  Eigen::VectorXd  X, mu;
  mu.setZero( Kobj->d);
  
  double lambda      = Rcpp::as<double>(theta_list["lambda"]);
  double tau         = Rcpp::as<double>(theta_list["tau"]);
  std::string noise = Rcpp::as<std::string>(theta_list["noise"]);
  int GAL = 0;
  if(noise == "GAL")
    GAL = 1;
  if(noise == "Normal")
    GAL = -1;
    
  int    RandomInter = Rcpp::as<double>(theta_list["randomIntercept"]);
  int sigma_rinter   = 0;
  if(RandomInter)
    sigma_rinter   = Rcpp::as<double>(theta_list["sigma_r"]); 
  //
  List out;
  List out_Y(obs_list.length());
  List out_X(obs_list.length());
  List out_V(obs_list.length());
  int counter = 0;
  double sigma_eps = Rcpp::as<double>(theta_list["sigma"]);
  Eigen::VectorXd beta      = Rcpp::as<Eigen::VectorXd>(theta_list["beta"]);
  
  Eigen::VectorXd V,iV;
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
     List obs_tmp = Rcpp::as<Rcpp::List>(*it);
      
    if(GAL >= 0){
      V = sampleV_pre(rgig, h, lambda ,GAL );  
      iV.resize(V.size());
      iV.array() = V.array().inverse();
    } else {
      V = h;
      iV.array() = h.array().inverse();
    }
    for(int i =0; i < Kobj->d; i++)
      z[i] =   normal(random_engine);
      
    Eigen::SparseMatrix<double,0,int> Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Q = Q  * iV.asDiagonal();
    Q = tau * Q  * Kobj->Q;     
    Solver.analyze(Q);
    Solver.compute(Q);
    X = Solver.rMVN(mu, z);
      
    Eigen::SparseMatrix<double,0,int> A  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    Eigen::MatrixXd B  = Rcpp::as<Eigen::MatrixXd >(obs_tmp["B"]);
    Eigen::VectorXd Y = A*X;
    Y += B * beta;
    if(RandomInter)
      Y.array() += sigma_rinter*normal(random_engine);
      
    for(int i =0; i < Y.size(); i++)
      Y[i] += normal(random_engine) * sigma_eps;
    
    out_X[counter]   = X;
    out_V[counter]   = V;
    out_Y[counter++] = Y;
  }
  out["Y"] = out_Y;
  out["X"] = out_X;
  out["V"] = out_V;
  return(out);
}

