#include <Rcpp.h>
#include <random>
#include <chrono>
#include <vector>
#include "Qmatrix.h"
#include "solver.h"
using namespace Rcpp;

// [[Rcpp::export]]
List estimateLong_cpp(List obs_list, List operator_list, List theta_list, double stepsize, int Niter, int nsim, int silent)
{
  

  //count number of patients  
  int nrep = obs_list.length();
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd> (theta_list["kappa"]);
  double tau = theta_list["tau"];
  Eigen::VectorXd  beta = Rcpp::as<Eigen::VectorXd>( theta_list["beta"]);
  double sigma = theta_list["sigma"];
  
  
  //Define operator
  Qmatrix* Kobj;   
  Kobj = new MaternMatrixOperator;
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  Kobj->vec_to_theta( kappa);
  
  //Prior solver
  cholesky_solver Qsolver;
  Qsolver.init(Kobj->d, 0, 0, 0);
  Qsolver.analyze(Kobj->Q);
  Qsolver.compute(Kobj->Q);
  
  Eigen::SparseMatrix<double,0,int> iC  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(operator_list["Ci"]);
  
  //Create solvers for each patient
  std::vector<  cholesky_solver >  Solver(nrep);
  int i = 0;
  Eigen::SparseMatrix<double,0,int> Q;
  std::vector< Eigen::SparseMatrix<double,0,int> > As(nrep);
  std::vector< Eigen::VectorXd > Ys(nrep);
  std::vector< Eigen::MatrixXd > Bs(nrep);
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    As[i] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    Ys[i] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    Bs[i] = Rcpp::as<Eigen::MatrixXd>(obs_tmp["B"]);
    
    Solver[i].init(Kobj->d, 0, 0, 0);
    Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Q = Q  * iC;
    Q = Q  * Kobj->Q;
    Q = Q + As[i].transpose()*As[i];
    Solver[i].analyze(Q);
    Solver[i].compute(Q);
    i++;
  }
  
  /*
  Simulation object
  */
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  Eigen::VectorXd  z;
  z.resize( Kobj->d);
  z.setZero(Kobj->d);
  
  Eigen::VectorXd  X, b, Ysim;
  X.resize( Kobj->d);
  X.setZero(Kobj->d);
  
  b.resize( Kobj->d);
  b.setZero(Kobj->d);
  
  
  
  Eigen::SparseMatrix<double,0,int> Qi;
  Eigen::VectorXd dbeta, dkappa;
  dbeta.resize( beta.size());
  dbeta.setZero(beta.size());
  dkappa.resize( kappa.size());
  dkappa.setZero( kappa.size());
  double dtau, dsigma;
  
  Eigen::MatrixXd B;
  Eigen::SparseMatrix<double,0,int> A;
  
  double d2sigma, d2tau, d2kappa;
  Eigen::MatrixXd d2beta;
  d2beta.setZero(beta.size(),beta.size());
  for(int iter=0;iter< Niter;iter++){
    if(silent == 0){
      Rcpp::Rcout << "i = " << iter << ": ";
      Rcpp::Rcout << "sigma = " << exp(0.5 * sigma);
      Rcpp::Rcout << ", beta = " << beta.transpose();
      Rcpp::Rcout << ", tau = " << exp(tau);
      Rcpp::Rcout << ", kappa = " << kappa.exp() << "\n";
    }
    
    //Update operator
    Kobj->vec_to_theta( kappa);
    Qsolver.compute(Kobj->Q);
    
    Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Q = Q  * iC;
    Q = exp(tau)*Q  * Kobj->Q;
    
    dtau = nsim*(nrep*Kobj->d)/2.0;
    
    dsigma = 0;
    dkappa[0] = nsim*nrep*Kobj->trace(0);
    dbeta.setZero(beta.size());
    
    d2sigma = 0;
    d2tau = 0;
    d2beta.setZero(beta.size(),beta.size());
    for( int i = 0; i < nrep; i++ ) {
      A = As[i];
      B  = Bs[i];
      Eigen::VectorXd Y = Ys[i];
      Qi = Q + (A.transpose()*A)/exp(sigma);
      
      Solver[i].compute(Qi);
      
      //compute mean
      Y -= B * beta;
      b = A.transpose()*Y/exp(sigma);
      
      for(int ii = 0; ii < nsim; ii ++){
        //Sample X|Y
        for(int j =0; j < Kobj->d; j++)
          z[j] =  normal(random_engine);
        
        X = Solver[i].rMVN(b, z);
        
        //Compute gradients
        Ysim = Y - A*X;
      
        dsigma += -Ysim.size()/2.0 + exp(-sigma)*Ysim.array().square().sum()/2.0;
        dbeta += exp(-sigma)* (B.transpose() * Ysim); 
        
        dkappa -= exp(tau)*X.transpose() * (Kobj->df(0) * iC * Kobj->Q) * X;
        Eigen::VectorXd vtmp = X.transpose()*Q;
        dtau -=  vtmp.dot(X)/2.0;
        
        //Hessian
        d2sigma -= exp(-sigma)*Ysim.array().square().sum()/2.0;
        d2beta  -= exp(-sigma)*B.transpose()*B;
        d2tau   -= vtmp.dot(X)/2.0;
      }
    } 
    //Rcpp::Rcout << d2sigma << " " << d2beta << " "<< d2tau << "\n";
    //take step in parameters
    sigma -= (dsigma/nsim)/d2sigma;
    beta  -= (d2beta.inverse()*dbeta)/nsim;
    tau   -= (dtau/nsim)/d2tau;
    kappa += dkappa * stepsize/nsim;
  }  
    
  Rcpp::List list;
  list["kappa"] = kappa.exp();
  list["sigma"] = exp(0.5*sigma);
  list["beta"]  = beta;
  list["tau"]   = exp(tau);
  
  return(list);
}

