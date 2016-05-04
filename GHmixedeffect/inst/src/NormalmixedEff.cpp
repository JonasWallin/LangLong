#include "MixedEffect.h"
#include "error_check.h"
#include  <cmath>
#include <unsupported/Eigen/KroneckerProduct>

#include <Eigen/LU>


Eigen::VectorXd sampleNormalCan(const Eigen::VectorXd & b,const Eigen::MatrixXd & Q)
{
  return b;
}



NormalMixedEffect::NormalMixedEffect(){
  counter = 0;
  //dlog_sigma2  = 0;
  //ddlog_sigma2 = 0;
} 
Rcpp::List NormalMixedEffect::toList()
{
  Rcpp::List out;
  out["B"]     = B;
  out["Sigma"] = Sigma;
  out["U"]     = U;
  return(out);
}
void NormalMixedEffect::initFromList(Rcpp::List const &init_list)
{
  
  std::vector<std::string> check_names =  {"B", "Sigma"};
  check_Rcpplist(init_list, check_names, "NormalmixedEff::initFromList");
  
  int count =0;
  Rcpp::List B_list = init_list["B"];
  B.resize(B_list.length());
  for( Rcpp::List::iterator it = B_list.begin(); it != B_list.end(); ++it ) {
    B[count] = Rcpp::as < Eigen::MatrixXd >( it[0]);
    count++;
  }
  D = duplicatematrix(B[0].cols());
  Dd = D.cast <double> (); 
  Sigma     =  Rcpp::as< Eigen::MatrixXd > (init_list["Sigma"]) ;
  Sigma_vech = vech(Sigma);
  UUt.setZero(Sigma.cols() * Sigma.rows());
  dSigma_vech.setZero(Sigma_vech.size());
  invSigma  = Sigma.inverse();
  if( init_list.containsElementNamed("U" )){
    U = Rcpp::as< Eigen::MatrixXd > (init_list["U"]) ;
  }else{
    U.setZero(B[0].cols(), B.size());
  }
  
}



void NormalMixedEffect::sampleU(const int i, 
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise)
{
  
          Eigen::VectorXd b   = exp( - log_sigma2_noise) * (B[i].transpose() * res);
          Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * B[i].transpose() * B[i];
          Q +=  invSigma;
          U.col(i) = sample_Nc(b, Q);
          
}

void NormalMixedEffect::gradient(const int i, const Eigen::VectorXd& res)
{
    counter++;
    Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
    
    UUt += vec( UUT);
}
void NormalMixedEffect::step_theta(double stepsize)
{
  double pos_def = 0;
  iSkroniS = kroneckerProduct(invSigma, invSigma); 
  UUt -= counter*vec(Sigma); 
  dSigma_vech = 0.5 * Dd.transpose() * iSkroniS * UUt;
  ddSigma = 0.5 * counter * Dd.transpose() * iSkroniS * Dd;
  dSigma_vech = ddSigma.ldlt().solve(dSigma_vech); 
  while(pos_def <= 0){
    Eigen::VectorXd Sigma_vech_temp = Sigma_vech;
    Sigma_vech_temp += stepsize * dSigma_vech;
    Eigen::VectorXd temp = Dd*Sigma_vech_temp;
    Sigma = veci(temp, Sigma.rows(), Sigma.cols());
    stepsize *= 0.5; 
    SelfAdjointEigenSolver<MatrixXd> eig(Sigma,EigenvaluesOnly);
    pos_def = eig.eigenvalues().minCoeff();
    if(stepsize <= 1e-16){
        Rcpp::Rcout << "Sigma = \n" << Sigma << "\n";
        Rcpp::Rcout << "pos_def = " << pos_def <<"\n"; 
        throw("in midexeffect not pos def \n");   
    }
    
   
  }
    
    dSigma_vech.setZero(Sigma_vech.size());
    UUt.setZero(Sigma.cols() * Sigma.rows());
    invSigma  = Sigma.inverse();
    Sigma_vech = vech(Sigma);
    counter = 0;
}

