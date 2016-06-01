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
  noise = "Normal";
  //dlog_sigma2  = 0;
  //ddlog_sigma2 = 0;
} 
Rcpp::List NormalMixedEffect::toList()
{
  Rcpp::List out;
  out["Bf"]          = Bf;
  out["Br"]          = Br;
  out["beta_random"] = beta_random;
  out["beta_fixed"]  = beta_fixed;
  out["Sigma"]       = Sigma;
  out["U"]           = U;
  out["noise"]       = noise;
  return(out);
}
void NormalMixedEffect::initFromList(Rcpp::List const &init_list)
{
  
 
  int count =0;
  if(init_list.containsElementNamed("B_fixed"))
  {
    Rcpp::List Bf_list = init_list["B_fixed"];
    Bf.resize(Bf_list.length());
    for( Rcpp::List::iterator it = Bf_list.begin(); it != Bf_list.end(); ++it ) {
      Bf[count++] = Rcpp::as < Eigen::MatrixXd >( it[0]);
    }
    grad_beta_f.setZero(Bf[0].cols());
    
    if(init_list.containsElementNamed("beta_fixed"))
      beta_fixed = Rcpp::as < Eigen::VectorXd >( init_list["beta_fixed"]);
    else
      beta_fixed.setZero(Bf[0].cols());
    
     H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());
  
    
  }else{ Bf.resize(0);}
  count = 0;
  
  
  if(init_list.containsElementNamed("B_random"))
  {
    Rcpp::List Br_list = init_list["B_random"];
    Br.resize(Br_list.length());
    for( Rcpp::List::iterator it = Br_list.begin(); it != Br_list.end(); ++it ) {
      Br[count++] = Rcpp::as < Eigen::MatrixXd >( it[0]);
    }
    grad_beta_r.setZero(Br[0].cols());
    grad_beta_r2.setZero(Br[0].cols());
    if(init_list.containsElementNamed("beta_random"))
      beta_random = Rcpp::as < Eigen::VectorXd >( init_list["beta_random"]);
    else
      beta_random.setZero(Br[0].cols());
  }else{ Br.resize(0);}
  
  H_beta_random.setZero(Br[0].cols(), Br[0].cols());
  if(Br[0].cols() > 0){
    D = duplicatematrix(Br[0].cols());
    Dd = D.cast <double> (); 
  }
  
  if(Br.size() > 0){
    if(init_list.containsElementNamed("Sigma"))
      Sigma     =  Rcpp::as< Eigen::MatrixXd > (init_list["Sigma"]) ;
    else
      Sigma.setIdentity(Br[0].cols(), Br[0].cols());
    
    Sigma_vech = vech(Sigma);
    UUt.setZero(Sigma.cols() * Sigma.rows());
    dSigma_vech.setZero(Sigma_vech.size());
    invSigma  = Sigma.inverse();
    if( init_list.containsElementNamed("U" ))
      U = Rcpp::as< Eigen::MatrixXd > (init_list["U"]);
    else
      U.setZero(Br[0].cols(), Br.size());
  }
  
}



void NormalMixedEffect::sampleU2(const int i, 
                                const Eigen::VectorXd& res,
                                const Eigen::VectorXd& iV,
                                const double log_sigma2_noise //= 0
                                )
{
    if(Br.size() == 0)
      return;

	
    Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res));
    Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * iV.asDiagonal() * Br[i];
    Q        +=  invSigma;
    U.col(i) =   sample_Nc(b, Q); 
}

void NormalMixedEffect::sampleU(const int i, 
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise )
{
    if(Br.size() == 0)
      return;

    Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * res);
    Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * Br[i];
    Q        +=  invSigma;
    U.col(i) =   sample_Nc(b, Q); 
}

void NormalMixedEffect::gradient2(const int i, 
                                 const Eigen::VectorXd& res,
                                 const Eigen::VectorXd& iV,
                                 const double log_sigma2_noise,  // = 0
                                 const double EiV // = 0
                                 )
{
    counter++;
    Eigen::VectorXd res_  = res;
    if(Br.size() > 0){
      res_ -= Br[i] * U.col(i);
      Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
      UUt += vec( UUT);
      grad_beta_r  += exp( - log_sigma2_noise) * (Br[i].transpose() * iV.cwiseProduct(res_));
      grad_beta_r2 += (invSigma * U.col(i));
      //H_beta_random +=   exp( - log_sigma2_noise) * (Br[i].transpose() * iV.asDiagonal()* Br[i]);
      H_beta_random +=  EiV *exp( - log_sigma2_noise) * (Br[i].transpose()* Br[i]);
    }
    if(Bf.size() > 0){
      grad_beta_f   +=  exp( - log_sigma2_noise) * (Bf[i].transpose() *  iV.cwiseProduct(res_));
      //H_beta_fixed  +=  exp( - log_sigma2_noise) * (Bf[i].transpose() *iV.asDiagonal()* Bf[i]);
      H_beta_fixed  +=  EiV * exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
    }
}
void NormalMixedEffect::gradient(const int i, 
                                 const Eigen::VectorXd& res,
                                 const double log_sigma2_noise)
{
    counter++;
    Eigen::VectorXd res_  = res;
    if(Br.size() > 0){
      res_ -= Br[i] * U.col(i);
      Eigen::MatrixXd UUT = U.col(i) * U.col(i).transpose();
      UUt += vec( UUT);
      grad_beta_r  += exp( - log_sigma2_noise) * (Br[i].transpose() * res_);
      grad_beta_r2 +=  (invSigma * U.col(i));
      H_beta_random +=  exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);
      
    }
    if(Bf.size() > 0){
      grad_beta_f   +=  exp( - log_sigma2_noise) * (Bf[i].transpose() * res_);
      H_beta_fixed  +=  exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
    }
}
void NormalMixedEffect::step_theta(double stepsize)
{
  if(Br.size() > 0){
    step_beta_random(stepsize);
    step_Sigma(stepsize);
  }
  if(Bf.size() > 0)
    step_beta_fixed(stepsize);
  
  
  counter = 0;
}
void NormalMixedEffect::step_beta_fixed(double stepsize)
{
    beta_fixed += stepsize *  H_beta_fixed.ldlt().solve(grad_beta_f);
    grad_beta_f.setZero(Bf[0].cols());
    H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());
  
}
void NormalMixedEffect::step_beta_random(double stepsize)
{
    beta_random += (0.5 * stepsize) *  H_beta_random.ldlt().solve(grad_beta_r);
    beta_random += (0.5 * stepsize) * (Sigma * grad_beta_r2)/ counter;
    grad_beta_r.setZero(Br[0].cols());
    grad_beta_r2.setZero(Br[0].cols());
    H_beta_random.setZero(Br[0].cols(), Br[0].cols());
}

void NormalMixedEffect::step_Sigma(double stepsize)
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
}

void NormalMixedEffect::remove_cov(const int i, Eigen::VectorXd & Y)
{
  if(Br.size() > 0 )
    Y -= Br[i] * beta_random;
  if(Bf.size() > 0)
    Y -= Bf[i] * beta_fixed;
}
void NormalMixedEffect::add_cov(const int    i, Eigen::VectorXd & Y)
{
  if(Br.size() > 0 )
    Y += Br[i] * beta_random;
  if(Bf.size() > 0)
    Y += Bf[i] * beta_fixed;
}