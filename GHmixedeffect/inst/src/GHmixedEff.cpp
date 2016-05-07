#include "MixedEffect.h"
#include "error_check.h"
#include <chrono>

NIGMixedEffect::NIGMixedEffect(){
  counter = 0;
  noise = "NIG";
} 
Rcpp::List NIGMixedEffect::toList()
{
  Rcpp::List out;
  out["Bf"]          = Bf;
  out["Br"]          = Br;
  out["beta_random"] = beta_random;
  out["beta_fixed"]  = beta_fixed;
  out["Sigma"]  = Sigma;
  out["U"]      = U;
  out["V"]      = V;
  out["nu"]     = nu;
  out["mu"]     = mu;
  return(out);
}
void NIGMixedEffect::initFromList(Rcpp::List const &init_list)
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
   
     H_beta_random.setZero(Br[0].cols(), Br[0].cols());
    if(Br[0].cols() > 0){
      D = duplicatematrix(Br[0].cols());
      Dd = D.cast <double> (); 
    }
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
    //
    // variance components
    if( init_list.containsElementNamed("V" ))
      V = Rcpp::as< Eigen::MatrixXd > (init_list["V"]) ;
    else
      V.setOnes(Br.size(), 1);
      
    if( init_list.containsElementNamed("mu" ))
      mu = Rcpp::as< Eigen::MatrixXd > (init_list["mu"]) ;
    else
      mu.setZero(Br[0].cols(), 1); 
      
    gradMu.setZero(Br[0].cols(), 1);
    gradMu_2.setZero(Br[0].cols(), 1);
  
    if( init_list.containsElementNamed("nu" ))
      nu = Rcpp::as< double > (init_list["nu"]) ;
    else
      nu = 1.;
      
    grad_nu = 0.;
    EV  = 1.; 
    EiV = 1. + 1./nu;  
    VV  = 1./nu;
      
    
    a_GIG = mu.transpose() * (invSigma *  mu);
    a_GIG += nu;
    rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      
  }else{ Br.resize(0);}
}


void NIGMixedEffect::sampleV(const int i)
{
  		  //  GIG	(p, a, b)
        double p  = 0.5 * (- 1 - Br[0].cols());
        Eigen::VectorXd U_ = U.col(i) + mu;
        double b  =  U_.transpose() * invSigma * U_;
        b += nu;
        V(i) = rgig.sample(p, a_GIG, b);
}

void NIGMixedEffect::sampleU(const int i, 
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise)
{
  		 if(Br.size() > 0){
          Eigen::VectorXd b   = exp( - log_sigma2_noise) * (Br[i].transpose() * res);
          Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * Br[i].transpose() * Br[i];
          Eigen::MatrixXd Qp  = invSigma / V(i);
          b += Qp * (- mu + V(i) * mu);
          Q += Qp;
          U.col(i) = sample_Nc(b, Q);
          sampleV(i);
  		 }
}

void NIGMixedEffect::gradient(const int i, 
                              const Eigen::VectorXd& res,
                              const double log_sigma2_noise)
{
    Eigen::VectorXd res_  = res;
    if(Br.size() > 0){
      
      Eigen::VectorXd U_ = U.col(i) - (-1 + V(i)) * mu; 
      gradient_sigma(i, U_);
      
      
     
      res_ -= Br[i] * U.col(i);
      grad_beta_r  += exp( - log_sigma2_noise) * (Br[i].transpose() * res_);
      grad_beta_r2 +=  (invSigma * U_)/V(i);
      H_beta_random +=  exp( - log_sigma2_noise) * (Br[i].transpose() * Br[i]);


      gradMu   += ((-1 + V(i) )/V(i) ) * (invSigma * U_); 
      
      gradMu_2 += (-1 + V(i) ) * exp( - log_sigma2_noise) * (Br[i].transpose() * res_);
      
      
      // dnu
      grad_nu += 0.5 * (1. / nu - V(i) - 1. / V(i) + 2. );
    }
    if(Bf.size() > 0){
      grad_beta_f   +=  exp( - log_sigma2_noise) * (Bf[i].transpose() * res_);
      H_beta_fixed  +=  exp( - log_sigma2_noise) * (Bf[i].transpose() * Bf[i]);
    }
    counter++;
}
void NIGMixedEffect::gradient_sigma(const int i, Eigen::VectorXd& U_ )
{
  Eigen::MatrixXd UUT =  (U_ * U_.transpose());
  UUT.array() /= V(i);
  UUt    += vec( UUT);
}
void NIGMixedEffect::step_theta(double stepsize)
{
  if(Br.size() > 0){
    step_beta_random(stepsize);
    step_mu(stepsize);
    step_Sigma(stepsize);
    step_nu(stepsize);
    a_GIG = mu.transpose() * (invSigma * mu);
    a_GIG += nu;
    H_beta_random.setZero(Br[0].cols(), Br[0].cols());
  }
  
  if(Bf.size() > 0)
    step_beta_fixed(stepsize);
    counter = 0;
}
void NIGMixedEffect::step_Sigma(double stepsize)
{
  double pos_def = 0;
  iSkroniS     = kroneckerProduct(invSigma, invSigma); 
  UUt         -= counter * vec(Sigma); 
  dSigma_vech  = 0.5 * Dd.transpose() * iSkroniS * UUt;
  ddSigma      = 0.5 * counter * Dd.transpose() * iSkroniS * Dd;
  dSigma_vech  = ddSigma.ldlt().solve(dSigma_vech); 
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
        throw("in NIGmidexeffect not pos def \n");   
    }
    
   
  }
    
    dSigma_vech.setZero(Sigma_vech.size());
    UUt.setZero(Sigma.cols() * Sigma.rows());
    invSigma  = Sigma.inverse();
    Sigma_vech = vech(Sigma);
}

void NIGMixedEffect::step_mu(double stepsize)
{
    
    mu += (0.5 * stepsize) *  H_beta_random.ldlt().solve(gradMu) / VV;
    // H_beta_random = H_mu_random
    mu += (0.5 * stepsize) * (Sigma * gradMu_2)/ (counter * (2*EiV - EV));
    gradMu.setZero(Br[0].cols(), 1);
    gradMu_2.setZero(Br[0].cols(), 1);
}

void NIGMixedEffect::step_nu(double stepsize)
{
   grad_nu  *=  (nu * nu) / (2. * counter); //hessian
  double nu_old = nu;
  nu += stepsize  * grad_nu;
  while(nu < 0){
    nu = nu_old;
    stepsize *= 0.5; 
    nu += stepsize  * grad_nu;
    if(stepsize <= 1e-16){
        Rcpp::Rcout << "nu = \n" << nu << "\n";
        Rcpp::Rcout << "grad_nu = " << grad_nu <<"\n"; 
        throw("in NIGmidexeffect nu is zero \n");   
    }
  }
  
  grad_nu = 0;
  EiV = 1. + 1./nu;
  VV = 1./nu;
}
void NIGMixedEffect::step_beta_fixed(double stepsize)
{
    beta_fixed += stepsize *  H_beta_fixed.ldlt().solve(grad_beta_f);
    grad_beta_f.setZero(Bf[0].cols());
    H_beta_fixed.setZero(Bf[0].cols(), Bf[0].cols());
  
}
void NIGMixedEffect::step_beta_random(double stepsize)
{
    beta_random += (0.5 * stepsize) *  H_beta_random.ldlt().solve(grad_beta_r);
    beta_random += (0.5 * stepsize) * (Sigma * grad_beta_r2)/ (counter * EiV);
    grad_beta_r.setZero(Br[0].cols());
    grad_beta_r2.setZero(Br[0].cols());   
}


void NIGMixedEffect::remove_cov(const int i, Eigen::VectorXd & Y)
{
  if(Br.size() > 0 )
    Y -= Br[i] * beta_random;
  if(Bf.size() > 0)
    Y -= Bf[i] * beta_fixed;
}
void NIGMixedEffect::add_cov(const int i, Eigen::VectorXd & Y)
{
  if(Br.size() > 0 )
    Y += Br[i] * beta_random;
  if(Bf.size() > 0)
    Y += Bf[i] * beta_fixed;
}