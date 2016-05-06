#include "MixedEffect.h"
#include "error_check.h"
#include <chrono>


Rcpp::List NIGMixedEffect::toList()
{
  Rcpp::List out;
  out["B"]      = B;
  out["Sigma"]  = Sigma;
  out["U"]      = U;
  out["V"]      = V;
  out["nu"]     = nu;
  out["mu"]     = mu;
  return(out);
}
void NIGMixedEffect::initFromList(Rcpp::List const &init_list)
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
  //
  // variance components
  if( init_list.containsElementNamed("V" ))
    V = Rcpp::as< Eigen::MatrixXd > (init_list["V"]) ;
  else
    V.setOnes(B.size(), 1);
    
  
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  //
  //  
  
  if( init_list.containsElementNamed("mu" ))
    mu = Rcpp::as< Eigen::MatrixXd > (init_list["mu"]) ;
  else
    mu.setZero(B[0].cols(), 1); 
  
  
  gradMu.setZero(B[0].cols(), 1);
  
  if( init_list.containsElementNamed("nu" ))
    nu = Rcpp::as< double > (init_list["nu"]) ;
  else
    nu = 1.;
    
  
  a_GIG = mu.transpose() * (invSigma *  mu);
  a_GIG += nu;
  
}


void NIGMixedEffect::sampleV(const int i)
{
  		  //  GIG	(p, a, b)
        double p  = 0.5 * (- 1 - B[0].cols());
        double b  =  U.col(i).transpose() *  U.col(i);
        b += nu;
        V(i) = rgig.sample(p, a_GIG, b);
}

void NIGMixedEffect::sampleU(const int i, 
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise)
{
  		 
          Eigen::VectorXd b   = exp( - log_sigma2_noise) * (B[i].transpose() * res);
          Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * B[i].transpose() * B[i];
          Eigen::MatrixXd Qp  = invSigma / V(i);
          b += Qp * (- mu + V(i) * mu);
          Q += Qp;
          U.col(i) = sample_Nc(b, Q);
          sampleV(i);
          
}

void NIGMixedEffect::gradient(const int i, 
                              const Eigen::VectorXd& res,
                              const double log_sigma2_noise)
{
    counter++;
    Eigen::VectorXd B_mu;
    B_mu.setOnes(U.col(i).size());
    B_mu         *= -1; 
    B_mu.array() += V(i); 
    Eigen::VectorXd U_ = U.col(i) - B_mu * mu; 
    Eigen::MatrixXd UUT =  (U_ * U_.transpose());
    UUT /= V(i);
    
    // V
    // X-V...
    UUt    += vec( UUT);
    gradMu += (-1 + V(i) ) * (invSigma * U_); 
    
    // dtau
}

void NIGMixedEffect::step_theta(double stepsize)
{
   
    a_GIG = mu.transpose() * (invSigma * mu);
    a_GIG += nu;
    gradMu.setZero(B[0].cols(), 1);
}

void NIGMixedEffect::step_Sigma(double stepsize)
{
  // I think I have done something smart here but I cant recall what?
  // is it the expected Fisher information?
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
    counter = 0;
}

void NIGMixedEffect::step_mu(const double stepsize)
{
  //expected fisher information should be
  // count*V[V] * \Sigma^{-1}
}