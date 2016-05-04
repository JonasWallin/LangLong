#include "MixedEffect.h"



Rcpp::List NIGMixedEffect::toList()
{
  Rcpp::List out;
  out["B"]     = B;
  out["Sigma"] = Sigma;
  out["U"]     = U;
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
  
  if( init_list.containsElementNamed("V" )){
    V = Rcpp::as< Eigen::MatrixXd > (init_list["V"]) ;
  }else{
    V.setZero(B.size(), 1);
  }
  


void NIGMixedEffect::sampleV(const int i)
{
  		  //  GIG	
          
}

void NIGMixedEffect::sampleU(const int i, 
                                const Eigen::VectorXd& res,
                                const double log_sigma2_noise)
{
  		  // sample V 	
          Eigen::VectorXd b   = exp( - log_sigma2_noise) * (B[i].transpose() * res);
          Eigen::MatrixXd Q   = exp( - log_sigma2_noise)  * B[i].transpose() * B[i];
          Q +=  invSigma;
          U.col(i) = sample_Nc(b, Q);
          
}

void NIGMixedEffect::gradient(const int i, const Eigen::VectorXd& res)
{
    counter++;
    Eigen::MatrixXd UUT = V(i) * (U.col(i) * U.col(i).transpose());
    // V
    // X-V...
    UUt += vec( UUT);
}