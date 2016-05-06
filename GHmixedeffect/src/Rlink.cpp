#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
#include "MixedEffect.h"

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
Eigen::VectorXd sample_Nccpp(Eigen::VectorXd b, Eigen::MatrixXd Q) {
  
  return sample_Nc(b, Q);
}


// [[Rcpp::export]]
Rcpp::List estimateME(Rcpp::List input)
{
  
  Rcpp::List Ys_list = Rcpp::as< Rcpp::List > (input["Y"]);
  int nobs = Ys_list.length();
  
  std::vector< Eigen::VectorXd > Ys(nobs);
  double log_sigma2_eps = 2 * log(Rcpp::as< double  > (input["sigma_eps"]));
  int Niter  = Rcpp::as< double  > (input["Niter"]);
  int count =0;
  
  double n = 0;
  // init data effect
  for( List::iterator it = Ys_list.begin(); it != Ys_list.end(); ++it ) {
    Ys[count] = Rcpp::as < Eigen::VectorXd >( it[0]);
    n += Ys[count].size();
    count++;
  }
  std::string noise;
  if( input.containsElementNamed("noise"))
      noise = Rcpp::as <std::string> (input["noise"]);
  else
     noise = "Normal";
  // init mixed effect
  MixedEffect *mixobj;
  
  if(noise == "Normal")
    mixobj = new NormalMixedEffect;
  else
    mixobj = new NormalMixedEffect;
   
  mixobj->initFromList(input);
  Eigen::VectorXd dbeta;
  Eigen::MatrixXd d2beta;
  
  
  
  for(int iter = 0; iter < Niter; iter++){
    double d2sigma  = 0;
    double dsigma   = 0;
    
    for(int i =0; i < nobs; i++)
    {
      Eigen::VectorXd  res = Ys[i];
      mixobj->remove_cov(i, res);
      mixobj->sampleU(i, res, log_sigma2_eps);
      mixobj->gradient(i, res, log_sigma2_eps);
      mixobj->remove_inter(i, res);
      d2sigma -=  0.5 * exp(-log_sigma2_eps)*res.array().square().sum();
      dsigma  += -0.5 * res.size();
      
    }
    
    mixobj->step_theta(0.33);
    dsigma -= d2sigma;
    log_sigma2_eps += dsigma / n; 
    
  }
  //mixobj.sampleU(0, Ys[0], 0);
  Rcpp::List output;
  output["beta"]        = mixobj->beta_random;
  output["sigma_eps"]   = exp(0.5 * log_sigma2_eps);
  output["mixedeffect"] = mixobj->toList();
  return(output);
}