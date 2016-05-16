#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
#include "measError.h"
#include "MixedEffect.h"


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
  GaussianMeasurementError *errObj;
  errObj = new GaussianMeasurementError;
  errObj->sigma = input["sigma_eps"];
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
    mixobj = new NIGMixedEffect;
   
  mixobj->initFromList(input);
  Eigen::VectorXd dbeta;
  Eigen::MatrixXd d2beta;
  
  Eigen::MatrixXd muVec;
  Eigen::MatrixXd betarVec;
  Eigen::VectorXd nuVec;
  Eigen::VectorXd sigmaVec;
  if(mixobj->Br.size() > 0){
    if(noise == "NIG"){
      muVec.resize(Niter, mixobj->Br[0].cols());
      nuVec.resize(Niter);
    }
    betarVec.resize(Niter, mixobj->Br[0].cols());
  }
  sigmaVec.resize(Niter);
  for(int iter = 0; iter < Niter; iter++){
    
    for(int i =0; i < nobs; i++)
    {
      Eigen::VectorXd  res = Ys[i];
      mixobj->remove_cov(i, res);
      mixobj->sampleU(i , res, 2 * log(errObj->sigma));
      mixobj->gradient(i, res, 2 * log(errObj->sigma));
      mixobj->remove_inter(i, res);
      errObj->gradient(i, res);
    }
    
    if(mixobj->Br.size() > 0){
      if(noise == "NIG"){
        muVec.row(iter) = ((NIGMixedEffect*) mixobj)->mu;
        nuVec(iter) = ((NIGMixedEffect*) mixobj)->nu;
      }
      betarVec.row(iter) = mixobj->beta_random;
    }
    mixobj->step_theta(0.33);
    errObj->step_theta(0.33);
    sigmaVec[iter] = errObj->sigma;
  }
  //mixobj.sampleU(0, Ys[0], 0);
  Rcpp::List output;
  output["beta"]        = mixobj->beta_random;
  output["sigma_eps"]   = errObj->sigma;
  output["mixedeffect"] = mixobj->toList();
  output["measerror"]   = errObj->toList();
  if(mixobj->Br.size() > 0){
    if(noise == "NIG"){
      output["mu"]      = muVec;
      output["nuVec"]   = nuVec;
      }
    output["betaVec"]   = betarVec;
    output["sigma_eps"]     = sigmaVec;
  }
  delete mixobj;
  delete errObj;
  return(output);
}


// [[Rcpp::export]]
Eigen::VectorXd samplePostV(Rcpp::List input)
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
 
  NIGMixedEffect mixobj;
   
  mixobj.initFromList(input);
  Eigen::VectorXd Vs;
  Vs.resize(Niter);
  for(int iter = 0; iter < Niter; iter++)
  {
    for(int i =0; i < nobs; i++)
    {
      Eigen::VectorXd  res = Ys[i];
      mixobj.remove_cov(i, res);
      mixobj.sampleU(i, res, log_sigma2_eps);
      
    }
    Vs[iter] = mixobj.V[0];
  } 
  return(Vs);
}