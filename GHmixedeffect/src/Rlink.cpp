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
  
  
  //        setting up measurement error            //
  //************************************************//
  //************************************************//
  std::string meas_noise = "Normal";
  if( input.containsElementNamed("meas_noise")){
  	meas_noise = Rcpp::as< std::string  > (input["meas_noise"]);
  }
  
  MeasurementError *errObj;
  if(meas_noise == "Normal"){
    errObj = new GaussianMeasurementError;
    errObj->sigma = input["sigma_eps"];
  }else{
    errObj = new NIGMeasurementError;
    if( input.containsElementNamed("meas_list") == 0)
      throw("in Rlink input list must contain list denoted meas_list! \n");  
    errObj->initFromList( Rcpp::as< Rcpp::List  >(input["meas_list"]));
  }
  //    end of measurement error setup              //
  //************************************************//
  //************************************************//
  
  
  
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
  Eigen::VectorXd nuVec_noise;
  Eigen::VectorXd sigmaVec;
  if(mixobj->Br.size() > 0){
    if(noise == "NIG"){
      muVec.resize(Niter, mixobj->Br[0].cols());
      nuVec.resize(Niter);
    }
    betarVec.resize(Niter, mixobj->Br[0].cols());
  }
  if(meas_noise == "NIG")
  	nuVec_noise.resize(Niter);
  sigmaVec.resize(Niter);
  for(int iter = 0; iter < Niter; iter++){
    
    for(int i =0; i < nobs; i++)
    {
      Eigen::VectorXd  res = Ys[i];
      mixobj->remove_cov(i, res);
      if(meas_noise == "Normal"){
        mixobj->sampleU( i, res, 2 * log(errObj->sigma));
        mixobj->gradient(i, res, 2 * log(errObj->sigma));
      }else{
        mixobj->sampleU2( i, res, errObj->Vs[i].cwiseInverse(), 2 * log(errObj->sigma));
        mixobj->gradient2(i,
                          res,
                          errObj->Vs[i].cwiseInverse(),
                          2 * log(errObj->sigma),
                          errObj->EiV);
      }
      mixobj->remove_inter(i, res);
      
      if(meas_noise == "NIG")
      	errObj->sampleV(i, res);
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
    if(meas_noise == "NIG")
    	nuVec_noise[iter] = ((NIGMeasurementError*) errObj)->nu;  
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
  
  if(meas_noise == "NIG")
  	output["nu_measerror"]     = nuVec_noise;
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

/*
	Function for estimating nu and sigma in measurement
	mainly a debug function but can be used for real
*/
// [[Rcpp::export]]
Rcpp::List estimateNIGnoise(Rcpp::List input)
{
	NIGMeasurementError *errObj;
	errObj = new NIGMeasurementError;
	errObj->initFromList(input);
	
	// read the data
	Rcpp::List Ys_list = Rcpp::as< Rcpp::List > (input["Y"]);
	int nobs = Ys_list.length();
	std::vector< Eigen::VectorXd > Ys(nobs);
	int count = 0;
    for( List::iterator it = Ys_list.begin(); it != Ys_list.end(); ++it ) 
    	Ys[count++] = Rcpp::as < Eigen::VectorXd >( it[0]);
  
  
  
    int Niter  = Rcpp::as< double  > (input["Niter"]); 
    Eigen::VectorXd sigmaVec, nuVec;
    sigmaVec.resize(Niter);
    nuVec.resize(Niter);
    for(int iter = 0; iter < Niter; iter++){
    
    	for(int i =  0; i < nobs; i++)
    	{
      		errObj->sampleV(i, Ys[i]);
      		errObj->gradient(i, Ys[i]);
    	}
    	errObj->step_theta(0.33);
    	sigmaVec[iter] = errObj->sigma;
    	nuVec[iter]    = ((NIGMeasurementError*) errObj)->nu;
    }
  
 	Rcpp::List output = errObj->toList();
 	output["sigmaVec"] = sigmaVec;
 	output["nuVec"]    = nuVec;
	delete errObj;
	return(output);
}


