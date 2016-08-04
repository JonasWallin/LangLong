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



/*
  Program for simulating from prior distribution of the of noise
*/
// [[Rcpp::export]]
Rcpp::List simulateNoise(Rcpp::List input)
{

  Rcpp::List Y_list = input["Y"];
  std::vector< Eigen::VectorXd > Y_in;
  Y_in.resize(Y_list.length());
  int i= 0 ;
  for( Rcpp::List::iterator it = Y_list.begin(); it != Y_list.end(); ++it ) {
    Y_in[i++] = Rcpp::as < Eigen::VectorXd >( it[0]);
  }


  MeasurementError *errObj;
  std::string noise = Rcpp::as <std::string> (input["noise"]);
  if(noise == "Normal")
    errObj = new GaussianMeasurementError;
  else
    errObj = new NIGMeasurementError;

  errObj->initFromList(input);
  std::vector< Eigen::VectorXd > Y = errObj->simulate( Y_in);

  Rcpp::List output;
  output["Y"]    = Y;
  delete errObj;
  return(output);
}
/*
  Program for simulating from prior distribution of the mixedeffect
*/
// [[Rcpp::export]]
Rcpp::List simulateMixed(Rcpp::List input)
{
  MixedEffect *mixobj;
  std::string noise = Rcpp::as <std::string> (input["noise"]);
  if(noise == "Normal")
    mixobj = new NormalMixedEffect;
  else
    mixobj = new NIGMixedEffect;

  mixobj->initFromList(input);
  mixobj->simulate();

  Rcpp::List output;
  output["U"]    = mixobj->U;
  output["beta_random"] = mixobj->beta_random;
  delete mixobj;
  return(output);
}

// [[Rcpp::export]]
Rcpp::List estimateME(Rcpp::List input)
{


  int Niter = Rcpp::as< int  > (input["Niter"]);
  int nSim  = Rcpp::as< int  > (input["nSim"]); 
  double alpha     = Rcpp::as< double    > (input["alpha"]);
  double step0     = Rcpp::as< double    > (input["step0"]);


  Rcpp::List Ys_list = Rcpp::as< Rcpp::List > (input["Y"]);
  int nobs = Ys_list.length();
  std::vector< Eigen::VectorXd > Ys(nobs);

  //        setting up measurement error            //
  //************************************************//
  //************************************************//
  std::string meas_noise = "Normal";
  
  Rcpp::List  meas_list = Rcpp::as< Rcpp::List  >(input["measurementError_list"]);
  if( meas_list.containsElementNamed("noise")){
  	meas_noise = Rcpp::as< std::string  > (meas_list["noise"]);
  }

  MeasurementError *errObj;
  if(meas_noise == "Normal"){
    errObj = new GaussianMeasurementError;
  }else{
    errObj = new NIGMeasurementError;
    
  }
  errObj->initFromList(meas_list );
  errObj->setupStoreTracj(Niter);
  //    end of measurement error setup              //
  //************************************************//
  //************************************************//



  int count = 0;

  double n = 0;
  // init data effect
  for( List::iterator it = Ys_list.begin(); it != Ys_list.end(); ++it ) {
    Ys[count] = Rcpp::as < Eigen::VectorXd >( it[0]);
    n += Ys[count].size();
    count++;
  }
  Rcpp::List  mixedEffect_list =  Rcpp::as< Rcpp::List  >(input["mixedEffect_list"]);

  std::string noise;
  if( mixedEffect_list.containsElementNamed("noise"))
      noise = Rcpp::as <std::string> (mixedEffect_list["noise"]);
  else
     noise = "Normal";
  // init mixed effect
  MixedEffect *mixobj;

  if(noise == "Normal")
    mixobj = new NormalMixedEffect;
  else
    mixobj = new NIGMixedEffect;
	
  mixobj->initFromList(mixedEffect_list);
  mixobj->setupStoreTracj(Niter);
  
  
  for(int iter = 0; iter < Niter; iter++){
	for(int ii = 0; ii < nSim; ii++){
    	for(int i =0; i < nobs; i++){
      	Eigen::VectorXd  res = Ys[i];
      	mixobj->remove_cov(i, res);
      	if(meas_noise == "Normal"){
       	 mixobj->sampleU( i, res, 2 * log(errObj->sigma));
       	 mixobj->gradient(i, res, 2 * log(errObj->sigma));
      	}else{
        	mixobj->sampleU2( i,
                         res,
                         errObj->Vs[i].cwiseInverse(),
                         2 * log(errObj->sigma));
        	mixobj->gradient2(i,
                          res,
                          errObj->Vs[i].cwiseInverse(),
                          2 * log(errObj->sigma),
                          errObj->EiV);
      	}
      	mixobj->remove_inter(i, res);

      	if(meas_noise != "Normal")
      		errObj->sampleV(i, res);
      	errObj->gradient(i, res);
    	}

    }
    	double stepsize = step0 / pow(iter  + 1, alpha);
    	mixobj->step_theta( stepsize);
    	errObj->step_theta( stepsize);
  
  }
  //mixobj.sampleU(0, Ys[0], 0);
  Rcpp::List output;
  output["Y_list"]      = Ys_list;
  output["mixedEffect_list"] = mixobj->toList();
  output["measurementError_list"]   = errObj->toList();
  
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


/*
	Function for estimating the Fisher information at \theta
*/
// [[Rcpp::export]]
Rcpp::List EstimateFisherInformation(Rcpp::List input)
{
  int count =0;
  Rcpp::List Ys_list = Rcpp::as< Rcpp::List > (input["Y"]);
  int nobs = Ys_list.length();
  std::vector< Eigen::VectorXd > Ys(nobs);


  //        setting up measurement error            //
  //************************************************//
  //************************************************//
  Rcpp::List meas_list = Rcpp::as< Rcpp::List  >(input["measurementError_list"]);
  std::string meas_noise = "Normal";
  if( meas_list.containsElementNamed("noise")){
  	meas_noise = Rcpp::as< std::string  > (meas_list["noise"]);
  }

  MeasurementError *errObj;
  if(meas_noise == "Normal"){
    errObj = new GaussianMeasurementError;
  }else{
    errObj = new NIGMeasurementError;
  }
  errObj->initFromList( meas_list);
  
  //    end of measurement error setup              //
  //************************************************//
  //************************************************//



  int Niter   = Rcpp::as< int  > (input["Niter"]);
  int nSim    = Rcpp::as< int  > (input["nSim"]);
  int nBurnin = Rcpp::as< int  > (input["nBurnin"]);
  double n = 0;
  // init data effect
  for( List::iterator it = Ys_list.begin(); it != Ys_list.end(); ++it ) {
    Ys[count] = Rcpp::as < Eigen::VectorXd >( it[0]);
    n += Ys[count].size();
    count++;
  }

  std::string noise;
  Rcpp::List mixedEffect_list = Rcpp::as< Rcpp::List  >(input["mixedEffect_list"]);
  if( mixedEffect_list.containsElementNamed("noise"))
      noise = Rcpp::as <std::string> (mixedEffect_list["noise"]);
  else
     noise = "Normal";
  // init mixed effect
  MixedEffect *mixobj;

  if(noise == "Normal")
    mixobj = new NormalMixedEffect;
  else
    mixobj = new NIGMixedEffect;


  mixobj->initFromList( mixedEffect_list);


  Eigen::MatrixXd Egrad2;
  Egrad2.setZero(mixobj->npars + errObj->npars, mixobj->npars + errObj->npars);
  Eigen::VectorXd Egrad;
  Egrad.setZero(mixobj->npars + errObj->npars);

  for(int iter = 0; iter < Niter; iter++){

	std::vector< Eigen::VectorXd > Y = errObj->simulate( Ys);
	mixobj->simulate(Y);
	for(int ii = 0; ii < nSim + nBurnin; ii++)
	{
    	for(int i =0; i < nobs; i++)
    	{
    
    		Eigen::VectorXd  res = Y[i];
      		mixobj->remove_cov(i, res);
      		if(meas_noise == "Normal"){
        		mixobj->sampleU( i, res, 2 * log(errObj->sigma));
        		if(ii >= nBurnin)
        			mixobj->gradient(i, res, 2 * log(errObj->sigma));
      		}else{
        		mixobj->sampleU2( i,
                		         res,
                        		 errObj->Vs[i].cwiseInverse(),
                         		2 * log(errObj->sigma));
                
        		if(ii >= nBurnin)
        			mixobj->gradient2(i,
                          		res,
                          		errObj->Vs[i].cwiseInverse(),
                          		2 * log(errObj->sigma),
                          		errObj->EiV);
      		}
      		mixobj->remove_inter(i, res);
      		
        	if(ii >= nBurnin){
      			errObj->gradient(i, res);
      			
          		
          	}
      	}
    
    }
    Eigen::VectorXd grad(mixobj->npars + errObj->npars);
    grad.segment(0, mixobj->npars) = mixobj->get_gradient();
    grad.segment(mixobj->npars, errObj->npars) = errObj->get_gradient();
    grad.array() /= nSim;
    Egrad2 += grad * grad.transpose();
    Egrad += grad;
    Rcpp::Rcout << "grad = " << grad.transpose() << "\n";
    
	mixobj->clear_gradient();
	errObj->clear_gradient();
  
  }
  Egrad2.array() /= Niter;
  Egrad.array() /= Niter; 
  Eigen::MatrixXd Vgrad = Egrad2 - Egrad * Egrad.transpose();
  Eigen::MatrixXd invVgrad  = Vgrad.inverse();
  mixobj->set_covariance(invVgrad.block(0, 0, mixobj->npars, mixobj->npars));
  errObj->set_covariance(invVgrad.block(mixobj->npars, mixobj->npars, errObj->npars, errObj->npars));
  Rcpp::Rcout << " std = " << invVgrad.diagonal().cwiseSqrt() << "\n";
  Rcpp::List output;
  output["Y_list"]      = Ys_list;
  
  output["mixedEffect_list"] = mixobj->toList();
  output["measurementError_list"]   = errObj->toList();
  
  
  delete mixobj;
  delete errObj;
  return(output);
}

