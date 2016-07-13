#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <math.h>
#include "Qmatrix.h"
#include "operator_helper.h"
#include "solver.h"
#include "GIG.h"
#include "MixedEffect.h"
#include "measError.h"
#include "process.h"
using namespace Rcpp;




double estDigamma(double x)
{
  return(R::digamma(x));
}





// [[Rcpp::export]]
List estimateLong_cpp(Rcpp::List in_list)
{

	
	//**********************************
	//      basic parameter
	//**********************************
	
	double pSubsample = Rcpp::as< double > (in_list["pSubsample"]);
	int nIter      = Rcpp::as< double > (in_list["nIter"]);
	int nSim       = Rcpp::as< double > (in_list["nSim"]);
  int nBurnin    = Rcpp::as< double > (in_list["nBurnin"] );
  int silent     = Rcpp::as< int    > (in_list["silent"]);
  double alpha     = Rcpp::as< double    > (in_list["alpha"]);
  double step0     = Rcpp::as< double    > (in_list["step0"]);
	//**********************************
	//     setting up the main data
	//**********************************
	Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
	int nindv = obs_list.length(); //count number of patients  
  int nSubsample = ceil(pSubsample * nindv);
	std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
	std::vector< Eigen::VectorXd > Ys( nindv);
	int count;
	count = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    	As[count] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    	Ys[count] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
      count++;
    }

	//**********************************
	//operator setup
	//***********************************
	Rcpp::List operator_list  = Rcpp::as<Rcpp::List> (in_list["operator_list"]);
	std::string type_operator = Rcpp::as<std::string>(operator_list["type"]); 
	Qmatrix* Kobj;   
	//Kobj = new MaternMatrixOperator;
	operator_select(type_operator, &Kobj);
	Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
	Eigen::VectorXd kappa;
  if(operator_list.containsElementNamed("kappa"))
    kappa = Rcpp::as<Eigen::VectorXd> ( operator_list["kappa"]);
  Eigen::VectorXd dkappa;
  dkappa.setZero(kappa.size());
  
  Eigen::VectorXd  tauVec;
  Eigen::MatrixXd kappaVec;
  if(kappa.size() > 0)
    kappaVec.resize(nIter, kappa.size());
  tauVec.resize(nIter,1);
	Kobj->vec_to_theta( kappa);
	Eigen::VectorXd h = Rcpp::as<Eigen::VectorXd>( operator_list["h"]);

	//Prior solver
	cholesky_solver Qsolver;
	Qsolver.init( Kobj->d, 0, 0, 0);
	Qsolver.analyze( Kobj->Q);
	Qsolver.compute( Kobj->Q);
	
	//Create solvers for each patient
	std::vector<  cholesky_solver >  Solver( nindv);
	Eigen::SparseMatrix<double, 0, int> Q;
	
	count = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>( *it);
    	Solver[count].init(Kobj->d, 0, 0, 0);
    	Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    	Q = Q  * Kobj->Q;
    	Q = Q + As[count].transpose()*As[count];
    	Solver[count].analyze(Q);
    	Solver[count].compute(Q);
    	count++;
  }	
	//**********************************
	// mixed effect setup
	//*********************************** 
	Rcpp::List mixedEffect_list  = Rcpp::as<Rcpp::List> (in_list["mixedEffect_list"]);
	std::string type_mixedEffect = Rcpp::as<std::string> (mixedEffect_list["noise"]);
	MixedEffect *mixobj;
	if(type_mixedEffect == "Normal")
    	mixobj = new NormalMixedEffect;
	else
		mixobj   = new NIGMixedEffect;
	
	mixobj->initFromList(mixedEffect_list);
	
  	Eigen::MatrixXd muVec_mixed;
  	Eigen::MatrixXd betarVec_mixed;
  	Eigen::MatrixXd betafVec_mixed;
  	Eigen::VectorXd nuVec_mixed;
  	Eigen::MatrixXd SigmaVec_mixed;
  	if(mixobj->Br.size() > 0){
    	if(type_mixedEffect == "NIG"){
      		muVec_mixed.resize(nIter, mixobj->Br[0].cols());
      		nuVec_mixed.resize(nIter);
    	}
      	betarVec_mixed.resize(nIter, mixobj->Br[0].cols());
      	SigmaVec_mixed.resize(nIter, mixobj->Br[0].cols() * mixobj->Br[0].cols());
    }
  	if(mixobj->Bf.size() > 0)
    	betafVec_mixed.resize(nIter, mixobj->Bf[0].cols());
  
  
  
  //**********************************
	// measurement setup
	//*********************************** 
  MeasurementError *errObj;
  Rcpp::List measurementError_list  = Rcpp::as<Rcpp::List> (in_list["measurementError_list"]);
  std::string type_MeasurementError= Rcpp::as <std::string> (measurementError_list["noise"]);
  if(type_MeasurementError == "Normal")
    errObj = new GaussianMeasurementError;
  else
    errObj = new NIGMeasurementError; 
    
  errObj->initFromList(measurementError_list);
  
  
  
  
  
	//**********************************
	// stochastic processes setup
	//*********************************** 
	Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
	Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V"]);
  
	std::string type_processes  = Rcpp::as<std::string> (processes_list["noise"]);
  
  
  Process *process;
  
  if (type_processes != "Normal"){
  	process  = new GHProcess;
  }else{ process  = new GaussianProcess;}
  
  process->initFromList(processes_list, h);
  	/*
  	Simulation objects
  	*/
  	std::mt19937 random_engine;
  	std::normal_distribution<double> normal;
  	std::default_random_engine gammagenerator;
  	random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  	gig rgig;
  	rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  	Eigen::VectorXd  z;
  	z.setZero(Kobj->d);
  	
  	Eigen::VectorXd b, Ysim;
  	b.setZero(Kobj->d);
  	
  	std::vector<int> longInd;
  	for (int i=0; i< nindv; i++) longInd.push_back(i);
  
  
    
  
  
  	for(int iter=0; iter < nIter + nBurnin; iter++){
    	/*
        printing output
      
      */
    	if(silent == 0){
      		Rcpp::Rcout << "i = " << iter << ": \n";
          if(mixobj->Br.size() > 0)
            Rcpp::Rcout << "mixObj->beta_r = " << mixobj->beta_random.transpose() << "\n";
          
          if(mixobj->Bf.size() > 0)
            Rcpp::Rcout << "mixObj->beta_f = " << mixobj->beta_fixed.transpose() << "\n";
            
          if(mixobj->Br.size() > 0)
            Rcpp::Rcout << "mixObj->Sigma = \n" << mixobj->Sigma << "\n"; 
            
            
            Rcpp::Rcout << "errObj->sigma = " << errObj->sigma << "\n"; 
            
            
          if(type_MeasurementError != "Normal")
            Rcpp::Rcout << "errObj->nu = " << ((NIGMeasurementError*) errObj)->nu << "\n";  
            
          Rcpp::Rcout << "tau = " << Kobj->tau << "\n"; 
          if(kappa.size() > 0 )
            Rcpp::Rcout << "kappa = " << kappa[0] << "\n"; 
    	}
      func_dkappa0(type_operator, dkappa, nSim, nSubsample, *Kobj);
      		
      		
      	Eigen::SparseMatrix<double,0,int> K = Eigen::SparseMatrix<double,0,int>(Kobj->Q);
      	K *= sqrt(Kobj->tau);
        
      	// subsampling
    	//sample Nlong values without replacement from 1:nrep
    	std::random_shuffle(longInd.begin(), longInd.end());
    
    	for(int ilong = 0; ilong < nSubsample; ilong++ ) {
      		int i = longInd[ilong];
          
      		Eigen::SparseMatrix<double,0,int> A = As[i];
      		Eigen::VectorXd  Y = Ys[i];	
      		
      		
      		for(int ii = 0; ii < nSim; ii ++){
      		
      			Eigen::VectorXd  res = Y;
      			//***************************************
      			//***************************************
      			//   building the residuals and sampling
      			//***************************************
      			//***************************************
      			
      			// removing fixed effect from Y
      			mixobj->remove_cov(i, res);
      			
      			res -= A * process->Xs[i]; 
      			//***********************************
      			// mixobj sampling
      			//***********************************
      			
            
      			if(type_MeasurementError == "Normal")
      				mixobj->sampleU( i, res, 2 * log(errObj->sigma));
      			else
      				mixobj->sampleU2( i,
                	         		  res, 
                    	     		  errObj->Vs[i].cwiseInverse(), 
                    	     		  2 * log(errObj->sigma));
                                 
            mixobj->remove_inter(i, res);
            
      

      			//***********************************
      			// sampling processes
      			//***********************************
      			
      			Eigen::SparseMatrix<double,0,int> Q = Eigen::SparseMatrix<double,0,int>(K.transpose());
      			Eigen::VectorXd iV(process->Vs[i].size());
      			iV.array() = process->Vs[i].array().inverse();
      			Q =  Q * iV.asDiagonal();
      			Q =  Q * K;   	
      			
        		
        		for(int j =0; j < Kobj->d; j++)
          			z[j] =  normal(random_engine);
      
      
      			res += A * process->Xs[i]; 
      			//Sample X|Y, V, sigma
            
            if(type_MeasurementError == "Normal")
      			  process->sample_X(i, 
                		z,
                		res,
                		Q,
                    	K,
                		A,
                		errObj->sigma,
                		Solver[i]);    
            else
              process->sample_Xv2( i, 
                  		z,
                		res,
                		Q,
                    	K,
                		A,
                		errObj->sigma,
                		Solver[i],
                        errObj->Vs[i].cwiseInverse());   
            res -= A * process->Xs[i];
            
            if(res.cwiseAbs().sum() > 1e16){
              throw("res outof bound\n");
            }
            
            // sample V| X
          	process->sample_V(i, rgig, K);
                                       
        		
           //***********************************
           // random variance noise sampling
      		 //***********************************
     			  if(type_MeasurementError != "Normal"){
      			  errObj->sampleV(i, res);                           
     			  }
      			//***************************************
      			//  computing gradients
      			//***************************************
      		  if(iter >= nBurnin){
             
      			  //***************************************
      			  // mixobj gradient
      			  //***************************************
              
      			  mixobj->add_inter(i, res);
      			  if(type_MeasurementError != "Normal")
        			  mixobj->gradient2(i,
        				  			          res,
                          	  		errObj->Vs[i].cwiseInverse(),
                              	  2 * log(errObj->sigma),
                                  errObj->EiV);
              else
                mixobj->gradient(i, res, 2 * log(errObj->sigma));
      		    
      			  mixobj->remove_inter(i, res);
      		
      		
      			  //***************************************
      			  // measurent error  gradient
      			  //***************************************
      			  errObj->gradient(i, res);
              
              
        		  //***************************************
      			  // operator gradient
      			  //***************************************
      			  Kobj->gradient( process->Xs[i],
      			  			      process->Vs[i].cwiseInverse());
      			  
      			   //***************************************
      			  // operator gradient
      			  //***************************************
      			  process->gradient(i);
              
      		  }  
      		}
      		
      		
		  	
      
    	}
      //**********************************
  	  	//  gradient step
		  	//*********************************** 
        if(iter >= nBurnin){
        	double stepsize = step0 / pow(iter - nBurnin + 1, alpha);
          mixobj->step_theta(stepsize);
          errObj->step_theta(stepsize);
          Kobj->step_theta(stepsize);
          process->step_theta(stepsize);
        }
        //**********************************
  	  	// storing the parameter traces
		//*********************************** 
      
        if(iter >= nBurnin)
        {
         
           if(mixobj->Br.size() > 0){
              if(type_mixedEffect != "Normal"){
                muVec_mixed.row(iter - nBurnin) = ((NIGMixedEffect*) mixobj)->mu;
                nuVec_mixed(iter - nBurnin)     = ((NIGMixedEffect*) mixobj)->nu;
      		    }
              betarVec_mixed.row(iter - nBurnin) = mixobj->beta_random;
              
          SigmaVec_mixed.row(iter - nBurnin) = Eigen::Map<Eigen::VectorXd>(mixobj->Sigma.data(), 
        														mixobj->Br.size()*mixobj->Br.size());
                                    
              
          }
    	    if(mixobj->Bf.size() > 0)
    		    betafVec_mixed.row(iter - nBurnin) = mixobj->beta_fixed;
        
          if(kappa.size() > 0)
            kappaVec.row(iter - nBurnin) = kappa;
        
          tauVec(iter - nBurnin) =  Kobj->tau;
        }
    }
    
  // storing the results  
  Rcpp::List out_list;
  out_list["pSubsample"]       = pSubsample;
  out_list["nIter"]            = nIter;
  out_list["nSim"]             = nSim;
  out_list["nBurnin"]          = nBurnin;
  out_list["silent"]           = silent;
  out_list["step0"]            = step0;
  out_list["alpha"]            = alpha;
  out_list["obs_list"]         = obs_list;
  out_list["Xs"]               = process->Xs;
  out_list["Vs"]               = process->Vs;
  out_list["tauVec"]           = tauVec;
  if(kappa.size() > 0)
    out_list["kappaVec"] = kappaVec;
  
  
  Rcpp::List mixobj_list       = mixobj->toList();
  if(mixobj->Br.size() > 0){
    if(type_mixedEffect != "Normal"){
      mixobj_list["muVec"]    = muVec_mixed;
      mixobj_list["nuVec"]    = nuVec_mixed;
    }
      mixobj_list["betarVec"] = betarVec_mixed;
      mixobj_list["SigmaVec"] = SigmaVec_mixed;
  }
  
  if(mixobj->Bf.size() > 0)
    	mixobj_list["betafVec"] = betafVec_mixed;	  
  
  
  out_list["mixedEffect_list"] = mixobj_list;
  
  Rcpp::List errobj_list           = errObj->toList();
  out_list["measurementError_list"] = errobj_list;
  return(out_list);
}