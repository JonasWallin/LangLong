
#include <Rcpp.h>
#include <random>
#include <chrono>
#include <string>
#include "Qmatrix.h"
#include "solver.h"
#include "operator_helper.h"
#include "GIG.h"
#include "MixedEffect.h"
#include "measError.h"
using namespace Rcpp;


List simulateLongV_cpp(List obs_list, List operator_list, List theta_list)
{
  /*
    CREATING THE OPERATOR AND SOLVER
  */
  std::string type_operator = Rcpp::as<std::string>(operator_list["type"]); 
  Qmatrix* Kobj;   
  //Kobj = new MaternMatrixOperator;
  operator_select(type_operator, &Kobj);
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  cholesky_solver Solver;
  Solver.init(Kobj->d, 0, 0, 0);
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd>(theta_list["kappa"]);
  Kobj->vec_to_theta(kappa.log()); 
  gig rgig;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  
  
  Eigen::VectorXd h  = Rcpp::as<Eigen::VectorXd>(operator_list["h"]);
  double lambda      = Rcpp::as<double>(theta_list["lambda"]);
  int    GAL         = Rcpp::as<double>(theta_list["GAL"]); 
  //
  List out;
  List out_V(obs_list.length());
  int counter = 0;
  
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
     List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    
    Eigen::VectorXd V = sampleV_pre(rgig, h, lambda ,GAL );

    out_V[counter]   = V;
  }
  return(out_V);
}

/*
	Simulating from the prior model
*/
// [[Rcpp::export]]
List simulateLongGH_cpp(Rcpp::List in_list)
{
 	//**********************************
	//setting up the main data
	//**********************************
	Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
	int nindv = obs_list.length(); //count number of patients  
	std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
	std::vector< Eigen::VectorXd > Ys( nindv);
	int counter = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    	As[counter] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    	Ys[counter++] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
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
	Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd> ( in_list["kappa"]);
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
	
  counter = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>( *it);
    	Solver[counter].init(Kobj->d, 0, 0, 0);
    	Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    	Q = Q  * Kobj->Q;
    	Q = Q + As[counter].transpose()*As[counter];
    	Solver[counter].analyze(Q);
    	Solver[counter++].compute(Q);
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
		
  
  //***********************************
	// measurement error setup
	//*********************************** 
  Rcpp::List MeasureError_list  = Rcpp::as<Rcpp::List> (in_list["MeasureError_list"]);
  MeasurementError *errObj;
  std::string MeasureNoise = Rcpp::as <std::string> (MeasureError_list["noise"]);
  if(MeasureNoise == "Normal")
    errObj = new GaussianMeasurementError;
  else
    errObj = new NIGMeasurementError;
  
  errObj->initFromList(MeasureError_list);
	//***********************************
	// stochastic processes setup
	//*********************************** 
	Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
	Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V_list"]);
	std::string type_processes = Rcpp::as<std::string> (processes_list["noise"]);
	std::vector< Eigen::VectorXd >   Vs( nindv);
  	std::vector< Eigen::VectorXd > Xs( nindv);
  	for(int i = 0; i < nindv; i++ ){ 
    	Xs[i].resize( Kobj->d );
    	if(type_processes == "Normal")
    		Vs[i] = h;
    	else
      	Vs[i] = Rcpp::as<Eigen::VectorXd>( V_list[i]);
  	}
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
  	
  	Eigen::VectorXd b;
  	b.setZero(Kobj->d);
  	
    //*********************************************
    //        simulating the measurement error
    //*********************************************
    std::vector< Eigen::VectorXd > Ysim = errObj->simulate( Ys);
  	mixobj->simulate();
   
    
    
    //*********************************************
    //        simulating the mixed effect
    //*********************************************
    for(int i = 0; i < Ysim.size(); i++)
    {
 		mixobj->add_inter(i, Ysim[i]);
		mixobj->add_cov(i, Ysim[i]);
 	  }
    
    Rcpp::List out_list;
    out_list["Y"]    = Ysim;
    out_list["U"]    = mixobj->U;
    
    
  return(out_list);

}

