#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <math.h>
#include "operatorMatrix.h"
#include "operator_helper.h"
#include "solver.h"
#include "GIG.h"
#include "MixedEffect.h"
#include "measError.h"
#include "process.h"
using namespace Rcpp;


// [[Rcpp::export]]
List estimateProcess_cpp(Rcpp::List in_list)
{


	//**********************************
	//      basic parameter
	//**********************************

	int nIter      = Rcpp::as< double > (in_list["nIter"]);
  int nBurnin    = Rcpp::as< double > (in_list["nBurnin"] );
  int silent     = Rcpp::as< int    > (in_list["silent"]);
  double alpha     = Rcpp::as< double    > (in_list["alpha"]);
  double step0     = Rcpp::as< double    > (in_list["step0"]);
  double sigma     = Rcpp::as< double    > (in_list["sigma"]);
  	//**********************************
	//     setting up the main data
	//**********************************
	Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
	int nindv = obs_list.length(); //count number of patients
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
	operator_list["nIter"] = nIter;
	std::string type_operator = Rcpp::as<std::string>(operator_list["type"]);
	operatorMatrix* Kobj;
	//Kobj = new MaternMatrixOperator;
	operator_select(type_operator, &Kobj);
	Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));


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


  	process->setupStoreTracj(nIter);
  	for(int iter=0; iter < nIter + nBurnin; iter++){

      		Rcpp::Rcout << "iter = " << iter << " ";
      		process->printIter() ;
      		Rcpp::Rcout << "\n";
      	Eigen::SparseMatrix<double,0,int> K = Eigen::SparseMatrix<double,0,int>(Kobj->Q);
//      	K *= sqrt(Kobj->tau);

    	for(int i = 0; i < nindv; i++ ) {

      		Eigen::SparseMatrix<double,0,int> A = As[i];
      		Eigen::VectorXd  Y = Ys[i];
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



             process->sample_X(i,
                		z,
                		Ys[i],
                		Q,
                    	K,
                		A,
                		sigma,
                		Solver[i]);


            // sample V| X
          	process->sample_V(i, rgig, K);



      		//***************************************
      		//  computing gradients
      		//***************************************
      		if(iter >= nBurnin){
      			   //***************************************
      			  // operator gradient
      			  //***************************************
      			  process->gradient(i,
                               K,
                               A,
                               Ys[i],
                               sigma,
                               Kobj->trace_variance(A));

      		  }
      		}
      //**********************************
  	  	//  gradient step
		  	//***********************************
        if(iter >= nBurnin){
        	double stepsize = step0 / pow(iter - nBurnin + 1, alpha);
          process->step_theta(stepsize);

        }

    }
  Rcpp::List process_list           = process->toList();
  return(process_list);
}