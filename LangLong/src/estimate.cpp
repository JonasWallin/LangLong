#include <Rcpp.h>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include "Qmatrix.h"
#include "operator_helper.h"
#include "solver.h"
#include "GIG.h"
#include "MixedEffect.h"
using namespace Rcpp;




double estDigamma(double x)
{
  return(R::digamma(x));
}

void sampleX( Eigen::VectorXd & X, 
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & V,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma2,
              cholesky_solver       & solver)
{
  
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()*Y/ sigma2;
  X = solver.rMVN(b, Z);
}


// [[Rcpp::export]]
List estimateLong_cpp(Rcpp::List in_list)
{
	//**********************************
	//setting up the main data
	//**********************************
	Rcpp::List obs_list  = Rcpp::as<Rcpp::List> (in_list["obs_list"]);
	int nindv = obs_list.length(); //count number of patients  
	std::vector< Eigen::SparseMatrix<double,0,int> > As( nindv);
	std::vector< Eigen::VectorXd > Ys( nindv);
	int i;
	i = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    	As[i] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    	Ys[i] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
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
	
	i = 0;
	for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    	List obs_tmp = Rcpp::as<Rcpp::List>( *it);
    	Solver[i].init(Kobj->d, 0, 0, 0);
    	Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    	Q = Q  * Kobj->Q;
    	Q = Q + As[i].transpose()*As[i];
    	Solver[i].analyze(Q);
    	Solver[i].compute(Q);
    	i++;
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
		
  
	//**********************************
	// stochastic processes setup
	//*********************************** 
	Rcpp::List processes_list   = Rcpp::as<Rcpp::List>  (in_list["processes_list"]);
	Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (processes_list["V_list"]);
	std::string type_processes = Rcpp::as<std::string> (processes_list["noise"]);
	std::vector< Eigen::VectorXd >   Vs( nindv);
  	std::vector< Eigen::VectorXd > Xs( nindv);
  	for( i = 0; i < nindv; i++ ){ 
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
  	
  	Eigen::VectorXd b, Ysim;
  	b.setZero(Kobj->d);
  	
  	std::vector<int> longInd;
  	for (i=0; i< nindv; i++) longInd.push_back(i);
  
  
 
  return(in_list);
}

// [[Rcpp::export]]
List samplePosteriorGH(List obs_list, 
                        List operator_list, 
                        List theta_list,
                        List mixed_list,
                        List V_list,
                        int nsim, 
                        int noise,
                        int commonsigma)
{
  int GAL = 0;
  if(noise == 1)
    GAL = 1; 

  //count number of patients  
  int nrep = obs_list.length();
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd> (theta_list["kappa"]);
  double tau = log(Rcpp::as< double> (theta_list["tau"]));
  Eigen::VectorXd  beta = Rcpp::as<Eigen::VectorXd>( theta_list["beta"]);
  double sigma = 0;
  Eigen::VectorXd sigma_v;
  if(commonsigma == 1){
    sigma  = 2 * log(Rcpp::as<double> (theta_list["sigma"]));  
  } else {
    sigma_v = Rcpp::as<Eigen::VectorXd>( theta_list["sigma"]);
    sigma_v.array() = 2 * sigma_v.array().log();
  }
  
  double asigma, bsigma;
  if(commonsigma == 0 ){
    asigma = theta_list["asigma"];  
    bsigma  = theta_list["bsigma"];
  }
  
  double lambda = log(Rcpp::as< double> (theta_list["lambda"]));
  double mu = theta_list["mu"];
  int    usingMixedEff = Rcpp::as<double>(mixed_list["on"]);
  NormalMixedEffect mixobj;
  if(usingMixedEff)
    mixobj.initFromList(mixed_list);
  
  //Define operator
  std::string type_operator = Rcpp::as<std::string>(operator_list["type"]); 
  Qmatrix* Kobj;   
  operator_select(type_operator, &Kobj);
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  Kobj->vec_to_theta( kappa);
  Eigen::VectorXd h = Rcpp::as<Eigen::VectorXd>(operator_list["h"]);
    
  //Prior solver
  cholesky_solver Qsolver;
  Qsolver.init(Kobj->d, 0, 0, 0);
  Qsolver.analyze(Kobj->Q);
  Qsolver.compute(Kobj->Q);
  
  
  //Create solvers for each patient
  std::vector<  cholesky_solver >  Solver(nrep);
  int i = 0;
  Eigen::SparseMatrix<double,0,int> Q;
  std::vector< Eigen::SparseMatrix<double,0,int> > As(nrep);
  std::vector< Eigen::VectorXd > Ys(nrep);
  std::vector< Eigen::MatrixXd > Bs(nrep);
  std::vector< Eigen::VectorXd > Xs(nrep);
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    As[i] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    Ys[i] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    Bs[i] = Rcpp::as<Eigen::MatrixXd>(obs_tmp["B"]);
    Xs[i].resize( Kobj->d );
    Solver[i].init(Kobj->d, 0, 0, 0);
    Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Q = Q  * Kobj->Q;
    Q = Q + As[i].transpose()*As[i];
    Solver[i].analyze(Q);
    Solver[i].compute(Q);
    i++;
  }
  
  std::vector< Eigen::VectorXd > Vs(nrep);
  for( i = 0; i < V_list.length(); i++ ) 
      Vs[i] = Rcpp::as<Eigen::VectorXd>(V_list[i]);
  /*
  Simulation object
  */
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  gig rgig;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  
  Eigen::VectorXd  z;
  z.setZero(Kobj->d);
  
  Eigen::VectorXd b, Ysim;
  b.setZero(Kobj->d);
  
  std::default_random_engine gammagenerator;
  
  Eigen::SparseMatrix<double,0,int> Qi;
 
  
  Eigen::MatrixXd B;
  Eigen::SparseMatrix<double,0,int> A;
  
 
  //Update operator
  Kobj->vec_to_theta( kappa);
  Qsolver.compute(Kobj->Q);
  std::vector< Eigen::MatrixXd > V_out(nrep), X_out(nrep);
  std::vector<Eigen::MatrixXd > U_out(nrep) ;  
  std::vector<Eigen::MatrixXd > sigma_out(nrep) ;  
    
  for( int i = 0; i < nrep; i++ ) {
    Rcpp::Rcout << "*";
    U_out[i].setZero(nsim, mixobj.U.rows());
    sigma_out[i].setZero(nsim,1);
    A = As[i];
    B  = Bs[i];
    Eigen::VectorXd Y = Ys[i];
    //compute mean
    Y -= B * beta;
      
    Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Eigen::VectorXd iV(Vs[i].size());
    iV.array() = Vs[i].array().inverse();
      
    Q = Q  * iV.asDiagonal();
    Q = exp(tau)*Q  * Kobj->Q;   
    V_out[i].setZero(nsim, Kobj->d);
    X_out[i].setZero(nsim, Kobj->d);
  
    if(usingMixedEff)
      mixobj.remove_inter(i, Y);
          
    if(commonsigma==0)
        sigma = sigma_v[i];
          
    for(int ii = 0; ii < nsim; ii ++){
      //Sample X|Y, V
        
      for(int j =0; j < Kobj->d; j++)
        z[j] =  normal(random_engine);
      
      sampleX(Xs[i], 
              z,
              Vs[i],
              Y,
              Q,
              A,
              exp(sigma),
              Solver[i]);
              
      X_out[i].row(ii) = Xs[i];     
      // sample V| X
      if(noise>=0){ 
        Vs[i] =   sampleV_post(rgig,
                             h, 
                             Kobj->Q*Xs[i],
                             exp(- 0.5 * tau),
                             mu,
                             exp(lambda),
                            GAL);
      }                      
      V_out[i].row(ii) = Vs[i]; 
      Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
      iV.array() = Vs[i].array().inverse();
      Q = Q  * iV.asDiagonal();
      Q =  exp(tau)* Q  * Kobj->Q;   
      
      if(usingMixedEff)
      {
          mixobj.add_inter(i, Y);
          mixobj.sampleU(i, Y - A * Xs[i], sigma);
          mixobj.remove_inter(i, Y);
          U_out[i].row(ii) = mixobj.U.col(i);
          //U_out(i, ii) = U[i];
      }
      Ysim = Y - A * Xs[i];
      // sample sigma | U,V,X  
      if(commonsigma==0){
        double alphai = exp(asigma) + Ysim.size()/2.0;
        double betai = exp(bsigma) + 0.5* Ysim.array().square().sum();
        std::gamma_distribution<double> distribution(alphai,1.0/betai);
        sigma = -log(distribution(gammagenerator)); 
        sigma_out[i](ii,1) = exp(0.5*sigma);
      }
    }
  }
  Rcpp::Rcout << "\n";
  Rcpp::List list;
  list["X"] = X_out;
  list["V"]  = V_out;
  list["U"]  = U_out;
  list["sigma"]  = sigma_out;
  delete Kobj;
  return(list);
}