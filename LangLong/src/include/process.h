#ifndef __PROCESS__
#define __PROCESS__


#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "Qmatrix.h"
#include "solver.h"
#include "GIG.h"
#include "rgig.h"

//object for storing the processes object for several classes
class Process {

  public:
  
  	int store_param;
  	int nindv; // number indiviuals 
  	std::vector< Eigen::VectorXd > Xs;
  	std::vector< Eigen::VectorXd >  Vs;
  	Eigen::SparseMatrix<double,0,int>  Q;
  	Eigen::VectorXd  h;
  	Eigen::VectorXd  iV;
  	std::string type_process;
    Process() {};
    
    // setups to store the tracjetory
    virtual void setupStoreTracj(const int Niter){};
    virtual ~Process(){};
    virtual Rcpp::List toList() {};
    virtual void gradient( const int i ,
    					   const Eigen::SparseMatrix<double,0,int> & K){};
    virtual void step_theta(const double step){};
    virtual void sample_V(const int, 
    					            gig &,
                          const Eigen::SparseMatrix<double,0,int> &){};
    
    
    virtual void initFromList(const Rcpp::List &, const Eigen::VectorXd  &){};
    // sampling where the measurement noise is normal
    virtual void sample_X(const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver){};
    // sampling where the measurement noise is non-normal
    virtual void sample_Xv2(  const int i, 
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise ){};
    

};

class GaussianProcess : public Process{

	void initFromList(const Rcpp::List  &,const  Eigen::VectorXd &);
	void sample_X( const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver);
    void sample_Xv2(  const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise);
	
    Rcpp::List toList();
};



class GHProcess : public Process{


	private:
	
		double dmu ;
		double dnu, ddnu ;
		Eigen::VectorXd EV;
		Eigen::VectorXd EiV;
		Eigen::VectorXd mu_vec;
		Eigen::VectorXd nu_vec;
		int vec_counter;
		double counter;
	public:
	
	
	double mu;
	double nu;
	void initFromList(const Rcpp::List  &, const Eigen::VectorXd &);
	void sample_X( const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver);
    
    
    void sample_Xv2(  const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise);

    void gradient( const int i ,
    					   const Eigen::SparseMatrix<double,0,int> & K);
    void step_theta(const double );	
    void step_mu(const double );	
    void step_nu(const double );
    Rcpp::List toList();
    void setupStoreTracj(const int);
    void sample_V(const int, 
    			  gig &,
                  const Eigen::SparseMatrix<double,0,int> &);
};



#endif