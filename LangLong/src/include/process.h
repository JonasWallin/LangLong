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
  
  	int nindv; // number indiviuals 
  	std::vector< Eigen::VectorXd > Xs;
  	std::vector< Eigen::VectorXd >  Vs;
  	Eigen::SparseMatrix<double,0,int>  Q;
  	Eigen::VectorXd  h;
  	Eigen::VectorXd  iV;
    Process() {};
    virtual ~Process(){};
    
    virtual void gradient( const int i ){};
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
};



class GHProcess : public Process{



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
	
    void sample_V(const int, 
    			  gig &,
                  const Eigen::SparseMatrix<double,0,int> &);
};



#endif