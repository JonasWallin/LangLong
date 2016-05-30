#ifndef __MEASE__
#define __MEASE__
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include "MatrixAlgebra.h"
#include "GIG.h"
class MeasurementError {
  protected:
  
  public:
  
   	double EV;  // if there the random variance in the Noise E[V]
    double EiV; // if there is random varoance in the noise E[V^-1]   
  	double sigma;
  	std::vector< Eigen::VectorXd > Vs;
    std::string noise;
    virtual void gradient(const int , const Eigen::VectorXd& ) = 0;
    virtual void step_theta(double stepsize) = 0;
    virtual void initFromList(Rcpp::List const &)=0;
    virtual Rcpp::List toList()=0;
    virtual void sampleV(const int , const Eigen::VectorXd& ) = 0;
  
  
};


class GaussianMeasurementError : public MeasurementError{
	private:
		double dsigma;
		double ddsigma;
    double counter;

	public:
		GaussianMeasurementError();
		void gradient(const int , const Eigen::VectorXd&);
		void step_theta(double stepsize);
		void initFromList(Rcpp::List const &);
		void sampleV(const int i, const Eigen::VectorXd& res) {};
		Rcpp::List toList();

};


class NIGMeasurementError : public MeasurementError{
  private:
		double dsigma;
		double ddsigma;
		double dnu;
		double ddnu;
		gig rgig;
    	double counter;

	public:
		double nu;
		NIGMeasurementError();
		void gradient(const int , const Eigen::VectorXd&);
		void step_theta(double stepsize);
		void step_sigma(double stepsize);
		void step_nu(double stepsize);
		void initFromList(Rcpp::List const &);
		void sampleV(const int , const Eigen::VectorXd& );
		Rcpp::List toList();

};

#endif