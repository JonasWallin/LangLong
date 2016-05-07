#ifndef __MEFF__
#define __MEFF__
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include "MatrixAlgebra.h"
#include "GIG.h"
//TODO: fixed and random effects should be updated jointly if the are very correclated!!
class MixedEffect {
  
  protected:
  
  public:
    Eigen::MatrixXd Sigma;
    std::string noise;
    Eigen::MatrixXd U;
    std::vector< Eigen::MatrixXd > Bf; // fixed covariates
    std::vector< Eigen::MatrixXd > Br; // mixed covariates
    Eigen::VectorXd beta_random;
    Eigen::VectorXd beta_fixed;
    virtual void initFromList(Rcpp::List const &)=0;
    virtual Rcpp::List toList()=0;
    virtual void sampleU(const int, const Eigen::VectorXd &,  const double ) = 0;
    virtual void remove_cov(const int , Eigen::VectorXd & )  = 0;
    virtual void add_cov(const int    , Eigen::VectorXd & )  = 0;
    virtual void add_inter(const int, Eigen::VectorXd &)     = 0;
    virtual void remove_inter(const int, Eigen::VectorXd &)  = 0;
    
    virtual void gradient(const int , const Eigen::VectorXd&, const double ) = 0;
    virtual void step_theta(double stepsize) = 0;
  
  
};

class NormalMixedEffect  : public MixedEffect{
  private:
    Eigen::MatrixXd invSigma;
    Eigen::MatrixXd iSkroniS; // helper matrix
    int counter;
    Eigen::VectorXd UUt;
    Eigen::VectorXd dSigma_vech;
    Eigen::MatrixXd ddSigma;
    Eigen::VectorXd Sigma_vech;
    
    Eigen::VectorXd grad_beta_r; // gradient for random intercept
    Eigen::VectorXd grad_beta_r2; //second gradient for random intercept
    Eigen::VectorXd grad_beta_f; // gradient for fixed intercept
    Eigen::MatrixXd H_beta_random; // obsereved fisher infromation for random effect
    Eigen::MatrixXd H_beta_fixed;// obsereved fisher infromation for fixed effect
  public:
  
    Eigen::MatrixXi D;
    Eigen::MatrixXd Dd;
  
  
    NormalMixedEffect();
    void initFromList(Rcpp::List const &);
    void sampleU(const int, const Eigen::VectorXd &, const double ) ;
    void remove_inter(const int i, Eigen::VectorXd & Y) {Y -= Br[i]*U.col(i);} ;
    void add_inter(const int i, Eigen::VectorXd & Y)    {Y += Br[i]*U.col(i);} ;
    void remove_cov(const int , Eigen::VectorXd & );
    void add_cov(const int    , Eigen::VectorXd & );
    void gradient(const int , const Eigen::VectorXd&, const double );
    
    void step_theta(double stepsize);
    void step_Sigma(double stepsize);
    void step_beta_fixed(double stepsize);
    void step_beta_random(double stepsize);
    Rcpp::List toList();
  
};



class NIGMixedEffect  : public MixedEffect{
  private:
    Eigen::MatrixXd invSigma;
    Eigen::MatrixXd iSkroniS; // helper matrix
    int counter;
    Eigen::VectorXd UUt;
    Eigen::VectorXd dSigma_vech;
    Eigen::MatrixXd ddSigma;
    Eigen::VectorXd Sigma_vech;
    
    
    double  grad_nu; // gradient for shape parameter
    Eigen::VectorXd gradMu; // gradient for skewness 
    Eigen::VectorXd gradMu_2;// second gradient for skewness 
    Eigen::VectorXd grad_beta_r; // gradient for random intercept
    Eigen::VectorXd grad_beta_r2; //second gradient for random intercept
    Eigen::VectorXd grad_beta_f; // gradient for fixed intercept
    Eigen::MatrixXd H_beta_random; // obsereved fisher infromation for random effect
    Eigen::MatrixXd H_beta_fixed;// obsereved fisher infromation for fixed effect
    double EV; //  prior expecation  of V, used for the Hessian of random effects
    double EiV; // prior expectation of 1/V used for the Hessian of random effects
    double VV; // prior variance of V used for the Hessian of random effects
    double a_GIG;
    gig rgig;
    
    void step_Sigma(double stepsize);
    void step_mu(double stepsize);
    void step_nu(double stepsize);
    
  public:
    Eigen::MatrixXi D;
    Eigen::MatrixXd Dd;
    Eigen::VectorXd V;
    Eigen::VectorXd mu;
    double          nu;
  
  
    NIGMixedEffect();
    void sampleV(const int);
    void initFromList(Rcpp::List const &);
    void sampleU(const int, const Eigen::VectorXd &, const double) ;
    void remove_inter(const int i, Eigen::VectorXd & Y) {Y -= Br[i]*U.col(i);} ;
    void add_inter(const int i, Eigen::VectorXd & Y)    {Y += Br[i]*U.col(i);} ;
    void remove_cov(const int , Eigen::VectorXd & );
    void add_cov(const int    , Eigen::VectorXd & );
    void gradient(const int , const Eigen::VectorXd& , const double );
    void gradient_sigma(const int , Eigen::VectorXd& );
    void step_theta(double stepsize);
    void step_beta_fixed(double stepsize);
    void step_beta_random(double stepsize);
    Rcpp::List toList();
  
  
  
};

// mixed effect util

// sampling NormalCanoncial N_c(b, Q^) \propto e^{-x'Qx + bx}
Eigen::VectorXd sample_Nc(const Eigen::VectorXd & ,const  Eigen::MatrixXd & );
#endif