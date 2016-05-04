#ifndef __MEFF__
#define __MEFF__
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include "MatrixAlgebra.h"
class MixedEffect {
  
  protected:
  
  public:
    Eigen::MatrixXd U;
    std::vector< Eigen::MatrixXd > B;
    virtual void initFromList(Rcpp::List const &)=0;
    virtual Rcpp::List toList()=0;
    virtual void sampleU(const int, const Eigen::VectorXd &,  const double ) = 0;
    //virtual void calcgrad() = 0;
    //virtual void init_grad() = 0;
    virtual void add_inter(const int, Eigen::VectorXd &)  = 0;
    virtual void remove_inter(const int, Eigen::VectorXd &)  = 0;
  
  
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
  public:
    Eigen::MatrixXi D;
    Eigen::MatrixXd Dd;
  
    Eigen::MatrixXd Sigma;
  
    NormalMixedEffect();
    void initFromList(Rcpp::List const &);
    void sampleU(const int, const Eigen::VectorXd &, const double sigma) ;
    void remove_inter(const int i, Eigen::VectorXd & Y) {Y -= B[i]*U.col(i);} ;
    void add_inter(const int i, Eigen::VectorXd & Y)    {Y += B[i]*U.col(i);} ;
    void gradient(const int , const Eigen::VectorXd& );
    void step_theta(double stepsize);
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
  public:
    Eigen::MatrixXi D;
    Eigen::MatrixXd Dd;
    Eigen::VectorXd V;
  
    Eigen::MatrixXd Sigma;
  
    NIGMixedEffect();
    void sampleV(const int);
    void initFromList(Rcpp::List const &);
    void sampleU(const int, const Eigen::VectorXd &, const double sigma) ;
    void remove_inter(const int i, Eigen::VectorXd & Y) {Y -= B[i]*U.col(i);} ;
    void add_inter(const int i, Eigen::VectorXd & Y)    {Y += B[i]*U.col(i);} ;
    void gradient(const int , const Eigen::VectorXd& );
    void step_theta(double stepsize);
    Rcpp::List toList();
  
};

// mixed effect util

// sampling NormalCanoncial N_c(b, Q^) \propto e^{-x'Qx + bx}
Eigen::VectorXd sample_Nc(const Eigen::VectorXd & ,const  Eigen::MatrixXd & );
#endif