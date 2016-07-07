#ifndef __Q__MATRIX__
#define __Q__MATRIX__
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "localprint.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "MatrixAlgebra.h"
#include "solver.h"
#ifdef TIME_IT
	#include "TimeIt.h"
#endif

#ifdef _OPENMP
	#include<omp.h>
#endif

class Qmatrix {
  protected:
    solver * Qsolver;
  public:
    Qmatrix() {Qsolver = NULL;};
    virtual ~Qmatrix(){delete Qsolver;};
    int d; //dimension
    int npars; //number of parameters
    Eigen::SparseMatrix<double,0,int> Q; // the generic matrix object
    Eigen::MatrixXd K;                   // the generic matrix object if Q is full!
    virtual void initFromList(Rcpp::List const &)=0;
    virtual void initFromList(Rcpp::List const &, Rcpp::List const &) {Rcpp::Rcout << "initFromList(list1,list2) not implimented in Qmatrix\n";};
    virtual Eigen::SparseMatrix<double,0,int> df(int)=0;
    virtual Eigen::SparseMatrix<double,0,int> d2f(int,int)=0;
    virtual double trace(int)=0;
    virtual double trace2(int,int)=0;

    virtual void vec_to_theta(const Eigen::VectorXd&)=0;
    virtual Eigen::VectorXd get_theta() = 0;
    virtual Rcpp::List output_list() = 0;
    virtual double f(const Eigen::VectorXd &)=0;
    virtual void show_results()=0;
    virtual double logdet()=0;
    // x  = A^-1b
    virtual Eigen::VectorXd solve(Eigen::VectorXd &b, Eigen::VectorXd &Guess) {return(Qsolver->solve(b, Guess));};
    virtual void update(){};
    virtual void update_param(){};
    
    virtual void gradient( const Eigen::VectorXd &, const Eigen::VectorXd & ){};
    virtual void step_theta(const double ){};
    
    
    double tau;
};

class constMatrix : public Qmatrix{
  protected:
    Eigen::SparseMatrix<double,0,int> m;
    Eigen::VectorXd v;
  public:
  
	void gradient(const Eigen::VectorXd &, const Eigen::VectorXd & );
  
	void step_theta(const double);
  	double  dtau;
  	double ddtau;
    void initFromList(Rcpp::List const &);
    Eigen::SparseMatrix<double,0,int> df(int){return m;};
    Eigen::SparseMatrix<double,0,int> d2f(int,int){return m;};
    void initFromList(Rcpp::List const &, Rcpp::List const &);
    double trace(int){return 0.0;};
    double trace2(int,int){return 0.0;};
    double f(const Eigen::VectorXd & theta_in){return 0.0;};

    void vec_to_theta(const Eigen::VectorXd&){};
    Eigen::VectorXd get_theta(){return v;};
    Rcpp::List output_list();
    void show_results(){};
    double logdet(){return 0;};
};


class MaternMatrixOperator : public Qmatrix{
  protected:
    double ldet;
    Eigen::VectorXd g,p;
    Eigen::SparseMatrix<double,0,int> G, C;
    Eigen::VectorXd kpv, phiv, dkappa, kappa, theta;
    Eigen::MatrixXd d2kappa,  Bkp;
    int nkp, nphi;
    Eigen::DiagonalMatrix<double,Dynamic> Kappa;
    Eigen::MatrixXd H;
    bool use_chol;
    double counter;
    int calc_det;
    SparseMatrix<double,0,int> * Phik, * Mk;
    VectorXd phi2;
    double dtau;
  	double ddtau;
    

  Eigen::VectorXd kappa_vec(Eigen::MatrixXd & B, Eigen::VectorXd & beta);
  Eigen::MatrixXd dkappa_mat(Eigen::MatrixXd & B, Eigen::VectorXd & beta); 

  public:
  
  	double tau;
  	int n;
    MaternMatrixOperator(){ counter = 0;};
    void initFromList(Rcpp::List const &);
    void initFromList(Rcpp::List const &, Rcpp::List const &);
    Eigen::SparseMatrix<double,0,int> df(int);
    Eigen::SparseMatrix<double,0,int> d2f(int,int);
    double trace(int);
    double trace2(int,int){ Rcpp::Rcout << "MaternMatrixOperator trace2 not defined\n";return(std::numeric_limits<double>::infinity());};
    double f(const Eigen::VectorXd &);

    void vec_to_theta(const Eigen::VectorXd&);
    Eigen::VectorXd get_theta();
    Rcpp::List output_list();
    void show_results();
    double logdet();
    
    void gradient( const Eigen::VectorXd &, const Eigen::VectorXd & );
    void step_theta(const double );
};


#endif