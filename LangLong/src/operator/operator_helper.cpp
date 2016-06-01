#include "operator_helper.h"

void func_dkappa(std::string type_operator,
                 Eigen::VectorXd &dkappa, 
                 const Eigen::VectorXd &X,
                 const Eigen::VectorXd &iV, 
                 Qmatrix &Kobj,
                 double tau)
{
  if(type_operator == "matern"){
    dkappa[0] -=  exp(tau)*X.transpose() * (Kobj.df(0) * iV.asDiagonal() * Kobj.Q) * X;
  }
}

void func_dkappa0(std::string type_operator,
                  Eigen::VectorXd &dkappa,
                  int nsim,
                  int nrep,
                  Qmatrix &Kobj)
{
  if(type_operator == "matern")
    dkappa[0] = nsim*nrep*Kobj.trace(0);
}

void operator_select(std::string type_operator, Qmatrix **Kobj)
{ 
  if(type_operator == "matern"){
        *Kobj = new MaternMatrixOperator;
  }else{
        *Kobj = new constMatrix;
  }
        
        
  
}

void print_kappa(std::string type_operator, Eigen::VectorXd &kappa)
{
  if(type_operator == "matern")
    Rcpp::Rcout << ", kappa = " << kappa.exp() << "\n";
}