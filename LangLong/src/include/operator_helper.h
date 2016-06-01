#ifndef __OPERATOR_HELP__
#define __OPERATOR_HELP__
#include <Rcpp.h>
#include <string>
#include "Qmatrix.h"

// computes the gradient for operator paramters
void func_dkappa(std::string type_operator,
                 Eigen::VectorXd &dkappa, 
                 const Eigen::VectorXd &X,
                 const Eigen::VectorXd &iV, 
                 Qmatrix &Kobj,
                 double tau);
                 
// inits the gradient for operator parameters
void func_dkappa0(std::string type_operator,
                  Eigen::VectorXd &dkappa,
                  int nsim,
                  int nrep,
                  Qmatrix &Kobj);
// starting up operator
void operator_select(std::string type_operator, Qmatrix **Kobj);

// printing parameter
void print_kappa(std::string type_operator, Eigen::VectorXd &kappa);
#endif