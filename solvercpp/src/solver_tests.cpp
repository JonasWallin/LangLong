#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include "MatrixAlgebra.h"
#include "solver.h"
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;
using namespace Eigen;


/*
// [[Rcpp::export]]
List chol_blockmatrix(List init_list_in)
{
   List init_list = clone(init_list_in);
   
   S = Rcpp::as<Eigen::MatrixXd>(init_list["S"]);
   nop = S.sum();
   Clist =      new SparseMatrix<double,0,int>[nop];
   for(int i=0;i<nop;i++){
    Clist[i] =  Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["C"]);
  }
   SparseMatrix<double,0,int> C;
   C.resize(d,d);
}
*/
// [[Rcpp::export]]
double traceWithLUSparse(Eigen::MatrixXd A, Eigen::SparseMatrix<double,0,int>  B)
{
  lu_solver SolverObj;

  SolverObj.init(A.cols(), 0, 0, 0);
  SolverObj.compute(A);
  return(SolverObj.trace(B));
}

// [[Rcpp::export]]
double traceWithLU(Eigen::MatrixXd A, Eigen::MatrixXd B)
{
  lu_solver SolverObj;

  SolverObj.init(A.cols(), 0, 0, 0);
  SolverObj.compute(A);
  return(SolverObj.trace(B));
}

// [[Rcpp::export]]
Eigen::VectorXd solveWithLU(Eigen::MatrixXd A, Eigen::VectorXd b)
{
  lu_solver SolverObj;

  SolverObj.init(A.cols(), 0, 0, 0);
  SolverObj.compute(A);
  return(SolverObj.solve(b,b));
}

// [[Rcpp::export]]
double logdetWithLU(Eigen::MatrixXd A)
{
  lu_solver SolverObj;

  SolverObj.init(A.cols(), 0, 0, 0);
  SolverObj.compute(A);
  return(SolverObj.logdet());
}

// [[Rcpp::export]]
double LU_sparse_logdet(Eigen::SparseMatrix<double,0,int> A)
{
  List  crap;
  solver * R;
  R = new lu_sparse_solver;
  R->initFromList(A.rows(),crap);
  R->analyze(A);
  R->compute(A);
  
  double d = R->logdet();
  delete R;
  
  return d;
}
// [[Rcpp::export]]
List LU_sparse_tests(List init_list_in)
{
  List init_list = clone(init_list_in);
  SparseMatrix<double,0,int>  A  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["A"]);
  
  solver * R;
  R = new lu_sparse_solver;
  R->initFromList(A.rows(),init_list);
  R->analyze(A);
  R->compute(A);
  Rcpp::List outList;
  if(init_list.containsElementNamed("B")==1){
    Eigen::SparseMatrix<double,0,int> B  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["B"]);
    double trace_AinvB = R->trace(B);
    outList["traceAinvB"] =trace_AinvB;
  }
  if(init_list.containsElementNamed("b")==1){
  
    VectorXd b  = Rcpp::as<VectorXd >(init_list["b"]);
    VectorXd x  = R->solve(b,b);
    outList["x"] =x;
  }
  return(outList);
}


// [[Rcpp::export]]
List solver_tests(List init_list_in)
{
  List init_list = clone(init_list_in);
  SparseMatrix<double,0,int> Q;
  solver * R;

  Q  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["Q"]);

  int solvertype = Rcpp::as<int>(init_list["type"]);
  int N = Rcpp::as<int>(init_list["trace.iter"]);
  double tol = Rcpp::as<double>(init_list["tol"]);
  int max_iter = Rcpp::as<int>(init_list["solver.max.iter"]);
  int n = Q.rows();

  if(solvertype==0) {
    R = new cholesky_solver;
  } else {
    R = new iterative_solver;
  }
  (*R).init(n,N,max_iter,tol);

  (*R).analyze(Q);
  (*R).compute(Q);

  int operation = Rcpp::as<int>(init_list["operation"]);

  Rcpp::List outList;
  if(operation==0){ //solve
    VectorXd v, X;
    v  = Rcpp::as<Eigen::VectorXd >(init_list["v"]);
    X = (*R).solve(v,v);
    outList["X"] = X;
  } else if(operation==1){ //trace
    SparseMatrix<double,0,int> M;
    M  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["M"]);
    double trace = (*R).trace(M);
    outList["trace"] = trace;
  } else if(operation==2){ //trace2
    SparseMatrix<double,0,int> M1, M2;
    M1 = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["M1"]);
    M2 = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["M2"]);
    double trace = (*R).trace2(M1,M2);
    outList["trace2"] = trace;
  } else if(operation==3){ //Qinv
    VectorXd vars = (*R).Qinv_diag();
    outList["vars"] = vars;
  } else if(operation==4){ //logdet
    double ld = (*R).logdet();
    outList["logdet"] = ld;
  } else if(operation == 5)
  {
    outList["Qinv"] = (*R).return_Qinv();
    
  }
  return(outList);
}


// [[Rcpp::export]]
List kronecker_tests(List init_list_in)
{
  List init_list = clone(init_list_in);

  SparseMatrix<double,0,int> M1, M2,Ma,Mb;
  M1 = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["M1"]);
  M2 = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["M2"]);

  const double ticks_per_ms = static_cast<double>(CLOCKS_PER_SEC);

  clock_t start = clock();
  Ma = kroneckerProduct(M1,M2);
  double time_Ma = static_cast<double>(clock()-start)  / ticks_per_ms;

  start = clock();
  Mb = kronecker(M1,M2);
  double time_Mb = static_cast<double>(clock()-start)  / ticks_per_ms;

  Rcpp::List outList;
  outList["Ma"] = Ma;
  outList["Mb"] = Mb;
  outList["time.Ma"] = time_Ma;
  outList["time.Mb"] = time_Mb;
  return(outList);
}
