#include "Qmatrix.h"
#include "error_check.h"

using namespace std;

void constMatrix::initFromList(Rcpp::List const & init_list)
{
 std::vector<std::string> check_names =  {"Q"};
  check_Rcpplist(init_list, check_names, "constMatrix::initFromList");
  Q  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["Q"]);
  d = Q.rows();
  npars = 0;
  v.setZero(1);
  m.resize(1,1);
}

void constMatrix::initFromList(Rcpp::List const & init_List, Rcpp::List const & solver_list) 
{
  this->initFromList(init_List);  
}

Rcpp::List constMatrix::output_list()
{
  Rcpp::List List;
  return(List);
}
