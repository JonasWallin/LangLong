#include "Qmatrix.h"
#include "error_check.h"

using namespace std;

void constMatrix::initFromList(Rcpp::List const & init_list)
{
 std::vector<std::string> check_names =  {"Q","loc", "h"};
  check_Rcpplist(init_list, check_names, "constMatrix::initFromList");
  Q  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["Q"]);
  d = Q.rows();
  npars = 0;
  v.setZero(1);
  m.resize(1,1);
  tau = 1.;
  if(init_list.containsElementNamed("tau"))
  	tau = Rcpp::as<double >( init_list["tau"]);
  
  dtau  = 0.;
  ddtau = 0.;
  
  loc  = Rcpp::as< Eigen::VectorXd >(init_list["loc"]);
  h  = Rcpp::as< Eigen::VectorXd >(init_list["h"]);
  h_average = h.sum() / h.size();
  m_loc = loc.minCoeff();
}

void constMatrix::initFromList(Rcpp::List const & init_List, Rcpp::List const & solver_list) 
{
  this->initFromList(init_List);  
}

void constMatrix::gradient( const Eigen::VectorXd & X, const Eigen::VectorXd & iV)
{
  Eigen::VectorXd vtmp = Q * X;
  
  double xtQx =  vtmp.dot(iV.asDiagonal() * vtmp); 
  dtau  	  +=  0.5 * d / tau;
  dtau        -=  0.5 * xtQx;
  ddtau       -=  0.5 * d / pow(tau, 2); 
}


double constMatrix::trace_variance( const Eigen::SparseMatrix<double,0,int>& A){
	Eigen::VectorXd  obs_loc = A * loc;
	obs_loc.array() -= m_loc;
	// covariance / h
	// covariance = t^3/3 (integrated brownian moition
	return(obs_loc.array().pow(3).sum() / (3 * h_average * tau)); 
	
}
void constMatrix::step_theta(const double stepsize)
{

	dtau  /= ddtau;
  dtau *= stepsize;
	double tau_temp = -1.;
    while(tau_temp < 0)
    {
    	dtau *= 0.5;
        tau_temp = tau - dtau;
    }
	tau = tau_temp;
	dtau   = 0;
	ddtau  = 0;
}

Rcpp::List constMatrix::output_list()
{
  Rcpp::List  List;
  List["tau"] = tau;
  return(List);
}
