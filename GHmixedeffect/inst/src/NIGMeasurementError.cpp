#include "measError.h"
#include "error_check.h"





NIGMeasurementError::NIGMeasurementError(){
  counter = 0;
  sigma   = 1;
  nu      = 1;
  dsigma  = 0;
  ddsigma = 0;
  dnu     = 0;
  ddnu    = 0;
  noise = "NIG";
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
} 
Rcpp::List NIGMeasurementError::toList()
{
  Rcpp::List out;
  out["sigma"]  = sigma;
  out["nu"]     = nu;
  out["Vs"]     = Vs;
  return(out);
}
void NIGMeasurementError::initFromList(Rcpp::List const &init_list)
{
  if(init_list.containsElementNamed("sigma"))
    sigma = Rcpp::as < double >( init_list["sigma"]);
    
  if(init_list.containsElementNamed("nu"))
    sigma = Rcpp::as < double >( init_list["nu"]);
 
    EV  = 1.; 
    EiV = 1. + 1./nu;  
 int i = 0;
 
 if( init_list.containsElementNamed("Vs" )){
 	Rcpp::List Vs_list = init_list["Vs"];
 	Vs.resize(Vs_list.length());
    for( Rcpp::List::iterator it = Vs_list.begin(); it != Vs_list.end(); ++it ) {
      Vs[i++] = Rcpp::as < Eigen::VectorXd >( it[0]);
    }
 }else
 	  throw("in NigMeasurementError::initFromList Vs must be set! \n");   
    
    
}
void NIGMeasurementError::sampleV(const int i, const Eigen::VectorXd& res){

        for(int j = 0; j < Vs[i].size(); j++)
        	Vs[i][j] = rgig.sample(-1., nu, pow(res[j]/sigma, 2) + nu);
};
void NIGMeasurementError::gradient(const int i, 
                                 const Eigen::VectorXd& res)
{
    counter++;
    Eigen::VectorXd res_ = res;
    Eigen::VectorXd iV = Vs[i].cwiseInverse();
    //res_.array() *= iV.array();
    dsigma += - res.size()/sigma + (res_.array().square()*iV.array()).sum() / pow(sigma, 3);
    // Expected fisher infromation
    // res.size()/pow(sigma, 2) - 3 * E[res.array().square().sum()] /pow(sigma, 4);
    ddsigma += - 2 * res.size()/pow(sigma, 2);
    
    dnu  += 0.5 * ( res.size() / nu + 2 * res.size() -   (Vs[i].array().sum() + iV.array().sum()) );
    ddnu += - res.size()/( nu * nu);
}
void NIGMeasurementError::step_theta(double stepsize)
{
  step_sigma(stepsize);
  step_nu(stepsize);
  counter = 0;
}

void NIGMeasurementError::step_sigma(double stepsize)
{
double sigma_temp = -1;
  dsigma /= ddsigma;
  while(sigma_temp < 0)
  {
    sigma_temp = sigma - stepsize * dsigma;
    stepsize *= 0.5;
    if(stepsize <= 1e-16)
        throw("in NIGMeasurementError:: can't make sigma it positive \n");   
  }
  sigma = sigma_temp;
  dsigma  = 0;
  ddsigma = 0;
}

void NIGMeasurementError::step_nu(double stepsize)
{
double nu_temp = -1;
  dnu /= ddnu;
  while(nu_temp < 0)
  {
    nu_temp = nu - stepsize * dnu;
    stepsize *= 0.5;
    if(stepsize <= 1e-16)
        throw("in NIGMeasurementError:: can't make nu it positive \n");   
  }
  nu = nu_temp;
  EV  = 1.; 
  EiV = 1. + 1./nu;  
  dnu  = 0;
  ddnu = 0;
  
}
