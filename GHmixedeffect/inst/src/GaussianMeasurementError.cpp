#include "measError.h"
#include "error_check.h"





GaussianMeasurementError::GaussianMeasurementError(){
  counter = 0;
  sigma   = 0;
  dsigma  = 0;
  ddsigma = 0;
  EV  = 1.;  // if there the random variance in the Noise E[V]
  EiV = 1.; 
  noise = "Normal";
} 
Rcpp::List GaussianMeasurementError::toList()
{
  Rcpp::List out;
  out["sigma"]  = sigma;
  out["noise"]  = noise;
  return(out);
}
void GaussianMeasurementError::initFromList(Rcpp::List const &init_list)
{
  if(init_list.containsElementNamed("sigma"))
    sigma = Rcpp::as < double >( init_list["sigma"]);
  else
    sigma = 1.;
}

void GaussianMeasurementError::gradient(const int i, 
                                 const Eigen::VectorXd& res)
{
    counter++;
    dsigma += - res.size()/sigma + res.array().square().sum() / pow(sigma, 3);
    // Expected fisher infromation
    // res.size()/pow(sigma, 2) - 3 * E[res.array().square().sum()] /pow(sigma, 4);
    ddsigma += - 2 * res.size()/pow(sigma, 2);
}
void GaussianMeasurementError::step_theta(double stepsize)
{
  double sigma_temp = -1;
  dsigma /= ddsigma;
  while(sigma_temp < 0)
  {
    sigma_temp = sigma - stepsize * dsigma;
    stepsize *= 0.5;
    if(stepsize <= 1e-16)
        throw("in GaussianMeasurementError:: can't make sigma it positive \n");   
  }
  sigma = sigma_temp;
  counter = 0;
  dsigma  = 0;
  ddsigma = 0;
}


std::vector< Eigen::VectorXd > GaussianMeasurementError::simulate(std::vector< Eigen::VectorXd > Y)
{
	std::vector< Eigen::VectorXd > residual( Y.size());
	for(int i = 0; i < Y.size(); i++)
		residual[i] =  sigma * (Rcpp::as< Eigen::VectorXd >(Rcpp::rnorm( Y[i].size()) ));
    
	return(residual);
}