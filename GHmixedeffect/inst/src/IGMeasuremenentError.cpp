#include "measError.h"
#include "error_check.h"



double Digamma(double x)
{
  return(R::digamma(x));
}

double Trigamma(double x)
{
  return(R::trigamma(x));
}

IGMeasurementError::IGMeasurementError() : NormalVarianceMixtureBaseError(){
  nu        = 1;
  dnu       = 0;
  ddnu      = 0;
  noise = "IG";
  
}



void IGMeasurementError::printIter() 
{
	NormalVarianceMixtureBaseError::printIter();
	Rcpp::Rcout << "\n nu = " << nu;

}
void IGMeasurementError::setupStoreTracj(const int Niter) // setups to store the tracjetory
{
	
	NormalVarianceMixtureBaseError::setupStoreTracj(Niter);
	nu_vec.resize(Niter);
}





Rcpp::List IGMeasurementError::toList()
{
  Rcpp::List out = NormalVarianceMixtureBaseError::toList();
  out["nu"]          = nu;
  
  if(store_param)
  	out["nu_vec"]    = nu_vec;
  
  return(out);
}
void IGMeasurementError::initFromList(Rcpp::List const &init_list)
{
  
  NormalVarianceMixtureBaseError::initFromList(init_list);
  
  if(init_list.containsElementNamed("nu"))
    nu = Rcpp::as < double >( init_list["nu"]);
  else
    nu = 1.;

    EV  = 1.; // not true it is the mode is alpha/(alpha - 1)
    EiV = 1.;

   npars += 1;
  digamma_nu  =  Digamma(nu);
  trigamma_nu =  Trigamma(nu);
  
 int i = 0;

 if( init_list.containsElementNamed("Vs" )){
 	Rcpp::List Vs_list = init_list["Vs"];
 	Vs.resize(Vs_list.length());
    for( Rcpp::List::iterator it = Vs_list.begin(); it != Vs_list.end(); ++it ) {
      Vs[i++] = Rcpp::as < Eigen::VectorXd >( it[0]);
    }
 }else
 	  throw("in IGMeasurementError::initFromList Vs must be set! \n");


}

double IGMeasurementError::simulate_V()
{
	return rgig.sample(-nu, 0, 2 * nu );
}

double IGMeasurementError::sample_V(const double res2_j, const int n_s)
{
	if(common_V == 0)
		return rgig.sample(-(nu + .5), 0 , res2_j + 2 * nu);
	
	return rgig.sample(-  (nu + .5 * n_s), 0, res2_j + 2 * nu );
}





void IGMeasurementError::gradient(const int i,
                                 const Eigen::VectorXd& res)
{
    NormalVarianceMixtureBaseError::gradient(i, res);
    Eigen::VectorXd iV = Vs[i].cwiseInverse();
    if(common_V == 0){
    	double logV = Vs[i].array().log().sum();
    	dnu  +=  res.size() *  (1. + log(nu)  - digamma_nu) - logV - iV.array().sum() ;
    	ddnu += res.size() * (1/ nu - trigamma_nu);
    }else{
    
    	double logV = log(Vs[i][0]);
    	dnu  +=   (1. + log(nu)  - digamma_nu) - logV - iV[0] ;
    	ddnu +=  (1/ nu - trigamma_nu);
    }
}

void IGMeasurementError::step_nu(double stepsize)
{
double nu_temp = -1;
  dnu /= ddnu;
  while(nu_temp < 0)
  {
    nu_temp = nu - stepsize * dnu;
    stepsize *= 0.5;
    if(stepsize <= 1e-16)
        throw("in IGMeasurementError:: can't make nu it positive \n");
  }
  nu = nu_temp;
  EV  = 1.;  // not true it is the mode that is 1.
  EiV = 1. ;
  ddnu = 0;
  digamma_nu  =  Digamma(nu);
  trigamma_nu =  Trigamma(nu);

}

void IGMeasurementError::step_theta(double stepsize)
{
  NormalVarianceMixtureBaseError::step_theta(stepsize);
  
  step_nu(stepsize);
  clear_gradient();
  
if(store_param)
  	nu_vec[vec_counter-1] = nu; // -1 since NormalVarianceMixtureBaseError increase vec_counter
  
}

void IGMeasurementError::clear_gradient()
{	
	NormalVarianceMixtureBaseError::clear_gradient();
	dnu    = 0;
}
Eigen::VectorXd IGMeasurementError::get_gradient()
{
	Eigen::VectorXd g = NormalVarianceMixtureBaseError::get_gradient();
	g[1] = dnu;
	return(g);
}