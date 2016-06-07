#include <Rcpp.h>
#include <random>
#include <chrono>
#include "GIG.h"
using namespace Rcpp;


double dlambda_V(const double loglambda,
                 const Eigen::VectorXd &V, 
                 const Eigen::VectorXd &h,
                 const int GAL)
{
  double dlambda = 0;
  double Vadj = 1e-10;
  if(GAL){
  for(int i=0; i < h.size(); i++){
    double h_lambda = exp(loglambda) * h[i];
    //digamma(0.1) = digamma(1.1) - 1/0.1;
    if(h_lambda > 1){
        dlambda -=  h_lambda * R::digamma(h_lambda);
      }else
      {
        dlambda -=  h_lambda * R::digamma(h_lambda + 1) - 1.;
      }
    dlambda += h_lambda *  log(V(i) - Vadj + 1e-10 ) ;
  }
  }else{
    double srqt_two = pow(2, 0.5);
    for(int i=0; i < h.size(); i++){
      dlambda +=  1 -  ( pow(h(i), 2) / V(i) ) * exp( 2 * loglambda);
      dlambda += srqt_two * h(i)  * exp(loglambda);  
    }
      
  }
  
  return(dlambda); 
}

Eigen::VectorXd sampleV_pre(gig &sampler,
                        const Eigen::VectorXd &h, 
                        const double tau,
                        const int GAL)
{
  Eigen::VectorXd V(h.size());
  double Vadj = 1e-10;
  if(GAL)
  {
    for(int i = 0; i < h.size(); i++)
      V[i] = sampler.sample( h[i] * tau , 2, 0) + Vadj; 
  }else{
    for(int i = 0; i < h.size(); i++)
      V[i] = sampler.sample(-0.5 , 2, pow(h[i] * tau, 2)); 
  }
  
  return(V);
}


Eigen::VectorXd sampleV_post(gig &sampler,
                        const Eigen::VectorXd &h, 
                        const Eigen::VectorXd &KX,
                        const double sigma,
                        const double mu,
                        const double tau,
                        const int GAL)
{
  Eigen::VectorXd  p, b;
  double Vadj = 1e-10;
  b = KX / sigma;
  b = b.array().square();
  double a  =  pow(mu / sigma, 2)  + 2.;
  if(GAL){
    p = h * tau;
    p.array() -= 0.5;
  }else{
    p.setOnes(h.size());
    p *= -1.;
    b.array() += (h * tau).array().square();
  }
  Eigen::VectorXd V(KX.size());
  double b_adj = 1e-14;
  for(int i = 0; i < KX.size(); i++)
      V[i] = sampler.sample( p[i], a, b[i] + b_adj) + Vadj; 
  
  return(V);
}
