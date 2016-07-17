#include "process.h"
#include "error_check.h"


void GaussianProcess::initFromList(const Rcpp::List & init_list,const Eigen::VectorXd & h_in)
{
  h = h_in;
  iV = h.cwiseInverse();
  std::vector<std::string> check_names =  {"X"};
  check_Rcpplist(init_list, check_names, "GaussianProcess::initFromList");
  Rcpp::List X_list = Rcpp::as<Rcpp::List>  (init_list["X"]);
  nindv = X_list.length();
  Xs.resize(nindv);
  Vs.resize(nindv);
  for(int i = 0; i < nindv; i++ ){
    	Xs[i] = Rcpp::as<Eigen::VectorXd>( X_list[i]);
    	Vs[i] = h;
  	}
  	
  
  
}

void GHProcess::initFromList(const Rcpp::List & init_list,const  Eigen::VectorXd & h_in)
{
  h = h_in;
  h2 = h.cwiseProduct(h); 
  h_sum = h.sum();
  std::vector<std::string> check_names =  {"X","V","mu","nu"};
  check_Rcpplist(init_list, check_names, "GHProcess::initFromList");
  Rcpp::List V_list           = Rcpp::as<Rcpp::List>  (init_list["V"]);
  Rcpp::List X_list = Rcpp::as<Rcpp::List>  (init_list["X"]);
  nindv = X_list.length();
  Xs.resize(nindv);
  Vs.resize(nindv);
  for(int i = 0; i < nindv; i++ ){ 
      
    	Xs[i] = Rcpp::as<Eigen::VectorXd>( X_list[i]);
    	//Vs[i] = h;
      	Vs[i] = Rcpp::as<Eigen::VectorXd>( V_list[i]);
  	}
  	
  	mu = Rcpp::as< double > (init_list["mu"]);
  	nu = Rcpp::as< double > (init_list["nu"]);
  	
  	type_process = Rcpp::as<std::string> (init_list["noise"]);
  	dmu  = 0;
  	dnu  = 0;
  	
  	if(type_process == "NIG")
  	{
  		EV = h;
  		EiV=  h2.cwiseInverse();
      	EiV *=  1./nu;
  		EiV += h.cwiseInverse();
  	}
  	
  	counter = 0;
  	store_param = 0;
  	
}


void GaussianProcess::sample_X(const int i, 
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver)
{
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()*Y/ sigma2;
  Xs[i] = solver.rMVN(b, Z);
}



void GHProcess::sample_X(const int i,  
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver)
{
  iV = Vs[i].cwiseInverse();
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()*Y/ sigma2;
  Eigen::VectorXd temp  =  - h;
  temp = temp.cwiseProduct(iV);
  temp.array() += 1.;
  temp *= mu;
  b +=  K.transpose() * temp;
  Xs[i] = solver.rMVN(b, Z);
}

void GaussianProcess::sample_Xv2( const int i, 
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise)
{
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*iV_noise.asDiagonal()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()* (iV_noise.cwiseProduct(Y) )/ sigma2;
  Xs[i] = solver.rMVN(b, Z);
}


void GHProcess::sample_Xv2(  const int i,
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & K,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma,
              cholesky_solver       & solver,
              const Eigen::VectorXd & iV_noise )
{
	iV = Vs[i].cwiseInverse();
  double sigma2  =pow(sigma, 2);
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*iV_noise.asDiagonal()*A)/ sigma2;
  

  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()* (iV_noise.cwiseProduct(Y) )/ sigma2;
  Eigen::VectorXd temp  =  - h;
  temp = temp.cwiseProduct(iV);
  temp.array() += 1.;
  temp *= mu;
  b +=  K.transpose() * temp;
   Xs[i] =solver.rMVN(b, Z);
}

void GHProcess::sample_V(const int i , 
    					              gig & rgig,
                            const Eigen::SparseMatrix<double,0,int> & K)
{
 	Vs[i] = sampleV_post(rgig,
                 h, 
                 K * Xs[i],
                 1.,
                 mu, 
                 nu, 
                 "NIG"); 

}


void GHProcess::gradient( const int i ,
			   			  const Eigen::SparseMatrix<double,0,int> & K)
{ 
	
  	counter++;
	
	iV = Vs[i].cwiseInverse();
	Eigen::VectorXd temp_1  =  - h;
  	temp_1 = temp_1.cwiseProduct(iV);
  	
  	temp_1.array() += 1.;
  	
	Eigen::VectorXd temp_2;
	Eigen::VectorXd temp_3 =  Vs[i] ;
  	temp_3 -= h;
	temp_3.array() *= mu;
	temp_2 = K * Xs[i];
	temp_2 -= temp_3;
	dmu += temp_1.dot(temp_2);
    
	
	// dnu
	if(type_process == "NIG")
		dnu  +=  iV.size() / nu - nu * (h2.dot(iV) + Vs[i].sum() - 2 * h_sum);
    	
    
};

void GHProcess::step_theta(const double stepsize)
{
	step_mu(stepsize);	
	step_nu(stepsize);
	counter = 0;
	
	if(store_param)
	{
		mu_vec[vec_counter] = mu; 
		nu_vec[vec_counter] = nu;
		vec_counter++;
	}
}

void GHProcess::step_mu(const double stepsize)
{
   Eigen::VectorXd temp = EV;
   temp.array()  = h2;
	double H_mu =  counter * (EV.sum() - EiV.dot(temp));
	mu -= (stepsize / H_mu ) * dmu; 
	dmu = 0;
}

void GHProcess::step_nu(const double stepsize)
{
  double nu_temp = -1;
  ddnu = - iV.size()/( nu * nu) - (h2.dot(EiV) + EV.sum()) + 2 * h_sum;
  ddnu *= counter;
  dnu /= ddnu;
  double stepsize_temp  = stepsize;
  while(nu_temp < 0)
  {
    nu_temp = nu - stepsize_temp * dnu;
    stepsize_temp *= 0.5;
    if(stepsize_temp <= 1e-16)
        throw("in GHProcess:: can't make nu it positive \n");   
  }
  nu = nu_temp;
  if(type_process == "NIG"){
  	EiV =  h2.cwiseInverse();
    EiV *=  1./nu;
  	EiV += h.cwiseInverse();
  }
  dnu  = 0;
  ddnu = 0;
}


void GHProcess::setupStoreTracj(const int Niter)
{

	mu_vec.resize(Niter);
	nu_vec.resize(Niter);
	vec_counter = 0;
	store_param = 1;
}

Rcpp::List GHProcess::toList()
{
  Rcpp::List out;
  out["nu"]     = nu;
  out["mu"]     = mu;
  out["Xs"]     = Xs;
  out["Vs"]     = Vs;
  
  if(store_param)
  {
  	out["mu_vec"]     = mu_vec;
  	out["nu_vec"]     = nu_vec;
  }
  
  return(out);
}

Rcpp::List GaussianProcess::toList()
{
  Rcpp::List out;
  out["Xs"]     = Xs;
  out["Vs"]     = Vs;
  return(out);
}