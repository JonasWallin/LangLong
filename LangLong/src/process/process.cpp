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
      	//Vs[i] = Rcpp::as<Eigen::VectorXd>( V_list[i]);
  	}
  
}

void GHProcess::initFromList(const Rcpp::List & init_list,const  Eigen::VectorXd & h_in)
{
  h = h_in;
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
  temp *= iV;
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
  temp *= iV;
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