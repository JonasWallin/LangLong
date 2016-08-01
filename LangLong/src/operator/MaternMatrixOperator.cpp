#include "Qmatrix.h"
#include "error_check.h"
#include "eigen_add_on.h"
#include "operator_helper.h"

void MaternMatrixOperator::initFromList(Rcpp::List const & init_list, Rcpp::List const & solver_list)
{
  
  std::vector<std::string> check_names =  {"C", "G", "kappa", "B.kappa", "h"};
  check_Rcpplist(init_list, check_names, "MaternMatrixOperator::initFromList");
  std::vector<std::string> check_names2 =  {"use.chol"};
  check_Rcpplist(solver_list, check_names2, "MaternMatrixOperator::initFromList");
  G  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["G"]);
  C  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["C"]);
  n = G.rows();
  d = n;
  kpv = Rcpp::as<Eigen::VectorXd>(init_list["kappa"]);
  nkp = kpv.size();	
  
  eigen_vector_from_list(kpv, init_list, "kappa");
  eigen_matrix_from_list(Bkp, init_list, "B.kappa") ;
  
  //kappa_vec  = &g_exp_transform;
  //dkappa_mat = &dg_exp_transform;
  kappa.resize(Bkp.rows());
  kappa = kappa_vec(Bkp,kpv);
  Q = G + kappa.asDiagonal()*C;
  dkappa.setZero(Bkp.rows());
  npars = nkp;

  theta = kpv;
  use_chol = Rcpp::as<int>(solver_list["use.chol"]);

  if(use_chol==1){
    Qsolver = new cholesky_solver;
  } else {
    Qsolver = new iterative_solver;
  }
  (*Qsolver).initFromList(n,solver_list);
  

  (*Qsolver).analyze(Q);
  (*Qsolver).compute(Q);
  Phik = new SparseMatrix<double,0,int>[npars];
  Eigen::MatrixXd dkappa = dkappa_mat(Bkp,kpv);


  for(int i = 0;i<nkp;i++){
    Phik[i].resize(n,n);
    for(int j=0;j<n;j++){
      Phik[i].insert(j,j) = dkappa(j,i);
    }
  }
  

  Mk = new SparseMatrix<double,0,int>[nkp];
  
  
  trace0 = this->trace(0);
    
}

void MaternMatrixOperator::initFromList(Rcpp::List const & init_list)
{
  G  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["G"]);
  C  = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(init_list["C"]);
  n = G.rows();
  d = n;
  kpv = Rcpp::as<Eigen::VectorXd>(init_list["kappa"]);

  nkp = kpv.size();

  Bkp = Rcpp::as<Eigen::MatrixXd>(init_list["B.kappa"]);
  //kappa_vec  = &g_exp_transform;
  //dkappa_mat = &dg_exp_transform;
  kappa = kappa_vec(Bkp,kpv);
  Q = G + kappa.asDiagonal()*C;

  npars = nkp;

  theta = kpv;

  use_chol = Rcpp::as<int>(init_list["use.chol"]);

  if(use_chol==1){
    Qsolver = new cholesky_solver;
  } else {
    Qsolver = new iterative_solver;
  }
  (*Qsolver).initFromList(n,init_list);
  (*Qsolver).analyze(Q);
  (*Qsolver).compute(Q);
  Phik = new SparseMatrix<double,0,int>[npars];
  
  
  Eigen::MatrixXd dkappa = dkappa_mat(Bkp,kpv);
  for(int i = 0;i<nkp;i++){
    Phik[i].resize(n,n);
    for(int j=0;j<n;j++){
      Phik[i].insert(j,j) = dkappa(j,i);
    }
  }
  

  Mk = new SparseMatrix<double,0,int>[nkp];
}

Eigen::SparseMatrix<double,0,int> MaternMatrixOperator::df(int i)
{
  if(i >= npars)
    throw("MaternMatrixOperator::df i must be less then npars");
  
  
  Mk[i] = C.transpose()*Phik[i];
  return Mk[i];
}

Eigen::SparseMatrix<double,0,int> MaternMatrixOperator::d2f(int i,int j)
{
  throw("not implimented MaternMatrixOperator::d2f\n");
}
double MaternMatrixOperator::trace(int i)
{
  if(i >= npars)
    throw("MaternMatrixOperator::df i must be less then npars");
  SparseMatrix<double,0,int> tmp1 = Phik[i]*C;
  return (*Qsolver).trace(tmp1);
}

double MaternMatrixOperator::logdet()
{
  if(calc_det == 0)
  {
    ldet = (*Qsolver).logdet();
    calc_det = 1;
  }
  return ldet;
  
};

void MaternMatrixOperator::vec_to_theta(const Eigen::VectorXd& theta_vec_in)
{
 
  int equal = 1;
  if(theta_vec_in.size() != Bkp.cols())
  {
    Rcpp::Rcout << "in MaternMatrixOperator:: theta wrong size\n";
    throw("error");
    
  }
  for(int i = 0; i < theta_vec_in.size(); i++)
  {
    
    if(theta[i] != theta_vec_in[i]){
      equal = 0;
      break;
    }
  }
  if(equal == 0){
     calc_det = 0 ;
     theta = theta_vec_in;
     kpv = theta;
     kappa = kappa_vec(Bkp,kpv);
     
      Q = G+kappa.asDiagonal()*C;
      
      (*Qsolver).compute(Q);
      Eigen::MatrixXd dkappa = dkappa_mat(Bkp,kpv);
      
      for(int i = 0;i<nkp;i++){
        Phik[i].resize(n,n);
        for(int j=0;j<n;j++){
          Phik[i].insert(j,j) = dkappa(j,i);
        }
      }
  }

}

double MaternMatrixOperator::f(const Eigen::VectorXd & theta_in)
{
  double v = 0;
  //if(kappa.minCoeff() < 0)
  //  v = std::numeric_limits<double>::infinity();
  return v;
}

void MaternMatrixOperator::show_results()
{
  Rcpp::Rcout << ", kappa2 = " << kpv.exp().transpose() <<"\n";
}
Eigen::VectorXd MaternMatrixOperator::get_theta()
{
  return theta;
}
Rcpp::List MaternMatrixOperator::output_list()
{
  Rcpp::List List;
  if(kpv.size() > 0)
    List["kappa"]     = kpv;
  
  if(G.size() > 0)
    List["G"]         = G;
  if(C.size() > 0)
    List["C"]         = C;
  
  if(Bkp.size() > 0)
    List["B.kappa"]   = Bkp;
    
  List["kappa_transformed"] = 1;
  
    List["tau"]     = tau;

  return(List);
}

Eigen::VectorXd MaternMatrixOperator::kappa_vec(Eigen::MatrixXd & B, Eigen::VectorXd & beta)
{
  Eigen::VectorXd temp =  B*beta; 
  temp.array() = temp.array().exp();
  return(temp);
}


Eigen::MatrixXd MaternMatrixOperator::dkappa_mat(Eigen::MatrixXd & B, Eigen::VectorXd & beta) 
{
  
  Eigen::VectorXd G = kappa_vec(B, beta); 
  Eigen::MatrixXd temp = G.asDiagonal() * B ;
  return(temp);
}

void MaternMatrixOperator::gradient( const Eigen::VectorXd & X, const Eigen::VectorXd & iV)
{
	counter++;
  Eigen::VectorXd vtmp = Q * X;
  
  double xtQx =  vtmp.dot(iV.asDiagonal() * vtmp); 
  dtau  	  +=  0.5 * d / tau;
  dtau        -=  0.5 * xtQx;
  ddtau       -=  0.5 * d / pow(tau, 2); 
  dkappa[0] += trace0;
  func_dkappa("matern",
             dkappa, 
             X,
             iV, 
             *this,
             log(tau));
  
}
void MaternMatrixOperator::step_theta(const double stepsize)
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
   	if(kappa.size() > 0){
    	dkappa[0]  /= (-counter * d);
        kappa[0] -= dkappa[0];
    }
   	this->vec_to_theta(kappa);
    dkappa.setZero(Bkp.rows());
    counter = 0;
    trace0 = this->trace(0);
}


