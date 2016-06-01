#include <Rcpp.h>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include "Qmatrix.h"
#include "operator_helper.h"
#include "solver.h"
#include "GIG.h"
#include "MixedEffect.h"
using namespace Rcpp;




double estDigamma(double x)
{
  return(R::digamma(x));
}

void sampleX( Eigen::VectorXd & X, 
              Eigen::VectorXd & Z,
              const Eigen::VectorXd & V,
              const Eigen::VectorXd & Y,
              const Eigen::SparseMatrix<double,0,int> & Q,
              const Eigen::SparseMatrix<double,0,int> & A,
              const double sigma2,
              cholesky_solver       & solver)
{
  
  Eigen::SparseMatrix<double,0,int> Qi = Q + (A.transpose()*A)/ sigma2;
  solver.compute(Qi);
  Eigen::VectorXd b = A.transpose()*Y/ sigma2;
  X = solver.rMVN(b, Z);
}


// [[Rcpp::export]]
List estimateLongGH_cpp(List obs_list, 
                        List operator_list, 
                        List theta_list, 
                        List mixed_list,
                        double stepsize, 
                        int Niter, 
                        int nsim, 
                        int burnin,
                        int noise, 
                        int commonsigma,
                        int silent,
                        List V_list,
                        Eigen::VectorXd  U,
                        int Nlong)
{
  int GAL = 0;
  if(noise == 1)
    GAL = 1;
  
  //count number of patients  
  int nrep = obs_list.length();
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd> (theta_list["kappa"]);
  double tau = theta_list["tau"];
  Eigen::VectorXd  beta = Rcpp::as<Eigen::VectorXd>( theta_list["beta"]);
  double sigma = theta_list["sigma"];  
  double asigma, bsigma; 
  
  if(commonsigma == 0 ){
    asigma = theta_list["asigma"];  
    bsigma  = theta_list["bsigma"];
  }
  
  double lambda = theta_list["lambda"];
  double mu = theta_list["mu"];
  //Define operator
  std::string type_operator = Rcpp::as<std::string>(operator_list["type"]); 
  
  
  // init mixed effect
  NormalMixedEffect mixobj;
  
  int    usingMixedEff = Rcpp::as<double>(mixed_list["on"]);
  
  if(usingMixedEff)
    mixobj.initFromList(mixed_list);
  Qmatrix* Kobj;   
  //Kobj = new MaternMatrixOperator;
  operator_select(type_operator, &Kobj);
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  Kobj->vec_to_theta( kappa);
  Eigen::VectorXd h = Rcpp::as<Eigen::VectorXd>(operator_list["h"]);
  
  //Prior solver
  cholesky_solver Qsolver;
  Qsolver.init(Kobj->d, 0, 0, 0);
  Qsolver.analyze(Kobj->Q);
  Qsolver.compute(Kobj->Q);
  
  //Create solvers for each patient
  std::vector<  cholesky_solver >  Solver(nrep);
  int i = 0;
  Eigen::SparseMatrix<double,0,int> Q;
  std::vector< Eigen::SparseMatrix<double,0,int> > As(nrep);
  std::vector< Eigen::VectorXd > Ys(nrep);
  std::vector< Eigen::MatrixXd > Bs(nrep);
  std::vector< Eigen::VectorXd > Xs(nrep);
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    As[i] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    Ys[i] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    Bs[i] = Rcpp::as<Eigen::MatrixXd>(obs_tmp["B"]);
    Xs[i].resize( Kobj->d );
    Solver[i].init(Kobj->d, 0, 0, 0);
    Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Q = Q  * Kobj->Q;
    Q = Q + As[i].transpose()*As[i];
    Solver[i].analyze(Q);
    Solver[i].compute(Q);
    i++;
  }
  
  std::vector< Eigen::VectorXd > Vs(nrep);
  for( i = 0; i < V_list.length(); i++ ){ 
    if(noise>=0){
      Vs[i] = Rcpp::as<Eigen::VectorXd>(V_list[i]);
    }else{
      Vs[i] = h;
    }
  }
      
  /*
  Simulation objects
  */
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  gig rgig;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  
  Eigen::VectorXd  z;
  z.setZero(Kobj->d);
  
  Eigen::VectorXd b, Ysim;
  b.setZero(Kobj->d);
  
  
  Eigen::VectorXd  sigma_v;
  sigma_v.setZero(nrep);
  //Vector for random longitudinal sampleing
  std::vector<int> longInd;
  for (i=0; i<nrep; i++) longInd.push_back(i);
  
  std::default_random_engine gammagenerator;
    
  Eigen::SparseMatrix<double,0,int> Qi;
  Eigen::VectorXd dbeta, dkappa;
  Eigen::MatrixXd beta_vec(Niter, beta.size()), kappa_vec(Niter, kappa.size());
  Eigen::VectorXd sigma_vec(Niter), tau_vec(Niter), lambda_vec(Niter), sigma_r_vec(Niter);
  Eigen::VectorXd asigma_vec(Niter),bsigma_vec(Niter);
  dbeta.resize( beta.size());
  dbeta.setZero(beta.size());
  dkappa.resize( kappa.size());
  dkappa.setZero( kappa.size());
  double dtau, dsigma, dlambda; 
  double d2sigma;
  double dasigma,dbsigma, d2asigma,d2bsigma; 
  double alphai,betai;
  
  
  Eigen::MatrixXd B;
  Eigen::SparseMatrix<double,0,int> A;
  
  double  d2tau;
  Eigen::MatrixXd d2beta;
  d2beta.setZero(beta.size(),beta.size());
  
  for(int iter=0;iter< Niter + burnin;iter++){
    if(silent == 0){
      Rcpp::Rcout << "i = " << iter << ": ";
      if(commonsigma==1){
        Rcpp::Rcout << "sigma = " << exp(0.5 * sigma);  
      }else{
        Rcpp::Rcout << "a = " << asigma;  
        Rcpp::Rcout << "b = " << bsigma;  
      }
      
      Rcpp::Rcout << ", beta = " << beta.transpose();
      Rcpp::Rcout << ", tau = " << exp(tau);
      print_kappa(type_operator, kappa);
      if(usingMixedEff)
        Rcpp::Rcout << ", Sigma = \n" << mixobj.Sigma ;
      Rcpp::Rcout << "\n, lambda = " << exp(lambda) << "\n";
    }
    
    //Update operator
    Kobj->vec_to_theta( kappa);
    Qsolver.compute(Kobj->Q);
    
   
    func_dkappa0(type_operator, dkappa, nsim, nrep, *Kobj);
    dtau = 0.5 * nsim * ( Nlong * Kobj->d);
    
    double dgamma, d2gamma;
    if(commonsigma==1){
      dsigma = 0;  
      d2sigma = 0;
      //start with prior
      //double asigma = 0.001; 
      //double bsigma = 0.05;
      //dsigma = nsim*(-(asigma+1) + bsigma*exp(-sigma));
      //d2sigma = -nsim*bsigma*exp(-sigma);
    }else{
      dasigma = 0;
      dbsigma = 0;
      d2asigma = 0;
      d2bsigma = 0;
      dgamma   =  R::digamma(exp(asigma)); 
      d2gamma  =  R::trigamma(exp(asigma));
    }
    
    dlambda = 0;
    
    dbeta.setZero(beta.size());
    
    
    d2tau = 0;
    d2beta.setZero(beta.size(),beta.size());
    
    //sample Nlong values without replacement from 1:nrep
    std::random_shuffle(longInd.begin(), longInd.end());
    
    for( int ilong = 0; ilong < Nlong; ilong++ ) {
      i = longInd[ilong];
      A = As[i];
      B  = Bs[i];
      Eigen::VectorXd Y = Ys[i];
      
      Y -= B * beta;
      
      Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
      Eigen::VectorXd iV(Vs[i].size());
      iV.array() = Vs[i].array().inverse();
      Q = Q  * iV.asDiagonal();
      Q = exp(tau) * Q * Kobj->Q;     
      
      if(usingMixedEff)
        mixobj.remove_inter(i, Y);
        
      if(commonsigma==0){
        
        sigma = sigma_v[i];
        dasigma += nsim*exp(asigma)*(bsigma - dgamma);
        dbsigma += nsim*exp(asigma);
        d2asigma += nsim*exp(asigma)*(bsigma - dgamma) - nsim*d2gamma*exp(2*asigma);
        
      }
      for(int ii = 0; ii < nsim; ii ++){
        //Sample X|Y, V, sigma
        for(int j =0; j < Kobj->d; j++)
          z[j] =  normal(random_engine);
      
        sampleX(Xs[i], 
                z,
                Vs[i],
                Y,
                Q,
                A,
                exp(sigma),
                Solver[i]);
        

        // sample V| X
        if(noise >=0){
          Vs[i] =   sampleV_post(rgig,
                               h, 
                               Kobj->Q * Xs[i],
                               exp(- 0.5 * tau),
                               mu,
                               exp(lambda),
                              GAL);
          Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
          iV.array() = Vs[i].array().inverse();
          Q = Q  * iV.asDiagonal();
          Q =  exp(tau) * Q  * Kobj->Q;             
        }                    
        // sample U | V, X, sigma
       
        if(usingMixedEff)
        { 
          mixobj.add_inter(i, Y);
          mixobj.sampleU(i, Y - A * Xs[i], sigma);
          mixobj.gradient(i,  Y - A * Xs[i], sigma);
          mixobj.remove_inter(i, Y);
        }
        
        Ysim = Y - A * Xs[i];
        // sample sigma | U,V,X  
        if(commonsigma==0){
          alphai = exp(asigma) + Ysim.size()/2.0;
          betai = exp(bsigma) + 0.5* Ysim.array().square().sum();
          std::gamma_distribution<double> distribution(alphai,1.0/betai);
          sigma = -log(distribution(gammagenerator));  
        }
  
        //Compute gradients
        if(commonsigma==1){
          dsigma += -Ysim.size()/2.0 + exp(-sigma)*Ysim.array().square().sum()/2.0;  
        }else{
          dasigma -= exp(asigma)*sigma;
          dbsigma -= exp(bsigma)*exp(-sigma); 
          d2asigma -= exp(asigma)*sigma;
          d2bsigma -= exp(bsigma)*exp(-sigma); 
        }
        
        dbeta += exp(-sigma)* (B.transpose() * Ysim); 
        
        func_dkappa(type_operator,
                    dkappa, 
                    Xs[i],
                    iV, 
                    *Kobj,
                    tau);
        Eigen::VectorXd vtmp = Xs[i].transpose()*Q;
        double xtQx =  vtmp.dot(Xs[i]); 
        dtau -=  0.5 * xtQx;
        
        //Hessian
        d2sigma -= exp(-sigma)*Ysim.array().square().sum()/2.0;
        d2beta  -= exp(-sigma)*B.transpose()*B;
        d2tau   -= 0.5 * xtQx;
        
        // GIG gradients
        if(noise>=0){
          dlambda += dlambda_V(lambda,
                             Vs[i], 
                             h,
                             GAL);
                             
        }
      }
      sigma_v[i] = sigma;
    } 
    if( iter >= burnin){
    //take step in parameters
      if(commonsigma==1){
        sigma -= (dsigma/nsim)/d2sigma; 
        //sigma = log(0.05*0.05);
      }else{
        asigma -= (dasigma/nsim)/d2asigma;  
        bsigma -= (dbsigma/nsim)/d2bsigma;  
      }
      
      //sigma += dsigma * stepsize/nsim;
      beta  -= (d2beta.inverse()*dbeta)/nsim;
      tau   -= (1./d2tau) * dtau; 
      //tau = log(Nlong*Kobj->d)-log(-2.0*exp(-tau)*d2tau/nsim);
      tau   += (dtau/nsim) *  stepsize;
      kappa  += dkappa * stepsize/nsim;
      lambda += dlambda * stepsize/nsim;
      kappa_vec.row(iter - burnin) = kappa;
      beta_vec.row(iter - burnin)  = beta;
      sigma_vec(iter - burnin)     = sigma;
      asigma_vec(iter - burnin)     = asigma;
      bsigma_vec(iter - burnin)     = bsigma;
      tau_vec(iter - burnin)       = tau;
      lambda_vec(iter - burnin)    = lambda;
      
      if(usingMixedEff)
        mixobj.step_theta(.3);
    }
  }  
    
  Rcpp::List list;
  if(usingMixedEff){
    list["Sigma"] = mixobj.Sigma;
  }
  list["kappa"] = kappa.exp();
  if(commonsigma==1){
    list["sigma"] = exp(0.5*sigma);  
    sigma_vec *= 0.5;
    sigma_vec.array() = sigma_vec.array().exp();
    list["sigma_vec"] = sigma_vec;
  }else{
    list["asigma"] = exp(asigma);
    list["bsigma"]= exp(bsigma);
    sigma_v *= 0.5;
    sigma_v.array() = sigma_v.array().exp();
    list["sigma"] = sigma_v;  
    asigma_vec.array() = asigma_vec.array().exp();
    bsigma_vec.array() = bsigma_vec.array().exp();
    list["asigma_vec"] = asigma_vec;
    list["bsigma_vec"] = bsigma_vec;
  }
  
  list["beta"]  = beta;
  list["tau"]   = exp(tau);
  list["lambda"]   = exp(lambda);
  list["mu"]   = 0.;
  kappa_vec.array() = kappa_vec.array().exp();
  list["kappa_vec"] = kappa_vec;
  
  list["beta_vec"]  = beta_vec;
  tau_vec.array() = tau_vec.array().exp();
  list["tau_vec"]   = tau_vec;
  lambda_vec.array() = lambda_vec.array().exp();
  list["lambda_vec"]   = lambda_vec;
  list["X"] = Xs;
  list["V"] = Vs;
  list["U"] = mixobj.U;
  list["mixedeffect"] = mixobj.toList();
  return(list);
}

// [[Rcpp::export]]
List samplePosteriorGH(List obs_list, 
                        List operator_list, 
                        List theta_list,
                        List mixed_list,
                        List V_list,
                        int nsim, 
                        int noise,
                        int commonsigma)
{
  int GAL = 0;
  if(noise == 1)
    GAL = 1; 

  //count number of patients  
  int nrep = obs_list.length();
  Eigen::VectorXd kappa = Rcpp::as<Eigen::VectorXd> (theta_list["kappa"]);
  double tau = log(Rcpp::as< double> (theta_list["tau"]));
  Eigen::VectorXd  beta = Rcpp::as<Eigen::VectorXd>( theta_list["beta"]);
  double sigma = 0;
  Eigen::VectorXd sigma_v;
  if(commonsigma == 1){
    sigma  = 2 * log(Rcpp::as<double> (theta_list["sigma"]));  
  } else {
    sigma_v = Rcpp::as<Eigen::VectorXd>( theta_list["sigma"]);
    sigma_v.array() = 2 * sigma_v.array().log();
  }
  
  double asigma, bsigma;
  if(commonsigma == 0 ){
    asigma = theta_list["asigma"];  
    bsigma  = theta_list["bsigma"];
  }
  
  double lambda = log(Rcpp::as< double> (theta_list["lambda"]));
  double mu = theta_list["mu"];
  int    usingMixedEff = Rcpp::as<double>(mixed_list["on"]);
  NormalMixedEffect mixobj;
  if(usingMixedEff)
    mixobj.initFromList(mixed_list);
  
  //Define operator
  std::string type_operator = Rcpp::as<std::string>(operator_list["type"]); 
  Qmatrix* Kobj;   
  operator_select(type_operator, &Kobj);
  Kobj->initFromList(operator_list, List::create(Rcpp::Named("use.chol") = 1));
  Kobj->vec_to_theta( kappa);
  Eigen::VectorXd h = Rcpp::as<Eigen::VectorXd>(operator_list["h"]);
    
  //Prior solver
  cholesky_solver Qsolver;
  Qsolver.init(Kobj->d, 0, 0, 0);
  Qsolver.analyze(Kobj->Q);
  Qsolver.compute(Kobj->Q);
  
  
  //Create solvers for each patient
  std::vector<  cholesky_solver >  Solver(nrep);
  int i = 0;
  Eigen::SparseMatrix<double,0,int> Q;
  std::vector< Eigen::SparseMatrix<double,0,int> > As(nrep);
  std::vector< Eigen::VectorXd > Ys(nrep);
  std::vector< Eigen::MatrixXd > Bs(nrep);
  std::vector< Eigen::VectorXd > Xs(nrep);
  for( List::iterator it = obs_list.begin(); it != obs_list.end(); ++it ) {
    List obs_tmp = Rcpp::as<Rcpp::List>(*it);
    As[i] = Rcpp::as<Eigen::SparseMatrix<double,0,int> >(obs_tmp["A"]);
    Ys[i] = Rcpp::as<Eigen::VectorXd>(obs_tmp["Y"]);
    Bs[i] = Rcpp::as<Eigen::MatrixXd>(obs_tmp["B"]);
    Xs[i].resize( Kobj->d );
    Solver[i].init(Kobj->d, 0, 0, 0);
    Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Q = Q  * Kobj->Q;
    Q = Q + As[i].transpose()*As[i];
    Solver[i].analyze(Q);
    Solver[i].compute(Q);
    i++;
  }
  
  std::vector< Eigen::VectorXd > Vs(nrep);
  for( i = 0; i < V_list.length(); i++ ) 
      Vs[i] = Rcpp::as<Eigen::VectorXd>(V_list[i]);
  /*
  Simulation object
  */
  std::mt19937 random_engine;
  std::normal_distribution<double> normal;
  random_engine.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count()); 
  gig rgig;
  rgig.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  
  Eigen::VectorXd  z;
  z.setZero(Kobj->d);
  
  Eigen::VectorXd b, Ysim;
  b.setZero(Kobj->d);
  
  std::default_random_engine gammagenerator;
  
  Eigen::SparseMatrix<double,0,int> Qi;
 
  
  Eigen::MatrixXd B;
  Eigen::SparseMatrix<double,0,int> A;
  
 
  //Update operator
  Kobj->vec_to_theta( kappa);
  Qsolver.compute(Kobj->Q);
  std::vector< Eigen::MatrixXd > V_out(nrep), X_out(nrep);
  std::vector<Eigen::MatrixXd > U_out(nrep) ;  
  std::vector<Eigen::MatrixXd > sigma_out(nrep) ;  
    
  for( int i = 0; i < nrep; i++ ) {
    Rcpp::Rcout << "*";
    U_out[i].setZero(nsim, mixobj.U.rows());
    sigma_out[i].setZero(nsim,1);
    A = As[i];
    B  = Bs[i];
    Eigen::VectorXd Y = Ys[i];
    //compute mean
    Y -= B * beta;
      
    Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
    Eigen::VectorXd iV(Vs[i].size());
    iV.array() = Vs[i].array().inverse();
      
    Q = Q  * iV.asDiagonal();
    Q = exp(tau)*Q  * Kobj->Q;   
    V_out[i].setZero(nsim, Kobj->d);
    X_out[i].setZero(nsim, Kobj->d);
  
    if(usingMixedEff)
      mixobj.remove_inter(i, Y);
          
    if(commonsigma==0)
        sigma = sigma_v[i];
          
    for(int ii = 0; ii < nsim; ii ++){
      //Sample X|Y, V
        
      for(int j =0; j < Kobj->d; j++)
        z[j] =  normal(random_engine);
      
      sampleX(Xs[i], 
              z,
              Vs[i],
              Y,
              Q,
              A,
              exp(sigma),
              Solver[i]);
              
      X_out[i].row(ii) = Xs[i];     
      // sample V| X
      if(noise>=0){ 
        Vs[i] =   sampleV_post(rgig,
                             h, 
                             Kobj->Q*Xs[i],
                             exp(- 0.5 * tau),
                             mu,
                             exp(lambda),
                            GAL);
      }                      
      V_out[i].row(ii) = Vs[i]; 
      Q = Eigen::SparseMatrix<double,0,int>(Kobj->Q.transpose());
      iV.array() = Vs[i].array().inverse();
      Q = Q  * iV.asDiagonal();
      Q =  exp(tau)* Q  * Kobj->Q;   
      
      if(usingMixedEff)
      {
          mixobj.add_inter(i, Y);
          mixobj.sampleU(i, Y - A * Xs[i], sigma);
          mixobj.remove_inter(i, Y);
          U_out[i].row(ii) = mixobj.U.col(i);
          //U_out(i, ii) = U[i];
      }
      Ysim = Y - A * Xs[i];
      // sample sigma | U,V,X  
      if(commonsigma==0){
        double alphai = exp(asigma) + Ysim.size()/2.0;
        double betai = exp(bsigma) + 0.5* Ysim.array().square().sum();
        std::gamma_distribution<double> distribution(alphai,1.0/betai);
        sigma = -log(distribution(gammagenerator)); 
        sigma_out[i](ii,1) = exp(0.5*sigma);
      }
    }
  }
  Rcpp::Rcout << "\n";
  Rcpp::List list;
  list["X"] = X_out;
  list["V"]  = V_out;
  list["U"]  = U_out;
  list["sigma"]  = sigma_out;
  delete Kobj;
  return(list);
}