// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
arma::mat tilde_u_generator_Gaussian_Copula_cpp_arma(const int& N){
  arma::mat tilde_u = zeros(N, 2);
  tilde_u.col(0) = arma::randn(N, distr_param(0,1));
  tilde_u.col(1) = arma::randn(N, distr_param(0,1));
  
  return tilde_u;
}

// [[Rcpp::export]]
arma::mat F_generator_Gaussian_Copula_cpp_arma(const arma::mat& tilde_u, const double& rho){
  int N = tilde_u.n_rows;
  arma::mat U = zeros(N, 2);
  U.col(0) = 0.5*(1.0 + erf(tilde_u.col(0)/pow(2, 0.5)));
  U.col(1) = 0.5*(1.0 + erf((rho*tilde_u.col(0) + sqrt(1-pow(rho, 2))*tilde_u.col(1))/pow(2, 0.5)));
  
  return U;
}


// [[Rcpp::export]]
arma::vec grad_F_generator_Gaussian_Copula_cpp_arma(const int& i_0, const arma::mat& tilde_u, const double& rho){
  int m = tilde_u.n_rows;
  //arma::mat out = zeros(m,2);
  arma::vec out = zeros(m);
  if(i_0 == 1){
    out = (tilde_u.col(0) - rho/(sqrt(1-pow(rho, 2)))*tilde_u.col(1))/(pow(2.0*M_PI, 0.5))%exp(-arma::pow((rho*tilde_u.col(0) + pow(1 - pow(rho, 2), 0.5)*tilde_u.col(1)), 2.0)/2);
  }
  
  return out;
}

// [[Rcpp::export]]
double grad_MMD_biRBF_grad_theta_Gaussian_Copula_cpp_arma(const int& i_0, const arma::mat& y, const arma::mat& tilde_u, const arma::vec& xi, const double& gamma){
  int n = y.n_rows;
  int m = tilde_u.n_rows;
  int i;
  int j;
  double rho = 2.0/(1+exp(-xi(0))) - 1.0;
  arma::mat U = zeros(m, 2);
  arma::mat grad_gen = zeros(m);
  double out = 0;
  if(i_0 == 1){
    U = F_generator_Gaussian_Copula_cpp_arma(tilde_u, rho);
    grad_gen = grad_F_generator_Gaussian_Copula_cpp_arma(1, tilde_u, rho)*2.0*exp(-xi(0))/pow((1+exp(-xi(0))), 2.0);// chain rule for the reparam of rho
    for(i = 0; i<n; i++){
      out -= 2.0/m*arma::sum((((y(i, 1)-U.col(1))%grad_gen)/
        (2.0*M_PI*pow(gamma, 2)))%exp(-arma::sum(((y.row(i) - U.each_row())%(y.row(i) - U.each_row())), 1)/(2.0*gamma)));
    }
    for(j = 0; j<m; j++){
      out -= arma::sum((((U(j, 1)-U.col(1))%(grad_gen(j) - grad_gen))/
        (2.0*M_PI*pow(gamma, 2)))%exp(-arma::sum(((U.row(j) - U.each_row())%(U.row(j) - U.each_row())), 1)/(2.0*gamma)))*1.0*n/(m*(m-1));
    }
    
  }
  
  return out;
}



// [[Rcpp::export]]
double grad_log_Beta_prior_cpp(const double& xi, const double& a, const double& b){
  return -b + (a + b)/(1+exp(xi));
}

// [[Rcpp::export]]
double grad_log_prior_biGaussian_copula_cpp_arma(const int& i_0, const arma::vec& xi, const double& a, const double& b){
  double out = 0;
  if(i_0 == 1){
    out = grad_log_Beta_prior_cpp(xi(i_0-1), a, b);
  }
  return out;
}


// [[Rcpp::export]]
double tilde_m_MMD_RBF_biGaussian_copula_cpp_arma(const int& i_0, const double& t, 
                                                  const arma::vec& theta_curr, 
                                                  const arma::vec& xi_curr, 
                                                  const arma::mat& y, const arma::mat& tilde_u,
                                                  const double& gamma, const double& w, const double& a, const double& b){

  double out = 0;
  double temp = -theta_curr(i_0-1)*grad_log_prior_biGaussian_copula_cpp_arma(i_0, xi_curr + theta_curr*t, a, b)
    + theta_curr(i_0-1)*w*grad_MMD_biRBF_grad_theta_Gaussian_Copula_cpp_arma(i_0, y, tilde_u, xi_curr + theta_curr*t, gamma);
  if(temp > 0){
    out = temp;
  }
  return out;
}

// [[Rcpp::export]]
double M_MMD_RBF_biGaussian_copula_cpp_arma(const int& i_0, const double& t, const arma::vec& theta, 
                                            const arma::vec& xi, const arma::mat& y,
                                            const double& gamma, const double& R, 
                                            const double& w, const double& a, const double& b){
  //int p = xi.n_elem - 1;
  int n = y.n_rows;
  double out = 0;
  if(i_0 == 1){
    out = max(a, b) + 4.0*R*n*w/(pow((2.0*M_PI*gamma),1.5))*exp(-0.5)*3.0;
  }
  return out;
  
}

// [[Rcpp::export]]
arma::vec sim_M_MMD_RBF_biGaussian_copula_cpp_arma(const arma::vec& theta, const arma::vec& xi, 
                                                   const arma::mat& y, const double& gamma, 
                                                   const double& R, const double& w, const double& a, const double& b){
  
  int n = y.n_rows;
  int p = xi.n_elem;
  double s_1 = 0;
  //double s_2 = 0;
  double tau_1 = 0;
  //double tau_2 = 0;
  arma::vec tau = zeros(1);
  int i;
  for(i=0; i<p; i++){
    
    // Cinlars Method
    s_1 = -log(R::runif(0, 1));
    tau_1 = s_1/(max(a, b) + 12.0*R*n*w/(pow((2.0*M_PI*gamma),1.5))*exp(-0.5));
    
    tau(i) = tau_1;
    
  }
  
  return tau;
}


// [[Rcpp::export]]
Rcpp::List ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(const double& T_end, const arma::vec& xi_0, 
                                                     const arma::vec& theta_0, const arma::mat& u_hat, const int& m, 
                                                     const double& gamma, const double& R, 
                                                     const double& w, const double& a, const double& b, const int& N_skeleton){
  // T_end is the length of time to ru this for 
  // xi_0 \in R is the initial value of the parameter
  // theta_0 \in \{-1, 1\} is the initial velocity
  
  int d;
  d = xi_0.n_elem;
  double T_current = 0.0;
  arma::vec skeleton_T(N_skeleton);
  arma::mat skeleton_Xi(N_skeleton, d);
  arma::mat skeleton_Theta(N_skeleton, d);
  double max_alpha = 0.0;
  double mean_alpha = 0.0;
  int k_alpha = 0;
  int k_skeleton = 0;
  int print_count = 0;
  
  arma::vec Xi_current = xi_0;
  arma::vec Theta_current = theta_0;
  
  skeleton_T(k_skeleton) = T_current;
  skeleton_Xi.row(k_skeleton) = Xi_current.t();
  skeleton_Theta.row(k_skeleton) = Theta_current.t();
  
  
  
  while(T_current < T_end){
    arma::vec tau_i = sim_M_MMD_RBF_biGaussian_copula_cpp_arma(Theta_current, Xi_current, u_hat, gamma, R, w, a, b);
    int i_0 = index_min(tau_i) + 1;// for R type indexing  
    double tau = tau_i(i_0 - 1);
    arma::mat u = tilde_u_generator_Gaussian_Copula_cpp_arma(m);
    double alpha = tilde_m_MMD_RBF_biGaussian_copula_cpp_arma(i_0, tau, Theta_current, Xi_current,  u_hat, u, gamma, w, a, b)/
      M_MMD_RBF_biGaussian_copula_cpp_arma(i_0, tau, Theta_current, Xi_current, u_hat, gamma, R, w, a, b);
    k_alpha += 1;
    max_alpha = max({max_alpha, alpha});
    mean_alpha = (mean_alpha*(k_alpha - 1) + alpha)/k_alpha;
    
    T_current += tau;
    Xi_current += Theta_current*tau;
    if(alpha > R::runif(0, 1)){
      double theta_new = -Theta_current(i_0-1);
      Theta_current(i_0-1) = theta_new;
      k_skeleton += 1;
      skeleton_T(k_skeleton) = T_current;
      skeleton_Xi.row(k_skeleton) = Xi_current.t();
      skeleton_Theta.row(k_skeleton) = Theta_current.t();
    }
    if(T_current > (print_count+1)*(T_end/10.0)){
      Rprintf("Time >  %f \n", T_current);
      // Rprintf("\r Time >  %f", T_current);
      print_count += 1;
    }
  }
  
  
  
  List out = List::create(Named("skeleton_T") = skeleton_T(span(0,k_skeleton)),
                          Named("skeleton_Xi") = skeleton_Xi.rows(span(0,k_skeleton)), 
                          Named("skeleton_Theta") = skeleton_Theta.rows(span(0,k_skeleton)),
                          Named("max_alpha") = max_alpha,
                          Named("mean_alpha") = mean_alpha);  
  
  
  return out;
  
}

// [[Rcpp::export]]
arma::mat skeleton_to_sample_cpp_arma(const arma::vec& skeleton_T, const arma::mat& skeleton_Xi, const arma::mat& skeleton_Theta, const int& N){
  int d = skeleton_Xi.n_cols;
  int k = skeleton_T.n_elem;
  double T_end = skeleton_T(k-1);
  arma::mat samples = zeros(N, d);
  arma::vec samplingTimes =  linspace(1, N, N)*T_end*1.0/(N);
  
  int i;
  
  for(i = 0; i<N; i++){
    
    int k_i = index_min(abs(skeleton_T - samplingTimes(i)));
    if(skeleton_T(k_i) <= samplingTimes(i)){
      samples.row(i) = skeleton_Xi.row(k_i) + skeleton_Theta.row(k_i)*(samplingTimes(i) - skeleton_T(k_i));
    } else{
      samples.row(i) = skeleton_Xi.row(k_i-1) + skeleton_Theta.row(k_i-1)*(samplingTimes(i) - skeleton_T(k_i-1));
    }
    
  }
  
  return samples;
}


