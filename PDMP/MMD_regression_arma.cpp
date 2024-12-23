// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double MMD_RBF_vec_cpp_arma(const arma::vec& y, const arma::vec& u, const double& gamma){
  int n = y.n_elem;
  int m = u.n_elem;
  double out1 = 0.0;
  double out2 = 0.0;
  int i;
  int k;
  for(i = 0; i<n; i++){
    out1 += - arma::sum(arma::normpdf(y(i) - u, 0.0, gamma).col(0));
  }
  out1 = out1*2.0/(m);
  for(k = 0; k<m; k++){
    arma::vec l = regspace(0,m-1);
    arma::uvec ind = arma::find(l != k);
    out2 += arma::sum(arma::normpdf(u(k) - u(ind), 0.0, gamma).col(0));
    }
  out2 = out2*1.0/(m*(m-1.0));
  return out1 + n*out2;
}

// [[Rcpp::export]]
arma::mat u_generator_Gaussian_regression_cpp_arma(const int& N){
  arma::mat u = zeros(N, 2);
  u.col(0) = arma::randu(N, distr_param(0,1));
  u.col(1) = arma::randu(N, distr_param(0,1));
  
  return u;
}

// [[Rcpp::export]]
arma::vec F_generator_Gaussian_regression_cpp_arma(const arma::mat& u, const arma::vec& beta, const double& sigma2, const arma::vec& X){
  arma::vec mu = (X.t()*beta);
  return mu(0) + (sqrt(sigma2)*sqrt(-2*log(u.col(0)))%cos(2*M_PI*u.col(1)));
}

// [[Rcpp::export]]
arma::vec grad_z_RBF_kernel_cpp_arma(const double& z, const arma::vec& x, const double& gamma){
  return (-(z-x)/(sqrt(2.0*M_PI)*pow(gamma, 1.5)))%exp(-(z-x)%(z-x)/(2.0*gamma));
}

// [[Rcpp::export]]
arma::vec grad_F_generator_Gaussian_regression_cpp_arma(const int& i_0, const arma::mat& u, const arma::vec& beta, const double& sigma2, const arma::vec& X){
  int p = beta.n_elem;
  int m = u.n_rows;
  arma::vec out = zeros(m);
  if(i_0 <= p){
    out = rep(X(i_0-1), m);
  }
  if(i_0 == (p+1)){
    out = (sqrt(sigma2)*sqrt(-2.0*log(u.col(0)))%cos(2.0*M_PI*u.col(1)));
  }
  
  return out;
}



// [[Rcpp::export]]
double grad_MMD_RBF_grad_theta_Gaussian_regression_cpp_arma(const int& i_0, const arma::vec& y, const arma::mat& u, const arma::vec& xi, const double& gamma, const arma::mat& X){
  int n = y.n_elem;
  int m = u.n_rows;
  int p = xi.n_elem - 1;
  int i;
  int j;
  arma::vec z = zeros(m);
  double out = 0;
  if(i_0 <= p){
    for(i = 0; i<n; i++){
      z = F_generator_Gaussian_regression_cpp_arma(u, xi.subvec(0, p-1), exp(2*xi(p)), X.row(i).t());
      out -= 2.0/m*arma::sum(-grad_z_RBF_kernel_cpp_arma(y(i), z, gamma)%grad_F_generator_Gaussian_regression_cpp_arma(i_0, u, xi.subvec(0, p-1), exp(2*xi(p)), X.row(i).t()));
    }
  }
  if(i_0 == (p+1)){
    for(i = 0; i<n; i++){
      z = F_generator_Gaussian_regression_cpp_arma(u, xi.subvec(0, p-1), exp(2*xi(p)), X.row(i).t());
      out -= 2.0/m*arma::sum(-grad_z_RBF_kernel_cpp_arma(y(i), z, gamma)%grad_F_generator_Gaussian_regression_cpp_arma(i_0, u, xi.subvec(0, p-1), exp(2*xi(p)), X.row(i).t()));
    }
    z = F_generator_Gaussian_regression_cpp_arma(u, xi.subvec(0, p-1), exp(2*xi(p)), zeros(p));
    for(j = 0; j<m; j++){
      arma::vec k = regspace(0,m-1);
      arma::uvec ind = arma::find(k != j);
      out += arma::sum(grad_z_RBF_kernel_cpp_arma(z(j), z(ind), gamma)%(grad_F_generator_Gaussian_regression_cpp_arma(i_0, u.row(j), xi.subvec(0, p-1), exp(2*xi(p)), zeros(p))(0,0) - 
        grad_F_generator_Gaussian_regression_cpp_arma(i_0, u.rows(ind), xi.subvec(0, p-1), exp(2*xi(p)), zeros(p))))*1.0*n/(m*(m-1));
    }
    
  }
  
  return out;
}



// [[Rcpp::export]]
double grad_log_Gaussian_prior_cpp(const double& xi, const double& m_0, const double& s_0){
  return -(xi - m_0)/pow(s_0, 2);
}

// [[Rcpp::export]]
double grad_log_IG_prior_cpp(const double& xi, const double& a_0, const double& b_0){
  return -2*a_0 + 2*b_0*exp(-2*xi);
}

// [[Rcpp::export]]
double grad_log_prior_Gaussian_regression_cpp_arma(const int& i_0, const arma::vec& xi,
                                                 const double& m_0, const double& s_0, 
                                                 const double& a_0, const double& b_0){
  int p = xi.n_elem - 1;
  double out = 0;
  if(i_0 <= p){
    out = grad_log_Gaussian_prior_cpp(xi(i_0-1), m_0, s_0);
  }
  if(i_0 == (p+1)){
    out = grad_log_IG_prior_cpp(xi(p), a_0, b_0);
  }
  return out;
}

// [[Rcpp::export]]
double tilde_m_MMD_RBF_Gaussian_regression_cpp_arma(const int& i_0, const double& t, 
                                                  const arma::vec& theta_curr, 
                                                  const arma::vec& xi_curr, 
                                                  const arma::vec& y, 
                                                  const arma::mat& X, const arma::mat& u,
                                                  const double& gamma, 
                                                  const double& m_0, const double& s_0, 
                                                  const double& a_0, const double& b_0, 
                                                  const double& w){
  
  double out = 0;
  double temp = -theta_curr(i_0-1)*grad_log_prior_Gaussian_regression_cpp_arma(i_0, xi_curr + theta_curr*t, m_0, s_0, a_0, b_0) + 
    theta_curr(i_0-1)*w*grad_MMD_RBF_grad_theta_Gaussian_regression_cpp_arma(i_0, y, u, xi_curr + theta_curr*t, gamma, X);
  if(temp > 0){
    out = temp;
  }
  return out;
}



// [[Rcpp::export]]
double M_MMD_RBF_Gaussian_regression_cpp_arma(const int& i_0, const double& t, const arma::vec& theta, 
                                            const arma::vec& xi, const arma::vec& y, 
                                            const arma::mat& X,
                                            const double& gamma, const double& m_0, 
                                            const double& s_0, const double& a_0, 
                                            const double& b_0, const double& R, 
                                            const double& w){
  int p = xi.n_elem - 1;
  int n = y.n_elem;
  double out = 0;
  if(i_0 <= p){
    out = abs(xi(i_0-1) - m_0)/pow(s_0, 2) + t/(pow(s_0, 2)) + 2.0*arma::sum(arma::abs(X.col(i_0-1)))*w/(sqrt(2.0*M_PI)*gamma)*exp(-0.5);
  }
  if(i_0 == (p+1)){
    out = a_0 + 2.0*b_0*exp(-2.0*(xi(i_0-1) + theta(i_0-1)*t)) + 2.0*n*w/(sqrt(2.0*M_PI)*gamma)*exp(-0.5)*exp(xi(i_0-1) + theta(i_0-1)*t)*R;
  }
  
  return out;
  
}



// [[Rcpp::export]]
arma::vec sim_M_MMD_RBF_Gaussian_regression_cpp_arma(const arma::vec& theta, const arma::vec& xi, 
                                                   const arma::vec& y, const arma::mat& X, const double& gamma, 
                                                   const double& m_0, const double& s_0, 
                                                   const double& a_0, const double& b_0, 
                                                   const double& R, const double& w){
  
  int n = y.n_elem;
  int p = xi.n_elem - 1;
  double a = 0;
  double b = 0;
  double s1 = 0;
  double s21 = 0;
  double s22 = 0;
  double s23 = 0;
  double tau_1 = 0;
  double tau_22 = 0;
  double tau_23 = 0;
  arma::vec tau = zeros(p+1);
  int i;
  for(i=0; i<p; i++){
    // Location
    a = abs(xi(i) - m_0)/pow(s_0, 2) + 2.0*arma::sum(arma::abs(X.col(i)))*w/(sqrt(2.0*M_PI)*gamma)*exp(-0.5);
    b = pow(s_0, -2);
  
    // Cinlars Method
    s1 = -log(R::runif(0, 1));
    // Case 1) - this case will always be activates 
    if((a > 0) & (b > 0)){
      tau_1 = (sqrt(pow(a, 2) + 2.0*b*s1) - a)/b;
    } 
    // Case 2)
    if((a < 0) & (b > 0)){
      tau_1 = (sqrt(2.0*b*s1) - a)/b;
    } 
    tau(i) = tau_1;
  }
  
  // Scale
  s21 = -log(R::runif(0, 1));
  s22 = -log(R::runif(0, 1));
  
  if((theta(p) > 0) & (b_0*exp(-2.0*xi(p)) <= s22)){
    tau_22 = R_PosInf;
  } else{
    tau_22 = -1.0/(2.0*theta(p))*log(1.0 - s22*theta(p)*exp(2.0*xi(p))/(b_0));
  }
  
  s23 = -log(R::runif(0, 1));
  
  if((theta(p) < 0) & ((2.0*n*w/(sqrt(2.0*M_PI)*gamma)*exp(-0.5)*R)*exp(xi(p)) <= s23)){
    tau_23 = R_PosInf;
  } else{
    tau_23 = 1.0/(theta(p))*log(1.0 + s23*theta(p)*exp(-xi(p))/(2.0*n*w/(sqrt(2.0*M_PI)*gamma)*exp(-0.5)*R));
  }
  tau(p) = arma::min(arma::vec({s21/a_0, tau_22, tau_23}));
  
  return tau;
}



// [[Rcpp::export]]
Rcpp::List ZigZag_MMD_RBF_Gaussian_regression_cpp_arma(const double& T_end, const arma::vec& xi_0, 
                                                     const arma::vec& theta_0, const arma::vec& y, 
                                                     const arma::mat& X, const int& m, 
                                                     const double& gamma, const double& m_0, 
                                                     const double& s_0, const double& a_0, 
                                                     const double& b_0, const double& R, 
                                                     const double& w, const int& N_skeleton){
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
    arma::vec tau_i = sim_M_MMD_RBF_Gaussian_regression_cpp_arma(Theta_current, Xi_current, y, X, gamma, m_0, s_0, a_0, b_0, R, w);
    int i_0 = index_min(tau_i) + 1;// for R type indexing  
    double tau = tau_i(i_0 - 1);
    arma::mat u = u_generator_Gaussian_regression_cpp_arma(m);
    double alpha = tilde_m_MMD_RBF_Gaussian_regression_cpp_arma(i_0, tau, Theta_current, Xi_current,  y, X, u, gamma, m_0, s_0, a_0, b_0, w)/
      M_MMD_RBF_Gaussian_regression_cpp_arma(i_0, tau, Theta_current, Xi_current, y, X, gamma, m_0, s_0, a_0, b_0, R, w);
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

