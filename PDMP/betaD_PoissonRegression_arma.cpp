// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;



// [[Rcpp::export]]
arma::vec dpoisRegression_vec_cpp_arma(const arma::vec& y, const arma::rowvec& X, const arma::vec& theta){
  int n = y.n_elem;
  return exp(y*(X*theta) - repelem(exp(X*theta), n, 1) - arma::lgamma(y+1));
}

// [[Rcpp::export]]
arma::vec dpoisRegression_cpp_arma(const double& y, const arma::rowvec& X, const arma::vec& theta){
  return exp(y*(X*theta) - exp(X*theta) - lgamma(y+1));
}

// [[Rcpp::export]]
arma::mat grad_log_dpoisRegression_cpp_arma(const arma::vec& y, const arma::mat& X, const arma::vec& theta){
  int n = y.n_elem;
  int p = X.n_cols;
  arma::mat out = zeros(n, p);
  int j;
  for(j = 0; j<p; j++){
    out.col(j) = X.col(j) % (y - exp(X*theta));
  }
  return out;
}

// [[Rcpp::export]]
arma::mat grad_dpoisRegression_cpp_arma(const arma::vec& y, const arma::mat& X, const arma::vec& theta){
  int n = y.n_elem;
  int p = X.n_cols;
  arma::mat out = zeros(n, p);
  int j;
  for(j = 0; j<p; j++){
    out.col(j) = X.col(j) % (y - exp(X*theta)) % dpoisRegression_vec_cpp_arma(y, X, theta);
  }
  return out;
}


// [[Rcpp::export]]
double grad_betaD_loss_poisRegressionOneObs_cpp_arma(const int& i_0, const double& y, const arma::vec& z, const arma::rowvec& X, const arma::vec& theta, const double& beta){
  int m = z.n_elem;
  double out;
  out = arma::sum(arma::pow(dpoisRegression_vec_cpp_arma(z, X, theta), beta - 1.0) * X(i_0-1) % (z - repelem(exp(X*theta), m, 1)))/m - 
      arma::sum(arma::pow(dpoisRegression_cpp_arma(y, X, theta), beta - 1) * X(i_0-1) * (y - exp(X*theta)));
  return out;
}

// [[Rcpp::export]]
double grad_betaD_loss_poisRegression_cpp_arma(const int& i_0, const arma::vec& y, const arma::mat& z, const arma::mat& X, const arma::vec& theta, const double& beta){
  int n = y.n_elem;
  //int m = z.n_cols;
  //int p = X.n_cols;
  double out = 0;
  int i;
  for(i = 0; i<n; i++){
    out +=  grad_betaD_loss_poisRegressionOneObs_cpp_arma(i_0, y(i), z.row(i).t(), X.row(i), theta, beta);
  }
  return out;
}


// [[Rcpp::export]]
double grad_log_Gaussian_prior_cpp(const double& xi, const double& m_0, const double& s_0){
  return -(xi - m_0)/pow(s_0, 2);
}

// [[Rcpp::export]]
double grad_log_prior_poisRegression_cpp_arma(const int& i_0, const arma::vec& xi,
                                                   const double& m_0, const double& s_0){
  //int p = xi.n_elem;
  double out = 0;
  out = grad_log_Gaussian_prior_cpp(xi(i_0-1), m_0, s_0);
  return out;
}


// [[Rcpp::export]]
double tilde_m_betaD_poisRegression_cpp_arma(const int& i_0, const double& t, 
                                                  const arma::vec& theta_curr, 
                                                  const arma::vec& xi_curr, 
                                                  const arma::vec& y, const arma::mat& z,
                                                  const arma::mat& X, const double& beta,
                                                  const double& w, const double& m_0, const double& s_0){
  //int n = y.n_rows;
  double out = 0;
  double temp = -theta_curr(i_0-1)*grad_log_prior_poisRegression_cpp_arma(i_0, xi_curr + theta_curr*t, m_0, s_0)
    +theta_curr(i_0-1)*w*grad_betaD_loss_poisRegression_cpp_arma(i_0, y, z, X, xi_curr + theta_curr*t, beta);
  if(temp > 0){
    out = temp;
  }
  return out;
}



// [[Rcpp::export]]
double M_betaD_poisRegression_cpp_arma(const int& i_0, const double& t, const arma::vec& theta, 
                                            const arma::vec& xi, const arma::vec& y,
                                            const arma::mat& X, const double& beta,
                                            const double& w, const double& m_0, 
                                            const double& s_0, const int& m, const double& sum_abs_X_max){
  //int p = xi.n_elem - 1;
  //int n = y.n_rows;
  double out = 0;

  //out = arma::sum(abs(X.col(i_0 - 1)) % (y + 1.0/(beta - 1))) + abs(xi(i_0-1) - m_0)/(pow(s_0,2.0)) + t/(pow(s_0,2.0)) + 
  // exp(sum_abs_X_max*t)*(beta)/(beta - 1)*sum(abs(X.col(i_0 - 1)) % exp(X*xi));
  out = (arma::sum(abs(X.col(i_0 - 1)) % (y + 1.0/(beta - 1))) + abs(xi(i_0-1) - m_0)/(pow(s_0,2.0)) + t/(pow(s_0,2.0)) + 
    exp(sum_abs_X_max*t)*(beta)/(beta - 1)*sum(abs(X.col(i_0 - 1)) % exp(X*xi)))/10.0;// improving efficiency 
  
  return out;
  
}


// [[Rcpp::export]]
arma::vec sim_M_betaD_poisRegression_cpp_arma(const arma::vec& theta, 
                                              const arma::vec& xi, const arma::vec& y,
                                              const arma::mat& X, const double& beta,
                                              const double& w, const double& m_0, 
                                              const double& s_0, const int& m, const double& sum_abs_X_max){
                                              
  
  //int n = y.n_elem;
  int p = xi.n_elem;
  double a = 0;
  double b = 0;
  double s1 = 0;
  double s2 = 0;
  arma::vec tau = zeros(p);
  double tau_1 = 0;
  double tau_2 = 0;
  int i;
  //double sum_abs_X_max = arma::max(X*theta);
  for(i=0; i<p; i++){
    // Location
    //a = arma::sum(abs(X.col(i)) % (y + 1.0/(beta - 1))) + abs(xi(i) - m_0)/(pow(s_0,2.0));
    //b = pow(s_0, -2);
    a = (arma::sum(abs(X.col(i)) % (y + 1.0/(beta - 1))) + abs(xi(i) - m_0)/(pow(s_0,2.0)))/10.0; // Improving efficiency
    b = pow(s_0, -2)/10.0; // Improving efficiency
    
    // Cinlars Method
    s1 = -log(R::runif(0, 1));
    tau_1 = (sqrt(pow(a, 2) + 2.0*b*s1) - a)/b;
    
    s2 = -log(R::runif(0, 1));
    //tau_2 = log(1.0 + s2*sum_abs_X_max/((beta)/(beta - 1)*sum(abs(X.col(i)) % exp(X*xi))))/sum_abs_X_max;
    tau_2 = log(1.0 + s2*sum_abs_X_max*10.0/((beta)/(beta - 1)*sum(abs(X.col(i)) % exp(X*xi))))/sum_abs_X_max; // Imprving efficiency
    

    tau(i) = arma::min(arma::vec({tau_1, tau_2}));
  }
  
  return tau;
}


// [[Rcpp::export]]
arma::mat z_generator_betaD_poisRegression_cpp_arma(const int& m, const arma::vec& theta,
                                                 const arma::mat& X){
  int n = X.n_rows;
  arma::mat z = zeros(n, m);
  int i;
  for(i=0; i<n; i++){
    arma::vec X_itheta = X.row(i)*theta;
    z.row(i) = as<arma::vec>(Rcpp::rpois(m, exp(X_itheta(0)))).t();
  }
  
  return z;
}



// [[Rcpp::export]]
Rcpp::List ZigZag_betaD_poisRegression_cpp_arma(const double& T_end, const arma::vec& xi_0, 
                                                     const arma::vec& theta_0, const arma::vec& y,
                                                     const arma::mat& X, const int& m, 
                                                     const double& beta, const double& sum_abs_X_max, 
                                                     const double& w, const double& m_0, const double& s_0, const int& N_skeleton){
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
    arma::vec tau_i = sim_M_betaD_poisRegression_cpp_arma(Theta_current, Xi_current, y, X, beta, w, m_0, s_0, m, sum_abs_X_max);
    int i_0 = index_min(tau_i) + 1;// for R type indexing  
    double tau = tau_i(i_0 - 1);
    arma::mat z = z_generator_betaD_poisRegression_cpp_arma(m, Xi_current + tau*Theta_current, X);
    double alpha = tilde_m_betaD_poisRegression_cpp_arma(i_0, tau, Theta_current, Xi_current, y, z, X, beta, w, m_0, s_0)/
      M_betaD_poisRegression_cpp_arma(i_0, tau, Theta_current, Xi_current, y, X, beta, w, m_0, s_0, m, sum_abs_X_max);
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


