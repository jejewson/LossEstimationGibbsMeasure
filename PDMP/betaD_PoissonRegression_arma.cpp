// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;


//arma::vec dpoisRegression_vec_cpp_arma(const arma::vec& y, const arma::mat& X, const arma::vec& theta){
//  return exp(y%(X*theta) - repexp(X*theta) - lgamma(y+1));
//}

// [[Rcpp::export]]
arma::vec dpoisRegression_vec_cpp_arma(const arma::vec& y, const arma::rowvec& X, const arma::vec& theta){
  int n = y.n_elem;
  //return exp(y%repelem((X*theta), n, 1) - repelem(exp(X*theta), n, 1) - arma::lgamma(y+1));
  return exp(y*(X*theta) - repelem(exp(X*theta), n, 1) - arma::lgamma(y+1));
}// this works for a vector y but with only one predictor i.e. what you need for the intrgal term

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
  //out = arma::sum(arma::pow(dpoisRegression_vec_cpp_arma(z, X, theta), beta - 1.0) * X(i_0-1) % (z - repelem(exp(X*theta), m, 1)))/m - 
  //    arma::sum(arma::pow(dpoisRegression_cpp_arma(y, X, theta), beta - 1) * X(i_0-1) * (y - exp(X*theta)));
  out = (beta)*(arma::sum(arma::pow(dpoisRegression_vec_cpp_arma(z, X, theta), beta - 1.0) * X(i_0-1) % (z - repelem(exp(X*theta), m, 1)))/m - 
      arma::sum(arma::pow(dpoisRegression_cpp_arma(y, X, theta), beta - 1) * X(i_0-1) * (y - exp(X*theta))));
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
  /*
  out = y_max*arma::sum(abs(X.col(i_0 - 1))) + abs(xi(i_0-1) - m_0)/(pow(s_0,2.0)) + t/(pow(s_0,2.0)) + 
    exp(sum_abs_X_max*t)*sum(abs(X.col(i_0 - 1)) % exp(X*theta));
    //exp(arma::max(X*theta)*t)*sum(abs(X.col(i_0 - 1)) % exp(X*theta));
  */
  /*
  out = arma::sum(abs(X.col(i_0 - 1)) % (y + 3*sqrt(y)/sqrt(m))) + abs(xi(i_0-1) - m_0)/(pow(s_0,2.0)) + t/(pow(s_0,2.0)) + 
    exp(sum_abs_X_max*t)*sum(abs(X.col(i_0 - 1)) % exp(X*theta));
  */
  // PREVIOUSLY HAD exp(X*theta) AT THE END WHICH I THINK WAS WRONG
  // PREVIOSULY HAD (2*BETA - 1) BUT NOTBOTH EXP(XTHETA) TERMS CAN BE ACTIVE
  //out = arma::sum(abs(X.col(i_0 - 1)) % (y + 1.0/(beta - 1))) + abs(xi(i_0-1) - m_0)/(pow(s_0,2.0)) + t/(pow(s_0,2.0)) + 
  // exp(sum_abs_X_max*t)*(beta)/(beta - 1)*sum(abs(X.col(i_0 - 1)) % exp(X*xi));
  // CORRECT
  //out = (arma::sum(abs(X.col(i_0 - 1)) % (y + 1.0/(beta - 1))) + abs(xi(i_0-1) - m_0)/(pow(s_0,2.0)) + t/(pow(s_0,2.0)) + 
  // exp(sum_abs_X_max*t)*(beta)/(beta - 1)*sum(abs(X.col(i_0 - 1)) % exp(X*xi)))/10.0;// improving efficiency 
  // Scaling in beta from the paper
  out = (beta*arma::sum(abs(X.col(i_0 - 1)) % (y + 1.0/(beta - 1))) + abs(xi(i_0-1) - m_0)/(pow(s_0,2.0)) + t/(pow(s_0,2.0)) + 
    exp(sum_abs_X_max*t)*(beta)*(beta)/(beta - 1)*sum(abs(X.col(i_0 - 1)) % exp(X*xi)))/10.0;// improving efficiency 
  
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
    //a = y_max*arma::sum(abs(X.col(i))) + abs(xi(i) - m_0)/(pow(s_0,2.0));
    //a = arma::sum(abs(X.col(i)) % (y + 3*sqrt(y)/sqrt(m))) + abs(xi(i) - m_0)/(pow(s_0,2.0));
    
    //a = arma::sum(abs(X.col(i)) % (y + 1.0/(beta - 1))) + abs(xi(i) - m_0)/(pow(s_0,2.0));
    //b = pow(s_0, -2);
    // CORRECT
    //a = (arma::sum(abs(X.col(i)) % (y + 1.0/(beta - 1))) + abs(xi(i) - m_0)/(pow(s_0,2.0)))/10.0; // Improving efficiency
    //b = pow(s_0, -2)/10.0; // Improving efficiency
    // Scaling in beta from the paper
    a = (beta*arma::sum(abs(X.col(i)) % (y + 1.0/(beta - 1))) + abs(xi(i) - m_0)/(pow(s_0,2.0)))/10.0; // Improving efficiency
    b = pow(s_0, -2)/10.0; // Improving efficiency
    
    // Cinlars Method
    s1 = -log(R::runif(0, 1));
    tau_1 = (sqrt(pow(a, 2) + 2.0*b*s1) - a)/b;
    
    s2 = -log(R::runif(0, 1));
    // PREVIOUSLY HAD exp(X*theta) AT THE END WHICH I THINK WAS WRONG
    //tau_2 = log(1.0 + s2*sum_abs_X_max/((beta)/(beta - 1)*sum(abs(X.col(i)) % exp(X*xi))))/sum_abs_X_max;
    // CORRECT
    //tau_2 = log(1.0 + s2*sum_abs_X_max*10.0/((beta)/(beta - 1)*sum(abs(X.col(i)) % exp(X*xi))))/sum_abs_X_max; // Improving efficiency
    // Scaling in beta from the paper
    tau_2 = log(1.0 + s2*sum_abs_X_max*10.0/((beta)*beta/(beta - 1)*sum(abs(X.col(i)) % exp(X*xi))))/sum_abs_X_max; // Improving efficiency

    tau(i) = arma::min(arma::vec({tau_1, tau_2}));
  }
  
  return tau;
}


// [[Rcpp::export]]
double M2_betaD_poisRegression_cpp_arma(const int& i_0, const double& t, const arma::vec& theta, 
                                       const arma::vec& xi, const arma::vec& y,
                                       const arma::mat& X, const double& beta,
                                       const double& w, const double& m_0, 
                                       const double& s_0, const int& m, const double& sum_abs_X_max){
  //int p = xi.n_elem - 1;
  //int n = y.n_rows;
  double out = 0;
  out = arma::sum(abs(X.col(i_0 - 1)) % (y + 1.0/(beta - 1))) + abs(xi(i_0-1) - m_0)/(pow(s_0,2.0)) + t/(pow(s_0,2.0)) + 
   (beta)/(beta - 1)*sum(abs(X.col(i_0 - 1)) % exp(X*theta*t) % exp(X*xi));
  
  return out;
  
}

// [[Rcpp::export]]
arma::vec sim_M2_betaD_poisRegression_cpp_arma(const arma::vec& theta, 
                                              const arma::vec& xi, const arma::vec& y,
                                              const arma::mat& X, const double& beta,
                                              const double& w, const double& m_0, 
                                              const double& s_0, const int& m, const double& sum_abs_X_max){
  
  
  int n = y.n_elem;
  int p = xi.n_elem;
  double a = 0;
  double b = 0;
  double s1 = 0;
  //double s2 = 0;
  arma::vec s2 = zeros(n);
  arma::vec tau = zeros(p);
  double tau_1 = 0;
  arma::vec tau_2 = zeros(n);
  //int i;
  int j;
  for(j=0; j<p; j++){
    a = (arma::sum(abs(X.col(j)) % (y + 1.0/(beta - 1))) + abs(xi(j) - m_0)/(pow(s_0,2.0)));
    b = pow(s_0, -2);
    
    // Cinlars Method
    s1 = -log(R::runif(0, 1));
    tau_1 = (sqrt(pow(a, 2) + 2.0*b*s1) - a)/b;
    
    //for(i=0; i<n; i++){
    //  s2 = -log(R::runif(0, 1));
    //  tau_2(i) = log(1.0 + s2*X.row(i)*theta/((beta)/(beta - 1)*(abs(X(i,j)) * exp(X.row(i)*xi))))/(X.row(i)*theta);
    //}
    s2 = -log(Rcpp::runif(n, 0, 1));
    tau_2 = log(1.0 + s2 % (X*theta)/((beta)/(beta - 1)*(abs(X.col(j)) % exp(X*xi))))/(X*theta);
    
    
    tau(j) = arma::min(join_cols(arma::vec({tau_1}), tau_2));
  }
  
  return tau;
}


// you could actually also sub-sample this bit
// [[Rcpp::export]]
arma::mat z_generator_betaD_poisRegression_cpp_arma(const int& m, const arma::vec& theta,
                                                 const arma::mat& X){
  int n = X.n_rows;
  arma::mat z = zeros(n, m);
  int i;
  for(i=0; i<n; i++){
    arma::vec X_itheta = X.row(i)*theta;
    //NumericVector X_itheta2 = NumericVector(X_itheta.begin(), X_itheta.end());
    //NumericVector z_temp = Rcpp::rpois(m, X_itheta(0));
    //z.col(i) = as<arma::vec>(z_temp);
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
    //arma::vec tau_i = sim_M_betaD_poisRegression_cpp_arma(Theta_current, Xi_current, y, X, beta, w, m_0, s_0, y_max, sum_abs_X_max);
    arma::vec tau_i = sim_M_betaD_poisRegression_cpp_arma(Theta_current, Xi_current, y, X, beta, w, m_0, s_0, m, sum_abs_X_max);
    int i_0 = index_min(tau_i) + 1;// for R type indexing  
    double tau = tau_i(i_0 - 1);
    arma::mat z = z_generator_betaD_poisRegression_cpp_arma(m, Xi_current + tau*Theta_current, X);
    //double alpha = tilde_m_betaD_poisRegression_cpp_arma(i_0, tau, Theta_current, Xi_current, y, z, X, beta, w, m_0, s_0)/
    //  M_betaD_poisRegression_cpp_arma(i_0, tau, Theta_current, Xi_current, y, X, beta, w, m_0, s_0, y_max, sum_abs_X_max);
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


/*

/// For weight calibration
 
// [[Rcpp::export]]
arma::vec grad_MMD_RBF_grad_theta_Gaussian_regression_weight_calib_cpp_arma(const double& y, const arma::mat& u, const arma::vec& xi, const double& gamma, const arma::vec& X){
 //int n = y.n_elem;// here y is a double
 int m = u.n_rows;
 int p = xi.n_elem - 1;
 int i;
 int j;
 arma::vec z = zeros(m);
 z = F_generator_Gaussian_regression_cpp_arma(u, xi.subvec(0, p-1), exp(2*xi(p)), X);// might need to transpose the X on its way in
 arma::vec grad = zeros(p+1);
 for(i = 0; i<p; i++){
  grad(i) -= 2.0/m*arma::sum(-grad_z_RBF_kernel_cpp_arma(y, z, gamma)%grad_F_generator_Gaussian_regression_cpp_arma(i+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X));
 }

 grad(p) -= 2.0/m*arma::sum(-grad_z_RBF_kernel_cpp_arma(y, z, gamma)%grad_F_generator_Gaussian_regression_cpp_arma(p+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X));
 for(j = 0; j<m; j++){
  arma::vec k = regspace(0,m-1);
  arma::uvec ind = arma::find(k != j);
  grad(p) += arma::sum(grad_z_RBF_kernel_cpp_arma(z(j), z(ind), gamma)%(grad_F_generator_Gaussian_regression_cpp_arma(p+1, u.row(j), xi.subvec(0, p-1), exp(2*xi(p)), X)(0,0) - 
 grad_F_generator_Gaussian_regression_cpp_arma(p+1, u.rows(ind), xi.subvec(0, p-1), exp(2*xi(p)), X)))*1.0/(m*(m-1));
 }

 
 return grad;
 }


// [[Rcpp::export]]
arma::vec Hess_F_generator_Gaussian_regression_cpp_arma(const int& i_0, const arma::mat& u, const arma::vec& beta, const double& sigma2, const arma::vec& X){
  int m = u.n_rows;
  int p = beta.n_elem;
  arma::vec out = zeros(m);
  if(i_0 <= p){
    out = zeros(m);
  }
  if(i_0 == (p+1)){
    out = (sqrt(sigma2)*sqrt(-2.0*log(u.col(0)))%cos(2.0*M_PI*u.col(1)));
  }
  
  return out;
}


// [[Rcpp::export]]
arma::vec Hess_z_RBF_kernel_cpp_arma(const double& z, const arma::vec& x, const double& gamma){
  return (1.0/(sqrt(2.0*M_PI)*pow(gamma, 1.5)))*exp(-(z-x)%(z-x)/(2.0*gamma))%((z-x)%(z-x)/gamma - 1.0);
}


// [[Rcpp::export]]
arma::mat Hess_MMD_RBF_grad_theta_Gaussian_regression_cpp_arma(const double& y, const arma::mat& u, const arma::vec& xi, const double& gamma, const arma::vec& X){
  //int n = y.n_elem; //only ever evaluated for one observation
  int m = u.n_rows;
  int p = xi.n_elem - 1;
  int j = 0;
  int i = 0;
  arma::vec z = zeros(m);
  z = F_generator_Gaussian_regression_cpp_arma(u, xi.subvec(0, p-1), exp(2*xi(p)), X);// might need to transpose the X on its way in
  arma::mat Hess_out = zeros(p+1,p+1);
  arma::vec Hess_z_RBF_kernel_eval = Hess_z_RBF_kernel_cpp_arma(y, z, gamma);// DO WE NEED ANEGATIVE
  for(i = 0; i<p; i++){
    for(j = 0; j<=i; j++){
      Hess_out(i, j) = - 2.0/m*arma::sum(Hess_z_RBF_kernel_eval%(grad_F_generator_Gaussian_regression_cpp_arma(i+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X)%grad_F_generator_Gaussian_regression_cpp_arma(j+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X)));
      if(i != j){
        Hess_out(j, i) = Hess_out(i, j);
      }
    }
    Hess_out(i, p) = - 2.0/(m)*arma::sum(Hess_z_RBF_kernel_eval%(grad_F_generator_Gaussian_regression_cpp_arma(i+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X)%grad_F_generator_Gaussian_regression_cpp_arma(p+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X)));
    Hess_out(p, i) = Hess_out(i, p);
  }
  
  Hess_out(p, p) = - 2.0/m*arma::sum(-grad_z_RBF_kernel_cpp_arma(y, z, gamma)%Hess_F_generator_Gaussian_regression_cpp_arma(p+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X) + 
    Hess_z_RBF_kernel_eval%grad_F_generator_Gaussian_regression_cpp_arma(p+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X)%grad_F_generator_Gaussian_regression_cpp_arma(p+1, u, xi.subvec(0, p-1), exp(2*xi(p)), X));
  // minus inside as we sawped y and z order
  for(j = 0; j<m; j++){
    arma::vec k = regspace(0,m-1);
    arma::uvec ind = arma::find(k != j);
    Hess_out(p, p) += (arma::sum(grad_z_RBF_kernel_cpp_arma(z(j), z(ind), gamma)%(Hess_F_generator_Gaussian_regression_cpp_arma(p+1, u.row(j), xi.subvec(0, p-1), exp(2*xi(p)), X)(0,0) - 
      Hess_F_generator_Gaussian_regression_cpp_arma(p+1, u.rows(ind), xi.subvec(0, p-1), exp(2*xi(p)), X))) + 
      arma::sum(Hess_z_RBF_kernel_cpp_arma(z(j), z(ind), gamma)%(grad_F_generator_Gaussian_regression_cpp_arma(p+1, u.row(j), xi.subvec(0, p-1), exp(2*xi(p)), X)(0,0) - 
      grad_F_generator_Gaussian_regression_cpp_arma(p+1, u.rows(ind), xi.subvec(0, p-1), exp(2*xi(p)), X))%(grad_F_generator_Gaussian_regression_cpp_arma(p+1, u.row(j), xi.subvec(0, p-1), exp(2*xi(p)), X)(0,0) - 
      grad_F_generator_Gaussian_regression_cpp_arma(p+1, u.rows(ind), xi.subvec(0, p-1), exp(2*xi(p)), X)))    
    )*1.0/(m*(m-1));
  }
  
  return Hess_out;
}



// [[Rcpp::export]]
double weight_calib_MMD_RBF_Gaussian_regression_cpp_arma(const arma::vec& y, const arma::mat u, const arma::vec& xi, const double& gamma, const arma::mat X){
  int d = xi.n_elem;
  int n = y.n_elem;
  int i = 0;
  
  
  //Hess_data <- array(NA, dim = c(n, p, p))
  arma::mat mean_grad2_data = zeros(d, d);
  arma::mat mean_Hess_data = zeros(d, d);
  for(i = 0; i<n; i++){
    arma::vec grad_data = grad_MMD_RBF_grad_theta_Gaussian_regression_weight_calib_cpp_arma(y(i), u, xi, gamma, X.row(i).t()); 
    mean_grad2_data += grad_data * grad_data.t();
    arma::mat Hess_data =  Hess_MMD_RBF_grad_theta_Gaussian_regression_cpp_arma(y(i), u, xi, gamma, X.row(i).t());
    mean_Hess_data += Hess_data;
  }
  //hat_I_theta_data <- mean_grad2_data/n
  //hat_J_theta_data <- mean_Hess_data/n

    
  //w_data <- sum(diag((hat_J_theta_data%*%solve(hat_I_theta_data)%*%t(hat_J_theta_data))))/sum(diag(hat_J_theta_data))
      
  return arma::trace(mean_Hess_data*arma::inv_sympd(mean_grad2_data)*mean_Hess_data.t())/arma::trace(mean_Hess_data);
}


*/

