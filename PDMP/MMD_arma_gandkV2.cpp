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

// Generator for g-and-k: uses standard normal random variables
// [[Rcpp::export]]
arma::vec u_generator_gandk_cpp_arma(const int& N){
  arma::vec u = arma::randn(N);
  return u;
}

// Transform from unconstrained ξ to constrained θ ∈ [0, 5]
// θ = 5 / (1 + exp(-ξ))
// [[Rcpp::export]]
double transform_to_bounded(const double& xi){
  return 5.0 / (1.0 + exp(-xi));
}

// Transform vector
// [[Rcpp::export]]
arma::vec transform_to_bounded_vec(const arma::vec& xi){
  return 5.0 / (1.0 + exp(-xi));
}

// g-and-k generator function
// Note: xi contains (xi_a, xi_b, xi_g, xi_k) - all on unconstrained scale
// [[Rcpp::export]]
arma::vec F_generator_gandk_cpp_arma(const arma::vec& zu, const double& xi_a, const double& xi_b, 
                                     const double& xi_g, const double& xi_k){
  int n = zu.n_elem;
  double c_val = 0.8;
  
  // Transform to constrained scale [0, 5]
  double a = transform_to_bounded(xi_a);
  double b = transform_to_bounded(xi_b);
  double g = transform_to_bounded(xi_g);
  double k = transform_to_bounded(xi_k);
  
  arma::vec exp_gz = exp(-g * zu);
  arma::vec t_term = (1.0 - exp_gz) / (1.0 + exp_gz);
  arma::vec inner_term = 1.0 + c_val * t_term;
  arma::vec zu_sq = zu % zu;
  arma::vec power_term = pow(1.0 + zu_sq, k);
  
  arma::vec samples = a + b * inner_term % power_term % zu;
  
  return samples;
}

// Gradient of g-and-k generator with respect to parameter i_0
// Chain rule: ∂F/∂ξ = ∂F/∂θ × ∂θ/∂ξ where ∂θ/∂ξ = θ(5-θ)/5
// [[Rcpp::export]]
arma::vec grad_F_generator_gandk_mat_cpp_arma(const int& i_0, const arma::vec& zu, 
                                              const double& xi_a, const double& xi_b, 
                                              const double& xi_g, const double& xi_k){
  int m = zu.n_elem;
  arma::vec out = zeros(m);
  double c_val = 0.8;
  
  // Transform to constrained scale
  double a = transform_to_bounded(xi_a);
  double b = transform_to_bounded(xi_b);
  double g = transform_to_bounded(xi_g);
  double k = transform_to_bounded(xi_k);
  
  arma::vec exp_gz = exp(-g * zu);
  arma::vec t_term = (1.0 - exp_gz) / (1.0 + exp_gz);
  arma::vec inner_term = 1.0 + c_val * t_term;
  arma::vec zu_sq = zu % zu;
  arma::vec power_term = pow(1.0 + zu_sq, k);
  
  if(i_0 == 1){  // Gradient with respect to xi_a
    // ∂F/∂a = 1, ∂a/∂ξ_a = a(5-a)/5
    double jacobian = a * (5.0 - a) / 5.0;
    out = ones(m) * jacobian;
  }
  if(i_0 == 2){  // Gradient with respect to xi_b
    // ∂F/∂b = inner_term * power_term * zu, ∂b/∂ξ_b = b(5-b)/5
    double jacobian = b * (5.0 - b) / 5.0;
    out = inner_term % power_term % zu * jacobian;
  }
  if(i_0 == 3){  // Gradient with respect to xi_g
    // ∂F/∂g = b * c * (∂t/∂g) * power_term * zu, ∂g/∂ξ_g = g(5-g)/5
    arma::vec dt_dg = zu % (1.0 - t_term % t_term) / 2.0;
    arma::vec dF_dg = b * c_val * dt_dg % power_term % zu;
    double jacobian = g * (5.0 - g) / 5.0;
    out = dF_dg * jacobian;
  }
  if(i_0 == 4){  // Gradient with respect to xi_k
    // ∂F/∂k = b * inner_term * zu * power_term * log(1+zu^2), ∂k/∂ξ_k = k(5-k)/5
    arma::vec log_term = log(1.0 + zu_sq);
    arma::vec dF_dk = b * inner_term % zu % power_term % log_term;
    double jacobian = k * (5.0 - k) / 5.0;
    out = dF_dk * jacobian;
  }
  
  return out;
}

// [[Rcpp::export]]
arma::vec grad_z_RBF_kernel_cpp_arma(const double& z, const arma::vec& x, const double& gamma){
  return (-(z-x)/(sqrt(2.0*M_PI)*pow(gamma, 1.5)))%exp(-(z-x)%(z-x)/(2.0*gamma));
}

// Gradient of MMD with respect to parameter i_0 for g-and-k
// [[Rcpp::export]]
double grad_MMD_RBF_grad_theta_gandk_cpp_arma(const int& i_0, const arma::vec& y, 
                                              const arma::vec& u, const arma::vec& xi, 
                                              const double& gamma){
  int n = y.n_elem;
  int m = u.n_elem;
  
  int j;
  arma::vec z = zeros(m);
  // xi contains (xi_a, xi_b, xi_g, xi_k) - all unconstrained
  z = F_generator_gandk_cpp_arma(u, xi(0), xi(1), xi(2), xi(3));
  double out = 0;
  
  // First term: gradient with respect to data
  arma::vec grad_all = grad_F_generator_gandk_mat_cpp_arma(i_0, u, xi(0), xi(1), xi(2), xi(3));
  
  for(j = 0; j<m; j++){
    arma::vec kernel_grad = grad_z_RBF_kernel_cpp_arma(z(j), y, gamma);
    out -= arma::sum(kernel_grad) * grad_all(j);
  }
  out = out*2.0*1.0/(n*m);
  
  // For all parameters, add the second term (since all affect the distribution)
  for(j = 0; j<m; j++){
    arma::vec k_vec = regspace(0, m-1);
    arma::uvec ind = arma::find(k_vec != j);
    
    if(ind.n_elem == 0){
      continue;
    }
    
    arma::vec z_ind = z(ind);
    arma::vec grad_ind = grad_all(ind);
    arma::vec kernel_grad_z = grad_z_RBF_kernel_cpp_arma(z(j), z_ind, gamma);
    
    out += arma::sum(kernel_grad_z % (grad_all(j) - grad_ind)) * 1.0/(m*(m-1));
  }
  
  return out;
}

// Gradient of log prior for uniform on [0, 5] with logit transformation
// For θ ~ Uniform(0, 5), with θ = 5/(1+exp(-ξ)), the gradient is: 1 - 2θ/5
// [[Rcpp::export]]
double grad_log_uniform_prior_logit_cpp(const double& xi){
  double theta = transform_to_bounded(xi);
  return 1.0 - 2.0 * theta / 5.0;
}

// Prior gradient for g-and-k parameters (all uniform on [0, 5])
// [[Rcpp::export]]
double grad_log_prior_gandk_cpp_arma(const int& i_0, const arma::vec& xi){
  return grad_log_uniform_prior_logit_cpp(xi(i_0 - 1));
}

// Thinning function for g-and-k
// [[Rcpp::export]]
double tilde_m_MMD_RBF_gandk_cpp_arma(const int& i_0, const double& t, 
                                      const arma::vec& theta_curr, 
                                      const arma::vec& xi_curr, 
                                      const arma::vec& y, const arma::vec& u,
                                      const double& gamma, const double& w){
  
  int n = y.n_elem;
  double out = 0;
  double temp = -theta_curr(i_0-1)*grad_log_prior_gandk_cpp_arma(i_0, xi_curr + theta_curr*t) + 
    theta_curr(i_0-1)*n*w*grad_MMD_RBF_grad_theta_gandk_cpp_arma(i_0, y, u, xi_curr + theta_curr*t, gamma);
  if(temp > 0){
    out = temp;
  }
  return out;
}

// Upper bound function for g-and-k
// [[Rcpp::export]]
double M_MMD_RBF_gandk_cpp_arma(const int& i_0, const double& t, const arma::vec& theta, 
                                const arma::vec& xi, const arma::vec& y, 
                                const double& gamma, const double& w){
  int n = y.n_elem;
  double kernel_bound = 2.0*n*w/(sqrt(2.0*M_PI)*gamma)*exp(-0.5);
  
  // For the prior gradient: |1 - 2θ/5| ≤ 1 for θ ∈ [0, 5]
  double prior_bound = 1.0;
  
  // For all parameters, use similar bounds
  // The gradient involves the Jacobian θ(5-θ)/5 which is maximized at θ=2.5
  // max(θ(5-θ)/5) = 2.5*2.5/5 = 1.25
  double jacobian_bound = 1.25;
  
  double out = abs(theta(i_0-1)) * prior_bound + kernel_bound * jacobian_bound;// * 10.0;
  
  return out;
}

// Simulation of event times for g-and-k
// [[Rcpp::export]]
arma::vec sim_M_MMD_RBF_gandk_cpp_arma(const arma::vec& theta, const arma::vec& xi, 
                                       const arma::vec& y, const double& gamma, 
                                       const double& w){
  
  int n = y.n_elem;
  int d = 4;
  arma::vec tau_vec = zeros(d);
  double kernel_bound = 2.0*n*w/(sqrt(2.0*M_PI)*gamma)*exp(-0.5);
  double prior_bound = 1.0;
  double jacobian_bound = 1.25;
  
  // All parameters have the same structure now
  // CHECK THESE AS V1 GOT IT WRONG!
  for(int i = 0; i < d; i++){
    double a_i = abs(theta(i)) * prior_bound + kernel_bound * jacobian_bound;// * 10.0;
    double s_i = -log(R::runif(0, 1));
    if(a_i > 0){
      tau_vec(i) = s_i / a_i;
    } else {
      tau_vec(i) = R_PosInf;
    }
  }
  
  return tau_vec;
}

// Main ZigZag sampler for g-and-k
// [[Rcpp::export]]
Rcpp::List ZigZag_MMD_RBF_gandk_cpp_arma(const double& T_end, const arma::vec& xi_0, 
                                         const arma::vec& theta_0, const arma::vec& y, 
                                         const int& m, const double& gamma, 
                                         const double& w, const int& N_skeleton){
  
  int d = 4; // g-and-k has 4 parameters
  double T_current = 0.0;
  arma::vec skeleton_T(N_skeleton);
  arma::mat skeleton_Xi(N_skeleton, d);
  arma::mat skeleton_Theta(N_skeleton, d);
  double max_alpha = 0.0;
  double mean_alpha = 0.0;
  int k_alpha = 0;
  int k_skeleton = 0;
  int print_count = 0;
  
  // Track which parameters are being updated
  arma::vec param_update_count = zeros(d);
  
  arma::vec Xi_current = xi_0;
  arma::vec Theta_current = theta_0;
  
  skeleton_T(k_skeleton) = T_current;
  skeleton_Xi.row(k_skeleton) = Xi_current.t();
  skeleton_Theta.row(k_skeleton) = Theta_current.t();
  
  while(T_current < T_end){
    if(k_skeleton >= N_skeleton - 1){
      Rcpp::warning("Skeleton storage full. Increase N_skeleton parameter.");
      break;
    }
    
    arma::vec tau_i = sim_M_MMD_RBF_gandk_cpp_arma(Theta_current, Xi_current, y, gamma, w);
    int i_0 = index_min(tau_i) + 1;
    double tau = tau_i(i_0 - 1);
    
    arma::vec u = u_generator_gandk_cpp_arma(m);
    double numerator = tilde_m_MMD_RBF_gandk_cpp_arma(i_0, tau, Theta_current, Xi_current, y, u, gamma, w);
    double denominator = M_MMD_RBF_gandk_cpp_arma(i_0, tau, Theta_current, Xi_current, y, gamma, w);
    double alpha = numerator / denominator;
    
    k_alpha += 1;
    max_alpha = max({max_alpha, alpha});
    mean_alpha = (mean_alpha*(k_alpha - 1) + alpha)/k_alpha;
    
    T_current += tau;
    Xi_current += Theta_current*tau;
    
    if(alpha > R::runif(0, 1)){
      double theta_new = -Theta_current(i_0-1);
      Theta_current(i_0-1) = theta_new;
      param_update_count(i_0-1) += 1;
      k_skeleton += 1;
      skeleton_T(k_skeleton) = T_current;
      skeleton_Xi.row(k_skeleton) = Xi_current.t();
      skeleton_Theta.row(k_skeleton) = Theta_current.t();
    }
    
    if(T_current > (print_count+1)*(T_end/10.0)){
      Rprintf("Time >  %f \n", T_current);
      print_count += 1;
    }
  }
  
  Rprintf("\nParameter update counts:\n");
  Rprintf("  a: %.0f\n", param_update_count(0));
  Rprintf("  b: %.0f\n", param_update_count(1));
  Rprintf("  g: %.0f\n", param_update_count(2));
  Rprintf("  k: %.0f\n", param_update_count(3));
  
  List out = List::create(Named("skeleton_T") = skeleton_T(span(0,k_skeleton)),
                          Named("skeleton_Xi") = skeleton_Xi.rows(span(0,k_skeleton)), 
                          Named("skeleton_Theta") = skeleton_Theta.rows(span(0,k_skeleton)),
                          Named("max_alpha") = max_alpha,
                          Named("mean_alpha") = mean_alpha,
                          Named("param_update_count") = param_update_count);  
  
  return out;
}

// [[Rcpp::export]]
arma::mat skeleton_to_sample_cpp_arma(const arma::vec& skeleton_T, const arma::mat& skeleton_Xi, 
                                      const arma::mat& skeleton_Theta, const int& N){
  int d = skeleton_Xi.n_cols;
  int k = skeleton_T.n_elem;
  
  if(k == 0){
    Rcpp::stop("Error: skeleton is empty");
  }
  
  double T_end = skeleton_T(k-1);
  arma::mat samples = zeros(N, d);
  arma::vec samplingTimes = linspace(1, N, N)*T_end*1.0/(N);
  
  int i;
  
  for(i = 0; i<N; i++){
    arma::vec diff = abs(skeleton_T - samplingTimes(i));
    uword k_i = index_min(diff);
    
    if(k_i >= k){
      Rcpp::stop("Error: k_i out of bounds in skeleton_to_sample");
    }
    
    if(skeleton_T(k_i) <= samplingTimes(i)){
      samples.row(i) = skeleton_Xi.row(k_i) + skeleton_Theta.row(k_i)*(samplingTimes(i) - skeleton_T(k_i));
    } else{
      if(k_i == 0){
        samples.row(i) = skeleton_Xi.row(0);
      } else {
        samples.row(i) = skeleton_Xi.row(k_i-1) + skeleton_Theta.row(k_i-1)*(samplingTimes(i) - skeleton_T(k_i-1));
      }
    }
  }
  
  return samples;
}

// Helper function to transform samples back to original scale [0, 5]
// [[Rcpp::export]]
arma::mat transform_samples_to_original_scale(const arma::mat& samples_unconstrained){
  arma::mat samples_original = samples_unconstrained;
  // All columns need to be transformed from unconstrained to [0, 5]
  for(int j = 0; j < 4; j++){
    samples_original.col(j) = transform_to_bounded_vec(samples_unconstrained.col(j));
  }
  return samples_original;
}


// From Gemini - need to test!
// [[Rcpp::export]]
double logSumExp_vec(arma::vec x) {
  // Find the maximum value in the vector
  double max_val = x.max();
  
  // Subtract the maximum value from all elements and then exponentiate
  arma::vec exp_shifted = arma::exp(x - max_val);
  
  // Sum the shifted exponentials and take the logarithm
  double sum_exp_shifted = arma::sum(exp_shifted);
  double log_sum_exp = std::log(sum_exp_shifted);
  
  // Add the maximum value back to the result
  return max_val + log_sum_exp;
}



// [[Rcpp::export]]
Rcpp::List K2_ABC_gandk_cpp_arma(const arma::vec& y, const double& epsilon, 
                              const double& N_MC, const int& m, const double& gamma){
  

  arma::vec a_sample = zeros(N_MC);
  arma::vec b_sample = zeros(N_MC);
  arma::vec g_sample = zeros(N_MC);
  arma::vec k_sample = zeros(N_MC);
  arma::vec log_tilde_weights = zeros(N_MC);
  arma::vec weights = zeros(N_MC);
  
  a_sample = arma::randu(N_MC, distr_param(0, 5));
  b_sample = arma::randu(N_MC, distr_param(0, 5));
  g_sample = arma::randu(N_MC, distr_param(0, 5));
  k_sample = arma::randu(N_MC, distr_param(0, 5));
  int print_count = 0;
  
  int j;
  for(j = 0; j<N_MC; j++){
    
    arma::vec tilde_u = u_generator_gandk_cpp_arma(m);
    arma::vec u = F_generator_gandk_cpp_arma(tilde_u, a_sample(j), b_sample(j), 
                                             -log(5/g_sample(j) - 1), -log(5/k_sample(j) - 1));
    
    double MMD = MMD_RBF_vec_cpp_arma(y, u, gamma);
    
    log_tilde_weights(j) = -MMD/epsilon;
    
    //if(j > (print_count+1)*(N_MC/10.0)){
    if(j > (print_count+0)*(N_MC/10.0)){
      Rprintf("Iteration  %i \n", j);
      print_count += 1;
    }
  }
  
  weights = exp(log_tilde_weights - logSumExp_vec(log_tilde_weights));
  
  List out = List::create(Named("a") = a_sample,
                          Named("b") = b_sample, 
                          Named("g") = g_sample,
                          Named("k") = k_sample,
                          Named("log_tilde_weights") = log_tilde_weights,
                          Named("weights") = weights);  
  
  
  return out;
  
}
