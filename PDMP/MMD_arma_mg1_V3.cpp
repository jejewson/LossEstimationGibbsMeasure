// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>

using namespace arma;
using namespace Rcpp;

// =============================================================================
// SHARED KERNEL FUNCTIONS
// =============================================================================


// Gradient of k_gamma(z, x) w.r.t. z.
// [[Rcpp::export]]
arma::vec grad_z_RBF_kernel_mg1_cpp_arma(const double& z, const arma::vec& x,
                                         const double& gamma) {
  return (-(z - x) / (std::sqrt(2.0 * M_PI) * std::pow(gamma, 1.5)))
  % arma::exp(-(z - x) % (z - x) / (2.0 * gamma));
}

// =============================================================================
// M/G/1 LATENT VARIABLE GENERATOR
//
// We generate m independent simulation runs, each requiring (n_obs + 1)
// customers.  The latent vector u has length  2 * m * (n_obs + 1):
//
//   u[ 0                    ... m*(n_obs+1) - 1 ]  =  Exp(1)      latents
//   u[ m*(n_obs+1)          ... 2*m*(n_obs+1) - 1 ]  =  Uniform(0,1) latents
//
// Within each half the m runs are stored consecutively:
//   run r  ->  Exp(1)  block:  u[ r*(n_obs+1) ... (r+1)*(n_obs+1) - 1 ]
//   run r  ->  Unif    block:  u[ m*(n_obs+1) + r*(n_obs+1) ... ]
//
// This mirrors the g-and-k convention where u has length m (one latent per
// simulated draw), but adapted for the M/G/1 case where each "draw" is a
// full queue run of (n_obs + 1) customers.
// [[Rcpp::export]]
arma::mat u_generator_mg1_cpp_arma(const int& m) {
  int T_queue = m + 1;            // customers per run
  arma::mat u(T_queue, 2);
  u.col(0) = -arma::log(arma::randu(T_queue)); // Exp(1)
  u.col(1) =  arma::randu(T_queue);            // Uniform(0,1)
  return u;
}



// =============================================================================
// INTERNAL HELPERS
// =============================================================================
// [[Rcpp::export]]
static arma::vec xi_to_theta_mg1(const arma::vec& xi,
                                 const arma::vec& p1,
                                 const arma::vec& p2) {
  arma::vec ex = arma::exp(xi);
  return (p2 % ex + p1) / (1.0 + ex);
}

// =============================================================================
// M/G/1 GENERATOR: m independent queue runs -> m * n_obs inter-departure times
//
// u      : latent vector of length 2*m*(n_obs+1) from u_generator_mg1_cpp_arma.
// xi     : unconstrained parameters.
// p1, p2 : prior bounds.
// n_obs  : number of observed inter-departure times (= number to simulate per run).
//
// Returns a vector of length m * n_obs:
//   [ x^{(1)}_1, ..., x^{(1)}_{n_obs},  x^{(2)}_1, ..., x^{(m)}_{n_obs} ]
// where x^{(r)} is the inter-departure times from the r-th independent run.
// [[Rcpp::export]]
arma::vec F_generator_mg1_cpp_arma(const arma::mat& u,
                                   const arma::vec& xi,
                                   const arma::vec& p1,
                                   const arma::vec& p2) {

  
  int m = u.n_rows - 1;
  int T_queue = m + 1;
  
  arma::vec theta = xi_to_theta_mg1(xi, p1, p2);
  double th1 = theta(0), th2 = theta(1), th3 = theta(2);
  
  arma::vec z(m);   // output
  
    
  arma::vec w_r = u.col(0) / th3;
  arma::vec s_r = th1 + (th2 - th1) * u.col(1);
    
  arma::vec D_r(T_queue);
  double d = 0.0, a = 0.0;
  for (int i = 0; i < T_queue; i++) {
    a    += w_r(i);
    d     = d + s_r(i) + std::max(0.0, a - d);
    D_r(i) = d;
  }
  z = arma::diff(D_r);

  return z;
}


// =============================================================================
// M/G/1 GENERATOR GRADIENT: gradients of each simulated x_i w.r.t. xi_{i_0}
//
// Returns a vector of length m * n_obs:
//   [ dx^{(1)}/dxi_{i_0}, ..., dx^{(m)}/dxi_{i_0} ]
// (same layout as F_generator_mg1_cpp_arma).
//
// Each run's gradients are computed independently via the queue backpropagation
// recursion.  The Exp(1) accumulator a_grad (for theta_3) is reset per run.
//
// Gradient recursion (see companion derivation):
//   Define  I_i = 1{A_i > D_{i-1}},   a_i^(j) = dA_i/dtheta_j
//
//   theta_1:  d_i^(1) = (1-I_i)*d_{i-1}^(1) + (1-v_i)
//   theta_2:  d_i^(2) = (1-I_i)*d_{i-1}^(2) + v_i
//   theta_3:  a_i^(3) = a_{i-1}^(3) - e_i/theta3^2
//             d_i^(3) = (1-I_i)*d_{i-1}^(3) + I_i * a_i^(3)
//
//   dx_i/dxi_{i_0} = [d_i^(j) - d_{i-1}^(j)] * dtheta_{i_0}/dxi_{i_0}
// [[Rcpp::export]]
arma::vec grad_F_generator_mg1_cpp_arma(const int& i_0,
                                        const arma::mat& u,
                                        const arma::vec& xi,
                                        const arma::vec& p1,
                                        const arma::vec& p2) {
  int m = u.n_rows - 1;
  int T_queue = m + 1;
  
  arma::vec theta = xi_to_theta_mg1(xi, p1, p2);
  double th1 = theta(0), th2 = theta(1), th3 = theta(2);
  
  // Jacobian scalar: dtheta_{i_0}/dxi_{i_0} = (p2-p1)*exp(xi)/(1+exp(xi))^2
  double ex_j       = std::exp(xi(i_0 - 1));
  double dtheta_dxi = (p2(i_0-1) - p1(i_0-1)) * ex_j /
    std::pow(1.0 + ex_j, 2.0);
  
  arma::vec grad_all(m, arma::fill::zeros);   // output
  
  arma::vec e_r = u.col(0);
  arma::vec v_r = u.col(1);
  
  arma::vec w_r = e_r / th3;
  arma::vec s_r = th1 + (th2 - th1) * v_r;
    
  // --- Forward queue pass: record A_i, D_i, I_i ---
  arma::vec A_r(T_queue), D_r(T_queue), I_r(T_queue);
  double d = 0.0, a = 0.0;
  for (int i = 0; i < T_queue; i++) {
    a      += w_r(i);
    A_r(i)  = a;
    I_r(i)  = (a > d) ? 1.0 : 0.0;
    d       = d + s_r(i) + std::max(0.0, a - d);
    D_r(i)  = d;
  }
    
  // --- Gradient recursion for departure times ---
  arma::vec d_grad(T_queue, arma::fill::zeros);
  double a_grad = 0.0;   // dA_i/dtheta_3; reset for each independent run
    
  for (int i = 0; i < T_queue; i++) {
    double prev_d = (i > 0) ? d_grad(i - 1) : 0.0;
    switch (i_0) {
    case 1:
      d_grad(i) = (1.0 - I_r(i)) * prev_d + (1.0 - v_r(i));
    break;
    case 2:
      d_grad(i) = (1.0 - I_r(i)) * prev_d + v_r(i);
      break;
    case 3:
      a_grad   -= e_r(i) / (th3 * th3);
      d_grad(i) = (1.0 - I_r(i)) * prev_d + I_r(i) * a_grad;
      break;
    }
  }
    
  // dx_i/dxi = diff(d_grad) * dtheta/dxi
  grad_all = arma::diff(d_grad) * dtheta_dxi;

  return grad_all;
}

// [[Rcpp::export]]
arma::mat grad_F_and_samples_generator_mg1_cpp_arma(const int& i_0,
                                                    const arma::mat& u,
                                                    const arma::vec& xi,
                                                    const arma::vec& p1,
                                                    const arma::vec& p2) {
  int m = u.n_rows - 1;
  int T_queue = m + 1;
  
  arma::vec theta = xi_to_theta_mg1(xi, p1, p2);
  double th1 = theta(0), th2 = theta(1), th3 = theta(2);
  
  // Jacobian scalar: dtheta_{i_0}/dxi_{i_0} = (p2-p1)*exp(xi)/(1+exp(xi))^2
  double ex_j       = std::exp(xi(i_0 - 1));
  double dtheta_dxi = (p2(i_0-1) - p1(i_0-1)) * ex_j /
    std::pow(1.0 + ex_j, 2.0);
  
  arma::vec grad_all(m, arma::fill::zeros);   // output
  arma::vec samples(m, arma::fill::zeros);   // output
  
  arma::vec e_r = u.col(0);
  arma::vec v_r = u.col(1);
    
  arma::vec w_r = e_r / th3;
  arma::vec s_r = th1 + (th2 - th1) * v_r;
    
  // --- Forward queue pass: record A_i, D_i, I_i ---
  arma::vec A_r(T_queue), D_r(T_queue), I_r(T_queue);
  double d = 0.0, a = 0.0;
  for (int i = 0; i < T_queue; i++) {
    a      += w_r(i);
    A_r(i)  = a;
    I_r(i)  = (a > d) ? 1.0 : 0.0;
    d       = d + s_r(i) + std::max(0.0, a - d);
    D_r(i)  = d;
  }
    
  samples = arma::diff(D_r);
    
  // --- Gradient recursion for departure times ---
  arma::vec d_grad(T_queue, arma::fill::zeros);
  double a_grad = 0.0;   // dA_i/dtheta_3; reset for each independent run
    
  for (int i = 0; i < T_queue; i++) {
    double prev_d = (i > 0) ? d_grad(i - 1) : 0.0;
    switch (i_0) {
    case 1:
      d_grad(i) = (1.0 - I_r(i)) * prev_d + (1.0 - v_r(i));
      break;
    case 2:
      d_grad(i) = (1.0 - I_r(i)) * prev_d + v_r(i);
      break;
    case 3:
      a_grad   -= e_r(i) / (th3 * th3);
      d_grad(i) = (1.0 - I_r(i)) * prev_d + I_r(i) * a_grad;
      break;
    }
  }
    
  // dx_i/dxi = diff(d_grad) * dtheta/dxi
  grad_all = arma::diff(d_grad) * dtheta_dxi;
  
  arma::mat out(m, 2);
  out.col(0) = samples;
  out.col(1) = grad_all;
  
  return out;
}

// Biased empirical MMD^2 between observed y (length n) and simulated z (length n_sim).
// [[Rcpp::export]]
double MMD_RBF_vec_mg1_cpp_arma(const arma::vec& y, const arma::mat& u,
                                const arma::vec& xi,
                                const arma::vec& p1,
                                const arma::vec& p2, const double& gamma) {
  
  int n     = y.n_elem;
  //int n_sim = u.n_elem;
  //int T_queue = n_obs + 1;
  //int m       = u.n_elem / (2 * T_queue);
  
  arma::vec z = F_generator_mg1_cpp_arma(u, xi, p1, p2);
  int n_sim = z.n_elem;
  
  //Rprintf("1");
  
  double out1 = 0.0, out2 = 0.0;
  
  for (int i = 0; i < n; i++){
    out1 -= arma::sum(arma::normpdf(y(i) - z, 0.0, std::sqrt(gamma)).col(0));
  }
  out1 *= 2.0 / n_sim;
  
  //Rprintf("2");
  
  for (int k = 0; k < n_sim; k++) {
    arma::vec l    = arma::regspace(0, n_sim - 1);
    arma::uvec ind = arma::find(l != k);
    out2 += arma::sum(arma::normpdf(z(k) - z(ind), 0.0, std::sqrt(gamma)).col(0));
  }
  out2 /= (n_sim * (n_sim - 1.0));
  
  return out1 + n * out2;
}


// =============================================================================
// MMD GRADIENT
// =============================================================================

// Gradient of the MMD objective w.r.t. xi_{i_0}.
//
// z        = F_generator(u, xi, ...)   has length m * n_obs
// grad_all = grad_F_generator(...)     has length m * n_obs
//
// The formula is identical to the g-and-k version; the only difference is
// that n_sim = m * n_obs instead of m.  The bound analysis shows that the
// resulting upper bound on n*w*|grad_MMD| is INDEPENDENT of m (the n_sim
// terms cancel), so the same constant M_j remains valid for all m.
// [[Rcpp::export]]
double grad_MMD_RBF_mg1_cpp_arma(const int& i_0,
                                 const arma::vec& y,
                                 const arma::mat& u,
                                 const arma::vec& xi,
                                 const arma::vec& p1,
                                 const arma::vec& p2,
                                 const double& gamma) {
  int n     = y.n_elem;                      // n_obs == length of y by construction
  arma::mat temp = grad_F_and_samples_generator_mg1_cpp_arma(i_0, u, xi, p1, p2);
  //arma::vec z        = F_generator_mg1_cpp_arma(u, xi, p1, p2);
  //int m          = z.n_elem;              // = m * n_obs
  //arma::vec grad_all = grad_F_generator_mg1_cpp_arma(i_0, u, xi, p1, p2);
  arma::vec z        = temp.col(0);
  int m          = z.n_elem;              // = m * n_obs
  arma::vec grad_all = temp.col(1);
  
  double out1 = 0.0, out2 = 0.0;
  
  // Term 1: gradient of cross-term  -2 E_{z,y}[k(z,y)]
  for (int j = 0; j < m; j++)
    out1 -= arma::sum(grad_z_RBF_kernel_mg1_cpp_arma(z(j), y, gamma)) * grad_all(j);
  //for (int j = 0; j < n; j++)
  //      out1 -= arma::sum(grad_z_RBF_kernel_mg1_cpp_arma(y(j), z, gamma) % grad_all);
  out1 *= 2.0 / (m);
  
  // Term 2: gradient of simulated self-term  E_{z,z'}[k(z,z')]
  for (int j = 0; j < m; j++) {
    arma::vec idx  = arma::regspace(0, m - 1);
    arma::uvec ind = arma::find(idx != j);
    out2 += arma::sum(
      grad_z_RBF_kernel_mg1_cpp_arma(z(j), z(ind), gamma) %
        (grad_all(j) - grad_all(ind))
    ) / (m * (m - 1.0));
  }
  return out1 + n * out2;
}

// =============================================================================
// PRIOR GRADIENT
// =============================================================================

// d/dxi_j  log pi(xi_j) = (1-exp(xi_j))/(1+exp(xi_j)) = -tanh(xi_j/2)
// Strictly bounded: |value| < 1  for all finite xi_j.
// [[Rcpp::export]]
double grad_log_prior_mg1_cpp_arma(const int& i_0, const arma::vec& xi) {
  double ex = std::exp(xi(i_0 - 1));
  return (1.0 - ex) / (1.0 + ex);
}


// =============================================================================
// ZIGZAG THINNING AND CONSTANT UPPER BOUND
// =============================================================================

// True ZigZag switching intensity (Poisson thinning numerator).
// u encodes m independent runs; its length implicitly carries m.
// [[Rcpp::export]]
double tilde_m_MMD_RBF_mg1_cpp_arma(const int& i_0,
                                    const double& t,
                                    const arma::vec& theta_curr,
                                    const arma::vec& xi_curr,
                                    const arma::vec& y,
                                    const arma::mat& u,
                                    const arma::vec& p1,
                                    const arma::vec& p2,
                                    const double& gamma,
                                    const double& w) {
  //int n = y.n_elem;
  arma::vec xi_t = xi_curr + theta_curr * t;
  
  double prior_grad = grad_log_prior_mg1_cpp_arma(i_0, xi_t);
  double mmd_grad   = grad_MMD_RBF_mg1_cpp_arma(i_0, y, u, xi_t, p1, p2, gamma);
  
  double rate = -theta_curr(i_0-1) * prior_grad
  + theta_curr(i_0-1) * w * mmd_grad;
  return (rate > 0.0) ? rate : 0.0;
}

// Constant upper bound M_j.
//
// Key insight: the bound on n*w*|grad_MMD| is INDEPENDENT of m because the
// m*n_obs factors in the cross- and self-terms both cancel.
// Specifically:
//   n*w * |Term 1| <= 2 * n*w * max|grad_K| * max|grad_x|
//                  = kernel_bound * max|grad_x|
// and similarly for Term 2, where
//   kernel_bound = 2 * n * w * exp(-0.5) / (sqrt(2*pi) * gamma).
//
// The generator gradient is further bounded via the Jacobian max:
//   max|dtheta_j/dxi_j| = (p2_j - p1_j) / 4
// The factor of 20 is a conservative multiplier for the queue gradient;
// reduce it if mean_alpha is very low, increase it if max_alpha > 1.
// [[Rcpp::export]]
double M_MMD_RBF_mg1_cpp_arma(const int& i_0,
                              const int& n,
                              const double& gamma,
                              const double& w,
                              const arma::vec& p1,
                              const arma::vec& p2) {
  double kernel_bound = 2.0 * n * w * std::exp(-0.5) /
    (std::sqrt(2.0 * M_PI) * gamma);
  double scale;
  switch (i_0) {
  case 1:  scale = (p2(0) - p1(0)) / 4.0 * 1;  break;
  case 2:  scale = (p2(1) - p1(1)) / 4.0 * 1;  break;
  default: scale = (p2(2) - p1(2)) / 4.0 * 40.0;  break;
  }
  return 1.0 + kernel_bound * scale;
}

// Simulate candidate event times for all three coordinates (constant-rate -> Exp).
// [[Rcpp::export]]
arma::vec sim_M_MMD_RBF_mg1_cpp_arma(const int& n,
                                     const double& gamma,
                                     const double& w,
                                     const arma::vec& p1,
                                     const arma::vec& p2) {
  arma::vec tau(3);
  for (int j = 0; j < 3; j++) {
    double M_j = M_MMD_RBF_mg1_cpp_arma(j + 1, n, gamma, w, p1, p2);
    tau(j) = -std::log(R::runif(0.0, 1.0)) / M_j;
  }
  return tau;
}

// =============================================================================
// MAIN ZIGZAG SAMPLER
// =============================================================================

// Run the ZigZag PDMP sampler for the M/G/1 queue posterior.
//
// m          : number of independent M/G/1 simulation runs per thinning step.
//              Increasing m reduces gradient variance (better mixing) at the
//              cost of m times more queue simulations per step.
//              Start with m = 5 (same default as the g-and-k example).
// [[Rcpp::export]]
Rcpp::List ZigZag_MMD_RBF_mg1_cpp_arma(const double& T_end,
                                       const arma::vec& xi_0,
                                       const arma::vec& theta_0,
                                       const arma::vec& y,
                                       const int& m,
                                       const arma::vec& p1,
                                       const arma::vec& p2,
                                       const double& gamma,
                                       const double& w,
                                       const int& N_skeleton) {
  int d = 3;
  int n = y.n_elem;
  double T_current = 0.0;
  
  arma::vec skeleton_T(N_skeleton);
  arma::mat skeleton_Xi(N_skeleton, d);
  arma::mat skeleton_Theta(N_skeleton, d);
  int k_skeleton = 0;
  double max_alpha = 0.0, mean_alpha = 0.0;
  int k_alpha = 0, print_count = 0;
  
  arma::vec Xi_current    = xi_0;
  arma::vec Theta_current = theta_0;
  
  skeleton_T(0)         = T_current;
  skeleton_Xi.row(0)    = Xi_current.t();
  skeleton_Theta.row(0) = Theta_current.t();
  
  while (T_current < T_end && k_skeleton < N_skeleton - 1) {
    
    // 1. Candidate event times from constant-rate Poisson processes.
    arma::vec tau_i = sim_M_MMD_RBF_mg1_cpp_arma(n, gamma, w, p1, p2);
    int i_0    = arma::index_min(tau_i) + 1;
    double tau = tau_i(i_0 - 1);
    
    // 2. Generate latent variables for m independent M/G/1 runs.
    arma::mat u = u_generator_mg1_cpp_arma(m);
    
    // 3. Thinning probability: true rate / constant upper bound.
    double tilde_m_val = tilde_m_MMD_RBF_mg1_cpp_arma(
      i_0, tau, Theta_current, Xi_current, y, u, p1, p2, gamma, w);
    double M_bound = M_MMD_RBF_mg1_cpp_arma(i_0, n, gamma, w, p1, p2);
    double alpha   = tilde_m_val / M_bound;
    
    k_alpha++;
    if (alpha > max_alpha) max_alpha = alpha;
    mean_alpha = (mean_alpha * (k_alpha - 1) + alpha) / k_alpha;
    
    // 4. Advance linearly.
    T_current   += tau;
    Xi_current  += Theta_current * tau;
    
    // 5. Accept flip.
    if (alpha > R::runif(0.0, 1.0)) {
      Theta_current(i_0 - 1) = -Theta_current(i_0 - 1);
      k_skeleton++;
      skeleton_T(k_skeleton)         = T_current;
      skeleton_Xi.row(k_skeleton)    = Xi_current.t();
      skeleton_Theta.row(k_skeleton) = Theta_current.t();
    }
    
    if (T_current > (print_count + 1) * (T_end / 10.0)) {
      Rprintf("Time > %.1f  |  skeleton: %d  |  mean alpha: %.4f  |  m: %d\n",
              T_current, k_skeleton, mean_alpha, m);
      print_count++;
    }
  }
  
  if (max_alpha > 1.0)
    Rprintf("WARNING: max_alpha = %.4f > 1. Increase scale factors in "
              "M_MMD_RBF_mg1_cpp_arma.\n", max_alpha);
  
  return Rcpp::List::create(
    Rcpp::Named("skeleton_T")     = skeleton_T.head(k_skeleton + 1),
    Rcpp::Named("skeleton_Xi")    = skeleton_Xi.rows(0, k_skeleton),
    Rcpp::Named("skeleton_Theta") = skeleton_Theta.rows(0, k_skeleton),
    Rcpp::Named("max_alpha")      = max_alpha,
    Rcpp::Named("mean_alpha")     = mean_alpha
  );
}

// =============================================================================
// POST-PROCESSING UTILITIES
// =============================================================================

// [[Rcpp::export]]
arma::mat skeleton_to_sample_mg1_cpp_arma(const arma::vec& skeleton_T,
                                          const arma::mat& skeleton_Xi,
                                          const arma::mat& skeleton_Theta,
                                          const int& N) {
  int d    = skeleton_Xi.n_cols;
  int k    = skeleton_T.n_elem;
  double T_end = skeleton_T(k - 1);
  arma::mat samples(N, d, arma::fill::zeros);
  arma::vec sampleTimes = arma::linspace(1, N, N) * T_end / N;
  
  for (int i = 0; i < N; i++) {
    int k_i = arma::index_min(arma::abs(skeleton_T - sampleTimes(i)));
    if (skeleton_T(k_i) <= sampleTimes(i)) {
      samples.row(i) = skeleton_Xi.row(k_i) +
        skeleton_Theta.row(k_i) *
        (sampleTimes(i) - skeleton_T(k_i));
    } else {
      samples.row(i) = skeleton_Xi.row(k_i - 1) +
        skeleton_Theta.row(k_i - 1) *
        (sampleTimes(i) - skeleton_T(k_i - 1));
    }
  }
  return samples;
}

// [[Rcpp::export]]
arma::mat transform_samples_mg1(const arma::mat& samples_xi,
                                const arma::vec& p1,
                                const arma::vec& p2) {
  int N = samples_xi.n_rows;
  arma::mat samples_theta(N, 3);
  for (int j = 0; j < 3; j++) {
    arma::vec ex = arma::exp(samples_xi.col(j));
    samples_theta.col(j) = (p2(j) * ex + p1(j)) / (1.0 + ex);
  }
  return samples_theta;
}
