simulate_gandk_with_gradient <- function(n, parms) {
  #' Simulate from the g-and-k distribution with gradient
  #'
  #' @param n Number of observations to simulate
  #' @param parms Vector of parameters [a, b, g, k] for the g-and-k distribution
  #' @return A list containing:
  #'   - samples: n draws from the g-and-k distribution
  #'   - gradient: an n x 4 matrix of gradients (∂q/∂a, ∂q/∂b, ∂q/∂g, ∂q/∂k)
  
  # Generate standard normal random variables
  zu <- rnorm(n)
  
  # Extract parameters
  a <- parms[1]
  b <- parms[2]
  g <- parms[3]
  k <- parms[4]
  c_val <- 0.8  # Fixed constant as in MATLAB code
  
  # Compute intermediate terms for gradient calculations
  exp_gz <- exp(-g * zu)
  t_term <- (1 - exp_gz) / (1 + exp_gz)  # The tanh-like term
  inner_term <- 1 + c_val * t_term
  zu_sq <- zu^2
  power_term <- (1 + zu_sq)^k
  
  # Compute the samples
  samples <- a + b * inner_term * power_term * zu
  
  # Compute gradients
  # Initialize gradient matrix (n x 4)
  gradient <- matrix(0, nrow = n, ncol = 4)
  colnames(gradient) <- c("da", "db", "dg", "dk")
  
  # 1. Gradient with respect to a (∂q/∂a)
  gradient[, 1] <- 1
  
  # 2. Gradient with respect to b (∂q/∂b)
  gradient[, 2] <- inner_term * power_term * zu
  
  # 3. Gradient with respect to g (∂q/∂g)
  # Using chain rule: ∂q/∂g = b * zu * (1 + zu^2)^k * ∂/∂g[1 + c*(1-exp(-g*zu))/(1+exp(-g*zu))]
  # Derivative of t_term with respect to g: ∂t/∂g = (2*c*zu*exp(-g*zu))/(1+exp(-g*zu))^2
  # But note: t = (1-exp(-gz))/(1+exp(-gz)) = tanh(g*zu/2)
  # So ∂t/∂g = zu * (1 - t^2) / 2
  dt_dg <- zu * (1 - t_term^2) / 2
  gradient[, 3] <- b * c_val * dt_dg * power_term * zu
  
  # 4. Gradient with respect to k (∂q/∂k)
  # Using chain rule: ∂q/∂k = b * (1 + c*t) * zu * ∂/∂k[(1 + zu^2)^k]
  # Derivative: ∂/∂k[(1 + zu^2)^k] = (1 + zu^2)^k * log(1 + zu^2)
  gradient[, 4] <- b * inner_term * zu * power_term * log(1 + zu_sq)
  
  return(list(samples = samples, gradient = gradient))
}


# Transform samples from log scale to original scale
# Input: samples_log_scale is an N x 4 matrix with columns (a, b, log_g, log_k)
# Output: samples_original is an N x 4 matrix with columns (a, b, g, k)

transform_samples_to_original_scale_R <- function(samples_log_scale) {
  samples_original <- samples_log_scale
  
  # Columns 1 and 2 (a and b) stay the same
  # Columns 3 and 4 (log_g and log_k) need to be exponentiated
  samples_original[, 3] <- exp(samples_log_scale[, 3])
  samples_original[, 4] <- exp(samples_log_scale[, 4])
  
  # Optionally add column names
  colnames(samples_original) <- c("a", "b", "g", "k")
  
  return(samples_original)
}

# Alternative: if you want to keep both versions
transform_samples_keep_both <- function(samples_log_scale) {
  samples_combined <- 5.0 / (1.0 + cbind(
    samples_log_scale,
    exp(samples_log_scale[, 3]),
    exp(samples_log_scale[, 4])))
  
  colnames(samples_combined) <- c("a", "b", "log_g", "log_k", "g", "k")
  
  return(samples_combined)
}


rgk <- function(n, parm1,parm2,parm3,parm4) {
  #' Simulate from the g-and-k distribution with gradient
  #'
  #' @param n Number of observations to simulate
  #' @param parms Vector of parameters [a, b, g, k] for the g-and-k distribution
  #' @return A list containing:
  #'   - samples: n draws from the g-and-k distribution
  #'   - gradient: an n x 4 matrix of gradients (∂q/∂a, ∂q/∂b, ∂q/∂g, ∂q/∂k)
  
  # Generate standard normal random variables
  zu <- rnorm(n)
  
  # Extract parameters
  a <- parm1
  b <- parm2
  g <- parm3
  k <- parm4
  c_val <- 0.8  # Fixed constant as in MATLAB code
  
  # Compute intermediate terms for gradient calculations
  exp_gz <- exp(-g * zu)
  t_term <- (1 - exp_gz) / (1 + exp_gz)  # The tanh-like term
  inner_term <- 1 + c_val * t_term
  zu_sq <- zu^2
  power_term <- (1 + zu_sq)^k
  
  # Compute the samples
  samples <- a + b * inner_term * power_term * zu
}