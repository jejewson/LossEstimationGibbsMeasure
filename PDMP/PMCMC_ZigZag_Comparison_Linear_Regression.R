
library(adaptMCMC)
library(actuar)
library(metRology)
library(LaplacesDemon)
library(invgamma)
library(MASS)

library(Rcpp)

sourceCpp("MMD_regression_arma.cpp")


# function for computing log loss
log_loss_mmd <- function(theta,y,X,m,w,gamma_mmd){

  beta = theta[1:p]
  sigma = exp(theta[p+1])
  
  mmd_all = rep(0,ndata)
  
  for (i in 1:ndata){
    x <- c(X[i,]%*%beta) + rnorm(m, 0, sigma)
    mmd_all[i] = MMD_RBF_vec_cpp_arma(y[i], x, gamma_mmd)
  }
  
  log_prior = sum(dnorm(beta, 0, s_0, log = TRUE)) + dinvgamma(sigma^2, shape = a_0, scale = b_0, log = TRUE) + 2*theta[p+1]
  
  return(-w*sum(mmd_all) + log_prior)

}

# function for computing log loss for a given set of random numbers
log_loss_mmd_block <- function(theta,z,y,X,m,w,gamma_mmd){
  
  beta = theta[1:p]
  sigma = exp(theta[p+1])
  
  mmd_all = rep(0,ndata)
  
  for (i in 1:ndata){
    x <- c(X[i,]%*%beta) + sigma*z[i,]
    mmd_all[i] = MMD_RBF_vec_cpp_arma(y[i], x, gamma_mmd)
  }
  
  log_prior = sum(dnorm(beta, 0, s_0, log = TRUE)) + dinvgamma(sigma^2, shape = a_0, scale = b_0, log = TRUE) + 2*theta[p+1]
  
  return(-w*sum(mmd_all) + log_prior)
  
}



# function for block correlated PMCMC algorithm
pmcmc_block <- function(theta_init,m,w,gamma_mmd,cov_rw,M){
  
  theta_curr = theta_init

  z_curr = matrix(rnorm(ndata*m), nrow = ndata, ncol = m)  # initialise random numbers for generating pseudo datasets
  
  
  logpost_curr = log_loss_mmd_block(theta_curr,z_curr,y,X,m,w,gamma_mmd)
  
  theta = matrix(0,M,p+1)
  
  for (i in 1:M){

    theta_prop = mvrnorm(n = 1, theta_curr, cov_rw);

    # update random numbers for block r
    r_prop = sample(ndata,1)
    z_prop = z_curr
    z_prop[r_prop,] = rnorm(m,0,1)
    
    logpost_prop = log_loss_mmd_block(theta_prop,z_prop,y,X,m,w,gamma_mmd)
    
    if (exp(logpost_prop - logpost_curr) > runif(1)){
      theta_curr = theta_prop
      logpost_curr = logpost_prop
      z_curr = z_prop
    }
    theta[i,] = theta_curr
    
  }
  
  return(theta)
  
}




# function for standard PMCMC algorithm
pmcmc <- function(theta_init,m,w,gamma_mmd,cov_rw,M){
  
  theta_curr = theta_init
  
  logpost_curr = log_loss_mmd(theta_curr,y,X,m,w,gamma_mmd)
  
  theta = matrix(0,M,p+1)
  
  for (i in 1:M){
    
    theta_prop = mvrnorm(n = 1, theta_curr, cov_rw);
    
    logpost_prop = log_loss_mmd(theta_prop,y,X,m,w,gamma_mmd)
    
    if (exp(logpost_prop - logpost_curr) > runif(1)){
      theta_curr = theta_prop
      logpost_curr = logpost_prop
    }
    theta[i,] = theta_curr
    
  }
  
  return(theta)
  
}




######## generate observed dataset of size n = 100

set.seed(1)  # set seed for reproducibility

N_rep <- 1

## "observed" data
beta_0 <- c(4, 4, 3, 3, 2, 2, 1, 1) 
p <- length(beta_0)
sigma_0 <- 1

ndata <- 100

X <- matrix(rnorm(ndata*p, 0, 1), nrow = ndata, ncol = p)

y <- X%*%beta_0 + rlaplace(ndata, 0, sigma_0)

OLS <- lm(y ~ X + 0)
OLS$coefficients


# prior info

m_0 <- 0
s_0 <- 5
a_0 <- 2
b_0 <- 0.5


w = sqrt(2*pi)   # learning rate (the learning rate is actually 1, but this takes into account a scalar factor of 1/sqrt(2*pi)) used for convenience in MMD calculation

gamma_mmd <- 1 # hyperparameter for MMD

sigma2 = 1



######### run standard pseudo-marginal sampler


load(file = "cov_rw_n100_zigzag.RData")   # load in tuned random walk proposal and initial values

m = 20
start.time <- Sys.time()
samples_pmcmc_m20 = pmcmc(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,500000)
end.time <- Sys.time()
time.taken.standard_m20 <- end.time - start.time


m = 10
start.time <- Sys.time()
samples_pmcmc_m10 = pmcmc(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,1000000)
end.time <- Sys.time()
time.taken.standard_m10 <- end.time - start.time


m = 5
start.time <- Sys.time()
samples_pmcmc_m5 = pmcmc(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,2000000)
end.time <- Sys.time()
time.taken.standard_m5 <- end.time - start.time

m = 2
start.time <- Sys.time()
samples_pmcmc_m2 = pmcmc(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,5000000)
end.time <- Sys.time()
time.taken.standard_m2 <- end.time - start.time


save(file = "ResultsPM.RData", samples_pmcmc_m20, time.taken.standard_m20, samples_pmcmc_m10, time.taken.standard_m10, samples_pmcmc_m5, time.taken.standard_m5, samples_pmcmc_m2, time.taken.standard_m2)





######### run correlated pseudo-marginal sampler

load(file = "cov_rw_n100_zigzag.RData")   # load in tuned random walk proposal and initial values


m = 100
start.time <- Sys.time()
samples_pmcmc_block_m100 = pmcmc_block(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,100000)
end.time <- Sys.time()
time.taken.m100 <- end.time - start.time


m = 20
start.time <- Sys.time()
samples_pmcmc_block_m20 = pmcmc_block(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,500000)
end.time <- Sys.time()
time.taken.m20 <- end.time - start.time

m = 10
start.time <- Sys.time()
samples_pmcmc_block_m10 = pmcmc_block(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,1000000)
end.time <- Sys.time()
time.taken.m10 <- end.time - start.time


m = 5
start.time <- Sys.time()
samples_pmcmc_block_m5 = pmcmc_block(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,2000000)
end.time <- Sys.time()
time.taken.m5 <- end.time - start.time


m = 2
start.time <- Sys.time()
samples_pmcmc_block_m2 = pmcmc_block(c(beta_0, log(1)),m,w,gamma_mmd,cov_rw_zigzag,5000000)
end.time <- Sys.time()
time.taken.m2 <- end.time - start.time

save(file = "ResultsPMBlock.RData", samples_pmcmc_block_m100, samples_pmcmc_block_m20, samples_pmcmc_block_m10, samples_pmcmc_block_m5, samples_pmcmc_block_m2, time.taken.m100, time.taken.m20, time.taken.m10, time.taken.m5, time.taken.m2)


########### compare to zig zag


library(actuar)
library(metRology)
library(LaplacesDemon)


library(Rcpp)

sourceCpp("MMD_regression_arma.cpp")


set.seed(1)

N_rep <- 1

## "observed" data
beta_0 <- c(4, 4, 3, 3, 2, 2, 1, 1) 
p <- length(beta_0)
sigma_0 <- 1

ndata <- 100

X <- matrix(rnorm(ndata*p, 0, 1), nrow = ndata, ncol = p)

y <- X%*%beta_0 + rlaplace(ndata, 0, sigma_0)

OLS <- lm(y ~ X + 0)
OLS$coefficients


# prior info

m_0 <- 0
s_0 <- 5
a_0 <- 2
b_0 <- 0.5




prob_all <- 10^{-6}

gamma <- 1


###### m = 2

m <- 2

prob_one <- prob_all^{1/m}
prob_one

sqrt(-2*log(punif(q=prob_one, 0, 1)))

R <- sqrt(-2*log(punif(q=prob_one, 0, 1)))

R

T_end <- 10000
xi_0 <- c(beta_0, log(1))
theta_0 <- rep(1, p+1)


timer_MMD_Gaussian_regression_m2.start <- Sys.time()
skeleton_Gaussian_regression_IntractLik_m2 <- ZigZag_MMD_RBF_Gaussian_regression_cpp_arma(T_end, xi_0, theta_0, y, X, m, gamma, m_0, s_0, a_0, b_0, R, w = sqrt(2*pi), N_skeleton = 2000000)
timer_MMD_Gaussian_regression_m2.end <- Sys.time()
timer_MMD_Gaussian_regression_m2 = timer_MMD_Gaussian_regression_m2.end - timer_MMD_Gaussian_regression_m2.start


sample_Gaussian_regression_IntractLik_m2 <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_Gaussian_regression_IntractLik_m2$skeleton_T, skeleton_Xi = skeleton_Gaussian_regression_IntractLik_m2$skeleton_Xi, skeleton_Theta = skeleton_Gaussian_regression_IntractLik_m2$skeleton_Theta, N = 50000)






###### m = 5

m <- 5

prob_one <- prob_all^{1/m}
prob_one

sqrt(-2*log(punif(q=prob_one, 0, 1)))

R <- sqrt(-2*log(punif(q=prob_one, 0, 1)))

R

T_end <- 10000
xi_0 <- c(beta_0, log(1))
theta_0 <- rep(1, p+1)


timer_MMD_Gaussian_regression_m5.start <- Sys.time()
skeleton_Gaussian_regression_IntractLik_m5 <- ZigZag_MMD_RBF_Gaussian_regression_cpp_arma(T_end, xi_0, theta_0, y, X, m, gamma, m_0, s_0, a_0, b_0, R, w = sqrt(2*pi), N_skeleton = 1000000)
timer_MMD_Gaussian_regression_m5.end <- Sys.time()
timer_MMD_Gaussian_regression_m5 = timer_MMD_Gaussian_regression_m5.end - timer_MMD_Gaussian_regression_m5.start


sample_Gaussian_regression_IntractLik_m5 <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_Gaussian_regression_IntractLik_m5$skeleton_T, skeleton_Xi = skeleton_Gaussian_regression_IntractLik_m5$skeleton_Xi, skeleton_Theta = skeleton_Gaussian_regression_IntractLik_m5$skeleton_Theta, N = 50000)




###### m = 10

m <- 10

prob_one <- prob_all^{1/m}
prob_one

sqrt(-2*log(punif(q=prob_one, 0, 1)))

R <- sqrt(-2*log(punif(q=prob_one, 0, 1)))

R

T_end <- 10000
xi_0 <- c(beta_0, log(1))
theta_0 <- rep(1, p+1)


timer_MMD_Gaussian_regression_m10.start <- Sys.time()
skeleton_Gaussian_regression_IntractLik_m10 <- ZigZag_MMD_RBF_Gaussian_regression_cpp_arma(T_end, xi_0, theta_0, y, X, m, gamma, m_0, s_0, a_0, b_0, R, w = sqrt(2*pi), N_skeleton = 1000000)
timer_MMD_Gaussian_regression_m10.end <- Sys.time()
timer_MMD_Gaussian_regression_m10 = timer_MMD_Gaussian_regression_m10.end - timer_MMD_Gaussian_regression_m10.start


sample_Gaussian_regression_IntractLik_m10 <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_Gaussian_regression_IntractLik_m10$skeleton_T, skeleton_Xi = skeleton_Gaussian_regression_IntractLik_m10$skeleton_Xi, skeleton_Theta = skeleton_Gaussian_regression_IntractLik_m10$skeleton_Theta, N = 50000)





###### m = 20

m <- 20

prob_one <- prob_all^{1/m}
prob_one

sqrt(-2*log(punif(q=prob_one, 0, 1)))

R <- sqrt(-2*log(punif(q=prob_one, 0, 1)))

R

T_end <- 10000
xi_0 <- c(beta_0, log(1))
theta_0 <- rep(1, p+1)


timer_MMD_Gaussian_regression_m20.start <- Sys.time()
skeleton_Gaussian_regression_IntractLik_m20 <- ZigZag_MMD_RBF_Gaussian_regression_cpp_arma(T_end, xi_0, theta_0, y, X, m, gamma, m_0, s_0, a_0, b_0, R, w = sqrt(2*pi), N_skeleton = 1000000)
timer_MMD_Gaussian_regression_m20.end <- Sys.time()
timer_MMD_Gaussian_regression_m20 = timer_MMD_Gaussian_regression_m20.end - timer_MMD_Gaussian_regression_m20.start

sample_Gaussian_regression_IntractLik_m20 <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_Gaussian_regression_IntractLik_m20$skeleton_T, skeleton_Xi = skeleton_Gaussian_regression_IntractLik_m20$skeleton_Xi, skeleton_Theta = skeleton_Gaussian_regression_IntractLik_m20$skeleton_Theta, N = 50000)




###### m = 50

m <- 50

prob_one <- prob_all^{1/m}
prob_one

sqrt(-2*log(punif(q=prob_one, 0, 1)))

R <- sqrt(-2*log(punif(q=prob_one, 0, 1)))

R

T_end <- 10000
xi_0 <- c(beta_0, log(1))
theta_0 <- rep(1, p+1)


timer_MMD_Gaussian_regression_m50.start <- Sys.time()
skeleton_Gaussian_regression_IntractLik_m50 <- ZigZag_MMD_RBF_Gaussian_regression_cpp_arma(T_end, xi_0, theta_0, y, X, m, gamma, m_0, s_0, a_0, b_0, R, w = sqrt(2*pi), N_skeleton = 1000000)
timer_MMD_Gaussian_regression_m50.end <- Sys.time()
timer_MMD_Gaussian_regression_m50 = timer_MMD_Gaussian_regression_m50.end - timer_MMD_Gaussian_regression_m50.start

sample_Gaussian_regression_IntractLik_m50 <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_Gaussian_regression_IntractLik_m50$skeleton_T, skeleton_Xi = skeleton_Gaussian_regression_IntractLik_m50$skeleton_Xi, skeleton_Theta = skeleton_Gaussian_regression_IntractLik_m50$skeleton_Theta, N = 50000)






save(file = "ResultsZZ.RData", skeleton_Gaussian_regression_IntractLik_m2, skeleton_Gaussian_regression_IntractLik_m5,skeleton_Gaussian_regression_IntractLik_m10,skeleton_Gaussian_regression_IntractLik_m20, skeleton_Gaussian_regression_IntractLik_m50, sample_Gaussian_regression_IntractLik_m50, timer_MMD_Gaussian_regression_m50, sample_Gaussian_regression_IntractLik_m20, timer_MMD_Gaussian_regression_m20, sample_Gaussian_regression_IntractLik_m10, timer_MMD_Gaussian_regression_m10, sample_Gaussian_regression_IntractLik_m5, timer_MMD_Gaussian_regression_m5, sample_Gaussian_regression_IntractLik_m2, timer_MMD_Gaussian_regression_m2)







###### m = 5 (Gold Standard Run)

m <- 5

prob_one <- prob_all^{1/m}
prob_one

sqrt(-2*log(punif(q=prob_one, 0, 1)))

R <- sqrt(-2*log(punif(q=prob_one, 0, 1)))

R

T_end <- 40000
xi_0 <- c(c(4, 4, 3, 3, 2, 2, 1, 1), log(1))
theta_0 <- rep(1, p+1)


timer_MMD_Gaussian_regression_gold.start <- Sys.time()
skeleton_Gaussian_regression_IntractLik_gold <- ZigZag_MMD_RBF_Gaussian_regression_cpp_arma(T_end, xi_0, theta_0, y, X, m, gamma, m_0, s_0, a_0, b_0, R, w = sqrt(2*pi), N_skeleton = 1000000)
timer_MMD_Gaussian_regression_gold.end <- Sys.time()
timer_MMD_Gaussian_regression_gold = timer_MMD_Gaussian_regression_gold.end - timer_MMD_Gaussian_regression_gold.start


sample_Gaussian_regression_IntractLik_gold <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_Gaussian_regression_IntractLik_gold$skeleton_T, skeleton_Xi = skeleton_Gaussian_regression_IntractLik_gold$skeleton_Xi, skeleton_Theta = skeleton_Gaussian_regression_IntractLik_gold$skeleton_Theta, N = 200000)



save(file = "ResultsZZ_n100_Gold.RData",skeleton_Gaussian_regression_IntractLik_gold, sample_Gaussian_regression_IntractLik_gold, timer_MMD_Gaussian_regression_gold)








