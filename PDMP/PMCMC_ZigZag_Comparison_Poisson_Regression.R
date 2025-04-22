

library(adaptMCMC)
library(actuar)
library(metRology)
library(LaplacesDemon)
library(invgamma)
library(MASS)
library(matrixStats)
library(pracma)

library(Rcpp)

sourceCpp("betaD_PoissonRegression_arma.cpp")
sourceCpp("rpois_inv.cpp")


# function for computing log loss
log_loss_betaD <- function(theta,y,X,m,w,beta){
  
  
  mu = exp(X%*%theta)
  dpois(y,mu,log=TRUE)
  
  
  second_term = (1 + 1/beta)*exp(logSumExp(beta*dpois(y,mu,log=TRUE)))
  
  
  betaD_all = rep(0,ndata)
  
  
  for (i in 1:ndata){
    z = rpois(m, mu[i])
    betaD_all[i] = exp(logSumExp(beta*dpois(z,mu[i],log=TRUE)))/m
  }
  
  log_prior = sum(dnorm(theta, m_0, s_0, log = TRUE)) 
  return(-w*(sum(betaD_all) - second_term) + log_prior)
  
}


# function for computing log loss for a given set of random numbers
log_loss_betaD_block <- function(theta,u,y,X,m,w,beta){
  
  
  mu = exp(X%*%theta)
  dpois(y,mu,log=TRUE)
  
  
  second_term = (1 + 1/beta)*exp(logSumExp(beta*dpois(y,mu,log=TRUE)))
  
  
  betaD_all = rep(0,ndata)
  
  
  for (i in 1:ndata){
    z = rpois_inverse_cpp(u[i,],mu[i])
    betaD_all[i] = exp(logSumExp(beta*dpois(z,mu[i],log=TRUE)))/m
  }
  
  log_prior = sum(dnorm(theta, m_0, s_0, log = TRUE)) 
  return(-w*(sum(betaD_all) - second_term) + log_prior)
  
}



# function for standard PMCMC algorithm
pmcmc <- function(theta_init,m,w,beta,cov_rw,M){
  
  theta_curr = theta_init
  
  logpost_curr = log_loss_betaD(theta_curr,y,X,m,w,beta)
  
  theta = matrix(0,M,p)
  
  for (i in 1:M){
    
    theta_prop = mvrnorm(n = 1, theta_curr, cov_rw);
    
    logpost_prop = log_loss_betaD(theta_prop,y,X,m,w,beta)
    
    if (exp(logpost_prop - logpost_curr) > runif(1)){
      theta_curr = theta_prop
      logpost_curr = logpost_prop
    }
    theta[i,] = theta_curr
    
  }
  
  return(theta)
  
}


#### this function now implemented in Rcpp
#rpois_inverse <- function(u, lambda){
  
#  n = length(u)
#  x = rep(0,n)
#  for (i in 1:n){
#    s = 0
#    while(TRUE){
#      if (u[i] <= ppois(s, lambda)){
#        x[i] = s
#        break
#      }else{
#        s = s+1
#      }
#    }
#  }
#  return(x)
#}


# function for block correlated PMCMC algorithm
pmcmc_block <- function(theta_init,m,w,beta,cov_rw,M){
  
  theta_curr = theta_init
  
  u_curr = matrix(runif(ndata*m), nrow = ndata, ncol = m)  # initialise random numbers for generating pseudo datasets
  
  
  logpost_curr = log_loss_betaD_block(theta_curr,u_curr,y,X,m,w,beta)
  
  theta = matrix(0,M,p)
  
  for (i in 1:M){
    
    theta_prop = mvrnorm(n = 1, theta_curr, cov_rw);
    
    # update random numbers for block r
    r_prop = sample(ndata,1)
    u_prop = u_curr
    u_prop[r_prop,] = runif(m)
    
    logpost_prop = log_loss_betaD_block(theta_prop,u_prop,y,X,m,w,beta)
    
    if (exp(logpost_prop - logpost_curr) > runif(1)){
      theta_curr = theta_prop
      u_curr = u_prop
      logpost_curr = logpost_prop
    }
    theta[i,] = theta_curr
    
  }
  
  return(theta)
  
}







########## Generate dataset of size n = 1000


set.seed(1)   # set seed for reproducibility



p <- 5

theta <- c(1, 0.5, 1.5, 0, 0)

ndata <- 1000

X <- cbind(1, matrix(rnorm(ndata*(p-1), 0, 0.25), nrow = ndata, ncol = p-1))

y <- rpois(ndata, exp(X%*%theta))




# prior hyperparameters
m_0 <- 0
s_0 <- 1






######### run standard pseudo-marginal sampler


load(file = "cov_rw_poisson_n1000_zigzag.RData")   # load in tuned random walk proposal and initial values



m = 10
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_pmcmc_m10 = pmcmc(theta,m,w,beta,cov_rw,500000)
end.time <- Sys.time()
time.taken.standard_m10 <- end.time - start.time


m = 20
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_pmcmc_m20 = pmcmc(theta,m,w,beta,cov_rw,250000)
end.time <- Sys.time()
time.taken.standard_m20 <- end.time - start.time


m = 50
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_pmcmc_m50 = pmcmc(theta,m,w,beta,cov_rw,100000)
end.time <- Sys.time()
time.taken.standard_m50 <- end.time - start.time



m = 100
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_pmcmc_m100 = pmcmc(theta,m,w,beta,cov_rw,50000)
end.time <- Sys.time()
time.taken.standard_m100 <- end.time - start.time




save(file = "ResultsPM_Poisson_n1000.RData", samples_pmcmc_m100, time.taken.standard_m100, samples_pmcmc_m50, time.taken.standard_m50, samples_pmcmc_m20, time.taken.standard_m20, samples_pmcmc_m10, time.taken.standard_m10)








######### run block pseudo-marginal sampler


load(file = "cov_rw_poisson_n1000_zigzag.RData")   # load in tuned random walk proposal and initial values



m = 2
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_block_pmcmc_m2 = pmcmc_block(theta,m,w,beta,cov_rw,2500000)
end.time <- Sys.time()
time.taken.block_m2 <- end.time - start.time


save(file = "ResultsPM_Poisson_block_n1000_m2.RData", samples_block_pmcmc_m2, time.taken.block_m2)




m = 5
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_block_pmcmc_m5 = pmcmc_block(theta,m,w,beta,cov_rw,1000000)
end.time <- Sys.time()
time.taken.block_m5 <- end.time - start.time


save(file = "ResultsPM_Poisson_block_n1000_m5.RData", samples_block_pmcmc_m5, time.taken.block_m5)


m = 10
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_block_pmcmc_m10 = pmcmc_block(theta,m,w,beta,cov_rw,500000)
end.time <- Sys.time()
time.taken.block_m10 <- end.time - start.time



save(file = "ResultsPM_Poisson_block_n1000_m10.RData", samples_block_pmcmc_m10, time.taken.block_m10)



m = 20
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_block_pmcmc_m20 = pmcmc_block(theta,m,w,beta,cov_rw,250000)
end.time <- Sys.time()
time.taken.block_m20 <- end.time - start.time


save(file = "ResultsPM_Poisson_block_n1000_m20.RData", samples_block_pmcmc_m20, time.taken.block_m20)




m = 50
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_block_pmcmc_m50 = pmcmc_block(theta,m,w,beta,cov_rw,100000)
end.time <- Sys.time()
time.taken.block_m50 <- end.time - start.time


save(file = "ResultsPM_Poisson_block_n1000_m50.RData", samples_block_pmcmc_m50, time.taken.block_m50)






m = 100
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_block_pmcmc_m100 = pmcmc_block(theta,m,w,beta,cov_rw,50000)
end.time <- Sys.time()
time.taken.block_m100 <- end.time - start.time


save(file = "ResultsPM_Poisson_block_n1000_m100.RData", samples_block_pmcmc_m100, time.taken.block_m100)





m = 200
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_block_pmcmc_m200 = pmcmc_block(theta,m,w,beta,cov_rw,25000)
end.time <- Sys.time()
time.taken.block_m200 <- end.time - start.time


save(file = "ResultsPM_Poisson_block_n1000_m200.RData", samples_block_pmcmc_m200, time.taken.block_m200)





m = 500
beta <- 0.5
w = 1
start.time <- Sys.time()
samples_block_pmcmc_m500 = pmcmc_block(theta,m,w,beta,cov_rw,10000)
end.time <- Sys.time()
time.taken.block_m500 <- end.time - start.time


save(file = "ResultsPM_Poisson_block_n1000_m500.RData", samples_block_pmcmc_m500, time.taken.block_m500)






#### zig zag sampling




m <- 2

sum_abs_X_max <- max(rowSums(abs(X)))

T_end <- 2000*sqrt(100)/sqrt(ndata)
#T_end <- 5000
xi_0 <- c(1, 0.5, 1.5, 0, 0)
theta_0 <- rep(1, p)

beta <- 1.5


timer_betaD_PoissonRegression_m2_rcpp.start <- Sys.time()
skeleton_betaD_PoissonRegression_IntractLik_m2_rcpp <- ZigZag_betaD_poisRegression_cpp_arma(T_end, xi_0, theta_0, y, X, m, beta, sum_abs_X_max, w = 1.0, m_0, s_0, N_skeleton = 200000)
timer_betaD_PoissonRegression_m2_rcpp.end <- Sys.time()
timer_betaD_PoissonRegression_m2_rcpp = timer_betaD_PoissonRegression_m2_rcpp.end - timer_betaD_PoissonRegression_m2_rcpp.start


sample_betaD_PoissonRegression_m2_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_betaD_PoissonRegression_IntractLik_m2_rcpp$skeleton_T, skeleton_Xi = skeleton_betaD_PoissonRegression_IntractLik_m2_rcpp$skeleton_Xi, skeleton_Theta = skeleton_betaD_PoissonRegression_IntractLik_m2_rcpp$skeleton_Theta, N = 2000)


save(file = "ResultsZZ_Poisson_n1000_m2.RData", sample_betaD_PoissonRegression_m2_rcpp, skeleton_betaD_PoissonRegression_IntractLik_m2_rcpp,timer_betaD_PoissonRegression_m2_rcpp)










m <- 5

sum_abs_X_max <- max(rowSums(abs(X)))

T_end <- 2000*sqrt(100)/sqrt(ndata)
#T_end <- 5000
xi_0 <- c(1, 0.5, 1.5, 0, 0)
theta_0 <- rep(1, p)

beta <- 1.5


timer_betaD_PoissonRegression_m5_rcpp.start <- Sys.time()
skeleton_betaD_PoissonRegression_IntractLik_m5_rcpp <- ZigZag_betaD_poisRegression_cpp_arma(T_end, xi_0, theta_0, y, X, m, beta, sum_abs_X_max, w = 1.0, m_0, s_0, N_skeleton = 200000)
timer_betaD_PoissonRegression_m5_rcpp.end <- Sys.time()
timer_betaD_PoissonRegression_m5_rcpp = timer_betaD_PoissonRegression_m5_rcpp.end - timer_betaD_PoissonRegression_m5_rcpp.start


sample_betaD_PoissonRegression_m5_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_betaD_PoissonRegression_IntractLik_m5_rcpp$skeleton_T, skeleton_Xi = skeleton_betaD_PoissonRegression_IntractLik_m5_rcpp$skeleton_Xi, skeleton_Theta = skeleton_betaD_PoissonRegression_IntractLik_m5_rcpp$skeleton_Theta, N = 2000)


save(file = "ResultsZZ_Poisson_n1000_m5.RData", sample_betaD_PoissonRegression_m5_rcpp, skeleton_betaD_PoissonRegression_IntractLik_m5_rcpp,timer_betaD_PoissonRegression_m5_rcpp)








m <- 10

sum_abs_X_max <- max(rowSums(abs(X)))

T_end <- 2000*sqrt(100)/sqrt(ndata)
#T_end <- 5000
xi_0 <- c(1, 0.5, 1.5, 0, 0)
theta_0 <- rep(1, p)

beta <- 1.5


timer_betaD_PoissonRegression_m10_rcpp.start <- Sys.time()
skeleton_betaD_PoissonRegression_IntractLik_m10_rcpp <- ZigZag_betaD_poisRegression_cpp_arma(T_end, xi_0, theta_0, y, X, m, beta, sum_abs_X_max, w = 1.0, m_0, s_0, N_skeleton = 200000)
timer_betaD_PoissonRegression_m10_rcpp.end <- Sys.time()
timer_betaD_PoissonRegression_m10_rcpp = timer_betaD_PoissonRegression_m10_rcpp.end - timer_betaD_PoissonRegression_m10_rcpp.start


sample_betaD_PoissonRegression_m10_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_betaD_PoissonRegression_IntractLik_m10_rcpp$skeleton_T, skeleton_Xi = skeleton_betaD_PoissonRegression_IntractLik_m10_rcpp$skeleton_Xi, skeleton_Theta = skeleton_betaD_PoissonRegression_IntractLik_m10_rcpp$skeleton_Theta, N = 2000)


save(file = "ResultsZZ_Poisson_n1000_m10.RData", sample_betaD_PoissonRegression_m10_rcpp, skeleton_betaD_PoissonRegression_IntractLik_m10_rcpp,timer_betaD_PoissonRegression_m10_rcpp)







m <- 20

sum_abs_X_max <- max(rowSums(abs(X)))

T_end <- 2000*sqrt(100)/sqrt(ndata)
#T_end <- 5000
xi_0 <- c(1, 0.5, 1.5, 0, 0)
theta_0 <- rep(1, p)

beta <- 1.5


timer_betaD_PoissonRegression_m20_rcpp.start <- Sys.time()
skeleton_betaD_PoissonRegression_IntractLik_m20_rcpp <- ZigZag_betaD_poisRegression_cpp_arma(T_end, xi_0, theta_0, y, X, m, beta, sum_abs_X_max, w = 1.0, m_0, s_0, N_skeleton = 200000)
timer_betaD_PoissonRegression_m20_rcpp.end <- Sys.time()
timer_betaD_PoissonRegression_m20_rcpp = timer_betaD_PoissonRegression_m20_rcpp.end - timer_betaD_PoissonRegression_m20_rcpp.start


sample_betaD_PoissonRegression_m20_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_betaD_PoissonRegression_IntractLik_m20_rcpp$skeleton_T, skeleton_Xi = skeleton_betaD_PoissonRegression_IntractLik_m20_rcpp$skeleton_Xi, skeleton_Theta = skeleton_betaD_PoissonRegression_IntractLik_m20_rcpp$skeleton_Theta, N = 2000)



save(file = "ResultsZZ_Poisson_n1000_m20.RData", sample_betaD_PoissonRegression_m20_rcpp, skeleton_betaD_PoissonRegression_IntractLik_m20_rcpp,timer_betaD_PoissonRegression_m20_rcpp)






m <- 50

sum_abs_X_max <- max(rowSums(abs(X)))

T_end <- 2000*sqrt(100)/sqrt(ndata)
#T_end <- 5000
xi_0 <- c(1, 0.5, 1.5, 0, 0)
theta_0 <- rep(1, p)

beta <- 1.5


timer_betaD_PoissonRegression_m50_rcpp.start <- Sys.time()
skeleton_betaD_PoissonRegression_IntractLik_m50_rcpp <- ZigZag_betaD_poisRegression_cpp_arma(T_end, xi_0, theta_0, y, X, m, beta, sum_abs_X_max, w = 1.0, m_0, s_0, N_skeleton = 200000)
timer_betaD_PoissonRegression_m50_rcpp.end <- Sys.time()
timer_betaD_PoissonRegression_m50_rcpp = timer_betaD_PoissonRegression_m50_rcpp.end - timer_betaD_PoissonRegression_m50_rcpp.start


sample_betaD_PoissonRegression_m50_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_betaD_PoissonRegression_IntractLik_m50_rcpp$skeleton_T, skeleton_Xi = skeleton_betaD_PoissonRegression_IntractLik_m50_rcpp$skeleton_Xi, skeleton_Theta = skeleton_betaD_PoissonRegression_IntractLik_m50_rcpp$skeleton_Theta, N = 2000)


save(file = "ResultsZZ_Poisson_n1000_m50.RData", sample_betaD_PoissonRegression_m50_rcpp, skeleton_betaD_PoissonRegression_IntractLik_m50_rcpp,timer_betaD_PoissonRegression_m50_rcpp)











#### Zig zag (Gold standard run)


m <- 5

sum_abs_X_max <- max(rowSums(abs(X)))

T_end <- 10000*sqrt(100)/sqrt(ndata)
#T_end <- 5000
xi_0 <- c(1, 0.5, 1.5, 0, 0)
theta_0 <- rep(1, p)

beta <- 1.5


timer_betaD_PoissonRegression_gold_rcpp.start <- Sys.time()
skeleton_betaD_PoissonRegression_IntractLik_gold_rcpp <- ZigZag_betaD_poisRegression_cpp_arma(T_end, xi_0, theta_0, y, X, m, beta, sum_abs_X_max, w = 1.0, m_0, s_0, N_skeleton = 1000000)
timer_betaD_PoissonRegression_gold_rcpp.end <- Sys.time()
timer_betaD_PoissonRegression_gold_rcpp = timer_betaD_PoissonRegression_gold_rcpp.end - timer_betaD_PoissonRegression_gold_rcpp.start


sample_betaD_PoissonRegression_gold_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_betaD_PoissonRegression_IntractLik_gold_rcpp$skeleton_T, skeleton_Xi = skeleton_betaD_PoissonRegression_IntractLik_gold_rcpp$skeleton_Xi, skeleton_Theta = skeleton_betaD_PoissonRegression_IntractLik_gold_rcpp$skeleton_Theta, N = 10000)


save(file = "ResultsZZ_Poisson_n1000_Gold.RData", sample_betaD_PoissonRegression_gold_rcpp, skeleton_betaD_PoissonRegression_IntractLik_gold_rcpp,timer_betaD_PoissonRegression_gold_rcpp)









