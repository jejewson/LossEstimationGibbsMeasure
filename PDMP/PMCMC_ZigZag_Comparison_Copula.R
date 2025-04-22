

library(LaplacesDemon)

library(mvtnorm)



library(Rcpp)


sourceCpp("MMD_copula_indirect_arma.cpp")






# function for computing log loss
log_loss_mmd <- function(theta,uhat,m,w,gamma_mmd){
  
  z = matrix(rnorm(2*m), nrow = m, ncol = 2)
  mmd = MMD_biRBF_invPhi_Gaussian_Copula_cpp_arma(uhat, z, theta, gamma_mmd)
  
  log_prior = -theta - 2*log(1 + exp(-theta))
  
  return(-w*mmd + log_prior)
  
}


# function for computing log loss for a given set of random numbers
log_loss_mmd_block <- function(theta,z,uhat,m,w,gamma_mmd){
  
  mmd = MMD_biRBF_invPhi_Gaussian_Copula_cpp_arma(uhat, z, theta, gamma_mmd)
  
  log_prior = -theta - 2*log(1 + exp(-theta))
  
  return(-w*mmd + log_prior)
  
}





# function for standard PMCMC algorithm
pmcmc <- function(theta_init,uhat,m,w,gamma_mmd,sd_rw,M){
  
  theta_curr = theta_init
  
  logpost_curr = log_loss_mmd(theta_curr,uhat,m,w,gamma_mmd)
  
  theta = matrix(0,M,1)
  
  for (i in 1:M){
    
    theta_prop = rnorm(n = 1, theta_curr, sd_rw);
    
    logpost_prop = log_loss_mmd(theta_prop,uhat,m,w,gamma_mmd)
    
    if (exp(logpost_prop - logpost_curr) > runif(1)){
      theta_curr = theta_prop
      logpost_curr = logpost_prop
    }
    theta[i,] = theta_curr
    
  }
  
  return(theta)
  
}


# function for block correlated PMCMC algorithm
pmcmc_block <- function(theta_init,uhat,m,w,gamma_mmd,sd_rw,M){
  
  theta_curr = theta_init
  
  z_curr = matrix(rnorm(2*m), nrow = m, ncol = 2)  # initialise random numbers for generating pseudo datasets
  
  logpost_curr = log_loss_mmd_block(theta_curr,z_curr,uhat,m,w,gamma_mmd)
  
  theta = matrix(0,M,1)
  
  for (i in 1:M){
    
    theta_prop = rnorm(n = 1, theta_curr, sd_rw);
    
    # update random numbers for block r
    r_prop = sample(m,1)
    z_prop = z_curr
    z_prop[r_prop,] = rnorm(2,0,1)
    
    logpost_prop = log_loss_mmd_block(theta_prop,z_prop,uhat,m,w,gamma_mmd)
    
    if (exp(logpost_prop - logpost_curr) > runif(1)){
      theta_curr = theta_prop
      logpost_curr = logpost_prop
      z_curr = z_prop
    }
    theta[i,] = theta_curr
    
  }
  
  return(theta)
  
}






####### Generate dataset of size n = 100





## Data


set.seed(1)    # set seed for reproducibility


N_rep <- 1


rho <- 0.5 
p <- 2


ndata <- 100


y <- rmvnorm(ndata, rep(0, 2), diag(1 - rho, p) + matrix(rho, nrow = p, ncol = p))








u_hat <- matrix(NA, nrow = ndata, ncol = p)
for(j in 1:p){
  for(i in 1:ndata){
    u_hat[i,j] <- length(which(y[,j] <= y[i,j]))/ndata
  }
  
}




plot(u_hat[,1], u_hat[,2])




w = 1


gamma_mmd <- 0.25 # hyperparameter for MMD






###### standard pseudo-marginal


load(file = "cov_rw_copula_n100_zigzag.RData") # load in posterior sd from zig zag to use in the proposal dist for pmcmc


m = 2
start.time <- Sys.time()
samples_pmcmc_m2 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,5000000)
end.time <- Sys.time()
time.taken.standard_m2 <- end.time - start.time




m = 5
start.time <- Sys.time()
samples_pmcmc_m5 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,2000000)
end.time <- Sys.time()
time.taken.standard_m5 <- end.time - start.time




m = 10
start.time <- Sys.time()
samples_pmcmc_m10 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,1000000)
end.time <- Sys.time()
time.taken.standard_m10 <- end.time - start.time




m = 20
start.time <- Sys.time()
samples_pmcmc_m20 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,500000)
end.time <- Sys.time()
time.taken.standard_m20 <- end.time - start.time




m = 50
start.time <- Sys.time()
samples_pmcmc_m50 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,200000)
end.time <- Sys.time()
time.taken.standard_m50 <- end.time - start.time






save(file = "ResultsPM_copula_n100.RData", samples_pmcmc_m50, time.taken.standard_m50, samples_pmcmc_m20, time.taken.standard_m20, samples_pmcmc_m10, time.taken.standard_m10, samples_pmcmc_m5, time.taken.standard_m5, samples_pmcmc_m2, time.taken.standard_m2)







###### block pseudo-marginal


load(file = "cov_rw_copula_n100_zigzag.RData") # load in posterior sd from zig zag to use in the proposal dist for pmcmc


m = 2
start.time <- Sys.time()
samples_pmcmc_block_m2 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,5000000)
end.time <- Sys.time()
time.taken.block_m2 <- end.time - start.time




m = 5
start.time <- Sys.time()
samples_pmcmc_block_m5 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,2000000)
end.time <- Sys.time()
time.taken.block_m5 <- end.time - start.time




m = 10
start.time <- Sys.time()
samples_pmcmc_block_m10 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,1000000)
end.time <- Sys.time()
time.taken.block_m10 <- end.time - start.time




m = 20
start.time <- Sys.time()
samples_pmcmc_block_m20 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,500000)
end.time <- Sys.time()
time.taken.block_m20 <- end.time - start.time




m = 50
start.time <- Sys.time()
samples_pmcmc_block_m50 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,200000)
end.time <- Sys.time()
time.taken.block_m50 <- end.time - start.time






save(file = "ResultsPM_copula_block_n100.RData", samples_pmcmc_block_m50, time.taken.block_m50, samples_pmcmc_block_m20, time.taken.block_m20, samples_pmcmc_block_m10, time.taken.block_m10, samples_pmcmc_block_m5, time.taken.block_m5, samples_pmcmc_block_m2, time.taken.block_m2)






###### zig zag sampling



m <- 2


#T_end <- 100
T_end <- 4000
xi_0 <- 0
theta_0 <- 1


gamma <- 0.25


a <- 1
b <- 1


R <- 3


timer_MMD_copula_m2_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m2_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 100000)
timer_MMD_copula_m2_rcpp.end <- Sys.time()
timer_MMD_copula_m2_rcpp = timer_MMD_copula_m2_rcpp.end - timer_MMD_copula_m2_rcpp.start


sample_copula_IntractLik_m2_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m2_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m2_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m2_rcpp$skeleton_Theta, N = 20000)


















m <- 5


#T_end <- 100
T_end <- 4000
xi_0 <- 0
theta_0 <- 1

gamma <- 0.25


a <- 1
b <- 1


R <- 3


timer_MMD_copula_m5_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m5_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 100000)
timer_MMD_copula_m5_rcpp.end <- Sys.time()
timer_MMD_copula_m5_rcpp = timer_MMD_copula_m5_rcpp.end - timer_MMD_copula_m5_rcpp.start


sample_copula_IntractLik_m5_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m5_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m5_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m5_rcpp$skeleton_Theta, N = 20000)












m <- 10


#T_end <- 100
T_end <- 4000
xi_0 <- 0
theta_0 <- 1


gamma <- 0.25


a <- 1
b <- 1


R <- 3


timer_MMD_copula_m10_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m10_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 100000)
timer_MMD_copula_m10_rcpp.end <- Sys.time()
timer_MMD_copula_m10_rcpp = timer_MMD_copula_m10_rcpp.end - timer_MMD_copula_m10_rcpp.start


sample_copula_IntractLik_m10_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m10_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m10_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m10_rcpp$skeleton_Theta, N = 20000)












m <- 20


#T_end <- 100
T_end <- 4000
xi_0 <- 0
theta_0 <- 1


gamma <- 0.25


a <- 1
b <- 1


R <- 3


timer_MMD_copula_m20_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m20_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 100000)
timer_MMD_copula_m20_rcpp.end <- Sys.time()
timer_MMD_copula_m20_rcpp = timer_MMD_copula_m20_rcpp.end - timer_MMD_copula_m20_rcpp.start


sample_copula_IntractLik_m20_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m20_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m20_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m20_rcpp$skeleton_Theta, N = 20000)








m <- 50


#T_end <- 100
T_end <- 4000
xi_0 <- 0
theta_0 <- 1


gamma <- 0.25


a <- 1
b <- 1


R <- 3


timer_MMD_copula_m50_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m50_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 100000)
timer_MMD_copula_m50_rcpp.end <- Sys.time()
timer_MMD_copula_m50_rcpp = timer_MMD_copula_m50_rcpp.end - timer_MMD_copula_m50_rcpp.start


sample_copula_IntractLik_m50_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m50_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m50_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m50_rcpp$skeleton_Theta, N = 20000)










save(file = "ResultsZZ_copula_n100.RData", sample_copula_IntractLik_m50_rcpp, skeleton_copula_IntractLik_m50_rcpp,timer_MMD_copula_m50_rcpp,
     sample_copula_IntractLik_m20_rcpp, skeleton_copula_IntractLik_m20_rcpp,timer_MMD_copula_m20_rcpp,
     sample_copula_IntractLik_m10_rcpp, skeleton_copula_IntractLik_m10_rcpp,timer_MMD_copula_m10_rcpp,
     sample_copula_IntractLik_m5_rcpp, skeleton_copula_IntractLik_m5_rcpp,timer_MMD_copula_m5_rcpp,
     sample_copula_IntractLik_m2_rcpp, skeleton_copula_IntractLik_m2_rcpp,timer_MMD_copula_m2_rcpp)










###### ZZ Gold Standard Run


m <- 5


#T_end <- 100
T_end <- 40000
rho = 0.5
xi_0 <- -log((1-rho)/(1+rho))
theta_0 <- 1


gamma <- 0.25


a <- 1
b <- 1


R <- 3


timer_MMD_copula_gold_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_gold_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 100000)
timer_MMD_copula_gold_rcpp.end <- Sys.time()
timer_MMD_copula_gold_rcpp = timer_MMD_copula_gold_rcpp.end - timer_MMD_copula_gold_rcpp.start


sample_copula_IntractLik_gold_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_gold_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_gold_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_gold_rcpp$skeleton_Theta, N = 200000)





save(file = "ResultsZZ_copula_n100_Gold.RData", sample_copula_IntractLik_gold_rcpp, skeleton_copula_IntractLik_gold_rcpp,timer_MMD_copula_gold_rcpp)















####### Generate dataset of size n = 1000




## Data

set.seed(1)

N_rep <- 1

rho <- 0.5 
p <- 2

ndata <- 1000

y <- rmvnorm(ndata, rep(0, 2), diag(1 - rho, p) + matrix(rho, nrow = p, ncol = p))




u_hat <- matrix(NA, nrow = ndata, ncol = p)
for(j in 1:p){
  for(i in 1:ndata){
    u_hat[i,j] <- length(which(y[,j] <= y[i,j]))/ndata
  }
  
}


plot(u_hat[,1], u_hat[,2])


w = 1

gamma_mmd <- 0.25 # hyperparameter for MMD





###### standard pseudo-marginal

load(file = "cov_rw_copula_n1000_zigzag.RData") # load in posterior sd from zig zag to use in the proposal dist for pmcmc

m = 2
start.time <- Sys.time()
samples_pmcmc_m2 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,500000)
end.time <- Sys.time()
time.taken.standard_m2 <- end.time - start.time


m = 5
start.time <- Sys.time()
samples_pmcmc_m5 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,200000)
end.time <- Sys.time()
time.taken.standard_m5 <- end.time - start.time


m = 10
start.time <- Sys.time()
samples_pmcmc_m10 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,100000)
end.time <- Sys.time()
time.taken.standard_m10 <- end.time - start.time


m = 20
start.time <- Sys.time()
samples_pmcmc_m20 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,50000)
end.time <- Sys.time()
time.taken.standard_m20 <- end.time - start.time


m = 50
start.time <- Sys.time()
samples_pmcmc_m50 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,20000)
end.time <- Sys.time()
time.taken.standard_m50 <- end.time - start.time



m = 100
start.time <- Sys.time()
samples_pmcmc_m100 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,10000)
end.time <- Sys.time()
time.taken.standard_m100 <- end.time - start.time



m = 200
start.time <- Sys.time()
samples_pmcmc_m200 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,5000)
end.time <- Sys.time()
time.taken.standard_m200 <- end.time - start.time


m = 500
start.time <- Sys.time()
samples_pmcmc_m500 = pmcmc(0.5,u_hat,m,w,gamma_mmd,sd_rw,2000)
end.time <- Sys.time()
time.taken.standard_m500 <- end.time - start.time



save(file = "ResultsPM_copula_n1000.RData",samples_pmcmc_m500, time.taken.standard_m500, samples_pmcmc_m200, time.taken.standard_m200, samples_pmcmc_m100, time.taken.standard_m100, samples_pmcmc_m50, time.taken.standard_m50, samples_pmcmc_m20, time.taken.standard_m20, samples_pmcmc_m10, time.taken.standard_m10, samples_pmcmc_m5, time.taken.standard_m5, samples_pmcmc_m2, time.taken.standard_m2)




# Block PMMH

load(file = "cov_rw_copula_n1000_zigzag.RData") # load in posterior sd from zig zag to use in the proposal dist for pmcmc



m = 5
start.time <- Sys.time()
samples_pmcmc_block_m5 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,4000000)
end.time <- Sys.time()
time.taken.block_m5 <- end.time - start.time


m = 10
start.time <- Sys.time()
samples_pmcmc_block_m10 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,2000000)
end.time <- Sys.time()
time.taken.block_m10 <- end.time - start.time



m = 20
start.time <- Sys.time()
samples_pmcmc_block_m20 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,1000000)
end.time <- Sys.time()
time.taken.block_m20 <- end.time - start.time




m = 50
start.time <- Sys.time()
samples_pmcmc_block_m50 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,400000)
end.time <- Sys.time()
time.taken.block_m50 <- end.time - start.time


m = 100
start.time <- Sys.time()
samples_pmcmc_block_m100 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,200000)
end.time <- Sys.time()
time.taken.block_m100 <- end.time - start.time



save(file = "ResultsPM_copula_block_n1000.RData", samples_pmcmc_block_m100, time.taken.block_m100, samples_pmcmc_block_m50, time.taken.block_m50, samples_pmcmc_block_m20, time.taken.block_m20, samples_pmcmc_block_m10, time.taken.block_m10, samples_pmcmc_block_m5, time.taken.block_m5)




m = 1000
start.time <- Sys.time()
samples_pmcmc_block_m1000 = pmcmc_block(0.5,u_hat,m,w,gamma_mmd,sd_rw,20000)
end.time <- Sys.time()
time.taken.block_m1000 <- end.time - start.time



save(file = "ResultsPM_copula_block_n1000_m1000.RData", samples_pmcmc_block_m1000, time.taken.block_m1000)



####### Zig Zag


###### m = 2


m <- 2

#T_end <- 100
T_end <- 4000
rho = 0.5
xi_0 <- -log((1-rho)/(1+rho))
theta_0 <- 1

gamma <- 0.25

a <- 1
b <- 1

R <- 3


timer_MMD_copula_m2_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m2_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 200000)
timer_MMD_copula_m2_rcpp.end <- Sys.time()
timer_MMD_copula_m2_rcpp = timer_MMD_copula_m2_rcpp.end - timer_MMD_copula_m2_rcpp.start

sample_copula_IntractLik_m2_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m2_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m2_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m2_rcpp$skeleton_Theta, N = 2000)









###### m = 5


m <- 5

#T_end <- 100
T_end <- 4000
rho = 0.5
xi_0 <- -log((1-rho)/(1+rho))
theta_0 <- 1

gamma <- 0.25

a <- 1
b <- 1

R <- 3

timer_MMD_copula_m5_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m5_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 200000)
timer_MMD_copula_m5_rcpp.end <- Sys.time()
timer_MMD_copula_m5_rcpp = timer_MMD_copula_m5_rcpp.end - timer_MMD_copula_m5_rcpp.start

sample_copula_IntractLik_m5_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m5_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m5_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m5_rcpp$skeleton_Theta, N = 2000)





###### m = 10


m <- 10

#T_end <- 100
T_end <- 4000
rho = 0.5
xi_0 <- -log((1-rho)/(1+rho))
theta_0 <- 1

gamma <- 0.25

a <- 1
b <- 1

R <- 3


timer_MMD_copula_m10_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m10_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 200000)
timer_MMD_copula_m10_rcpp.end <- Sys.time()
timer_MMD_copula_m10_rcpp = timer_MMD_copula_m10_rcpp.end - timer_MMD_copula_m10_rcpp.start

sample_copula_IntractLik_m10_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m10_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m10_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m10_rcpp$skeleton_Theta, N = 2000)







###### m = 20


m <- 20

#T_end <- 100
T_end <- 4000
rho = 0.5
xi_0 <- -log((1-rho)/(1+rho))
theta_0 <- 1

gamma <- 0.25

a <- 1
b <- 1

R <- 3


timer_MMD_copula_m20_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m20_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 200000)
timer_MMD_copula_m20_rcpp.end <- Sys.time()
timer_MMD_copula_m20_rcpp = timer_MMD_copula_m20_rcpp.end - timer_MMD_copula_m20_rcpp.start

sample_copula_IntractLik_m20_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m20_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m20_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m20_rcpp$skeleton_Theta, N = 2000)






###### m = 50


m <- 50

#T_end <- 100
T_end <- 4000
rho = 0.5
xi_0 <- -log((1-rho)/(1+rho))
theta_0 <- 1

gamma <- 0.25

a <- 1
b <- 1

R <- 3

timer_MMD_copula_m50_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_m50_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 200000)
timer_MMD_copula_m50_rcpp.end <- Sys.time()
timer_MMD_copula_m50_rcpp = timer_MMD_copula_m50_rcpp.end - timer_MMD_copula_m50_rcpp.start

sample_copula_IntractLik_m50_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_m50_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_m50_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_m50_rcpp$skeleton_Theta, N = 2000)






save(file = "ResultsZZ_copula_n1000.RData", sample_copula_IntractLik_m50_rcpp, skeleton_copula_IntractLik_m50_rcpp,timer_MMD_copula_m50_rcpp,
     sample_copula_IntractLik_m20_rcpp, skeleton_copula_IntractLik_m20_rcpp,timer_MMD_copula_m20_rcpp,
     sample_copula_IntractLik_m10_rcpp, skeleton_copula_IntractLik_m10_rcpp,timer_MMD_copula_m10_rcpp,
     sample_copula_IntractLik_m5_rcpp, skeleton_copula_IntractLik_m5_rcpp,timer_MMD_copula_m5_rcpp,
     sample_copula_IntractLik_m2_rcpp, skeleton_copula_IntractLik_m2_rcpp,timer_MMD_copula_m2_rcpp)






###### m = 5 LONG RUN


m <- 5

#T_end <- 100
T_end <- 20000
rho = 0.5
xi_0 <- -log((1-rho)/(1+rho))
theta_0 <- 1

gamma <- 0.25

a <- 1
b <- 1

R <- 3

timer_MMD_copula_gold_rcpp.start <- Sys.time()
skeleton_copula_IntractLik_gold_rcpp <- ZigZag_MMD_RBF_biGaussian_copula_cpp_arma(T_end, array(xi_0, dim = c(1)), array(theta_0, dim = c(1)), u_hat, m, gamma, R, w = 1, a, b, N_skeleton = 1000000)
timer_MMD_copula_gold_rcpp.end <- Sys.time()
timer_MMD_copula_gold_rcpp = timer_MMD_copula_gold_rcpp.end - timer_MMD_copula_gold_rcpp.start

sample_copula_IntractLik_gold_rcpp <- skeleton_to_sample_cpp_arma(skeleton_T = skeleton_copula_IntractLik_gold_rcpp$skeleton_T, skeleton_Xi = skeleton_copula_IntractLik_gold_rcpp$skeleton_Xi, skeleton_Theta = skeleton_copula_IntractLik_gold_rcpp$skeleton_Theta, N = 200000)




save(file = "ResultsZZ_copula_n1000_Gold.RData", sample_copula_IntractLik_gold_rcpp, skeleton_copula_IntractLik_gold_rcpp, timer_MMD_copula_gold_rcpp)






