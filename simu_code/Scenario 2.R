

# In this script, we consider the influence of model paramters on estimation.

rm(list=ls())
source("util_funcs2.R")
library(SpaCOAP)
width <- 20; height <- 30
n <- width*height
p=500
q = 5; d <- 40; k <- 3; r <- 3
N <- 100
methodNames <- c("SpaCOAP", "COAP", "PLNPCA", "MRRR", "FAST")
n_methods <- length(methodNames)
metricList <- list(F_tr = matrix(NA,N, n_methods), B_tr = matrix(NA,N, n_methods),
                   F_CCor=matrix(NA,N, n_methods), B_CCor=matrix(NA, N, n_methods),
                   alpha_norm2=matrix(NA,N, n_methods), beta_norm2 = matrix(NA,N, n_methods),
                   alpha_norm1 = matrix(NA, N, n_methods), beta_norm1 = matrix(NA, N, n_methods),
                   timeMat = matrix(NA, N, n_methods))

for(ii in seq_along(metricList)) colnames(metricList[[ii]]) <- methodNames
bandwidth <- 1
rho<-   c(8,0.6) #varying rho over c(8,0.6) ,c(12,0.6), c(12,1.5)#   
sigma2_eps=1
eta<-0.5
for(i in 1: N)
{ 
  datlist <- gendata_simu2(seed=i, width=width, height = height,
                           p=p, q=q, d=d, k=k, rank0 = r, bandwidth=1,
                           eta0 = eta, rho=rho, sigma2_eps=sigma2_eps)
  X_count <- datlist$X; H <- datlist$H; Z <- datlist$Z
  F0 <- datlist$F0; B0 <- datlist$B0
  bbeta0 <- datlist$bbeta0; alpha0 <- datlist$alpha0
  Adj_sp <- getneighbor_weightmat(datlist$pos, 1.1,  bandwidth)
  tic <- proc.time()
  reslist <- SpaCOAP(X_count,Adj_sp, H, Z = Z, rank_use = 3, q=q, epsELBO = 1e-20, algo = 1)
  toc <- proc.time()
  time_use <- toc[3] - tic[3]
  metricList$F_tr[i,1] <- trace_statistic_fun(reslist$F, F0)
  metricList$F_CCor[i,1] <- MCCor(reslist$F, F0)
  metricList$B_tr[i,1] <- trace_statistic_fun(reslist$B, B0)
  metricList$B_CCor[i,1] <- MCCor(reslist$B, B0)
  metricList$alpha_norm1[i,1] <- norm1_vec(reslist$alpha- alpha0)/mean(abs(alpha0))  
  metricList$alpha_norm2[i,1] <- norm2_vec(reslist$alpha- alpha0)/sqrt(mean(alpha0^2))
  metricList$beta_norm1[i,1] <- norm1_vec(reslist$bbeta- bbeta0)/mean(abs(bbeta0))
  metricList$beta_norm2[i,1] <- norm2_vec(reslist$bbeta- bbeta0)/sqrt(mean(bbeta0^2))
  metricList$timeMat[i,1] <- time_use
  
  
  library(COAP)
  tic <- proc.time()
  res_coap <- RR_COAP(X_count, Z = cbind(Z, H), rank_use= k+r, q=5, epsELBO = 1e-20)
  toc <- proc.time()
  time_coap <- toc[3] - tic[3]
  metricList$F_tr[i,2] <- trace_statistic_fun(res_coap$H, F0)
  metricList$F_CCor[i,2] <- MCCor(res_coap$H, F0)
  metricList$B_tr[i,2] <- trace_statistic_fun(res_coap$B, B0)
  metricList$B_CCor[i,2] <- MCCor(res_coap$B, B0)
  alpha_coap <- res_coap$bbeta[,1:k]
  beta_coap <- res_coap$bbeta[,(k+1):(k+d)]
  metricList$alpha_norm1[i,2] <- norm1_vec(alpha_coap- alpha0)/mean(abs(alpha0))
  metricList$alpha_norm2[i,2] <- norm2_vec(alpha_coap- alpha0)/sqrt(mean(alpha0^2))
  metricList$beta_norm1[i,2] <- norm1_vec(beta_coap- bbeta0)/mean(abs(bbeta0))
  metricList$beta_norm2[i,2] <- norm2_vec(beta_coap- bbeta0)/sqrt(mean(bbeta0^2))
  metricList$timeMat[i,2] <- time_coap
  
  
  tic <- proc.time()
  res_plnpca <- PLNPCA_run(X_count, cbind(Z[,-1],H), q=q)
  toc <- proc.time()
  time_plnpca <- toc[3] - tic[3]
  
  metricList$F_tr[i,3] <- trace_statistic_fun(res_plnpca$PCs, F0)
  metricList$F_CCor[i,3] <- MCCor(res_plnpca$PCs, F0)
  metricList$B_tr[i,3] <- trace_statistic_fun(res_plnpca$loadings, B0)
  metricList$B_CCor[i,3] <- MCCor(res_plnpca$loadings, B0)
  alpha_plnpca <- t(res_plnpca$bbeta[1:k,])
  beta_plnpca <- t(res_plnpca$bbeta[(k+1):(k+d),])
  metricList$alpha_norm1[i,3] <- norm1_vec(alpha_plnpca- alpha0)/mean(abs(alpha0))
  metricList$alpha_norm2[i,3] <- norm2_vec(alpha_plnpca- alpha0)/sqrt(mean(alpha0^2))
  metricList$beta_norm1[i,3] <- norm1_vec(beta_plnpca- bbeta0)/mean(abs(bbeta0))
  metricList$beta_norm2[i,3] <- norm2_vec(beta_plnpca- bbeta0)/sqrt(mean(bbeta0^2))
  metricList$timeMat[i,3] <- time_plnpca
  
  ## MRRR
  tic <- proc.time()
  res_mrrr <- mrrr_run(X_count, cbind(Z,H), r+ncol(Z), q=q, truncflag=T , trunc = 1000)
  str(res_mrrr)
  toc <- proc.time()
  time_mrrr <- toc[3] - tic[3]
  hbbeta_mrrr <-t(res_mrrr$coef[1:ncol(cbind(Z,H)), ])
  Theta_hb <- (res_mrrr$coef[(ncol(cbind(Z,H))+1): (nrow(cbind(Z,H))+ncol(cbind(Z,H))), ])
  svdTheta <- svd(Theta_hb, nu=q, nv=q)
  metricList$F_tr[i,4] <- trace_statistic_fun(svdTheta$u, F0)
  metricList$F_CCor[i,4] <- MCCor(svdTheta$u, F0)
  metricList$B_tr[i,4] <- trace_statistic_fun(svdTheta$v, B0)
  metricList$B_CCor[i,4] <- MCCor(svdTheta$v, B0)
  alpha_mrrr <- hbbeta_mrrr[,1:k]
  beta_mrrr <- hbbeta_mrrr[,(k+1):(k+d)]
  metricList$alpha_norm1[i,4] <- norm1_vec(alpha_mrrr- alpha0)/mean(abs(alpha0))
  metricList$alpha_norm2[i,4] <- norm2_vec(alpha_mrrr- alpha0)/sqrt(mean(alpha0^2))
  metricList$beta_norm1[i,4] <- norm1_vec(beta_mrrr- bbeta0)/mean(abs(bbeta0))
  metricList$beta_norm2[i,4] <- norm2_vec(beta_mrrr- bbeta0)/sqrt(mean(bbeta0^2))
  metricList$timeMat[i,4] <- time_mrrr
  
  ## FAST
  tic <- proc.time()
  res_fast <- fast_run(X_count, Adj_sp, q=q, verbose=TRUE, epsELBO=1e-8)
  toc <- proc.time()
  time_fast <- toc[3] - tic[3]
  metricList$F_tr[i,5] <- trace_statistic_fun(res_fast$hV, F0)
  metricList$F_CCor[i,5] <- MCCor(res_fast$hV, F0)
  metricList$B_tr[i,5] <- trace_statistic_fun(res_fast$W, B0)
  metricList$B_CCor[i,5] <- MCCor(res_fast$W, B0)
  metricList$timeMat[i,5] <- time_fast
}

sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD)



