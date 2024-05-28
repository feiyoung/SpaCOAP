

# Generate data -----------------------------------------------------------------

gendata_simu2 <- function (seed = 1, width=20, height=30, p = 500, d=40, k=3, q = 5,
          rank0=3, eta0 = 0.5, bandwidth = 1,
          rho = c(10, 1), sigma2_eps=1, seed.beta=1){
  require(MASS)
  library(LaplacesDemon)
  if(rank0<1) stop("gendata_simu1: rank0 must be greater than 0!")
  if(sigma2_eps<=0) stop("gendata_simu1: sigma2_eps must be a postive real!")
  cor.mat<-function (p, rho, type = "toeplitz") {
    if (p == 1) 
      return(matrix(1, 1, 1))
    mat <- diag(p)
    if (type == "toeplitz") {
      for (i in 2:p) {
        for (j in 1:i) {
          mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
        }
      }
    }
    if (type == "identity") {
      mat[mat == 0] <- rho
    }
    return(mat)
  }
  Diag<-function (vec){
    q <- length(vec)
    if (q > 1) {
      y <- diag(vec)
    }
    else {
      y <- matrix(vec, 1, 1)
    }
    return(y)
  }
  normalizeMatrix <- function(adj, alpha0){
    require(Matrix)
    s <- Matrix::colSums(adj)
    n <- nrow(adj)
    adj_norm <- adj
    for(i in 1:n){
      adj_norm[i, i:n] <- alpha0*adj[i, i:n]/ s[i]
      adj_norm[i:n, i] <-  adj_norm[i, i:n]
    }
    return(adj_norm)
  }
  
  if(length(rho)<2) stop("rho must be a numeric vector of length 2!")
  
  
  fac_B <- rho[1]
  fac_z <- rho[2]
  set.seed(seed.beta) # Fixed bbeta0
  alpha0 <- matrix(rnorm(p* k), p, k) * fac_z
  rank_true <- rank0
  bbeta0 <- t(matrix(rnorm(d*rank_true), d, rank_true) %*% matrix(rnorm(rank_true* p), rank_true, p))*fac_z
  
  B = matrix(rnorm(p*q), p, q)
  qrB <- qr(B)
  B0 <- qr.Q(qrB) %*% diag(sqrt(seq(q, 1, by=-1))) *fac_B
  B0 <- B0 %*% GFM:::Diag(sign(B0[1,]))
  
  set.seed(seed)
  if(k<2) stop("k must be greater than 1!")
  n <- width*height
  
  set.seed(seed)
  ## generate spatially dependent factors: F0
  pos <- cbind(rep(seq_len(width), each=height), rep(seq_len(height), width))
  # Adj_sp <- PRECAST::getAdj_reg(pos, platform='ST')
  Adj_sp <- getneighbor_weightmat(pos, 1.1,  bandwidth)
  M <- matrix(0, nrow=n, ncol=q)
  W <- normalizeMatrix(Adj_sp, alpha0=eta0)
  U <- as(diag(rep(1, n)), "sparseMatrix") -  W
  U2 <- qr.solve(U)
  V <-  diag(q)
  U1 <- round(as.matrix(U2),4)
  F0 <- rmatrixnorm(M, U1, V)
  
  Z <- cbind(1, matrix(rnorm(n*(k-1)), n, k-1))
  H <- MASS::mvrnorm(n, mu=rep(2, d), Sigma=cor.mat(p=d, rho=0.5)) *0.1
  epsilon = matrix(rnorm(n*p, sd= sqrt(sigma2_eps)), n, p)
  Xi = Z %*% t(alpha0) + H %*% t(bbeta0) + F0 %*%t(B0) + epsilon
  
  X <- matrix(rpois(n = n*p, lambda =exp(Xi)) ,n, p)
  
  return(list(X = X, Z=Z, H= H, Adj_sp = Adj_sp, alpha0=alpha0, bbeta0=bbeta0, B0 = B0, F0 = F0, 
              rank=rank0, q=q, eta0=eta0, pos=pos))
}



# Metrics -----------------------------------------------------------------


moranI <- function(x, w, sample_rate=1, seed=1){
  require(ape)
  set.seed(seed)
  idx <- sample(length(x), sample_rate*length(x))
  Moran.I(x[idx], as.matrix(w[idx,idx]))$observed
}

colSD <- function(mat, na.rm=TRUE){
  apply(mat, 2, sd, na.rm=na.rm)
}

norm2_vec <- function(x){
  if(is.matrix(x)) x <- as.vector(x)
  sqrt(sum(x^2/ length(x)))
} 

norm1_vec <- function(x) mean(abs(x))
MCCor <- function(H, H0){
  mean(cancor(H, H0)$cor)
}

trace_statistic_fun <- function(H, H0){
  
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% qr.solve(t(H) %*% H) %*% t(H) %*% H0
  
  tr_fun(mat1) / tr_fun(t(H0) %*% H0)
  
}



# Other mehthods ----------------------------------------------------------


factorm <- function(X, q=NULL){
  
  signrevise <- GFM:::signrevise
  if ((!is.null(q)) && (q < 1)) 
    stop("q must be NULL or other positive integer!")
  if (!is.matrix(X)) 
    stop("X must be a matrix.")
  mu <- colMeans(X)
  X <- scale(X, scale = FALSE)
  
  n <- nrow(X)
  p <- ncol(X)
  if (p > n) {
    svdX <- eigen(X %*% t(X))
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hatF <- as.matrix(svdX$vector[, 1:q] * sqrt(n))
    B2 <- n^(-1) * t(X) %*% hatF
    sB <- sign(B2[1, ])
    hB <- B2 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
    hH <- sapply(1:q, function(k) hatF[, k] * sign(B2[1, 
                                                      ])[k])
  }
  else {
    svdX <- eigen(t(X) %*% X)
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hB1 <- as.matrix(svdX$vector[, 1:q])
    hH1 <- n^(-1) * X %*% hB1
    svdH <- svd(hH1)
    hH2 <- signrevise(svdH$u * sqrt(n), hH1)
    if (q == 1) {
      hB1 <- hB1 %*% svdH$d[1:q] * sqrt(n)
    }
    else {
      hB1 <- hB1 %*% diag(svdH$d[1:q]) * sqrt(n)
    }
    sB <- sign(hB1[1, ])
    hB <- hB1 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
    hH <- sapply(1:q, function(j) hH2[, j] * sB[j])
  }
  sigma2vec <- colMeans((X - hH %*% t(hB))^2)
  res <- list()
  res$hH <- hH
  res$hB <- hB
  res$mu <- mu
  res$q <- q
  res$sigma2vec <- sigma2vec
  res$propvar <- sum(evalues[1:q])/sum(evalues)
  res$egvalues <- evalues
  attr(res, "class") <- "fac"
  return(res)
}


PLNPCA_run <- function(X_count, covariates, q,  Offset=rep(1, nrow(X_count)), workers=NULL,
                       maxIter=10000,ftol_rel=1e-8, xtol_rel= 1e-4){
  require(PLNmodels)
  if(!is.null(workers)){
    future::plan("multisession", workers = workers)
  }
  if(!is.character(Offset)){
    dat_plnpca <- prepare_data(X_count, covariates)
    dat_plnpca$Offset <- Offset
  }else{
    dat_plnpca <- prepare_data(X_count, covariates, offset = Offset)
  }
  
  d <- ncol(covariates)
  #  offset(log(Offset))+
  formu <- paste0("Abundance ~ 1 + offset(log(Offset))+",paste(paste0("V",1:d), collapse = '+'))
  control_use  <- list(maxeval=maxIter, ftol_rel=ftol_rel, xtol_rel= ftol_rel)
  control_main <- PLNPCA_param(
    backend = "nlopt",
    trace = 1,
    config_optim = control_use,
    inception = NULL
  )
  
  myPCA <- PLNPCA(as.formula(formu), data = dat_plnpca, ranks = q,  control = control_main)
  
  myPCA1 <- getBestModel(myPCA)
  myPCA1$scores
  
  res_plnpca <- list(PCs= myPCA1$scores, bbeta= myPCA1$model_par$B, 
                     loadings=myPCA1$model_par$C)
  
  return(res_plnpca)
}

## Compare with MRRR
mrrr_run <- function(Y, X, rank0, q=NULL, family=list(poisson()),
                     familygroup=rep(1,ncol(Y)), epsilon = 1e-4, sv.tol = 1e-2,
                     maxIter = 2000, trace=TRUE, truncflag=FALSE, trunc=500){
  # epsilon = 1e-4; sv.tol = 1e-2; maxIter = 30; trace=TRUE
  # Y <- X_count; X <- cbind(Z, H); rank0 = r + ncol(Z)
  
  require(rrpack)
  
  n <- nrow(Y); p <- ncol(Y)
  
  if(!is.null(q)){
    rank0 <- rank0+q
    X <- cbind(X, diag(n))
  }
  if(truncflag){
    ## Trunction
    Y[Y>trunc] <- trunc
    
  }
  
  svdX0d1 <- svd(X)$d[1]
  init1 = list(kappaC0 = svdX0d1 * 5)
  offset = NULL
  control = list(epsilon = epsilon, sv.tol = sv.tol, maxit = maxIter,
                 trace = trace, gammaC0 = 1.1, plot.cv = TRUE,
                 conv.obj = TRUE)
  fit.mrrr <- mrrr(Y=Y, X=X[,-1], family = family, familygroup = familygroup,
                   penstr = list(penaltySVD = "rankCon", lambdaSVD = 1),
                   control = control, init = init1, maxrank = rank0)
  
  return(fit.mrrr)
}

gllvm_run <- function(Y, X, q){
  library(gllvm)
  nc <- ncol(X) -1 
  colnames(X) <- c("Intercept", paste0("V",1:nc))
  tic <- proc.time()
  res_gllvm <- gllvm(y=Y, X=X, family=poisson(), num.lv= q, control = list(trace=T))
  # hH_gllvm <- predictLVs.gllvm(res_gllvm)
  toc <- proc.time()
  time_gllvm <- toc[3] - tic[3]
  
}


fast_run <- function(X_count, Adj_sp, q, verbose=TRUE, epsELBO=1e-8){
  require(ProFAST)

  reslist <- FAST_run(XList = list(X_count), 
                      AdjList = list(Adj_sp), q = q, fit.model = 'poisson', 
                      verbose=verbose, epsLogLik=epsELBO)
  reslist$hV <- reslist$hV[[1]]
  return(reslist)
}

## SpatialPCA

SpatialPCA_run <- function(X_count, covariates,pos,  SpatialPCnum=15, fast=TRUE, bandwidthtype = c("Silverman","SJ")){
  
  bandwidthtype <- match.arg(bandwidthtype)
 
  Xmat <- t(X_count)
  colnames(Xmat) <- paste0("spot",1:ncol(Xmat))
  row.names(Xmat) <- paste0("gene",1:nrow(Xmat))
  require(SpatialPCA)
  location <- pos ## location should plus a big number for each sample.
  row.names(location) <- colnames(Xmat)
 
  stereo_seq = CreateSpatialPCAObject(counts=Xmat, 
                                        location=location, project = "SpatialPCA",
                                        customGenelist=row.names(Xmat),
                                        min.loctions = 0, min.features=0)

  
  
  stereo_seq = SpatialPCA_buildKernel(stereo_seq, kerneltype="gaussian", bandwidthtype=bandwidthtype,
                                      bandwidth.set.by.user=NULL,sparseKernel=TRUE,
                                      sparseKernel_tol=1e-20,sparseKernel_ncore=1)
  stereo_seq = SpatialPCA_EstimateLoading(stereo_seq,fast=fast,SpatialPCnum=SpatialPCnum)
  stereo_seq = SpatialPCA_SpatialPCs(stereo_seq, fast=fast)
bbeta<-solve(t(covariates)%*%covariates)%*%t(covariates)%*%(X_count-t(as.matrix(stereo_seq@SpatialPCs))%*%t(stereo_seq@W))
  res_spaPCA <- list(loading=stereo_seq@W, PCs=t(as.matrix(stereo_seq@SpatialPCs)),bbeta=t(bbeta))

  return(res_spaPCA)
}


