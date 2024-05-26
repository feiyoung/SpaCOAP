# generate man files
# devtools::document()
# R CMD check --as-cran SpaCOAP_1.2.tar.gz
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_reference()
# pkgdown::build_article("simu")
# pkgdown::build_article("mouseSpleen")

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('html_document'))

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('pdf_document'), clean = F)



#' Fit the SpaCOAP model
#' @description Fit the spatial covariate-augmented overdispersed Poisson factor model
#' @param X_count a count matrix, the observed count matrix with shape n-by-p.
#' @param Adj_sp a sparse matrix, the weighted adjacency matrix;
#' @param H a n-by-d matrix, the covariate matrix with low-rank regression coefficient matrix;
#' @param Z an optional matrix, the fixed-dimensional covariate matrix with control variables; default as a full-one column vector if there is no additional covariates.
#' @param offset an optional vector, the offset for each unit; default as full-zero vector.
#' @param rank_use an optional integer, specify the rank of the regression coefficient matrix; default as 5.
#' @param q an optional string, specify the number of factors; default as 15.
#' @param epsELBO  an optional positive vlaue, tolerance of relative variation rate of the envidence lower bound value, defualt as '1e-8'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param add_IC_inter a logical value, add the identifiability condition in iterative algorithm or add it after algorithm converges; default as FALSE.
#' @param seed an integer, set the random seed in initialization, default as 1;
#' @param algo an optional integer taking value 1 0r 2, select the algorithm used, default as 1, representing variational EM algorithm.
#' @return return a list including the following components:
#' \itemize{
#'   \item \code{F} - the predicted factor matrix; 
#'   \item \code{B} - the estimated loading matrix; 
#'   \item \code{bbeta} - the estimated low-rank large coefficient matrix; 
#'   \item \code{alpha0} - the estimated regression coefficient matrix corresponing to Z;
#'   \item \code{invLambda} - the inverse of the estimated variances of error;
#'   \item \code{eta} - the estimated spatial autocorrelation parameter;
#'   \item \code{S} - the approximated posterior covariance for each row of F;
#'   \item \code{ELBO} -  the ELBO value when algorithm stops;
#'   \item \code{ELBO_seq} - the sequence of ELBO values.
#'   \item \code{time_use} - the running time in model fitting of SpaCOAP;
#' }
#' @details None
#' @seealso None
#' @references Liu W, Zhong Q. High-dimensional covariate-augmented overdispersed poisson factor model. Biometrics. 2024 Jun;80(2):ujae031.
#' @export
#' @useDynLib SpaCOAP, .registration = TRUE
#' @importFrom  irlba irlba
#' @importFrom  Rcpp evalCpp
#'
#'
#' @examples
#' width <- 20; height <- 15; p <- 100
#' d <- 20; k <- 3; q <- 6; r <- 3
#' datlist <- gendata_spacoap(width=width, height=height, p=p, d=20, k=k, q=q, rank0=r)
#' fitlist <- SpaCOAP(X_count=datlist$X, Adj_sp = datlist$Adj_sp, H= datlist$H, Z = datlist$Z, q=6, rank_use=3)
#' str(fitlist)

SpaCOAP <- function(X_count, Adj_sp, H, Z=matrix(1, nrow(X_count),1),offset=rep(0, nrow(X_count)), rank_use=5,
                    q=15, epsELBO=1e-8, 
                    maxIter=30,  verbose=TRUE,  add_IC_inter=FALSE,
                    seed=1, algo=1){
  # offset=rep(0, nrow(X_count));rank_use=5;q=15; algo=1; seed=1;add_IC_inter=FALSE
  # q=5; epsELBO=1e-6;  maxIter=10; verbose=TRUE
  
  if(ncol(H) < 2) stop("SpaCOAP: ncol(H) must be greater than 1!")
  Diag <- function(vec){
    q <- length(vec)
    if(q > 1){
      y <- diag(vec)
    }else{
      y <- matrix(vec, 1,1)
    }
    return(y)
  }
  scaleMatrix_byrow <- function(adj){
    require(Matrix)
    s <- Matrix::rowSums(adj)+1e-8
    n <- length(s)
    adj_norm <- adj / matrix(s, nrow=n, ncol=n, byrow = FALSE)
    #adj_norm <- adj
    # for(i in 1:n){
    #   adj_norm[i, ] <- adj[i,]/ s[i]
    # }
    
    return(adj_norm)
  }
  get_initials <- function(X, q){
    require(irlba)
    n <- nrow(X); p <- ncol(X)
    mu <- colMeans(X)
    X <- X - matrix(mu, nrow=n, ncol=p, byrow=TRUE)
    svdX  <- irlba(A =X, nv = q)
    PCs <- sqrt(n) * svdX$u
    # loadings <- svdX$v %*% diag(svdX$d[1:q]) / sqrt(n)
    loadings <- svdX$v %*% Diag(svdX$d[1:q]) / sqrt(n)
    dX <- PCs %*% t(loadings) - X
    Lam_vec <- colSums(dX^2)/n
    return(list(hH = PCs, hB = loadings, hmu=mu,sigma2vec = Lam_vec))
    
  }
  message("Calculate initial values...")
  n <- nrow(X_count); p <- ncol(X_count); 
  if(any(Z[,1]!=1)) warning("The first column of covariates Z is not a full-one column vector, so it will fit a model without intercept!")
  Mu_y_int = log(1+ X_count)
  S_y_int = matrix(1, n, p);
  a <- offset
  lm1 <- lm(Mu_y_int~0+cbind(Z, H))
  coefM <- t(coef(lm1))
  resi <- resid(lm1)
  set.seed(seed)
  fit_approxPCA <- get_initials(resi, q=q)
  B_int <- fit_approxPCA$hB
  M_int <-  fit_approxPCA$hH
  w_plus_vec = Matrix::rowSums(Adj_sp)
  w_plus_vec[w_plus_vec<1e-20] <- 1e-5
  wt_sq_sum <- sum(Matrix::rowSums(scaleMatrix_byrow(Adj_sp)^2) * w_plus_vec)
  S_int <- diag(rep(1, q))
  invLambda_int = rep(1, p); ## the default p is 100, the default q is 15
  k <- ncol(Z)
  d <- ncol(H)
  rank_use <- min(rank_use, d)
  alpha_int <- coefM[,1:k]
  if(k==1){
    alpha_int <- matrix(alpha_int, ncol=1)
  }
  bbeta_int <- coefM[,(k+1):(k+d)]
  message("Model fitting...")
  tic <- proc.time()
  reslist <- spacoap_cpp(X_count, a, Z, H, Adj_sp, rank_use, Mu_y_int, S_y_int, 
                         invLambda_int, B_int, alpha_int, bbeta_int, M_int, S_int, 
                         epsELBO, maxIter, verbose, TRUE, w_plus_vec, wt_sq_sum, algo=1L, 
                         fast_svd = TRUE, add_IC_inter) 
  toc <- proc.time()
  reslist$time_use <- toc[3] - tic[3]
  
  return(reslist)
  
}

#' Select the parameters in COAP models
#' @description Select the number of factors and the rank of coefficient matrix in the covariate-augmented overdispersed Poisson factor model
#' @param X_count a count matrix, the observed count matrix with shape n-by-p.
#' @param Adj_sp a sparse matrix, the weighted adjacency matrix;
#' @param H a n-by-d matrix, the covariate matrix with low-rank regression coefficient matrix;
#' @param Z an optional matrix, the fixed-dimensional covariate matrix with control variables; default as a full-one column vector if there is no additional covariates.
#' @param offset an optional vector, the offset for each unit; default as full-zero vector.
#' @param q_max an optional string, specify the upper bound for the number of factors; default as 15.
#' @param r_max an optional integer, specify the upper bound for the rank of the regression coefficient matrix; default as 24.
#' @param threshold  an optional 2-dimensional positive vector, specify the the thresholds that filters the singular values of beta and B, respectively.
#' @param verbose a logical value, whether output the information in iteration.
#' @param ..., other arguments passed to the function \code{\link{SpaCOAP}}.
#' @return return a named vector with names `hr` and `hq`, the estimated rank and number of factors.
#' @details The threshold is to filter the singular values with  low signal, to assist the identification of underlying model structure.
#' @seealso \code{\link{SpaCOAP}}
#' @references None
#' @export
#'
#'
#'
#' @examples
#' width <- 20; height <- 15; p <- 300
#' d <- 20; k <- 3; q <- 6; r <- 3
#' datlist <- gendata_spacoap(width=width, height=height, p=p, d=d, k=k, q=q, rank0=r)
#' set.seed(1)
#' para_vec <- chooseParams(X_count=datlist$X, Adj_sp=datlist$Adj_sp, H= datlist$H, Z = datlist$Z, r_max=6)
#' print(para_vec)


chooseParams <- function(X_count, Adj_sp, H, Z=matrix(1, nrow(X_count),1), offset=rep(0, nrow(X_count)), 
         q_max=15, r_max=24,
         threshold=c(1e-1, 1e-2),verbose=TRUE, ...){
  
  reslist <- SpaCOAP(X_count, Adj_sp, H, Z = Z, rank_use = r_max, q= q_max,
                     verbose=verbose,...)
 
  thre1 <- threshold[1]
  beta_svalues <- svd(reslist$bbeta)$d
  beta_svalues <- beta_svalues[beta_svalues>thre1]
  ratio1 <- beta_svalues[-length(beta_svalues)] / beta_svalues[-1]
  hr <- which.max(ratio1[-length(ratio1)])
  
  
  thre2 <- threshold[2]
  B_svalues <- svd(reslist$B)$d
  B_svalues <- B_svalues[B_svalues>thre2]
  ratio_fac <- B_svalues[-length(B_svalues)] / B_svalues[-1]
  hq <- which.max(ratio_fac)
  
  return(c(hr=hr, hq=hq))
}



