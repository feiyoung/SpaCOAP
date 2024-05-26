
#' Generate simulated data
#' @description Generate simulated data from spaital covariate-augmented Poisson factor models
#' @param seed a postive integer, the random seed for reproducibility of data generation process.
#' @param width a postive integer, specify the width of the spatial grid. 
#' @param height a postive integer, specify the height of the spatial grid. 
#' @param p a postive integer, specify the dimension of count variables.
#' @param d a postive integer,  specify the dimension of covariate matrix with low-rank regression coefficient matrix.
#' @param r a postive integer,  specify the dimension of  covariate matrix as control variables.
#' @param q a postive integer,  specify the number of factors.
#' @param rank0 a postive integer, specify the rank of the coefficient matrix.
#' @param eta0 a real between 0 and 1, specify the spatial autocorrelation parameter.
#' @param rho a numeric vector with length 2 and positive elements, specify the signal strength of loading matrix and regression coefficient, respectively. 
#' @param sigma2_eps a positive real, the variance of overdispersion error.
#' @param seed.beta a postive integer, the random seed for reproducibility of data generation process by fixing the regression coefficient matrix beta.
#' @return return a list including the following components:
#' \itemize{
#'   \item \code{X} - the high-dimensional count matrix; 
#'   \item \code{Z} - the low-dimensional covariate matrix with control variables.
#'   \item \code{H} - the high-dimensional covariate matrix;
#'   \item \code{Adj_sp} - the weighted adjacence matrix;
#'   \item \code{alpha0} - the regression coefficient matrix corresponing to Z;
#'   \item \code{bbeta0} - the low-rank large regression coefficient matrix corresponing to H; 
#'   \item \code{B0} - the loading matrix;
#'   \item \code{F0} -  the laten factor matrix;
#'   \item \code{rank0} - the true rank of bbeta0;
#'   \item \code{q} - the true number of factors;
#'   \item \code{eta0} - spatial autocorrelation parameter;
#'   \item \code{pos} - spatial coordinates for each observation.
#' }
#' @details None
#' @seealso \code{\link{SpaCOAP}}
#' @references None
#' @export
#' @useDynLib SpaCOAP, .registration = TRUE
#' @importFrom  MASS mvrnorm
#' @importFrom  stats cov lm residuals rnorm rpois
#' @importFrom  LaplacesDemon rmatrixnorm 
#' @importFrom  Rcpp evalCpp
#' @examples
#' width <- 20; height <- 15; p <- 100
#' d <- 20; k <- 3; q <- 6; r <- 3
#' datlist <- gendata_spacoap(width=width, height=height, p=p, d=20, k=k, q=q, rank0=r)
#' str(datlist)

gendata_spacoap <- function (seed = 1, width=20, height=30, p = 500, d=40, k=3, q = 5,
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


