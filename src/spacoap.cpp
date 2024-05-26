// This script implement spatial covariate-augumented Poisson factor model.
// Date: 2022-12-27


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//' 
// Calculate the weighted matrix based on distance.
// [[Rcpp::export]]
arma::sp_mat getneighbor_weightmat(const arma::mat x, const double& radius, const double& width)	{
  int N = x.n_rows;
  arma::sp_mat D(N, N);
  double dis;
  uvec idx, idx2;
  for (int j = 0; j < N-1; ++j)
  {    
    idx = find(abs(x(j,0) - x.col(0))<radius); 
    idx2 = find(idx>j);
    int p = idx2.n_elem;
    for (int i = 0; i < p; ++i)
    {
      dis = norm(x.row(idx(idx2(i))) - x.row(j), 2);
      if (dis < radius){
        D(idx(idx2(i)),j) = exp(-dis/width);
        D(j,idx(idx2(i))) = exp(-dis/width);
      }
    }
  }
  return D;
}



// diag(W0^t* Cki * W0)
vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}  



List irlbaCpp(const mat& X, const int& q){
  Rcpp::Environment irlba("package:irlba");
  Rcpp::Function f = irlba["irlba"];
  
  return f(Named("A") = X, Named("nv") = q);
  
}  
// search the best alpha from grid;
double obj_fun_alpha(const mat& M, const mat& Muf, const double& alpha){
  
  mat mat_tmp = (M - alpha* Muf) % (M - alpha* Muf);
  return accu(mat_tmp);
}
double obj_fun_alpha2(const mat& M, const mat& Muf, const double& alpha, const vec& w_plus_vec){
  
  mat mat_tmp = (M - alpha* Muf) % (M - alpha* Muf);
  mat_tmp.each_col() %= w_plus_vec;
  return accu(mat_tmp);
}
double update_eta(const mat& M, const mat& Muf, const vec& alpha_grid, const vec& w_plus_vec){
  
  int ig, ngrid= alpha_grid.n_elem;
  vec objvec(ngrid);
  for(ig=0; ig< ngrid; ++ig){
    objvec(ig) = obj_fun_alpha2(M, Muf, alpha_grid(ig), w_plus_vec);
  }
  
  return alpha_grid(index_min(objvec));
}


double update_eta2(const mat& M, const mat& Muf, const vec& alpha_grid, 
                   const double& S_tr, const vec& w_plus_vec){
  
  int ig, ngrid= alpha_grid.n_elem;
  vec objvec(ngrid);
  for(ig=0; ig< ngrid; ++ig){
    objvec(ig) = obj_fun_alpha2(M, Muf, alpha_grid(ig), w_plus_vec) +alpha_grid(ig)*alpha_grid(ig)*S_tr;
  }
  
  return alpha_grid(index_min(objvec));
}


double calELBO(const mat& X_count, const vec& a, const mat& Z, const mat& H, const mat& Mu_y,
               const mat& S_y, const mat& M, const cube& S, const mat& Muf,
               const vec& invLambda, const mat& B, const mat& bbeta,
               const mat& alpha, const double& eta, const vec& w_plus_vec, const double& wt_sq_sum, const int& algo=1){
  
  int i, n = X_count.n_rows, p = X_count.n_cols, q= B.n_cols;
  double pois_term1, pois_term2, logS, entropy=0.0, ELBO=0.0;
  pois_term1 = accu(X_count % Mu_y); 
  pois_term2 = -accu(exp(Mu_y + S_y/2));
  
  // Rprintf("pois_term2= %4f \n", pois_term2);
  // log P(Y)
  mat dX = (Mu_y -repmat(a, 1, p) - Z*alpha.t() -H*bbeta.t() - M * B.t()) % repmat(sqrt(invLambda.t()), n, 1);
  mat LamB = B % repmat(sqrt(invLambda), 1, q); 
  mat S_sum = sum(S,2);
  double dimreduction_term1 = -0.5*(accu(dX % dX)+ trace(LamB * S_sum * LamB.t())+
                                    accu(S_y % repmat(invLambda.t(), n, 1)) - n* accu(log(invLambda)) );
  
  // Rprintf("dimreduction_term1= %4f \n", dimreduction_term1);
  // log P(F)
  double F_term = 0.0;
  vec Strvec(n, fill::zeros);
  logS = 0;
  for(i=0; i<n; i++){
    Strvec(i) = trace(S.slice(i));
    logS += log(det(S.slice(i)));
  }
  F_term =  accu((M-eta*Muf)%(M-eta*Muf) % repmat(w_plus_vec, 1, q)) + accu(w_plus_vec % Strvec);
  // Rprintf("F_term1= %4f \n", F_term);
  if(algo==1){
    F_term += eta*eta * wt_sq_sum * trace(S_sum);
  }
  // Rprintf("F_term2= %4f \n", F_term);
  F_term = -0.5* F_term;
  
  
  
  entropy = 0.5*accu(log(S_y)) + 0.5 * logS;
  
  // Rprintf("entropy= %4f \n", entropy);
  
  ELBO = pois_term1 + pois_term2 + dimreduction_term1 + F_term + entropy;
  return ELBO;
}

// update bbeta
// separate method
mat update_bbeta_sep(const mat& Z, const mat&Y, const vec& invLambda, const int& rank_use, const bool& fast_svd=true){ // can be speeded up !!!!
  
  // Perform singular value decomposition of X
  int  d = Z.n_cols;
  mat C_ls = pinv(Z.t() * Z) * Z.t()*Y; // d*p
  rowvec sqrt_sigma_inv = sqrt(invLambda.t());
  //Rprintf("good2\n");
  mat ZC = Z * (C_ls % repmat(sqrt_sigma_inv, d, 1)); // n*p
  mat V;
  if(fast_svd){
    Rcpp::List svdX = irlbaCpp(ZC, rank_use);
    mat V1 = svdX["v"];
    V = V1;
    V1.reset();
  }else{
    mat U, V1;
    vec s;
    svd(U, s, V1, ZC);  // how to speed up using approximate SVD
    V = V1;
    U.reset();
    V1.reset();
  }
  ZC.reset();
  mat C_rr = (C_ls % repmat(sqrt_sigma_inv, d, 1))* (V.cols(0, rank_use-1)) * (trans(V.cols(0, rank_use-1)) % repmat(1/sqrt_sigma_inv, rank_use, 1)) ; // d*p
  //Rprintf("good2\n");
  
  return C_rr.t();
  
}

void add_IC_Orth(mat& B){
  // Try the orthogonal matrix method
  int qs1 = B.n_cols;
  mat U1, V1;
  vec s1;
  svd(U1, s1, V1, B);
  vec signU1 = sign(U1.row(0).t());
  B = U1.cols(0,qs1-1) * diagmat(s1 % signU1.subvec(0,qs1-1));
}


arma::mat get_Vmean(const arma::mat& V, const arma::sp_mat& Adj){
  int i, n = V.n_rows, q= V.n_cols;
  vec m(n);
  mat Uv(n, q, fill::zeros);
  for (i = 0; i < n; i++)
  {
    arma::vec col(Adj.col(i)); // the class label of neighbors of i-th sample.
    uvec q1 = find(col > 0);
    // cout<<q1.n_rows<<endl;
    if( q1.n_rows>0){
      Uv.row(i) = mean(V.rows(q1));
    }
    
  }
  return Uv;
}

arma::mat get_weightedmean(const arma::mat& V, const arma::sp_mat& weight_Adj_mat){
  int i, n = V.n_rows, q= V.n_cols;
  vec m(n);
  mat Uv(n, q, fill::zeros);
  for (i = 0; i < n; i++)
  {
    arma::vec col(weight_Adj_mat.col(i)); // the class label of neighbors of i-th sample.
    uvec q1 = find(col > 0);
    // cout<<q1.n_rows<<endl;
    if( q1.n_rows>0){
      Uv.row(i) = sum(V.rows(q1) % repmat(col(q1),1, q)) / accu(col(q1));
    }
    
  }
  return Uv;
}


//  VB_PEstep(X_count, a, Z, H, Mu_y, S_y, invLambda, B, bbeta, alpha, Muf, M, S);
void VB_PEstep(const mat& X_count, const vec& a, const mat& Z, const mat& H,
              const sp_mat& weight_Adj_sp,  mat& Mu_y,  mat& S_y,
              const vec& invLambda, const mat& B, const mat& bbeta, const mat& alpha,
              const vec& alpha_grid,
              const bool& up_eta,const vec& w_plus_vec, const double& wt_sq_sum,  const int& algo, mat& Muf, mat& M,
              cube& S, double& eta){
  
  // algo specify the algorithm used.
  int  i, n = X_count.n_rows, p = X_count.n_cols, q= B.n_cols;
  //  ## VB P-step and E-step
  // update posterior variance of y: S_y
  // update posterior mean of y: Mu_y
  // double elbo2 =  calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
  // Laplace approximation plus Taylor approximation.
  // Matrix operation for this method:
  // Rprintf("Update Mu_y and S_y\n");
  mat tmp_mat1 = repmat(a, 1, p)+ Z* alpha.t() + H * bbeta.t();
  mat tmp_mat = tmp_mat1 + M * B.t();
  Mu_y = (X_count -  exp(Mu_y) % (1- Mu_y) + repmat(invLambda.t(), n,1) % tmp_mat) / 
    (exp(Mu_y) + repmat(invLambda.t(), n,1) );
  S_y = 1.0 / (exp(Mu_y) + repmat(invLambda.t(), n,1) );
  
  // double elbo3 =  calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
  // Rprintf("dMu_S_y= %4f \n", elbo3 - elbo2);
  
  // ## update posterior variance of f_i, m_i and S_i
  // Rprintf("Update M and S\n");
  mat S_tmp = B.t() *  (B % repmat(invLambda, 1, q));
  mat M_tmp = ((Mu_y - tmp_mat1) * (B % repmat(invLambda, 1, q))+ eta* Muf);
  double Str_tmp = 0.0;
  for(i=0; i<n; ++i){
    S.slice(i) = inv(S_tmp + w_plus_vec(i) * eye(q,q));
    M.row(i) = M_tmp.row(i) * S.slice(i);
    Str_tmp += trace(S.slice(i));
  }
  
  // double elbo4 =  calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
  // Rprintf("dM_S= %4f \n", elbo4 - elbo3);
  
  // update eta
  //If using VPEM algorithm:
  // update eta
  if(up_eta && (algo==2) ){
    // Rprintf("update eta using VPEM\n");
    eta = update_eta(M, Muf, alpha_grid, w_plus_vec); // use the previous M.
  }
  
  // P-step
  Muf = get_weightedmean(M, weight_Adj_sp);
  if(algo==1){
    double S_tr = Str_tmp* wt_sq_sum / n; // approximate the term
    // If using VEM algorithm
    if(up_eta){
      // Rprintf("update eta using VEM\n");
      eta = update_eta2(M, Muf, alpha_grid, S_tr, w_plus_vec); // use the previous M.
    }
  }
  // double elbo5 =  calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
  // Rprintf("deta= %4f \n", elbo5 - elbo4);
  
}




// Augumented Poisson factor model
// [[Rcpp::export]]
Rcpp::List spacoap_cpp(const arma::mat& X_count, const arma::vec& a, const arma::mat& Z, 
                       const arma::mat& H, const arma::sp_mat& Adj_sp, const int& rank_use,
                              const arma::mat& Mu_y_int, 
                              const arma::mat& S_y_int,
                              const arma::vec& invLambda_int, const arma::mat& B_int, 
                              const arma::mat& alpha_int, const arma::mat& bbeta_int,
                              const arma::mat& M_int, const arma::mat& S_int,
                              const double& epsELBO, const int& maxIter, const bool& verbose, 
                              const bool& up_eta, const arma::vec& w_plus_vec,const double& wt_sq_sum, 
                              const int& algo, const bool& fast_svd=true,
                              const bool& add_IC_inter=false){
  
  int  n = X_count.n_rows, p = X_count.n_cols, q= B_int.n_cols;
  
  //bool sep_opt_beta = true;
  
  // Initialize
  mat Mu_y(Mu_y_int), S_y(S_y_int), B(B_int), alpha(alpha_int), bbeta(bbeta_int), M(M_int);
  cube S(q,q,n, fill::zeros);
  for(int i=0; i<n; ++i){
    S.slice(i) = S_int;
  }
  vec invLambda(invLambda_int);
  vec alpha_grid = linspace(0.05, 0.95, 30);
  double eta;
  if(up_eta){
    eta = 0.5;
  }else{
    eta = 0.0;
  }
  if(add_IC_inter){
    add_IC_Orth(B);
  }
  
  mat Muf(n,q, fill::zeros);
  vec ELBO_vec(maxIter), Lambda;
  ELBO_vec(0) = INT_MIN;
  mat S_bar, dX, tY;
  int iter;
  
  // Rprintf("Good initialization!\n");
  
  for(iter = 1; iter < maxIter; ++iter){
    
    
    // Rprintf("E step starting!\n");
    // VB P-step and E-step
    VB_PEstep(X_count, a, Z, H, Adj_sp, Mu_y, S_y, invLambda, B, bbeta,
              alpha, alpha_grid, up_eta, w_plus_vec, wt_sq_sum, algo, Muf, M, S, eta);
    // Rprintf("Finish E step!\n");
    //VB M-step
    // double elbo1 = calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
    //update B
    // Rprintf("Update B!\n");
    mat S_sum = sum(S,2);
    B = trans(Mu_y - repmat(a, 1, p)- Z* alpha.t() - H * bbeta.t()) * M * inv(M.t() * M + S_sum);
    // double elbo2 = calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
    // Rprintf("dB= %4f \n", elbo2 - elbo1);
    if(add_IC_inter){
      add_IC_Orth(B);
    }
    
    // update alpha
    // Rprintf("Update alpha!\n");
    alpha = trans(Mu_y - repmat(a, 1, p)- H * bbeta.t()- M*B.t()) * Z * inv_sympd(Z.t()*Z);
    // double elbo3 = calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
    // Rprintf("dalpha= %4f \n", elbo3 - elbo2);
    
    // update Lambda
    // Rprintf("Update Lambda!\n");
    dX = Mu_y -repmat(a, 1, p)- Z*alpha.t()- H * bbeta.t()- M * B.t();
    Lambda =  trans(mean(dX % dX + S_y)) + decomp(mean(S, 2), B);
    Lambda.elem( find(Lambda < 1e-4) ).fill(1e-4); // increase the numerical stability!
    invLambda = 1.0 / Lambda;
    // double elbo4 = calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
    // Rprintf("dLambda= %4f \n", elbo4 - elbo3);
    
    
    // update bbeta
    // Rprintf("Update bbeta!\n");
    tY = Mu_y - repmat(a, 1, p)- Z*alpha.t()- M * B.t();
    bbeta = update_bbeta_sep(H, tY, invLambda, rank_use, fast_svd); // Fixed-rank regression.
    // double elbo5 = calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
    // Rprintf("dbeta= %4f \n", elbo5 - elbo4);
    
    
    
   
    
    ELBO_vec(iter) = calELBO( X_count, a, Z, H, Mu_y, S_y, M, S, Muf, invLambda, B, bbeta, alpha, eta, w_plus_vec, wt_sq_sum,algo);
    
    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n", 
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;
  }
  
  // Add identifiability condition for F and B;
  if(!add_IC_inter){
    mat T = join_rows(Z, H);
    // Add identifiability condition for H and B;
    mat alpha_pri = pinv(T.t() * T) * T.t()* M; // d*q
    mat M_new = M - T * alpha_pri;
    // vec mu = alpha.col(0) + B * alpha_pri.row(0).t();
    Rcpp::List svdX = irlbaCpp(M_new*B.t(), q);
    mat V2 = svdX["v"];
    mat U2 = svdX["u"];
    vec s1 = svdX["d"];
    rowvec signB1 = sign(V2.row(0));
    M = sqrt(n) * U2 * diagmat(signB1.t());
    B = V2 * diagmat(s1.subvec(0, q-1) % signB1.t()) / sqrt(n);
    alpha = trans(Mu_y - repmat(a, 1, p)- H * bbeta.t()- M*B.t()) * Z * inv_sympd(Z.t()*Z);
    tY = Mu_y - repmat(a, 1, p)- Z*alpha.t()- M * B.t();
    bbeta = update_bbeta_sep(H, tY, invLambda, rank_use, fast_svd); // Fixed-rank regression.
  }
  // output return value
  List resList = List::create(
    Rcpp::Named("F") = M,
    Rcpp::Named("B") = B,
    Rcpp::Named("bbeta") = bbeta,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("invLambda") = invLambda,
    Rcpp::Named("eta") = eta,
    //Rcpp::Named("Mu_y") = Mu_y,
    //Rcpp::Named("S_y") = S_y,
    Rcpp::Named("S") = S,
    //Rcpp::Named("Muf") = Muf,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    //Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);
  
}
