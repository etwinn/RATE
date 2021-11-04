// Emily Winn
// We are running a sampler for multivariate normal using code written by Ahmadou Dicko (https://gallery.rcpp.org/articles/simulate-multivariate-normal/) make sure
// to cite in the paper.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat GaussKernel(arma::mat X, double h = 1){
    int i,j;
    double n = X.n_cols;
    double p = X.n_rows;
    mat K = zeros<mat>(n,n);
    for (i = 0; i<n; i++){
        for(j = 0; j<n; j++){
            if(i==j){
                break;
            }
            else{
                K(i,j) = exp(-h/(p)*sum(pow(X.col(i)-X.col(j),2)));
            }
        }
    }
    return K + trans(K);
}

// [[Rcpp::export]]
arma::mat GaussCoKernel(arma::mat X, arma::mat Y, double h = 1){
    int i,j;
    double n = X.n_cols;
    double p = X.n_rows;
	if (Y.n_cols != n || Y.n_rows != p) {
		throw std::invalid_argument( "Dimensions of matrices don't match.");
	}
    mat K = zeros<mat>(n,n);
    for (i = 0; i<n; i++){
        for(j = 0; j<n; j++){
            K(i,j) = exp(-h/(p)*sum(pow(X.col(i)-Y.col(j),2)));
        }
    }
    return K;
}