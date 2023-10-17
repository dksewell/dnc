#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat drawU(const IntegerVector & dims,
                 const IntegerVector & TAU,
                 const IntegerVector & RR,
                 const NumericVector & Xvec,
                 const IntegerMatrix & Z){
  const double PI = 3.14159265358979323846;
int n=dims(0),p=dims(1), TT=dims(2), M=dims(3);
arma::cube X(Xvec.begin(),p,TT,n);
arma::colvec tau=Rcpp::as<arma::colvec>(TAU);
arma::colvec rr=Rcpp::as<arma::colvec>(RR);

arma::mat nu = arma::zeros(p,M);

for(int tt=0;tt<TT;tt++){
for(int i=0;i<n;i++){
nu.col(Z(i,tt)) = nu.col(Z(i,tt)) + tau(i)*rr(i)*X.slice(i).col(tt);
}
}


return nu;
}
