#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec drawTau(const IntegerVector & dims,
                     const NumericVector & Xvec,
                     const NumericVector & RR,
                     const NumericMatrix & UU,
                     const IntegerMatrix & Z){
  const double PI = 3.14159265358979323846;
                       
int n=dims(0),p=dims(1), TT=dims(2);
arma::cube X(Xvec.begin(),p,TT,n);
arma::colvec rr=Rcpp::as<arma::colvec>(RR);
arma::mat uu = Rcpp::as<arma::mat>(UU);

arma::colvec bi2Inv = arma::zeros(n,1);

for(int i=0;i<n;i++){

for(int tt=0;tt<TT;tt++){
bi2Inv(i) = bi2Inv(i) + 
pow(arma::norm(X.slice(i).col(tt)-rr(i)*uu.col(Z(i,tt)),2),2);
}
bi2Inv(i) = 0.5*bi2Inv(i);

}

return bi2Inv;
}
