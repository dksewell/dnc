#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List drawS(const IntegerVector & dims,
                 const NumericVector & Yvec,
                 const NumericVector & OmVec,
                 const NumericVector & XVec,
                 const NumericVector & ALPHA){
  const double PI = 3.14159265358979323846;
int n=dims(0),p=dims(1), TT=dims(2);

arma::cube X(XVec.begin(),p,TT,n);
arma::cube Y(Yvec.begin(),n,n,TT);
arma::cube Om(OmVec.begin(),n,n,TT);
double alpha=Rcpp::as<double>(ALPHA);

arma::colvec bj4 = arma::zeros(n,1);
arma::colvec aj4 = arma::zeros(n,1);

arma::colvec mat1by1 = arma::zeros(1,1);

for(int j=0;j<n;j++){

for(int tt=0;tt<TT;tt++){
for(int i=0;i<n;i++){
if(i!=j){
mat1by1 = trans(X.slice(i).col(tt))*X.slice(j).col(tt);
bj4(j) = bj4(j)+Om(i,j,tt)*mat1by1(0)*mat1by1(0);
aj4(j) = aj4(j)+(Y(i,j,tt)-0.5-Om(i,j,tt)*alpha)*mat1by1(0);
}
}
}
bj4(j)=1/bj4(j);
aj4(j)=bj4(j)*(aj4(j)-1);

}


return Rcpp::List::create(aj4,bj4);
}