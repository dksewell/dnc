#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List drawMiss(const IntegerVector & dims,
                    const NumericVector & Yvec,
                    const NumericVector & SS1,
                    const IntegerMatrix & MISS,
                    const NumericVector & XVec,
                    const NumericVector & ALPHA){
  const double PI = 3.14159265358979323846;
    
int n=dims(0),p=dims(1), TT=dims(2);
arma::cube X(XVec.begin(),p,TT,n);
arma::cube Y(Yvec.begin(),n,n,TT);
double alpha=Rcpp::as<double>(ALPHA);
arma::colvec ss1 = Rcpp::as<arma::colvec>(SS1);
Rcpp::IntegerMatrix Miss(MISS);//Miss is a 4xM matrix, colnames=row,col,time,count

double YijProb=0;
arma::colvec uu = arma::zeros(1,1);
int MM = Miss.nrow();

for(int mm=0;mm<MM;mm++){
YijProb = alpha + ss1(Miss(mm,1)-1)*arma::dot(X.slice(Miss(mm,0)-1).col(Miss(mm,2)-1),
X.slice(Miss(mm,1)-1).col(Miss(mm,2)-1));
YijProb = 1/(1+exp(-YijProb));
uu = arma::randu(1);
if(uu(0)<YijProb){
Miss(mm,3) += 1;
Y(Miss(mm,0)-1,Miss(mm,1)-1,Miss(mm,2)-1)=1;
}
}

return Rcpp::List::create(Y,Miss);
}                        
