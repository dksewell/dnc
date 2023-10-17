#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double cLogLikelihood(const NumericVector & Yvec, 
                      const NumericVector & Xvec, 
                      const IntegerVector & DIMS, 
                      const NumericVector & ALPHA,
                      const NumericVector & SS1){
  const double PI = 3.14159265358979323846;


/*Dims is c(n,p,TT)
Z nxT; X pxTxn; Y pxpxT; 
*/
Rcpp::IntegerVector dims(DIMS);
arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
arma::cube Y(Yvec.begin(),dims[0],dims[0],dims[2]);
double alpha = Rcpp::as<double>(ALPHA);
arma::colvec ss1 = Rcpp::as<arma::colvec>(SS1);

double ret =0, eta=0;

for(int tt=0;tt<dims(2);tt++)
{
for(int i=0;i<dims(0);i++)
{
for(int j=0;j<dims(0);j++)
{
if(j!=i)
{
eta = alpha + ss1(j)*arma::dot(X.slice(i).col(tt),X.slice(j).col(tt));
ret += Y(i,j,tt)*eta-log(1+exp(eta));
}
}
}
}

return ret;

}