#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube csimulateData(const NumericVector & Xvec,
                         const IntegerVector & dims,
                         const NumericVector & ALPHA,
                         const NumericVector & SS1){
  const double PI = 3.14159265358979323846;
    
  /*Dims is c(n,p,TT)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  double alpha = Rcpp::as<double>(ALPHA);
  arma::colvec ss1 = Rcpp::as<arma::colvec>(SS1);
  arma::cube Y= arma::zeros(dims(0),dims(0),dims(2));  
  
  
  double uu=0, eta=0, ProbY=0;
  arma::mat insides = arma::zeros(1,1);
  
  
  for(int tt=0;tt<dims(2);tt++)
{
  for(int i=0;i<dims(0);i++)
{
  for(int j=0;j<dims(0);j++)
{
  if(j!=i)
{
  eta = alpha + ss1(j)*arma::dot(X.slice(i).col(tt),X.slice(j).col(tt));
  ProbY = 1/(1+exp(-eta));
  uu = arma::randu();
  if(uu<ProbY)
{
  Y(i,j,tt)=1;
}
  
  
}
}
}
}
  
  return Y;
}

