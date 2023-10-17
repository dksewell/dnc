#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double cposteriorNoOmega(const IntegerVector & dims,
                         const NumericVector & Yvec,
                         const NumericVector & XVec,
                         const NumericVector & RR,
                         const NumericVector & TAU,
                         const NumericVector & SS1,
                         const NumericVector & ALPHA,
                         const NumericMatrix & BETA,
                         const IntegerMatrix & Z,
                         const NumericMatrix & UU,
                         const NumericVector & CC,
                         const NumericMatrix & GAMSTAR,
                         const NumericVector & A2STAR,
                         const NumericVector & B2STAR,
                         const NumericVector & B3STAR){
  const double PI = 3.14159265358979323846;
  //dims = c(n,p,TT,MM)  
  //Variables to be read in:
  arma::cube Y(Yvec.begin(),dims[0],dims[0],dims[2]);
  double cc = Rcpp::as<double>(CC);
  double a2Star = Rcpp::as<double>(A2STAR);
  double b2Star = Rcpp::as<double>(B2STAR);
  double b3Star = Rcpp::as<double>(B3STAR);
  
  int n=dims(0),p=dims(1), TT=dims(2), M=dims(3);
  arma::cube X(XVec.begin(),p,TT,n);
  arma::colvec rr= Rcpp::as<arma::colvec>(RR);
  arma::colvec tau= Rcpp::as<arma::colvec>(TAU);
  arma::mat uu = Rcpp::as<arma::mat>(UU);
  arma::colvec ss1=Rcpp::as<arma::colvec>(SS1);
  double alpha = Rcpp::as<double>(ALPHA);
  arma::mat Beta = Rcpp::as<arma::mat>(BETA);
  arma::mat GamStar = Rcpp::as<arma::mat>(GAMSTAR);
  
  
  //Nuisance Variables
  double ret =0, eta=0;
  arma::colvec mat1by1 =arma::zeros(1,1);
  
  // Likelihood
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
  
  //Priors sans Omega
  ret -= 0.5*alpha*alpha/b3Star;
  
  for(int i=0;i<n;i++){
  ret -= ss1(i) - log(tau(i)) +rr(i)*tau(i)/cc - (a2Star-1)*log(tau(i)) +tau(i)/b2Star;
  }
  
  for(int i=0;i<n;i++){
  mat1by1 = (trans(X.slice(i).col(0))-rr(i)*uu.row(Z(i,0)))*(X.slice(i).col(0)-rr(i)*trans(uu.row(Z(i,0))));
  ret += log(Beta(0,Z(i,0)))+0.5*TT*p*log(tau(i)) - 0.5*tau(i)*mat1by1(0);
  for(int tt=1;tt<TT;tt++){
  mat1by1 = (trans(X.slice(i).col(tt))-rr(i)*uu.row(Z(i,tt)))*(X.slice(i).col(tt)-rr(i)*trans(uu.row(Z(i,tt))));
  ret += log(Beta(Z(i,tt-1)+1,Z(i,tt))) - 0.5*tau(i)*mat1by1(0);
  }
  }
  
  ret -= M*log(2*pow(PI,p/2.0)/tgamma(p/2.0));
  
  for(int hh=0; hh<M+1;hh++){
  ret += lgamma(sum(GamStar.row(hh)));
  for(int gg=0; gg<M;gg++){
  ret += (GamStar(hh,gg)-1.0)*log(Beta(hh,gg)) - lgamma(GamStar(hh,gg));
  }
  }
  
  
  return ret;
}