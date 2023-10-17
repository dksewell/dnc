#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double cllikProj(const NumericVector & x,
                   const IntegerVector & dims,
                   const NumericVector & ALPHA0,
                   const NumericMatrix & BETA0,
                   const NumericMatrix & UU,
                   const NumericVector & RR,
                   const NumericVector & TAU){
  const double PI = 3.14159265358979323846;

// ###Main objects
//DIMS = c(n,p,K,TT);
int n = dims[0], p=dims[1], K=dims[2], TT=dims[3];
arma::cube X(x.begin(),dims[1],dims[3],dims[0]);

arma::mat uu = Rcpp::as<arma::mat>(UU);
arma::colvec rr = Rcpp::as<arma::colvec>(RR);
arma::colvec tau = Rcpp::as<arma::colvec>(TAU);


arma::mat calpha= Rcpp::as<arma::mat>(BETA0);
arma::rowvec cbeta = Rcpp::as<arma::rowvec>(ALPHA0);


//###Nuisance Objects
arma::mat xiSlice =arma::zeros(p,TT);
arma::mat mvnInsides = arma::zeros(1,1);
arma::colvec Ztm1 = arma::zeros(K,1);
arma::colvec Zt = arma::zeros(K,1);
double ret =0, tempSum =0, tempSum1 =0,JDt=0, JDtm1=0;
//double sumQ =0, normSum =0, numer =0, denom =0;
arma::colvec mui = arma::zeros(p,1);
arma::colvec Xcond = arma::zeros(K,1);

for(int i =0; i<n; i++)
{
xiSlice = X.slice(i);

// Compute unnormalized P(Z_1|X_1);
for(int k=0; k<K;k++)
{
mvnInsides = trans(xiSlice.col(0)- rr(i)*trans(uu.row(k)))*(xiSlice.col(0)- rr(i)*trans(uu.row(k)));
Ztm1(k) = cbeta(k)*pow(tau(i),p*0.5)*pow(2*PI,-0.5*p)*exp(-0.5*mvnInsides(0,0));
}
// Compute marginal P(X_1);
JDtm1 = sum(Ztm1);
// Normalize P(Z_1|X_1);
for(int k=0; k<K;k++)
{
Ztm1(k) = Ztm1(k)/JDtm1;
}

for(int tt=2; tt<TT+1; tt++)
{
// Compute P(X_t|Z_t);
for(int k=0;k<K;k++)
{
mui = rr(i)*trans(uu.row(k));
mvnInsides = tau(i)*trans(xiSlice.col(tt-1)-mui)*(xiSlice.col(tt-1)-mui);
Xcond(k) = pow(tau(i),0.5*p)*pow(2*PI,-0.5*p)*exp(-0.5*mvnInsides(0,0));
}

if(tt<TT)
{
// Compute unnormalized P(Z_t|X_1,...,X_t);
for(int k=0;k<K;k++)
{
Zt(k)=0;
for(int ell=0;ell<K;ell++)
{
Zt(k) = Zt(k) + calpha(ell,k)*Ztm1(ell);
}
Zt(k) = Zt(k)*Xcond(k);
// Zt(k) = (trans(Ztm1)*calpha.col(k))*Xcond(k);
}
// Normalize P(Z_t|X_1,...,X_t);
tempSum = sum(Zt);
for(int k=0;k<K;k++)
{
Zt(k) = Zt(k)/tempSum;
}
}

// Compute P(X_1,...,X_t);
tempSum=0;
for(int elltm1=0;elltm1<K;elltm1++)
{
tempSum1=0;
for(int ellt=0;ellt<K;ellt++)
{
tempSum1 = tempSum1 + calpha(elltm1,ellt)*Xcond(ellt);
}
tempSum = tempSum + Ztm1(elltm1)*tempSum1;
}
JDt = JDtm1*tempSum;

Ztm1 = Zt;
JDtm1 = JDt;

}
ret += log(JDt);


}

return ret;
}


