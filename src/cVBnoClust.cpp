#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cVBnoClust(const IntegerVector & dims,
                      const NumericVector & Yvec,
                      const NumericVector & muVec,
                      const NumericVector & SigVec,
                      const NumericVector & A2STAR,
                      const NumericVector & B2STAR,
                      const NumericVector & B3STAR,
                      const NumericVector & EOmVec,
                      const NumericVector & A0,
                      const NumericVector & B0,
                      const NumericVector & A0STAR,
                      const NumericVector & B0STAR, 
                      const NumericVector & AI2,
                      const NumericVector & BI2,
                      const NumericVector & A3,
                      const NumericVector & B3,
                      const NumericVector & ES1,
                      const NumericVector & ES2){
  const double PI = 3.14159265358979323846;

RNGScope scope;
Environment base("package:base");
Function digamma = base["digamma"];


//dims = c(n,p,TT,MM)

//Variables to be read in:
int n=dims(0),p=dims(1), TT=dims(2);
arma::cube Y(Yvec.begin(),dims[0],dims[0],dims[2]);
arma::cube mu(muVec.begin(),p,TT,n);
arma::cube Sig(SigVec.begin(),TT*p,p,n);
double a0Star = Rcpp::as<double>(A0STAR);
double b0Star = Rcpp::as<double>(B0STAR);
double a2Star = Rcpp::as<double>(A2STAR);
double b2Star = Rcpp::as<double>(B2STAR);
double b3Star = Rcpp::as<double>(B3STAR);

double a0=Rcpp::as<double>(A0);
double b0=Rcpp::as<double>(B0);
arma::cube EOm(EOmVec.begin(),n,n,TT);  
arma::colvec ai2= Rcpp::as<arma::colvec>(AI2);
arma::colvec bi2= Rcpp::as<arma::colvec>(BI2);
double a3=Rcpp::as<double>(A3);
double b3=Rcpp::as<double>(B3);
arma::colvec Es1=Rcpp::as<arma::colvec>(ES1);
arma::colvec Es2=Rcpp::as<arma::colvec>(ES2);
arma::colvec ai4= arma::zeros(n,1);
arma::colvec bi4= arma::zeros(n,1);

//Nuisance Variables
double const1=0, const2=0, const3=0;
arma::colvec mat1by1 =arma::zeros(1,1);
arma::colvec cvecp1 = arma::zeros(p,1);
arma::colvec cvecp2 = arma::zeros(p,1);
arma::mat matpbyp1 = arma::zeros(p,p);

// omega_{ijt} ----------------------------------------

for(int tt=0;tt<TT;tt++){
for(int i=0;i<(n-1);i++){
for(int j=i+1;j<n;j++){

mat1by1 = trans(mu.slice(i).col(tt))*mu.slice(j).col(tt);
const1 = b3+a3*a3+2*a3*Es1(j)*mat1by1(0);
const2 = b3+a3*a3+2*a3*Es1(i)*mat1by1(0);
mat1by1 = mat1by1*mat1by1 + 
trans(mu.slice(j).col(tt))*Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1)*mu.slice(j).col(tt) +
trans(mu.slice(i).col(tt))*Sig.slice(j).submat(tt*p,0,(tt+1)*p-1,p-1)*mu.slice(i).col(tt);
const3 = arma::trace(Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1)*Sig.slice(j).submat(tt*p,0,(tt+1)*p-1,p-1)) +
mat1by1(0);
const1 += Es2(j)*const3;
const2 += Es2(i)*const3;
const1 = sqrt(const1);
const2 = sqrt(const2);
EOm(i,j,tt) = 0.5/const1*(1-exp(-const1))/(1+exp(-const1));
EOm(j,i,tt) = 0.5/const2*(1-exp(-const2))/(1+exp(-const2));

}
}
}



// X_{it} ----------------------------------------


for(int i=0; i<n;i++){
//t=1
cvecp1 = arma::zeros(p,1);
for(int j=0; j<n;j++){
if(j != i){
cvecp1 = cvecp1 + ( (Y(i,j,0)-0.5)*Es1(j)+(Y(j,i,0)-0.5)*Es1(i)- 
a3*(EOm(i,j,0)*Es1(j)+EOm(j,i,0)*Es1(i)))*mu.slice(j).col(0);
}
}
cvecp1 = cvecp1 + ai2(i)*bi2(i)*mu.slice(i).col(1);

Sig.slice(i).submat(0,0,p-1,p-1) = arma::zeros(p,p);
for(int j=0;j<n;j++){
if(j != i){
Sig.slice(i).submat(0,0,p-1,p-1) = Sig.slice(i).submat(0,0,p-1,p-1) +
(EOm(i,j,0)*Es2(j)+EOm(j,i,0)*Es2(i))*(Sig.slice(j).submat(0,0,p-1,p-1)+
mu.slice(j).col(0)*trans(mu.slice(j).col(0)));
}
}
const1=a0/b0+ai2(i)*bi2(i);
for(int pp=0;pp<p;pp++){
Sig(pp,pp,i) = Sig(pp,pp,i)+ const1;
}
Sig.slice(i).submat(0,0,p-1,p-1) = inv(Sig.slice(i).submat(0,0,p-1,p-1));

mu.slice(i).col(0) = Sig.slice(i).submat(0,0,p-1,p-1)*cvecp1;



//2 >= t < T
for(int tt=1;tt<TT-1;tt++){

cvecp1 = arma::zeros(p,1);
for(int j=0; j<n;j++){
if(j != i){
cvecp1 = cvecp1 + ( (Y(i,j,tt)-0.5)*Es1(j)+(Y(j,i,tt)-0.5)*Es1(i)- 
a3*(EOm(i,j,tt)*Es1(j)+EOm(j,i,tt)*Es1(i)))*mu.slice(j).col(tt);
}
}
cvecp1 = cvecp1 + ai2(i)*bi2(i)*(mu.slice(i).col(tt-1)+mu.slice(i).col(tt+1));

Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1) = arma::zeros(p,p);
for(int j=0;j<n;j++){
if(j != i){
Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1) = Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1) +
(EOm(i,j,tt)*Es2(j)+EOm(j,i,tt)*Es2(i))*(Sig.slice(j).submat(tt*p,0,(tt+1)*p-1,p-1)+
mu.slice(j).col(tt)*trans(mu.slice(j).col(tt)));
}
}
const1=2*ai2(i)*bi2(i);

for(int pp=0;pp<p;pp++){
Sig(tt*p + pp,pp,i) = Sig(tt*p+pp,pp,i)+ const1;
}
Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1) = inv(Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1));
mu.slice(i).col(tt) = Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1)*cvecp1;

}



//t=T
cvecp1 = arma::zeros(p,1);
for(int j=0; j<n;j++){
if(j != i){
cvecp1 = cvecp1 + ( (Y(i,j,TT-1)-0.5)*Es1(j)+(Y(j,i,TT-1)-0.5)*Es1(i)- 
a3*(EOm(i,j,TT-1)*Es1(j)+EOm(j,i,TT-1)*Es1(i)))*mu.slice(j).col(TT-1);
}
}
cvecp1 = cvecp1 + ai2(i)*bi2(i)*mu.slice(i).col(TT-2);

Sig.slice(i).submat((TT-1)*p,0,TT*p-1,p-1) = arma::zeros(p,p);
for(int j=0;j<n;j++){
if(j != i){
Sig.slice(i).submat((TT-1)*p,0,TT*p-1,p-1) = Sig.slice(i).submat((TT-1)*p,0,TT*p-1,p-1) +
(EOm(i,j,TT-1)*Es2(j)+EOm(j,i,TT-1)*Es2(i))*(Sig.slice(j).submat((TT-1)*p,0,TT*p-1,p-1)+
mu.slice(j).col(TT-1)*trans(mu.slice(j).col(TT-1)));
}
}
const1=ai2(i)*bi2(i);
for(int pp=0;pp<p;pp++){
Sig((TT-1)*p+pp,pp,i) = Sig((TT-1)*p+pp,pp,i)+ const1;
}
Sig.slice(i).submat((TT-1)*p,0,TT*p-1,p-1) = inv(Sig.slice(i).submat((TT-1)*p,0,TT*p-1,p-1));

mu.slice(i).col(TT-1) = Sig.slice(i).submat((TT-1)*p,0,TT*p-1,p-1)*cvecp1;

}


// tau_i ----------------------------------------

for(int i=0;i<n;i++){

ai2(i) = a2Star +p*(TT-1)/2.0;
bi2(i) = 1/b2Star;

const1=0;
for(int tt=1;tt<TT;tt++){
mat1by1 = trans(mu.slice(i).col(tt))*mu.slice(i).col(tt)+
trans(mu.slice(i).col(tt-1))*mu.slice(i).col(tt-1)-
2*trans(mu.slice(i).col(tt))*mu.slice(i).col(tt-1);
bi2(i) = bi2(i) + 0.5*(
arma::trace(Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1))+
arma::trace(Sig.slice(i).submat((tt-1)*p,0,tt*p-1,p-1))+
mat1by1(0) );
}

bi2(i) = 1/bi2(i);

}


// alpha ----------------------------------------

a3=0;
b3=0;

for(int tt=0;tt<TT;tt++){
for(int i=0;i<n-1;i++){
for(int j=i+1;j<n;j++){
if(j != i){
mat1by1 = (EOm(i,j,tt)*Es1(j) + EOm(j,i,tt)*Es1(i))*
trans(mu.slice(i).col(tt))*mu.slice(j).col(tt);
a3 = a3 + Y(i,j,tt)+Y(j,i,tt)-1.0-mat1by1(0);
b3= b3 + EOm(i,j,tt) + EOm(j,i,tt);
}
}
}
}
b3 = 1/(b3 +1/b3Star);
a3 = a3*b3;


// s_j ----------------------------------------

for(int j=0;j<n;j++){

const1 = 0;
const2 = 0;
for(int tt=0;tt<TT;tt++){
cvecp1 = arma::zeros(p,1);
for(int i=0;i<n;i++){
if(i != j){
cvecp1 = cvecp1 + (Y(i,j,tt)-0.5-a3*EOm(i,j,tt))*mu.slice(i).col(tt);

mat1by1 = trans(mu.slice(j).col(tt))*mu.slice(i).col(tt);
mat1by1 = mat1by1*mat1by1 +
trans(mu.slice(j).col(tt))*Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1)*mu.slice(j).col(tt)+
trans(mu.slice(i).col(tt))*Sig.slice(j).submat(tt*p,0,(tt+1)*p-1,p-1)*mu.slice(i).col(tt);
const2 += EOm(i,j,tt)*(
arma::trace(Sig.slice(j).submat(tt*p,0,(tt+1)*p-1,p-1)*Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1)) +
mat1by1(0) );
}
}
mat1by1 = trans(mu.slice(j).col(tt))*cvecp1;
const1 += mat1by1(0);
}
const2= 1/const2;
const1 = const2*(const1-1);

ai4(j) = const1;
bi4(j) = const2;
Es1(j) = const1+sqrt(const2)*R::dnorm(const1/sqrt(const2),0.0,1.0,0)/
R::pnorm(const1/sqrt(const2),0.0,1.0,1,0);
Es2(j) = const1*const1 + const2 +const1*sqrt(const2)*R::dnorm(const1/sqrt(const2),0.0,1.0,0)/
R::pnorm(const1/sqrt(const2),0.0,1.0,1,0);

}

//Recompute with constraint in place
/*
const3=n/sum(Es1);
Es1 = Es1*const3;
Es2 = Es2*const3*const3;
*/

// sigma^2 --------------------------------

a0= a0Star+n*p*0.5;
b0=b0Star;
for(int i=0;i<n;i++){
mat1by1 = trans(mu.slice(i).col(0))*mu.slice(i).col(0);
b0 = b0 + 0.5*(mat1by1(0) + 
arma::trace(Sig.slice(i).submat(0,0,p-1,p-1)));
}

return Rcpp::List::create(EOm,mu,Sig,ai2,bi2,a0,b0,a3,b3,Es1,Es2,ai4,bi4);

}