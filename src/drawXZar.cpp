#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List drawXZar(const IntegerVector & dims,
                    const NumericVector & XVec,
                    const NumericVector & Yvec,
                    const NumericVector & OmVec,
                    const NumericVector & SS1,
                    const NumericVector & ALPHA,
                    const NumericVector & TAU,
                    const NumericVector & RR,
                    const NumericMatrix & UU,
                    const IntegerMatrix & ZZ,
                    const NumericVector & rnVec,
                    const NumericMatrix & BB,
                    const NumericVector & RNa,
                    const NumericVector & B3Star,
                    const NumericVector & CC){
  const double PI = 3.14159265358979323846;
                        
Environment base("package:base");
Function sampleInt = base["sample"];

int n=dims(0),p=dims(1), TT=dims(2), M=dims(3);

arma::cube X(XVec.begin(),p,TT,n);
arma::cube Y(Yvec.begin(),n,n,TT);
arma::cube Om(OmVec.begin(),n,n,TT);
arma::colvec ss1=Rcpp::as<arma::colvec>(SS1);
double alpha=Rcpp::as<double>(ALPHA);
arma::colvec tau=Rcpp::as<arma::colvec>(TAU);
arma::colvec rr=Rcpp::as<arma::colvec>(RR);
arma::mat uu = Rcpp::as<arma::mat>(UU);
arma::cube rn(rnVec.begin(),p,TT,n);
arma::mat bb = Rcpp::as<arma::mat>(BB);
double rna=Rcpp::as<double>(RNa);
double b3Star=Rcpp::as<double>(B3Star);
double cc= Rcpp::as<double>(CC);
arma::colvec ai1 = arma::zeros(n,1);
Rcpp::IntegerMatrix Z(ZZ);

//Nuisance variables
//for X
arma::mat Sigit = arma::zeros(p,p);
arma::vec eigval;
arma::mat eigvec;
arma::colvec mat1by1 =arma::zeros(1,1);
arma::colvec cvecp1 = arma::zeros(p,1);
arma::mat B = arma::zeros(p,p);
//For Z
Rcpp::NumericVector probs(M);
Rcpp::IntegerVector Ztemp(1);
//For alpha
double a3=0, b3=0;


// X_{it} ----------------------------------------
for(int i=0; i<n;i++){
for(int tt=0;tt<TT;tt++){

cvecp1 = rr(i)*tau(i)*uu.col(Z(i,tt));
Sigit = arma::zeros(p,p);
Sigit.diag() += tau(i);
for(int j=0; j<n;j++){
if(j != i){
cvecp1 = cvecp1 + 
( (Y(i,j,tt)-0.5-alpha*Om(i,j,tt))*ss1(j)+(Y(j,i,tt)-0.5-alpha*Om(j,i,tt))*ss1(i) )*X.slice(j).col(tt);
Sigit = Sigit +
(Om(i,j,tt)*ss1(j)*ss1(j) + Om(j,i,tt)*ss1(i)*ss1(i))*X.slice(j).col(tt)*trans(X.slice(j).col(tt));
}
}
Sigit = inv_sympd(Sigit);
cvecp1 = Sigit*cvecp1;
eig_sym(eigval,eigvec,Sigit);
B = eigvec*sqrt(diagmat(eigval));
X.slice(i).col(tt) = cvecp1 + B*rn.slice(i).col(tt);

}
}

// Z_{it} ----------------------------------------
for(int i=0;i<n;i++){

for(int g=0;g<M;g++){
probs(g) = bb(0,g)*pow(tau(i),0.5*p)*
exp(-0.5*tau(i)*pow(arma::norm(X.slice(i).col(0)-rr(i)*uu.col(g),2),2));
}
Ztemp = sampleInt(M,Named("size",1),Named("prob",probs));
Z(i,0)=Ztemp(0)-1;

for(int tt=1;tt<TT;tt++){
for(int g=0;g<M;g++){
probs(g) = bb(Z(i,tt-1)+1,g)*pow(tau(i),0.5*p)*
exp(-0.5*tau(i)*pow(arma::norm(X.slice(i).col(tt)-rr(i)*uu.col(g),2),2));
}
Ztemp = sampleInt(M,Named("size",1),Named("prob",probs));
Z(i,tt)=Ztemp(0)-1;
}

}


// alpha ----------------------------------------
for(int tt=0;tt<TT;tt++){
for(int i=0;i<n-1;i++){
for(int j=i+1;j<n;j++){
if(j != i){
mat1by1 = (Om(i,j,tt)*ss1(j) + Om(j,i,tt)*ss1(i))*
trans(X.slice(i).col(tt))*X.slice(j).col(tt);
a3 = a3 + Y(i,j,tt)+Y(j,i,tt)-1.0-mat1by1(0);
b3= b3 + Om(i,j,tt) + Om(j,i,tt);
}
}
}
}
b3 = 1/(b3 +1/b3Star);
a3 = a3*b3;
alpha = a3 + sqrt(b3)*rna;


// rr ----------------------------------------
for(int i=0;i<n;i++){
for(int tt=0;tt<TT;tt++){
mat1by1 =trans(X.slice(i).col(tt))*uu.col(Z(i,tt));
ai1(i) = ai1(i)+ mat1by1(0);
}
ai1(i)=1.0*(ai1(i)-1.0/cc)/TT;
}

return Rcpp::List::create(X,Z,alpha,ai1);

}
