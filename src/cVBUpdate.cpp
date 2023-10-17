#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cVBUpdate(const IntegerVector & DIMS,
            const NumericVector & Yvec,
            const NumericVector & CC,
            const NumericMatrix & GAMSTAR,
            const NumericVector & A2STAR,
            const NumericVector & B2STAR,
            const NumericVector & B3STAR,
            const NumericVector & EOmVec,
            const NumericVector & muVec,
            const NumericVector & SigVec,
            const NumericMatrix & BI0G,
            const NumericVector & BitbarVec,
            const NumericVector & BithkVec,
            const NumericVector & ER,
            const NumericVector & ER2,
            const NumericVector & AI2,
            const NumericVector & BI2,
            const NumericMatrix & NU,
            const NumericVector & A3,
            const NumericVector & B3,
            const NumericVector & ES1,
            const NumericVector & ES2,
            const NumericMatrix & GAM){
  const double PI = 3.14159265358979323846;
  
  RNGScope scope;
  Environment base("package:base");
  Function digamma = base["digamma"];
  
  
  //dims = c(n,p,TT,MM)
  
  //Variables to be read in:
  Rcpp::IntegerVector dims(DIMS);
  int n=dims(0),p=dims(1), TT=dims(2), M=dims(3);
  arma::cube Y(Yvec.begin(),dims[0],dims[0],dims[2]);
  double cc = Rcpp::as<double>(CC);
  double a2Star = Rcpp::as<double>(A2STAR);
  double b2Star = Rcpp::as<double>(B2STAR);
  double b3Star = Rcpp::as<double>(B3STAR);
  
  arma::cube EOm(EOmVec.begin(),n,n,TT);  
  arma::cube mu(muVec.begin(),p,TT,n);
  arma::cube Sig(SigVec.begin(),TT*p,p,n);
  arma::mat Bi0g = Rcpp::as<arma::mat>(BI0G);
  arma::cube Bithk(BithkVec.begin(),TT*M,M,n);
  arma::colvec Er= Rcpp::as<arma::colvec>(ER);
  arma::colvec Er2= Rcpp::as<arma::colvec>(ER2);
  arma::colvec ai1 = arma::zeros(n,1);
  arma::colvec bi1 = arma::zeros(n,1);
  arma::colvec ai2= Rcpp::as<arma::colvec>(AI2);
  arma::colvec bi2= Rcpp::as<arma::colvec>(BI2);
  arma::colvec ai4 = arma::zeros(n,1);
  arma::colvec bi4 = arma::zeros(n,1);
  arma::mat nu = Rcpp::as<arma::mat>(NU);
  arma::colvec kappas = arma::zeros(M,1);
  double a3=Rcpp::as<double>(A3);
  double b3=Rcpp::as<double>(B3);
  arma::colvec Es1=Rcpp::as<arma::colvec>(ES1);
  arma::colvec Es2=Rcpp::as<arma::colvec>(ES2);
  arma::cube Bitbar(BitbarVec.begin(),TT,M,n);
  arma::mat GamStar = Rcpp::as<arma::mat>(GAMSTAR);
  arma::mat Gam = Rcpp::as<arma::mat>(GAM);
  arma::cube cijt = arma::zeros(n,n,TT);
  
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
  cijt(i,j,tt) = const1;
  cijt(j,i,tt) = const2;
  
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
  cvecp2 = arma::zeros(p,1);
  for(int gg=0;gg<M;gg++){
  cvecp2 = cvecp2 + Bi0g(i,gg)*trans(nu.row(gg));
  }
  cvecp1 = cvecp1 + Er(i)*ai2(i)*bi2(i)*cvecp2;
  
  Sig.slice(i).submat(0,0,p-1,p-1) = arma::zeros(p,p);
  for(int j=0;j<n;j++){
  if(j != i){
  Sig.slice(i).submat(0,0,p-1,p-1) = Sig.slice(i).submat(0,0,p-1,p-1) +
  (EOm(i,j,0)*Es2(j)+EOm(j,i,0)*Es2(i))*(Sig.slice(j).submat(0,0,p-1,p-1)+
  mu.slice(j).col(0)*trans(mu.slice(j).col(0)));
  }
  }
  const1=0;
  for(int gg=0;gg<M;gg++){
  const1 += Bi0g(i,gg);
  }
  for(int pp=0;pp<p;pp++){
  Sig(pp,pp,i) = Sig(pp,pp,i)+ const1*ai2(i)*bi2(i);
  }
  Sig.slice(i).submat(0,0,p-1,p-1) = inv(Sig.slice(i).submat(0,0,p-1,p-1));
  
  mu.slice(i).col(0) = Sig.slice(i).submat(0,0,p-1,p-1)*cvecp1;
  
  
  //t>=2
  for(int tt=1;tt<TT;tt++){
  
  cvecp1 = arma::zeros(p,1);
  for(int j=0; j<n;j++){
  if(j != i){
  cvecp1 = cvecp1 + ( (Y(i,j,tt)-0.5)*Es1(j)+(Y(j,i,tt)-0.5)*Es1(i)- 
  a3*(EOm(i,j,tt)*Es1(j)+EOm(j,i,tt)*Es1(i)))*mu.slice(j).col(tt);
  }
  }
  cvecp2 = arma::zeros(p,1);
  for(int hh=0;hh<M;hh++){
  for(int kk=0;kk<M;kk++){
  cvecp2 = cvecp2 + Bitbar(tt-1,hh,i)*
  Bithk(tt*M+hh,kk,i)*trans(nu.row(kk));
  }
  }
  cvecp1 = cvecp1 + Er(i)*ai2(i)*bi2(i)*cvecp2;
  
  Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1) = arma::zeros(p,p);
  for(int j=0;j<n;j++){
  if(j != i){
  Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1) = Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1) +
  (EOm(i,j,tt)*Es2(j)+EOm(j,i,tt)*Es2(i))*(Sig.slice(j).submat(tt*p,0,(tt+1)*p-1,p-1)+
  mu.slice(j).col(tt)*trans(mu.slice(j).col(tt)));
  }
  }
  const1=0;
  for(int hh=0;hh<M;hh++){
  for(int kk=0;kk<M;kk++){
  const1 += Bitbar(tt-1,hh,i)*Bithk(tt*M+hh,kk,i);
  }
  }
  for(int pp=0;pp<p;pp++){
  Sig(tt*p+pp,pp,i) = Sig(tt*p+pp,pp,i)+ const1*ai2(i)*bi2(i);
  }
  Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1) = inv(Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1));
  
  mu.slice(i).col(tt) = Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1)*cvecp1;
  
  }
  }
  
  
  // Z_{it} ----------------------------------------
  
  for(int i=0;i<n;i++){
  
  //Beta_{i0g}
  for(int gg=0;gg<M;gg++){
  mat1by1 = Er(i)*ai2(i)*bi2(i)*nu.row(gg)*mu.slice(i).col(0);
  Bi0g(i,gg) = Rcpp::as<double>(digamma(Gam(0,gg))) +  mat1by1(0);
  }
  Bi0g.row(i) = exp(Bi0g.row(i) - max(Bi0g.row(i)));
  const1 = sum(Bi0g.row(i));
  for(int gg=0;gg<M;gg++){
  Bi0g(i,gg) = Bi0g(i,gg)/const1;
  Bitbar(0,gg,i) = Bi0g(i,gg);
  }
  
  //Beta_{ithk}
  for(int tt=1;tt<TT;tt++){
  for(int hh=0;hh<M;hh++){
  for(int kk=0;kk<M;kk++){
  mat1by1 = Er(i)*ai2(i)*bi2(i)*nu.row(kk)*mu.slice(i).col(tt);
  Bithk(tt*M+hh,kk,i) = Rcpp::as<double>(digamma(Gam(hh+1,kk))) +  mat1by1(0) ;
  }
  Bithk.slice(i).row(tt*M+hh) = exp(Bithk.slice(i).row(tt*M+hh)-max(Bithk.slice(i).row(tt*M+hh)));
  
  const1 = sum(Bithk.slice(i).row(tt*M+hh));
  for(int kk=0;kk<M;kk++){
  Bithk(tt*M+hh,kk,i) = Bithk(tt*M+hh,kk,i)/const1;
  }
  }
  Bitbar.slice(i).submat(tt,0,tt,M-1) =
  Bitbar.slice(i).submat(tt-1,0,tt-1,M-1)*
  Bithk.slice(i).submat(tt*M,0,tt*M+M-1,M-1);
  }
  
  }
  
  
  
  // r_i ----------------------------------------
  for(int i=0;i<n;i++){
  
  mat1by1(0)=0.0;
  for(int gg=0;gg<M;gg++){
  mat1by1 = mat1by1+ 
  Bi0g(i,gg)*nu.row(gg)*mu.slice(i).col(0);
  }
  for(int tt=1;tt<TT;tt++){
  for(int hh=0;hh<M;hh++){
  for(int kk=0;kk<M;kk++){
  mat1by1 = mat1by1 +
  Bitbar(tt-1,hh,i)*Bithk(tt*M+hh,kk,i)*
  nu.row(kk)*mu.slice(i).col(tt);
  }
  }
  }
  ai1(i) = 1.0*(mat1by1(0) -1.0/cc)/TT;
  
  const2=1.0/(TT*ai2(i)*bi2(i));
  bi1(i) = const2;  

  Er(i) = ai1(i)+sqrt(const2)*R::dnorm(ai1(i)/sqrt(const2),0.0,1.0,0)/
  R::pnorm(ai1(i)/sqrt(const2),0.0,1.0,1,0);
  Er2(i) = pow(ai1(i),2)+const2 +ai1(i)*sqrt(const2)*
  R::dnorm(ai1(i)/sqrt(const2),0.0,1.0,0)/
  R::pnorm(ai1(i)/sqrt(const2),0.0,1.0,1,0);
  
  }
  
  
  // tau_i ----------------------------------------
  
  for(int i=0;i<n;i++){
  
  ai2(i) = a2Star +p*TT/2.0+1.0;
  bi2(i) = Er(i)/cc+1/b2Star;
  
  for(int gg=0;gg<M;gg++){
  mat1by1 = trans(mu.slice(i).col(0))*mu.slice(i).col(0)-
  2*Er(i)*nu.row(gg)*mu.slice(i).col(0);
  bi2(i) = bi2(i) + 0.5*Bi0g(i,gg)*(
  arma::trace(Sig.slice(i).submat(0,0,p-1,p-1)) +
  mat1by1(0)+Er2(i) );
  }
  for(int tt=1;tt<TT;tt++){
  for(int hh=0;hh<M;hh++){
  for(int kk=0;kk<M;kk++){
  mat1by1 = trans(mu.slice(i).col(tt))*mu.slice(i).col(tt)-
  2*Er(i)*nu.row(kk)*mu.slice(i).col(tt);
  bi2(i) = bi2(i) + 0.5*Bitbar(tt-1,hh,i)*Bithk(tt*M+hh,kk,i)*(
  arma::trace(Sig.slice(i).submat(tt*p,0,(tt+1)*p-1,p-1))+
  mat1by1(0)+Er2(i) );
  }
  }
  }
  
  
  bi2(i) = 1/bi2(i);
  
  
  }
  
  
  // u_g ----------------------------------------
  
  nu = arma::zeros(M,p);
  
  for(int gg=0;gg<M;gg++){
  
  for(int i=0;i<n;i++){
  cvecp1 = arma::zeros(p,1);
  for(int tt=1;tt<TT;tt++){
  for(int hh=0;hh<M;hh++){
  cvecp1 = cvecp1+
  Bitbar(tt-1,hh,i)*Bithk(tt*M+hh,gg,i)*mu.slice(i).col(tt);
  }
  }
  
  nu.row(gg) = nu.row(gg) +
  ai2(i)*bi2(i)*Er(i)*(
  Bi0g(i,gg)*trans(mu.slice(i).col(0)) + trans(cvecp1) );
  }
  const1 = arma::norm(nu.row(gg),2);
  nu.row(gg) = nu.row(gg)/const1;
  kappas(gg) = const1;  

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
  
  ai4(j) = 0;
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
  ai4(j) += mat1by1(0);
  }
  const2= 1/const2;
  ai4(j) = const2*(ai4(j)-1);
  bi4(j) = const2;
  
  Es1(j) = ai4(j)+sqrt(const2)*R::dnorm(ai4(j)/sqrt(const2),0.0,1.0,0)/
  R::pnorm(ai4(j)/sqrt(const2),0.0,1.0,1,0);
  Es2(j) = ai4(j)*ai4(j) + const2 +ai4(j)*sqrt(const2)*R::dnorm(ai4(j)/sqrt(const2),0.0,1.0,0)/
  R::pnorm(ai4(j)/sqrt(const2),0.0,1.0,1,0);
  
  }
  
  
  //Recompute with constraint in place
  /*
  const3=n/sum(Es1);
  Es1 = Es1*const3;
  Es2 = Es2*const3*const3;
  */
  
  // beta ----------------------------------------
  
  for(int gg=0;gg<M;gg++){
  Gam(0,gg)= GamStar(0,gg);
  for(int i=0;i<n;i++){
  Gam(0,gg) = Gam(0,gg) + Bi0g(i,gg);
  }
  }
  
  for(int hh=0;hh<M;hh++){
  for(int kk=0;kk<M;kk++){
  Gam(hh+1,kk)= GamStar(hh+1,kk);
  for(int tt=1;tt<TT;tt++){
  for(int i=0;i<n;i++){
  Gam(hh+1,kk)=Gam(hh+1,kk)+
  Bitbar(tt-1,hh,i)*Bithk(tt*M+hh,kk,i);
  }
  }
  }
  }
  
/*  
  return Rcpp::List::create(EOm,mu,Sig,Bi0g,Bitbar,Bithk,Er,Er2,ai2,bi2,nu,
              a3,b3,Es1,Es2,Gam,cijt,ai1,ai4,bi4,bi1,kappas);
  */
  for(int hh=0;hh<M;hh++){
    nu.row(hh) = nu.row(hh)*kappas(hh);
  }

  return Rcpp::List::create(EOm,mu,Sig,Bi0g,Bitbar,Bithk,Er,Er2,ai2,bi2,nu,
              a3,b3,Es1,Es2,Gam,ai1,ai4,bi4,bi1);

}  

