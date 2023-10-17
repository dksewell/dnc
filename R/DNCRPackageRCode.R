require("skmeans")
require("plot3D")
require("plot3Drgl")
require("BayesLogit")
require("movMF")
require("tcltk")

.VBCluster <- function(DIMS,YY,CC,A2STAR,B2STAR,B3STAR,EOM,MU,SIG,
                       BI0G,BITBAR,BITHK,ER,ER2,AI2,BI2,NU,A3,B3,ES,ES2,
                       GAMSTAR,GAM,MaxIt=100,epsilon=1e-5){
  
  Crit <- numeric(MaxIt+1)
  Crit[1] =  cLogLikelihood(YY,MU,DIMS,A3,ES)
  
  if(MaxIt>1) pb <- tkProgressBar(min=1,max=MaxIt,title="VB Algorithm")
  TIME=system.time({
    for(iter in 1:MaxIt){
      Update <- cVBUpdate(DIMS,YY,CC,GAMSTAR,A2STAR,B2STAR,B3STAR,EOM,MU,SIG,
                          BI0G,BITBAR,BITHK,ER,ER2,AI2,BI2,NU,A3,B3,ES,ES2,GAM)
      
      non0ind <- which(abs(EOM) > 1e-10)
      non0ind1 <- which(abs(BI0G) > 1e-10)
      non0ind2 <- which(abs(BITBAR) > 1e-10)
      non0ind3 <- which(abs(BITHK) > 1e-10)
      non0ind4 <- which(abs(NU) > 1e-10)
      
      
      
      EOM <- Update[[1]]
      MU <- Update[[2]]
      SIG <- Update[[3]]
      BI0G <- Update[[4]]
      BITBAR <- Update[[5]]
      BITHK <- Update[[6]]
      ER <- drop(Update[[7]])
      ER2 <- drop(Update[[8]])
      AI2 <- drop(Update[[9]])
      BI2 <- drop(Update[[10]])
      NU <- Update[[11]]
      A3 <- Update[[12]]
      B3 <- Update[[13]]
      ES <- drop(Update[[14]])
      ES2 <- drop(Update[[15]])
      GAM <- Update[[16]]
      AI1 <- Update[[17]]
      AI4 <- Update[[18]]
      BI4 <- Update[[19]]
      BI1 <- Update[[20]]
      
      KAPPAS <- apply(NU,1,function(x)sqrt(x%*%x))
      NU = t(apply(NU,1,function(x)x/sqrt(x%*%x)))
      
      
      Crit[iter+1] =  cLogLikelihood(YY,MU,DIMS,A3,ES)
      
      if(MaxIt>1) setTkProgressBar(pb,iter)
      try({if( abs(Crit[iter+1]-Crit[iter])/abs(Crit[iter]) < epsilon ) break },silent=TRUE)
    }
    close(pb)
  })
  
  return(list(EOm=EOM,mu=MU,Sig=SIG,Bi0g=BI0G,Bitbar=BITBAR,Bithk=BITHK,
              Er=ER,Er2=ER2,ai1=AI1,bi1=BI1,ai2=AI2,ai4=AI4,bi4=BI4,
              bi2=BI2,nu=NU,kappa=KAPPAS,a3=A3,b3=B3,Es=ES,
              Es2=ES2,Gam=GAM,Crit=Crit[1:(iter+1)],Time=TIME))#cijt=CIJT,
}


# VB No clustering --------------------------------------------------------

.VBnoClustering <- function(DIMS,YY,MU,SIG,A2STAR,B2STAR,B3STAR,
                            EOM,A0,B0,A0STAR,B0STAR,AI2,BI2,A3,
                            B3,ES,ES2,MaxIt=1000,epsilon=1e-5){
  Crit1 <- numeric(MaxIt+1)
  Crit1[1] <- cLogLikelihood(YY,MU,DIMS,A3,ES)
  
  pb <- tkProgressBar(min=1,max=MaxIt,title="VB without clustering (stage 2 initialization if M>0)")
  TIME=system.time({
    for(iter in 1:MaxIt){
      Update <- cVBnoClust(DIMS,YY,MU,SIG,A2STAR,B2STAR,B3STAR,
                           EOM,A0,B0,A0STAR,B0STAR,AI2,BI2,A3,
                           B3,ES,ES2)
      
      if(sum(is.na(unlist(Update)))>0){
        print("NAs in Update")
        break
      }else{
        EOM <- Update[[1]]
        MU <- Update[[2]]
        SIG <- Update[[3]]
        AI2 <- drop(Update[[4]])
        BI2 <- drop(Update[[5]])
        A0 <- Update[[6]]
        B0 <- Update[[7]]
        A3 <- Update[[8]]
        B3 <- Update[[9]]
        ES <- drop(Update[[10]])
        ES2 <- drop(Update[[11]])
        AI4 <-drop(Update[[12]])
        BI4 <-drop(Update[[13]])
      }
      setTkProgressBar(pb,iter)
      Crit1[iter+1] <- cLogLikelihood(YY,MU,DIMS,A3,ES)
      if(abs(Crit1[iter+1]-Crit1[iter])/abs(Crit1[iter]) < epsilon) break  
    }
  })
  close(pb)
  
  return(list(EOm=EOM,mu=MU,Sig=SIG,ai2=AI2,bi2=BI2,a0=A0,b0=B0,a3=A3,
              b3=B3,Es=ES,Es2=ES2,Crit=Crit1[1:(iter+1)],Time=TIME,
              ai4=AI4,bi4=BI4))
  
}


# Projection Clustering (main) ---------------------------------------------

.clustProjectionInitialize = function(Y,M,p=3,hyperparms=NULL,
                                      MaxIt=1000,epsilon=1e-15){
  n = dim(Y)[1]
  TT = dim(Y)[3]
  if(is.null(hyperparms)){
    cc=NULL;a2Star=NULL;b2Star=NULL;b3Star=NULL;
    GamStar=NULL;a0Star=NULL;b0Star=NULL
  }else{
    cc=hyperparms$cc;a2Star=hyperparms$a2Star;
    b2STAR=hyperparms$b2STAR;b3Star=hyperparms$b3Star;
    GamStar=hyperparms$GamStar;a0Star=hyperparms$a0Star;b0Star=hyperparms$b0Star
  }
  
  mu = array(0,c(p,TT,n))
  Sig = array(0,c(TT*p,p,n))
  for(i in 1:n){
    for(tt in 1:TT){
      Sig[((tt-1)*p+1):(tt*p),,i] = diag(p)
    }
  }
  EOm = array(0,c(n,n,TT))
  a3 = sum(c(Y))/(TT*n*(n-1))
  a3 = log(a3/(1-a3))
  if(is.null(b3Star)){
    b3=b3Star=100
  }else{
    b3=b3Star
  }
  
  Er2 = Er = tau = Es2 = Es = numeric(n)
  
  ###Find s and rtilde
  for(tt in 1:TT){
    Er <- Er + apply(Y[,,tt],1,sum)+1/TT
    Es <- Es + apply(Y[,,tt],2,sum)+1/TT
  }
  Es <- Es/Er
  Es2 = Es^2
  Er <- Er/sum(Er)
  
  ###Find direction for X
  dissim <- array(0,dim=dim(Y))
  for(tt in 1:TT) dissim[,,tt] <- igraph::shortest.paths(graph=igraph::graph.adjacency(Y[,,tt]),mode="all")
  dissim[which(dissim==Inf,arr.ind=TRUE)] <- max(c(dissim[which(dissim!=Inf,arr.ind=TRUE)]))
  VV <- array(0,c(p,TT,n))
  VV[,1,] <- t(cmdscale(d=dissim[,,1],k=p))
  temp.lambda <- 10
  H <- matrix(-1/n,n,n)+diag(n)
  
  for(tt in 2:TT){
    temp <- 1/(1+temp.lambda)*H%*%(-1/2*dissim[,,tt]^2)%*%H +
      temp.lambda/(1+temp.lambda)*t(VV[,tt-1,])%*%VV[,tt-1,]
    temp <- eigen(temp)
    VV[,tt,] <- t(temp$vectors[,1:p]%*%diag(temp$values[1:p])^(1/2))
    VV[,tt,] <- t(vegan::procrustes(X=t(VV[,tt-1,]),Y=t(VV[,tt,]),scale=FALSE)$Yrot)
  }
  dissim[which(is.na(dissim))] <- 2
  
  for(tt in 1:TT){
    VV[,tt,] <- scale(VV[,tt,],center=FALSE,
                      scale=apply(VV[,tt,],2,function(x) sqrt(x%*%x)))
  }
  
  ###Find r and X
  for(tt in 1:TT){
    mu[,tt,] <- matrix(Er,p,n,byrow=TRUE)*VV[,tt,]
  }
  .llik.wrap <- function(x){
    cLogLikelihood(Y,x*mu,c(n,p,TT),a3,Es)
  }
  Optim <- optimize(f=.llik.wrap,interval=c(1e-20,200),maximum=TRUE)
  Er <- Er*Optim$max
  mu <- mu*Optim$max
  rm(Optim,.llik.wrap,VV)
  for(tt in 2:TT){
    mu[,tt,] <- t(vegan::procrustes(X=t(mu[,tt-1,]),Y=t(mu[,tt,]),scale=FALSE)$Yrot)
  }
  
  
  ###Find tau
  if(is.null(cc))cc=mean(16/Er)
  if(is.null(b2Star)){
    b2Star=1
  }
  if(is.null(a2Star)){
    a2Star=mean((4/Er)^2)/b2Star#1
  }
  ai2 <- rep(a2Star,n)
  bi2= rep(b2Star,n)
  
  
  ###Find Omega
  for(tt in 1:TT){
    for(i in 1:n){
      for(j in c(1:n)[-i]){
        eta <- a3 + Es[j]*mu[,tt,i]%*%mu[,tt,j]
        EOm[i,j,tt] <- 1/(2*eta)*(exp(eta)-1)/(1+exp(eta))
      }
    }
  }
  
  
  ###Hyperparms for initial variance
  if(is.null(a0Star)){
    a0Star <- 2+0.01
  }
  if(is.null(b0Star)){
    b0Star <- cov(t(mu[,1,]))
    b0Star <- mean(diag(b0Star))*(a0Star-1)
  }
  
  ###VB no clustering
  init1 <- .VBnoClustering(DIMS=c(n,p,TT),YY=Y,MU=mu,SIG=Sig,A2STAR=a2Star,B2STAR=b2Star,
                           B3STAR=b3Star,EOM=EOm,A0=a0Star,B0=b0Star,A0STAR=a0Star,B0STAR=b0Star,
                           AI2=ai2,BI2=bi2,A3=a3,B3=b3,ES=Es,ES2=Es2,epsilon=epsilon,MaxIt=MaxIt)
  
  EOm <- init1$EOm
  mu <- init1$mu
  Sig <- init1$Sig
  ai2 <- init1$ai2
  bi2 <- init1$bi2
  a0 <- init1$a0
  b0 <- init1$b0
  a3 <- init1$a3
  b3 <- init1$b3
  Es <- init1$Es
  Es2 <- init1$Es2
  
  ###Assign clusters
  VV <- mu
  for(tt in 1:TT){
    VV[,tt,] <- t(t(VV[,tt,])/apply(t(VV[,tt,]),1,function(x) sqrt(x%*%x)))
  }
  
  VAggreg <- matrix(0.0,TT*n,p)
  for(tt in 1:TT){
    VAggreg[(tt-1)*n+1:n,] <- t(VV[,tt,])
  }
  
  Z <- matrix(0,n,TT)
  km <- skmeans(x=VAggreg,k=M)
  for(tt in 1:TT){
    Z[,tt] <- km$cluster[(n*(tt-1)+1):(n*tt)]
  }
  nu <- km$proto
  
  if(is.null(GamStar)){
    GamStar = matrix(0,M+1,M)
    GamStar[1,] <- 1.5
    for(hh in 1:M){
      GamStar[hh+1,hh] <- 4*(M-1)*1.5
      GamStar[hh+1,-hh] <- 1.5
    }
  }
  Gam = GamStar
  
  Bi0g <- matrix(0.2/(M-1),n,M)
  for(i in 1:n)Bi0g[i,Z[i,1]] <- 0.8
  Bithk <- array(0.2/(M-1),c(TT*M,M,n))
  for(i in 1:n){
    for(tt in 1:TT){
      for(h in 1:M){
        Bithk[(tt-1)*M+h,Z[i,tt],i] <- 0.8
      }
    }
  }
  Bitbar <- array(0,c(TT,M,n))
  for(i in 1:n){
    Bitbar[1,1:M,i] <- Bi0g[i,]
    for(tt in 2:TT){
      Bitbar[tt,,i] <-
        t(Bithk[(tt-1)*M+1:M,,i])%*%
        Bitbar[tt-1,,i]
    }
  }
  
  
  return(list(cc=cc,a2Star=a2Star,b2Star=b2Star,b3Star=b3Star,EOm=EOm,mu=mu,
              Sig=Sig,Bi0g=Bi0g,Bitbar=Bitbar,Bithk=Bithk,Er=Er,
              Er2=Er2,ai2=ai2,bi2=bi2,nu=nu,a3=a3,b3=b3,Es=Es,Es2=Es2,
              GamStar=GamStar,Gam=Gam))
  
}


.clustProjectionInitializeNoClust = function(Y,p=3,hyperparms=NULL,
                                             MaxIt=5000,epsilon=1e-5){
  n = dim(Y)[1]
  TT = dim(Y)[3]
  if(is.null(hyperparms)){
    cc=NULL;a2Star=NULL;b2Star=NULL;b3Star=NULL;
    GamStar=NULL;a0Star=NULL;b0Star=NULL
  }else{
    cc=hyperparms$cc;a2Star=hyperparms$a2Star;
    b2STAR=hyperparms$b2STAR;b3Star=hyperparms$b3Star;
    a0Star=hyperparms$a0Star;b0Star=hyperparms$b0Star
  }
  
  mu = array(0,c(p,TT,n))
  Sig = array(0,c(TT*p,p,n))
  for(i in 1:n){
    for(tt in 1:TT){
      Sig[((tt-1)*p+1):(tt*p),,i] = diag(p)
    }
  }
  EOm = array(0,c(n,n,TT))
  a3 = sum(c(Y))/(TT*n*(n-1))
  a3 = log(a3/(1-a3))
  if(is.null(b3Star)){
    b3=b3Star=100
  }else{
    b3=b3Star
  }
  
  Er2 = Er = tau = Es2 = Es = numeric(n)
  
  ###Find s and rtilde
  for(tt in 1:TT){
    Er <- Er + apply(Y[,,tt],1,sum)+1/TT
    Es <- Es + apply(Y[,,tt],2,sum)+1/TT
  }
  Es <- Es/Er
  Es2 = Es^2
  Er <- Er/sum(Er)
  
  ###Find direction for X
  dissim <- array(0,dim=dim(Y))
  for(tt in 1:TT) dissim[,,tt] <- igraph::shortest.paths(graph=igraph::graph.adjacency(Y[,,tt]),mode="all")
  dissim[which(dissim==Inf,arr.ind=TRUE)] <- max(c(dissim[which(dissim!=Inf,arr.ind=TRUE)]))
  VV <- array(0,c(p,TT,n))
  VV[,1,] <- t(cmdscale(d=dissim[,,1],k=p))
  temp.lambda <- 10
  H <- matrix(-1/n,n,n)+diag(n)
  
  for(tt in 2:TT){
    temp <- 1/(1+temp.lambda)*H%*%(-1/2*dissim[,,tt]^2)%*%H +
      temp.lambda/(1+temp.lambda)*t(VV[,tt-1,])%*%VV[,tt-1,]
    temp <- eigen(temp)
    VV[,tt,] <- t(temp$vectors[,1:p]%*%diag(temp$values[1:p])^(1/2))
    VV[,tt,] <- t(vegan::procrustes(X=t(VV[,tt-1,]),Y=t(VV[,tt,]),scale=FALSE)$Yrot)
  }
  dissim[which(is.na(dissim))] <- 2
  
  for(tt in 1:TT){
    VV[,tt,] <- scale(VV[,tt,],center=FALSE,
                      scale=apply(VV[,tt,],2,function(x) sqrt(x%*%x)))
  }
  
  ###Find r and X
  for(tt in 1:TT){
    mu[,tt,] <- matrix(Er,p,n,byrow=TRUE)*VV[,tt,]
  }
  .llik.wrap <- function(x){
    cLogLikelihood(Y,x*mu,c(n,p,TT),a3,Es)
  }
  Optim <- optimize(f=.llik.wrap,interval=c(1e-20,200),maximum=TRUE)
  Er <- Er*Optim$max
  mu <- mu*Optim$max
  rm(Optim,.llik.wrap,VV)
  for(tt in 2:TT){
    mu[,tt,] <- t(vegan::procrustes(X=t(mu[,tt-1,]),Y=t(mu[,tt,]),scale=FALSE)$Yrot)
  }
  
  
  ###Find tau
  if(is.null(cc))cc=mean(16/Er)
  if(is.null(b2Star)){
    b2Star=1
  }
  if(is.null(a2Star)){
    a2Star=mean((4/Er)^2)/b2Star#1
  }
  ai2 <- rep(a2Star,n)
  bi2= rep(b2Star,n)
  
  
  ###Find Omega
  for(tt in 1:TT){
    for(i in 1:n){
      for(j in c(1:n)[-i]){
        eta <- a3 + Es[j]*mu[,tt,i]%*%mu[,tt,j]
        EOm[i,j,tt] <- 1/(2*eta)*(exp(eta)-1)/(1+exp(eta))
      }
    }
  }
  
  
  ###Hyperparms for initial variance
  if(is.null(a0Star)){
    a0Star <- 2+0.01
  }
  if(is.null(b0Star)){
    b0Star <- cov(t(mu[,1,]))
    b0Star <- mean(diag(b0Star))*(a0Star-1)
  }
  
  ###VB no clustering
  init1 <- .VBnoClustering(DIMS=c(n,p,TT),YY=Y,MU=mu,SIG=Sig,A2STAR=a2Star,B2STAR=b2Star,
                           B3STAR=b3Star,EOM=EOm,A0=a0Star,B0=b0Star,A0STAR=a0Star,B0STAR=b0Star,
                           AI2=ai2,BI2=bi2,A3=a3,B3=b3,ES=Es,ES2=Es2,epsilon=epsilon,MaxIt=MaxIt)
  
  EOm <- init1$EOm
  mu <- init1$mu
  Sig <- init1$Sig
  ai2 <- init1$ai2
  bi2 <- init1$bi2
  a0 <- init1$a0
  b0 <- init1$b0
  a3 <- init1$a3
  b3 <- init1$b3
  ai4 <- init1$ai4
  bi4 <- init1$bi4
  Es <- init1$Es
  Es2 <- init1$Es2
  
  return(list(cc=cc,a2Star=a2Star,b2Star=b2Star,b3Star=b3Star,EOm=EOm,mu=mu,
              Sig=Sig,a0=a0,b0=b0,ai2=ai2,bi2=bi2,a3=a3,b3=b3,ai4=ai4,bi4=bi4,
              Es=Es,Es2=Es2))
  
}



# dnc ---------------------------------------------------------------------

dnc = function(Y,M,p=3,method="VB",init=NULL,hyperparms=NULL,Missing=NULL,
               controls=list(MaxIt=500,epsilon=1e-5,MaxItStg2=100,epsilonStg2=1e-15,
                             nDraws=10000,burnin=1000)){
  
  if(!is.null(init) & is.null(hyperparms))stop("Hyperparameters must be specified")
  if(!method%in%c("VB","Gibbs"))stop("Invalid method")
  if(class(Y)!="array")stop("Y must be an array")
  M <- as.integer(M)
  if(M < 0) stop("M must be a non-negative integer")
  Missing1 = which(is.na(Y),arr.ind=TRUE)
  if(length(Missing1)>0){
    dens = sum(na.omit(c(Y)))/(nrow(Y)*(nrow(Y)-1)*dim(Y)[3]-nrow(Missing1))
    Y[which(is.na(Y))] <- sample(1:0,size=nrow(Missing1),replace=T,
                                 prob=c(dens,1-dens))
    if(!is.null(Missing)){
      Missing = rbind(Missing,Missing1)
      Missing = Missing[!duplicated(Missing),]
    }else{
      Missing = Missing1
    }
  }
  if(!is.null(Missing)){
    if(method=="VB"){
      method="Gibbs"
      warning("Missing data is only incorporated into the Gibbs sampler. Switching methods from 'VB' to 'Gibbs'",
              immediate.=TRUE)
    }
    if(M==0){
      print("M=0 only currently compatible with VB algorithm.  Please either enter  a new M value or 'exit' to exit.")
      temp <- readline()
      if(temp=="exit"){ 
        break
      }else{
        M <- as.integer(temp)
      }
      
    }
  }
  
  
  if(is.null(controls$MaxIt))controls$MaxIt=5000
  if(is.null(controls$MaxItStg2))controls$MaxItStg2=1000
  if(is.null(controls$epsilon))controls$epsilon=1e-5
  if(is.null(controls$epsilonStg2))controls$epsilonStg2=1e-15
  if(is.null(controls$nDraws))controls$nDraws=10000
  if(is.null(controls$burnin))controls$burnin=1000
  
  
  ###VB Approach
  if(method=="VB"){
    n <- dim(Y)[1]
    TT <- dim(Y)[3]
    
    if(M==0){
      
      clust = .clustProjectionInitializeNoClust(Y=Y,p=p,hyperparms=hyperparms,
                                                MaxIt=controls$MaxIt,epsilon=controls$epsilon)
      clust$method = "VB"
      clust$Y = Y
      clust <- clust[c("method","Y","mu","Sig","a0","b0","ai2","bi2",
                       "a3","b3","ai4","bi4","Es","Es2")]
      class(clust) <- "dnc"
      clust$pm = dncMAP(clust)
      return(clust)
      
    }else{
      
      if(is.null(init)){
        init= .clustProjectionInitialize(Y=Y,M=M,p=p,hyperparms=hyperparms,
                                         MaxIt=controls$MaxItStg2,
                                         epsilon=controls$epsilonStg2)
        clust=.VBCluster(DIMS=c(n,p,TT,M),YY=Y,CC=init$cc,A2STAR=init$a2Star,B2STAR=init$b2Star,
                         B3STAR=init$b3Star,EOM=init$EOm,MU=init$mu,SIG=init$Sig,BI0G=init$Bi0g,
                         BITBAR=init$Bitbar,BITHK=init$Bithk,ER=init$Er,
                         ER2=init$Er2,AI2=init$ai2,BI2=init$bi2,NU=init$nu,A3=init$a3,
                         B3=init$b3,ES=init$Es,ES2=init$Es2,
                         GAMSTAR=init$GamStar,GAM=init$Gam,MaxIt=controls$MaxIt,
                         epsilon=controls$epsilon)
        clust$a2Star=init$a2Star
        clust$b2Star=init$b2Star
        clust$b3Star=init$b3Star
        clust$cc=init$cc
        clust$GamStar=init$GamStar
      }else{
        clust=.VBCluster(DIMS=c(n,p,TT,M),YY=Y,CC=hyperparms$cc,A2STAR=hyperparms$a2Star,B2STAR=hyperparms$b2Star,
                         B3STAR=hyperparms$b3Star,EOM=init$EOm,MU=init$mu,SIG=init$Sig,BI0G=init$Bi0g,
                         BITBAR=init$Bitbar,BITHK=init$Bithk,ER=init$Er,
                         ER2=init$Er2,AI2=init$ai2,BI2=init$bi2,NU=init$nu,A3=init$a3,
                         B3=init$b3,ES=init$Es,ES2=init$Es2,
                         GAMSTAR=hyperparms$GamStar,GAM=init$Gam,MaxIt=controls$MaxIt,epsilon=controls$epsilon)
        clust$a2Star=hyperparms$a2Star
        clust$b2Star=hyperparms$b2Star
        clust$b3Star=hyperparms$b3Star
        clust$cc=hyperparms$cc
        clust$GamStar=hyperparms$GamStar
      }
      
      
      clust$Z=matrix(NA,n,TT)
      for(tt in 1:TT){
        clust$Z[,tt] <- unlist(apply(as.matrix(clust$Bitbar[tt,,]),2,which.max))
      }
      clust$method = "VB"
      clust$Y = Y
      clust <- clust[c("method","Y","mu","Sig","ai1","bi1","Er","Er2",
                       "ai2","bi2","a3","b3","ai4","bi4","Es","Es2",
                       "nu","kappa","Z","Bi0g","Bithk","Bitbar","Gam")]
      class(clust) <- "dnc"
      clust$pm = dncMAP(clust)
      
      return(clust)
    }
  }else{###End of if(VB)
    n=dim(Y)[1]
    TT=dim(Y)[3]
    Om = array(0,dim(Y))
    X = array(0,c(p,TT,n,controls$nDraws))
    Z = array(0,c(n,TT,controls$nDraws))
    rr = matrix(0,n,controls$nDraws)
    tau = matrix(0,n,controls$nDraws)
    uu = array(0,c(M,p,controls$nDraws))
    alpha = numeric(controls$nDraws)
    ss = matrix(0,n,controls$nDraws)
    bb = array(0,c(M+1,M,controls$nDraws))
    PostVal=numeric(controls$nDraws+controls$burnin)
    
    if(M==0){
      
      clust = .clustProjectionInitializeNoClust(Y=Y,p=p,hyperparms=hyperparms,
                                                MaxIt=controls$MaxIt,epsilon=controls$epsilon)
      clust$method = "VB"
      clust$Y = Y
      clust <- clust[c("method","Y","mu","Sig","a0","b0","ai2","bi2",
                       "a3","b3","ai4","bi4","Es","Es2")]
      class(clust) <- "dnc"
      clust$pm = dncMAP(clust)
      return(clust)
      
    }else{
      if(is.null(init)){
        init= .clustProjectionInitialize(Y=Y,M=M,p=p,hyperparms=hyperparms,
                                         MaxIt=controls$MaxItStg2,
                                         epsilon=controls$epsilonStg2)
        cc=init$cc
        a2Star=init$a2Star
        b2Star=init$b2Star
        b3Star=init$b3Star
        GamStar=init$GamStar
        for(tt in 1:TT){
          Z[,tt,1] <- unlist(apply(as.matrix(init$Bitbar[tt,,]),2,which.max))
        }
        bb[,,1]=t(scale(t(init$Gam),center=FALSE,scale=rowSums(init$Gam)))
        
      }else{
        cc=hyperparms$cc
        a2Star=hyperparms$a2Star
        b2Star=hyperparms$b2Star
        b3Star=hyperparms$b3Star
        GamStar=hyperparms$GamStar    
        
        Z[,,1] = init$Z
        bb[,,1] = init$beta
        
      }
      
      Om = init$EOm
      X[,,,1] = init$mu
      rr[,1] = init$Er
      tau[,1] = init$ai2*init$bi2
      uu[,,1] = init$nu
      alpha[1] = init$a3
      ss[,1] = init$Es
      
      PostVal[1] = cposteriorNoOmega(c(n,p,TT,M),Y,X[,,,1],
                                       rr[,1],tau[,1],ss[,1],alpha[1],
                                       bb[,,1],Z[,,1]-1,uu[,,1],cc,
                                       GamStar,a2Star,b2Star,b3Star)
      
      if(!is.null(Missing)) Missing = cbind(Missing,rep(0,nrow(Missing)))
      pb=tkProgressBar(min=1,max=controls$burnin,title="Burn-in period")
      for(it in 1:controls$burnin){
        
        ###Missing data
        if(!is.null(Missing)){
          temp = drawMiss(c(n,p,TT),Y,ss[,1],Missing,X[,,,1],alpha[1])
          Y = temp[[1]]
          Missing = temp[[2]]
        }
        
        ###Omega: Auxiliary variables
        for(tt in 1:TT){
          temp = c(t(X[,tt,,1])%*%t(apply(X[,tt,,1],1,function(x) x*ss[,1])))
          Om[,,tt] = matrix(rpg.devroye(num=n^2,1,z=alpha[1]+temp),n,n)
        }
        
        ###X,Z, alpha and rr
        RN = array(rnorm(p*n*TT),c(p,TT,n))
        RNa=rnorm(1)
        draw1 = drawXZar(c(n,p,TT,M),X[,,,1],Y,Om,ss[,1],
                          alpha[1],tau[,1],rr[,1],
                          t(uu[,,1]),Z[,,1]-1,RN,bb[,,1],
                          RNa,b3Star,cc)
        X[,,,1] = draw1[[1]]
        Z[,,1] = draw1[[2]]+1
        alpha[1] = draw1[[3]]
        rr[,1] = .rtrunc(n,"norm",0,Inf,mean=drop(draw1[[4]]),sd=1/sqrt(tau[,1]*TT))
        
        ###tau
        draw2 = drop(drawTau(c(n,p,TT,M),X[,,,1],rr[,1],t(uu[,,1]),Z[,,1]-1))
        draw2 =  1/b2Star + rr[,1]/cc + draw2 
        tau[,1] = rgamma(n,shape=a2Star+p*TT/2+1,rate=draw2)
        
        ###uu
        draw3=t(drawU(c(n,p,TT,M),tau[,1],rr[,1],X[,,,1],Z[,,1]-1))
        for(g in 1:M){
          uu[g,,1] = movMF::rmovMF(n=1,theta=draw3[g,])
        }
        
        ###beta
        JhgtCard=matrix(0,M,M)
        for(tt in 1:TT){
          if(tt>1){
            temp = table(Z[,tt-1,1],Z[,tt,1])
            ind = which(temp!=0,arr.ind=TRUE)
            rnames=as.integer(rownames(temp))
            cnames=as.integer(colnames(temp))
            for(k in 1:dim(ind)[1]){
              JhgtCard[rnames[ind[k,1]],cnames[ind[k,2]]]=
                JhgtCard[rnames[ind[k,1]],cnames[ind[k,2]]] + 
                temp[ind[k,1],ind[k,2]]
            }
          }
        }
        temp=table(Z[,1,1])
        Ig1=numeric(M)
        for(k in 1:dim(temp)){
          Ig1[as.integer(names(temp)[k])]=temp[k]
        }
        bb[1,,1] = MCMCpack::rdirichlet(1,GamStar[1,]+Ig1)
        for(h in 1:M){
          bb[1+h,,1] = MCMCpack::rdirichlet(1,GamStar[1+h,]+JhgtCard[h,])
        }
        
        ###Draw ss
        draw4 = drawS(c(n,p,TT),Y,Om,X[,,,1],alpha[1])
        ss[,1]= .rtrunc(n,"norm",0,Inf,mean=drop(draw4[[1]]),sd=sqrt(drop(draw4[[2]])))
        ss[,1]=ss[,1]/mean(ss[,1])
        
        PostVal[1+it] = cposteriorNoOmega(c(n,p,TT,M),Y,X[,,,1],rr[,1],
                                          tau[,1],ss[,1],alpha[1],bb[,,1],
                                          Z[,,1]-1,uu[,,1],cc,GamStar,
                                          a2Star,b2Star,b3Star)
        
        setTkProgressBar(pb,it)
      }
      close(pb)
      
      if(!is.null(Missing)) Missing[,4] = 0L
      pb=tkProgressBar(min=1,max=controls$nDraws,title="Gibbs sampler")
      for(it in 2:controls$nDraws){
        
        ###Missing data
        if(!is.null(Missing)){
          temp = drawMiss(c(n,p,TT),Y,ss[,it-1],Missing,X[,,,it-1],alpha[it-1])
          Y = temp[[1]]
          Missing = temp[[2]]
        }
        
        ###Omega: Auxiliary variables
        for(tt in 1:TT){
          temp = c(t(X[,tt,,it-1])%*%t(apply(X[,tt,,it-1],1,function(x) x*ss[,it-1])))
          Om[,,tt] = matrix(rpg.devroye(n^2,1,alpha[it-1]+temp),n,n)
        }
        
        ###X,Z, alpha and rr
        RN = array(rnorm(p*n*TT),c(p,TT,n))
        RNa=rnorm(1)
        draw1 = drawXZar(c(n,p,TT,M),X[,,,it-1],Y,Om,ss[,it-1],
                          alpha[it-1],tau[,it-1],rr[,it-1],
                          t(uu[,,it-1]),Z[,,it-1]-1,RN,bb[,,it-1],
                          RNa,b3Star,cc)
        X[,,,it] = draw1[[1]]
        Z[,,it] = draw1[[2]]+1
        alpha[it] = draw1[[3]]
        rr[,it] = .rtrunc(n,"norm",0,Inf,mean=drop(draw1[[4]]),sd=1/sqrt(tau[,it-1]*TT))
        
        ###tau
        draw2 = drop(drawTau(c(n,p,TT,M),X[,,,it],rr[,it],t(uu[,,it-1]),Z[,,it]-1))
        draw2 =  1/b2Star + rr[,it]/cc + draw2 
        tau[,it] = rgamma(n,shape=a2Star+p*TT/2+1,rate=draw2)
        
        ###uu
        draw3=t(drawU(c(n,p,TT,M),tau[,it],rr[,it],X[,,,it],Z[,,it]-1))
        for(g in 1:M){
          uu[g,,it] = movMF::rmovMF(n=1,theta=draw3[g,])
        }
        
        ###beta
        JhgtCard=matrix(0,M,M)
        for(tt in 1:TT){
          if(tt>1){
            temp = table(Z[,tt-1,it],Z[,tt,it])
            ind = which(temp!=0,arr.ind=TRUE)
            rnames=as.integer(rownames(temp))
            cnames=as.integer(colnames(temp))
            for(k in 1:dim(ind)[1]){
              JhgtCard[rnames[ind[k,1]],cnames[ind[k,2]]]=
                JhgtCard[rnames[ind[k,1]],cnames[ind[k,2]]] + 
                temp[ind[k,1],ind[k,2]]
            }
          }
        }
        temp=table(Z[,1,1])
        Ig1=numeric(M)
        for(k in 1:dim(temp)){
          Ig1[as.integer(names(temp)[k])]=temp[k]
        }
        bb[1,,it] = MCMCpack::rdirichlet(1,GamStar[1,]+Ig1)
        for(h in 1:M){
          bb[1+h,,it] = MCMCpack::rdirichlet(1,GamStar[1+h,]+JhgtCard[h,])
        }
        
        ###Draw ss
        draw4 = drawS(c(n,p,TT),Y,Om,X[,,,it],alpha[it])
        ss[,it]= .rtrunc(n,"norm",0,Inf,mean=drop(draw4[[1]]),sd=sqrt(drop(draw4[[2]])))
        ss[,it]=ss[,it]/mean(ss[,it])
        
        ###Evaluate posterior to find mode
        PostVal[controls$burnin+it] = 
          cposteriorNoOmega(c(n,p,TT,M),Y,X[,,,it],rr[,it],tau[,it],ss[,it],
                            alpha[it],bb[,,it],Z[,,it]-1,
                            uu[,,it],cc,GamStar,a2Star,b2Star,b3Star)
        
        
        setTkProgressBar(pb,it)
      }
      close(pb)
      
      if(is.null(Missing)){
        clust = list(X=X,Z=Z,r=rr,tau=tau,u=uu,alpha=alpha,s=ss,
                     beta=bb,posterior=PostVal)
        clust$method = "Gibbs"
        clust$Y = Y
        clust <- clust[c("method","Y","X","r","tau","alpha","s","u","Z","beta","posterior")]
        class(clust) <- "dnc"
        clust$pm = dncMAP(clust)
        return(clust)
        
      }else{
        Missing[,4] = Missing[,4]/controls$nDraws
        colnames(Missing) = c("Row","Col","Time","Prob(Y_{ij}=1)")
        clust = list(X=X,Z=Z,r=rr,tau=tau,u=uu,alpha=alpha,s=ss,
                     beta=bb,posterior=PostVal,Missing=Missing)
        clust$method = "Gibbs"
        clust$Y = Y
        clust <- clust[c("method","Y","X","r","tau","alpha","s","u","Z","beta","posterior","Missing")]
        class(clust) <- "dnc"
        clust$pm = dncMAP(clust)
        return(clust)
        
      }
    }
  }
}


# Stable trunc norm functions ---------------------------------------------


.qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
}

.rtrunc <- function(n, spec, a = -Inf, b = Inf, ...)
{
  x <- u <- runif(n, min = 0, max = 1)
  x <- .qtrunc(u, spec, a = a, b = b,...)
  return(x)
}


# Simulate data -----------------------------------------------------------

simulate.dnc <- function(object,nsim=1,seed=NULL, ...){
  if(!is.null(seed))set.seed(seed)
  n = dim(object$Y)[1]
  TT = dim(object$Y)[3]
  pp = dim(object$pm$X)[1]
  Ysims = array(0L,c(n,n,TT,nsim))
  if(nsim>1) pb = tkProgressBar(min=1,max=nsim,title="Simulating data")
  for(ss in 1:nsim){
    Ysims[,,,ss] =
      csimulateData(object$pm$X,c(n,pp,TT),object$pm$alpha,object$pm$s)
    if(nsim>1) setTkProgressBar(pb,ss)
  }
  if(nsim>1) close(pb)
  return(Ysims)
}


# BIC for model selection -------------------------------------------------

BIC.dnc <- function(object,...){
  
  if(is.null(object$pm$u))stop("Current functionality for BIC is for G>1")
  
  n = dim(object$Y)[1]
  pp = dim(object$pm$X)[1]
  TT = dim(object$Y)[3]
  M = dim(object$pm$u)[1]
  
  ret = -2*cLogLikelihood(object$Y,object$pm$X,DIMS=c(n,pp,TT),
                          object$pm$alpha,object$pm$s)-
    2*cllikProj(object$pm$X,c(n,pp,M,TT),object$pm$beta[1,],
                object$pm$beta[-1,],object$pm$u,object$pm$r,
                object$pm$tau) +
    (n+1)*log(sum(object$Y)) + (2*n+(pp-1)*M+(M+1)*(M-1))*log(n*TT)
  
  return(ret)
}



# Posterior Mode ----------------------------------------------------------

dncMAP =function(dncObj){
  
  if(class(dncObj)!="dnc")stop("Not a 'dnc' object")
  
  if(dncObj$method=="VB"){
    if(is.null(dncObj$nu)){
      p=dim(dncObj$mu)[1]
      TT=dim(dncObj$mu)[2]
      n=dim(dncObj$mu)[3]
      tauhat <- (dncObj$ai2-1)*dncObj$bi2
      s2hat = dncObj$b0/(dncObj$a0+1)
      return(list(alpha=dncObj$a3,
                  X=dncObj$mu,
                  s=sapply(dncObj$ai4,function(x)max(x,0)),
                  tau=tauhat,
                  sigma2= s2hat))
      
    }else{
      M = dim(dncObj$nu)[1]
      p=dim(dncObj$mu)[1]
      TT=dim(dncObj$mu)[2]
      n=dim(dncObj$mu)[3]
      tauhat <- (dncObj$ai2-1)*dncObj$bi2
      betahat <- matrix(0,M+1,M)
      zhat <- matrix(0,n,TT)
      for(hh in 1:(M+1)){
        betahat[hh,] <- (dncObj$Gam[hh,]-1)/(sum(dncObj$Gam[hh,])-M)
      }
      for(i in 1:n){
        Z2 <- numeric(TT)
        Z1 <- apply(dncObj$Bitbar[,,i],1,which.max)
        for(count in 1:1000){
          Z1 -> Z2
          Z1[1] <- which.max(log(dncObj$Bi0g[i,])+log(dncObj$Bithk[M+1:M,Z1[2],i]))
          for(tt in 2:(TT-1)){
            Z1[tt] <- which.max(log(dncObj$Bithk[(tt-1)*M+Z1[tt-1],,i])+log(dncObj$Bithk[tt*M+1:M,Z1[tt+1],i]))
          }
          Z1[TT] <- which.max(log(dncObj$Bithk[(TT-1)*M+Z1[TT-1],,i]))
          if(sum(Z1 != Z2) == 0) break
        }
        zhat[i,] <- Z1
      }
      
      
      return(list(alpha=dncObj$a3,
                  X=dncObj$mu,
                  s=sapply(dncObj$ai4,function(x)max(x,0)),
                  tau=tauhat,
                  r=sapply(dncObj$ai1,function(x)max(x,0)),
                  u=dncObj$nu,
                  Z=zhat,
                  beta=betahat))
    }
  }
  if(dncObj$method=="Gibbs"){
    burnin = length(dncObj$posterior) - dim(dncObj$X)[4]
    ind= which.max(dncObj$posterior[-c(1:burnin)])
    return(list(alpha=dncObj$alpha[ind],
                X=dncObj$X[,,,ind],
                s=dncObj$s[,ind],
                tau=dncObj$tau[,ind],
                r=dncObj$r[,ind],
                u=dncObj$u[,,ind],
                Z=dncObj$Z[,,ind],
                beta=dncObj$beta[,,ind]))
    
  }
}



# print dnc ---------------------------------------------------------------

print.dnc = function(x,printDens=FALSE,...){
  cat("Dynamic Network Clustering Object\n")
  cat("Number of actors = ",dim(x$Y)[1],"\n")
  cat("Number of time points = ",dim(x$Y)[3],"\n")
  cat("Estimation method: ",x$method,"\n")
  if(is.null(x$pm$r)){
    M = dim(x$pm$X)[1]
    if(M==2){
      cat("Latent space = circle \n")
    }
    if(M==3){
      cat("Latent space = sphere \n")
    }
    if(M>3){
      cat("Latent space = ",dim(x$pm$X)[1],"- dimensional hypersphere \n")
    }
    cat("No clustering performed \n")
  }else{
    M = dim(x$pm$X)[1]
    if(M==2){
      cat("Latent space = circle \n")
    }
    if(M==3){
      cat("Latent space = sphere \n")
    }
    if(M>3){
      cat("Latent space = ",dim(x$pm$X)[1],"- dimensional hypersphere \n")
    }    
    cat("Number of clusters: ",dim(x$pm$beta)[2],"\n")
  }
  if(printDens){
    temp = numeric(dim(x$Y)[3])
    for(tt in 1:dim(x$Y)[3]){
      temp[tt] <- sum(x$Y[,,tt])/dim(x$Y)[1]/(dim(x$Y)[1]-1)
    }
    temp = cbind(1:length(temp),temp);
    colnames(temp)=c("Time","Network Density")
    print(as.data.frame(temp),row.names=FALSE)
  }
}



# plot.dnc ----------------------------------------------------------------

plot.dnc = function(x,aggregated=TRUE,plotRGL=TRUE,
                    Lines=TRUE,colByComm=TRUE,
                    INDEX=1:min(dim(x$pm$X)[1],3),...){
  if(!length(INDEX)%in%c(2,3))stop("Choose either 2 or 3 latent dimensions.")
  if(max(INDEX)>dim(x$pm$X)[1])stop("Maximum index must not be greater than p.")
  
  if(is.null(x$pm$u)){Lines=FALSE;colByComm=FALSE}
  
  if(length(INDEX)==3){
    M <- mesh(seq(0, 2*pi, length.out = 20),
              seq(0, pi, length.out = 20))
    u <- M$x ; v <- M$y
    xx <- cos(u)*sin(v)
    yy <- sin(u)*sin(v)
    zz <- cos(v)
    VV <- x$pm$X[INDEX,,]
    TT <- dim(x$Y)[3]
    for(tt in 1:TT){
      VV[,tt,] <- t(t(VV[,tt,])/apply(t(VV[,tt,]),1,function(x) sqrt(x%*%x)))
    }
    if(aggregated){
      surf3D(x = xx, y = yy, z = zz,
             colvar = NULL, lighting = TRUE, #plot = FALSE,
             facets = NA, col = gray(0.8), lwd = 0.75,...)
      if(colByComm){
        Cols = rainbow(dim(x$pm$u)[1],s=0.5)
        if(Lines){
          for(i in 1:dim(x$pm$u)[1]){
            lines3D(c(0,x$pm$u[i,INDEX[1]]),c(0,x$pm$u[i,INDEX[2]]),c(0,x$pm$u[i,INDEX[3]]),
                    add=TRUE,col=Cols[i],colorkey=FALSE,...)
          }
        }
        for(tt in 1:TT){
          for(mm in 1:dim(x$pm$u)[1]){
            ind = which(x$pm$Z[,tt]==mm)
            if(length(ind)>0){
              points3D(VV[1,tt,ind],VV[2,tt,ind],VV[3,tt,ind],add=TRUE,
                       col=Cols[mm],colkey=FALSE,...)
            }
          }
        }
      }else{
        if(Lines){
          Cols = rainbow(dim(x$pm$u)[1],s=0.5)
          for(i in 1:dim(x$pm$u)[1]){
            lines3D(c(0,x$pm$u[i,INDEX[1]]),c(0,x$pm$u[i,INDEX[2]]),c(0,x$pm$u[i,INDEX[3]]),
                    add=TRUE,col=Cols[i],colorkey=FALSE,...)
          }
        }
        for(tt in 1:TT){
          points3D(VV[1,tt,],VV[2,tt,],VV[3,tt,],add=TRUE,col="black",colkey=FALSE,...)
        }
      }
      if(plotRGL) plotrgl()
    }else{
      for(tt in 1:TT){
        surf3D(x = xx, y = yy, z = zz,
               colvar = NULL, lighting = TRUE, #plot = FALSE,
               facets = NA, col = gray(0.8), lwd = 0.75,...)
        if(colByComm){
          Cols = rainbow(dim(x$pm$u)[1],s=0.5)
          for(mm in 1:dim(x$pm$u)[1]){
            ind = which(x$pm$Z[,tt]==mm)
            if(length(ind)>0){
              points3D(VV[1,tt,ind],VV[2,tt,ind],VV[3,tt,ind],add=TRUE,
                       col=Cols[mm],colkey=FALSE,...)
            }
          }
          if(Lines){
            for(i in 1:dim(x$pm$u)[1]){
              lines3D(c(0,x$pm$u[i,INDEX[1]]),c(0,x$pm$u[i,INDEX[2]]),c(0,x$pm$u[i,INDEX[3]]),
                      add=TRUE,col=Cols[i],colorkey=FALSE,...)
            }
          }
        }else{
          Cols = rainbow(dim(x$pm$u)[1],s=0.5)
          Cols = Cols[unique(c(x$pm$Z))]
          points3D(VV[1,tt,],VV[1,tt,],VV[3,tt,],add=TRUE,col="black",colkey=FALSE,...)
          if(Lines){
            for(i in 1:dim(x$pm$u)[1]){
              lines3D(c(0,x$pm$u[i,INDEX[1]]),c(0,x$pm$u[i,INDEX[2]]),c(0,x$pm$u[i,INDEX[3]]),
                      add=TRUE,col=Cols[i],colorkey=FALSE,...)
            }
          }
        }
        if(plotRGL) plotrgl()
        print("Hit 'Enter' to continue")
        readline()
      }
    }
  }
  if(length(INDEX)==2){
    xyCirc = seq(0,2*pi,length.out = 500)
    xyCirc = cbind(cos(xyCirc),sin(xyCirc))
    
    VV <- x$pm$X[INDEX,,]
    TT <- dim(x$Y)[3]
    for(tt in 1:TT){
      VV[,tt,] <- t(t(VV[,tt,])/apply(t(VV[,tt,]),1,function(x) sqrt(x%*%x)))
    }
    if(Lines){
      UU = x$pm$u[,INDEX]
      UU <- t(apply(UU,1,function(x)x/sqrt(x%*%x)))
    }
    if(aggregated){
      plot(xyCirc,type="l",xaxt="n",yaxt="n",xlab="",ylab="",
           bty="n",col=gray(0.8))
      if(colByComm){
        Cols = rainbow(dim(x$pm$u)[1],s=0.5)
        if(Lines){
          for(i in unique(c(x$pm$Z))){
            lines(c(0,UU[i,1]),c(0,UU[i,2]),col=Cols[i])
          }
        }
        for(tt in 1:TT){
          points(VV[1,tt,],VV[2,tt,],col=Cols[x$pm$Z[,tt]],...)
        }
      }else{
        if(Lines){
          Cols = rainbow(dim(x$pm$u)[1],s=0.5)
          for(i in unique(c(x$pm$Z))){
            lines(c(0,UU[i,1]),c(0,UU[i,2]),col=Cols[i])
          }
        }
        for(tt in 1:TT){
          points(VV[1,tt,],VV[2,tt,],...)
        }
      }
    }else{
      for(tt in 1:TT){
        plot(xyCirc,type="l",xaxt="n",yaxt="n",xlab="",ylab="",
             bty="n",col=gray(0.8))
        if(colByComm){
          Cols = rainbow(dim(x$pm$u)[1],s=0.5)
          if(Lines){
            for(i in unique(c(x$pm$Z))){
              lines(c(0,UU[i,1]),c(0,UU[i,2]),
                    col=Cols[i])
            }
          }
          points(VV[1,tt,],VV[2,tt,],col=Cols[x$pm$Z[,tt]],...)
        }else{
          if(Lines){
            Cols = rainbow(dim(x$pm$u)[1],s=0.5)
            for(i in unique(c(x$pm$Z))){
              lines(c(0,UU[i,1]),c(0,UU[i,2]),
                    col=Cols[i])
            }
          }
          points(VV[1,tt,],VV[2,tt,],...)
        }
        print("Hit 'Enter' to continue")
        readline()
      }
    }
  }
}
