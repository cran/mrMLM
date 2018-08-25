FASTmrEMMA<-function(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,Likelihood,CLO){


if(is.null(CLO)==FALSE){
  OUT<-list(gen=gen,phe=phe)
  return(OUT)
}else{
if(Likelihood=="REML"){
  flagREMLE<-1
}else if(Likelihood=="ML"){
  flagREMLE<-0 
}

inputform<-Genformat

if(is.null(kk)){

  emma.kinship <- function(snps, method="additive", use="all") {
    n0 <- sum(snps==0,na.rm=TRUE)
    nh <- sum(snps==0.5,na.rm=TRUE)
    n1 <- sum(snps==1,na.rm=TRUE)
    nNA <- sum(is.na(snps))
    stopifnot(n0+nh+n1+nNA == length(snps))
    if ( method == "dominant" ) {
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
      snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
    }
    else if ( method == "recessive" ) {
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
      snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
    }
    else if ( ( method == "additive" ) && ( nh > 0 ) ) {
      dsnps <- snps
      rsnps <- snps
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
      los<-intersect(which(!is.na(snps)),which(snps==0.5))
      dsnps[los] <- flags[los]
      rm(flags)
      gc()
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
      rsnps[los] <- flags[los]
      rm(flags,snps)
      gc()
      snps <- rbind(dsnps,rsnps)
      rm(dsnps,rsnps)
      gc()
    }
    if ( use == "all" ) {
      mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
      losna<-which(is.na(snps))
      snps[losna] <- mafs[losna]
      rm(mafs)
      gc()
    }
    else if ( use == "complete.obs" ) {
      snps <- snps[rowSums(is.na(snps))==0,]
    }
    n <- ncol(snps)
    K<-(t(snps)%*%snps+t(1-snps)%*%(1-snps))/nrow(snps)
    diag(K) <- 1
    return(K)
  }




if(is.null(gen)==TRUE)
{
  warning("Please input correct genotype dataset !")
}else{
  snp8<-deepcopy(gen,3:ncol(gen))
  kk<-emma.kinship(snp8[,]) 
  rm(snp8)
  gc()
 }
}

if(is.null(psmatrix)){
  flagps<-1
}else{
  
  flagps<-0
}


if(is.null(svpal)==TRUE||is.null(svmlod)==TRUE){
  warning("Please set parameter!")
}

if((svpal<0)||(svpal>1))
{
  warning("Please input critical P-value between 0 and 1!")
}

if(svmlod<0)
{
  warning("Please input critical LOD score: >0!")
}
if(exists("gen")==FALSE)
{
  warning("Please input correct genotype dataset !")
}
if(exists("phe")==FALSE)
{
  warning("Please input correct phenotype dataset !")
}
if(exists("kk")==FALSE)
{
  warning("Please input correct kinship (K) dataset !")
}
if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(ncol(gen)!=(nrow(phe)+2)))
{
  warning("Sample sizes between genotypic and phenotypic datasets do not equal !")
}

if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2)))&&(svpal>=0)&&(svpal<=1)&&(svmlod>=0))
{

parmsShow=NULL
wan=NULL
parms=NULL
ress1=NULL
mannewp=NULL

multinormal<-function(y,mean,sigma)
{
  pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
  return (pdf_value)
}

ebayes_EM<-function(x,z,y)
{
  n<-nrow(z);k<-ncol(z)
  
  if(abs(min(eigen(crossprod(x,x))$values))<1e-6){
    b<-solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y)
  }else{
    b<-solve(crossprod(x,x))%*%(crossprod(x,y))
  }
  
  v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
  u<-matrix(rep(0,k),k,1)
  v<-matrix(rep(0,k),k,1)
  s<-matrix(rep(0,k),k,1)
  for(i in 1:k)
  {
    zz<-z[,i]
    s[i]<-((crossprod(zz,zz)+1e-100)^(-1))*v0
    u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
    v[i]<-u[i]^2+s[i]
  }
  
  vv<-matrix(rep(0,n*n),n,n);
  for(i in 1:k)
  {
    zz<-z[,i]
    vv=vv+tcrossprod(zz,zz)*v[i]
  }
  vv<-vv+diag(n)*v0
  
  iter<-0;err<-1000;iter_max<-500;err_max<-1e-8
  tau<-0;omega<-0
  while((iter<iter_max)&&(err>err_max))
  {
    iter<-iter+1
    v01<-v0
    v1<-v
    b1<-b
    vi<-solve(vv)
    xtv<-crossprod(x,vi)
    
    if(ncol(x)==1)
    {
      b<-((xtv%*%x)^(-1))*(xtv%*%y)
    }else{
      if(abs(min(eigen(xtv%*%x)$values))<1e-6){
        b<-solve((xtv%*%x)+diag(ncol(x))*1e-8)%*%(xtv%*%y)
      }else{
        b<-solve(xtv%*%x)%*%(xtv%*%y)
      }
    }
    r<-y-x%*%b
    ss<-matrix(rep(0,n),n,1)
    for(i in 1:k)
    {
      zz<-z[,i]
      zztvi<-crossprod(zz,vi)
      u[i]<-v[i]*zztvi%*%r
      s[i]<-v[i]*(1-zztvi%*%zz*v[i])
      v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
      ss<-ss+zz*u[i]
    }
    v0<-as.numeric(crossprod(r,(r-ss))/n)
    
    vv<-matrix(rep(0,n*n),n,n)
    for(i in 1:k)
    {
      zz<-z[,i]
      vv<-vv+tcrossprod(zz,zz)*v[i]
    }
    vv<-vv+diag(n)*v0
    
    err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
    beta<-t(b)
    sigma2<-v0
  }
  
  wang<-matrix(rep(0,k),k,1)
  for (i in 1:k){
    stderr<-sqrt(s[i]+1e-20)
    t<-abs(u[i])/stderr
    f<-t*t
    p<-1-pchisq(f,1)
    wang[i]<-p
  }
  
  return(list(u=u,sigma2=sigma2,wang=wang))
}

likelihood<-function(xxn,xxx,yn,bbo)
{
  nq<-ncol(xxx)
  ns<-nrow(yn)
  at1<-0
  
  if(is.null(bbo)==TRUE){
    ww1<-1:ncol(xxx)
    ww1<-as.matrix(ww1)
  }else{
    ww1<-as.matrix(which(abs(bbo)>1e-5))
  }
  at1<-dim(ww1)[1]
  lod<-matrix(rep(0,nq),nq,1)
  if(at1>0.5)
    ad<-cbind(xxn,xxx[,ww1])
  else
    ad<-xxn
  if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6)
    bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
  else
    bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
  vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns);
  ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
  sub<-1:ncol(ad);
  if(at1>0.5)
  {
    for(i in 1:at1)
    {
      ij<-which(sub!=sub[i+ncol(xxn)])
      ad1<-ad[,ij]
      if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6)
        bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
      else
        bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn) 
      vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
      ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
      lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
    }
  }
  return (lod)
}
emma.eigen.L <- function(Z,K,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.L.wo.Z(K))
  }
  else {
    return(emma.eigen.L.w.Z(Z,K,complete))
  }
}
#likelihood
emma.eigen.L.wo.Z <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}
#likelihood
emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
  return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
}


emma.eigen.R <- function(Z,K,X,complete=TRUE) {
  if ( ncol(X) == 0 ) {
    return(emma.eigen.L(Z,K))
  }
  else if ( is.null(Z) ) {
    return(emma.eigen.R.wo.Z(K,X))
  }
  else {
    return(emma.eigen.R.w.Z(Z,K,X,complete))
  }
}
#restricted likelihood
emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)
  
  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)    
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1:(t-q)],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                           complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
}


emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(delta*lambda+1))))-sum(log(delta*xi+1))) )  
}

emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*xi.1+1)) ))
  
}

emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- delta*lambda+1
  return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(xi/(delta*xi+1))) )
}

emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- delta*lambda+1
  return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(xi.1/(delta*xi.1+1))) )
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(delta*lambda+1))))-sum(log(delta*lambda+1))) )
}

emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- delta*lambda+1
  return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(lambda/ldelta)) )
}

emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
  t <- t1
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- delta*lambda+1
  return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
}


emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                     esp=1e-10, eig.L = NULL, eig.R = NULL)
{
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(ML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z(K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1<-matrix(eig.R$values,n-q,m)    
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
    Xis.1<-matrix(eig.L$values,n,m)
    Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1
    Etasq <- matrix(etas*etas,n-q,m)
    dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
      }
    }
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z(Z,K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1<-matrix(eig.R$values,t-q,m)
    Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
    Xis.1<-matrix(eig.L$values,t,m)
    Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
      }
    }
  }
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  optLL=replaceNaN(optLL)
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
  }
  maxvg <- maxve*maxdelta
  
  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
}

emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1<-matrix(eig.R$values,n-q,m)
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
    Etasq <- matrix(etas*etas,n-q,m)
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
      }
    }
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1 <- matrix(eig.R$values,t-q,m) 
    Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
      }
    }
  }  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  optLL=replaceNaN(optLL)
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
  }
  maxvg <- maxve*maxdelta
  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
}

emma.maineffects.B<-function(Z=NULL,K,deltahat.g,complete=TRUE){
  if( is.null(Z) ){
    return(emma.maineffects.B.Zo(K,deltahat.g))
  }
  else{
    return(emma.maineffects.B.Z(Z,K,deltahat.g,complete))
  }
}


emma.maineffects.B.Zo <-function(K,deltahat.g){
  t <- nrow(K)
  stopifnot(ncol(K) == t)
  B<-deltahat.g*K+diag(1,t)
  eig<-eigen(B,symmetric=TRUE)
  qr.B<-qr(B)
  q<-qr.B$rank
  stopifnot(!is.complex(eig$values))
  A<-diag(1/sqrt(eig$values[1:q]))
  Q<-eig$vectors[,1:q]
  C<-Q%*%A%*%t(Q)
  return(list(mC=C,Q=Q,A=A))
}

emma.maineffects.B.Z <- function(Z,K,deltahat.g,complete=TRUE){
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)  
  B <- deltahat.g*Z%*%K%*%t(Z)+diag(1,n)
  eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
  qr.B<-qr(B)
  q<-qr.B$rank
  stopifnot(!is.complex(eig$values))
  A<-diag(1/sqrt(eig$values[1:q]))
  Q<-eig$vectors[,1:q]
  C<-Q%*%A%*%t(Q)
  return(list(mC=C,Q=Q,A=A,complete=TRUE))
}
emma.MLE0.c <- function(Y_c,W_c){
  n <- length(Y_c)
  stopifnot(nrow(W_c)==n)
  M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
  etas<-crossprod(M_c,Y_c)
  LL <- 0.5*n*(log(n/(2*pi))-1-log(sum(etas*etas)))
  return(list(ML=LL))
}

emma.REMLE0.c <- function(Y_c,W_c){
  n <- length(Y_c)
  stopifnot(nrow(W_c)==n)
  M_c <-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
  eig <-eigen(M_c)
  t <-qr(W_c)$rank
  v <-n-t
  U_R <-eig$vector[,1:v]
  etas<-crossprod(U_R,Y_c)
  LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
  return(list(REML=LL))
}

replaceNaN<-  function(LL) {
  index=(LL=="NaN")
  if(length(index)>0) theMin=min(LL[!index])
  if(length(index)<1) theMin="NaN"
  LL[index]=theMin
  return(LL)    
}

emma.eigen.L.c <- function(Z,K,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.L.wo.Z.c(K))
  }
  else {
    return(emma.eigen.L.w.Z.c(Z,K,complete))
  }
}
#likelihood
emma.eigen.L.wo.Z.c <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}#don't run

#likelihood
emma.eigen.L.w.Z.c <- function(Z,K,complete=TRUE) {
  if ( complete == FALSE ) {#don't run
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }#don't run
  eig <- eigen(K%*%crossprod(Z,Z),symmetric=TRUE)
  #symmetric=FALSE change to TRUE,EISPACK=TRUE  20140827
  return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
}

#restricted likelihood
emma.eigen.R.c <- function(Z,K,X,complete=TRUE) {
  if ( ncol(X) == 0 ) {
    return(emma.eigen.L.c(Z,K))
  }
  else if ( is.null(Z) ) {
    return(emma.eigen.R.wo.Z.c(K,X))
  }
  else {
    return(emma.eigen.R.w.Z.c(Z,K,X,complete))
  }
}
#restricted likelihood
emma.eigen.R.wo.Z.c <- function(K, X) {
  
  if(is.matrix(X)) {    
    n<-nrow(X)   
    q<-ncol(X) 
  } 
  else{   
    n<-length(X) 
    q<-1 
  }
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

#restricted likelihood
emma.eigen.R.w.Z.c <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  if(!is.matrix(Z)) n<-length(Z) 
  t <-1
  if(is.matrix(X)) {
    q<-ncol(X)
  }    
  else  {
    q<-1
  }
  
  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)    
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:t],qr.X)),
                           complete=TRUE)[,c(1:t,(t+q+1):n)]))
  
}

emma.delta.REML.LL.w.Z.c <- function(logdelta, lambda, etas.1, n, q, etas.2.sq ) {
  nq<-n-q
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
}

emma.delta.REML.dLL.w.Z.c <- function(logdelta, lambda, etas.1, n, q1, etas.2.sq ) {
  
  q<-q1
  nq<-n-q
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- delta*lambda+1
  return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
}

emma.MLE.c <- function(y, X, K=1, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL)
{
  if(is.matrix(y)){
    n<-nrow(y)
  }
  else{
    n <- length(y)
  }
  
  t <- 1
  if ( is.matrix(X) ){
    q <- ncol(X)
    stopifnot(nrow(X)==n)
  }
  else{
    q <- 1
    stopifnot(length(X)==n)
  }
  
  stopifnot(K == 1)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(ML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z.c(K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z.c(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,n-q,m)    
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
    Xis.1<-matrix(eig.L$values,n,m)
    Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1
    Etasq <- matrix(etas*etas,n-q,m)
    dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
      }
    }
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z.c(Z,K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z.c(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1<-etas[1:t]
    etas.2<-etas[(t+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1<-matrix(eig.R$values,t,m)
    Lambdas <- Lambdas.1 * matrix(delta,t,m,byrow=TRUE) + 1
    Xis.1<-matrix(eig.L$values,t,m)
    Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t,m)
    dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
      }
    }
  }
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  optLL=replaceNaN(optLL)
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
  }
  maxvg <- maxve*maxdelta
  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,U_R=eig.R$vectors,etas.1=etas.1,etas=etas,lambda=eig.R$values))
}

emma.REMLE.c <- function(y, X, K=1, Z=NULL, ngrids=100, llim=-10, ulim=10,
                         esp=1e-10, eig.L = NULL, eig.R = NULL)
{
  if(is.matrix(y)){
    n<-nrow(y)
  }
  else{
    n <- length(y)
  }
  
  t <- 1
  if ( is.matrix(X) ){
    q <- ncol(X)
  }
  else{
    q <- 1
  }
  
  stopifnot(K == 1)
  stopifnot(nrow(X) == n)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }
  
  
  if ( is.null(Z) ) {
    
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z.c(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1<-matrix(eig.R$values,n-q,m)
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
    Etasq <- matrix(etas*etas,n-q,m)
    
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
      }
    }
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z.c(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:t]
    etas.2 <- etas[(t+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1 <- matrix(eig.R$values,t,m) 
    Lambdas <- Lambdas.1 * matrix(delta,t,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t,m)
    
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(llim,eig.R$values,etas.1,n,q,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(ulim,eig.R$values,etas.1,n,q,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.w.Z.c, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, q1=q, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(r$root,eig.R$values, etas.1, n, q, etas.2.sq ))
      }
    }
    
  }  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  optLL=replaceNaN(optLL)
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
  }
  maxvg <- maxve*maxdelta
  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,U_R=eig.R$vectors,etas.1=etas.1,etas=etas,lambda=eig.R$values))
}



emma.ML.LRT.c.noalpha <- function(ys, xs, K=1, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
  stopifnot(K == 1)
  ys <- Z%*%ys
  xs <- Z%*%xs
  X0 <- Z%*%X0
  
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  n <- nrow(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)
  MLE0<-emma.MLE0.c(ys,X0)
  
  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  vgs <- vector(length=t)
  ves <- vector(length=t)
  deltas<-vector(length=t)
  bhats<-vector(length=t)
  var.bhats.ratio <- vector(length=t)
  d <- vector(length=t)
  
  stats <- vector(length=t)
  ps <- vector(length=t)
  
  for (i in 1:t){
    
    vids <- !is.na(xs[,i])
    xv <- xs[vids,i]
    yv <- ys[vids]
    x0v<-X0[vids,]
    MLE1 <- emma.MLE.c (yv, x0v, K=1, xv, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
    
    if(length(MLE1$vg)!=0){
      ML1s[i]<-MLE1$ML
      ML0s[i]<-MLE0$ML
      vgs[i]<-MLE1$vg
      ves[i]<-MLE1$ve
      deltas[i]<-MLE1$delta
      nv<-length(MLE1$etas)
      Lam<-diag(c(1/(MLE1$delta*MLE1$lambda+1),rep(1,nv-1)))
      Lam1 <- diag(c(1/(MLE1$delta*MLE1$lambda+1)^2,rep(1,nv-1)))
      temp <- crossprod(xv,MLE1$U_R)
      bhats[i] <- MLE1$delta*temp%*%Lam%*%MLE1$etas
      var.bhats.ratio[i] <- MLE1$delta^2*temp%*%Lam%*%t(temp)%*%temp%*%Lam%*%t(temp)+MLE1$delta*temp%*%Lam1%*%t(temp)
      d[i] <- (1-var.bhats.ratio[i])*((1-var.bhats.ratio[i])>=0)
      stats[i]<- 2*(MLE1$ML-MLE0$ML)
      ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2
    }else{
      ps[i]<-1
    }
  }
  return(list(ps=ps,bhats=bhats,deltas=deltas,d=d,ML1s=ML1s,ML0s=ML0s,stats=stats,vgs=vgs,ves=ves))
} 


yraw<-matrix(phe[,1],,1)
xnames<-gen[,1:2]
snp1<-deepcopy(gen,3:ncol(gen))
mydataq<-snp1[,]
mydataq2<-t(mydataq)
rm(mydataq)
gc()
mydata<-big.matrix(nrow(mydataq2),ncol(mydataq2),type="double",shared = FALSE)
mydata[,]<-mydataq2[,]
rm(mydataq2)
gc()
m<-dim(mydata)[2]
n<-dim(mydata)[1]
Y<-yraw
K<-matrix(kk,nrow=dim(kk)[1])
W0<-matrix(1,n,1)

if(is.null(psmatrix)==FALSE){
  W1<-psmatrix
  W<-cbind(W0,W1)
  
}
if(is.null(psmatrix)==TRUE){
  W<-W0 
  
}

rm(kk)
gc()

maf.fun<-function(snp){
  leng<-length(snp)
  id.1<-length(which(snp==1))
  id.0<-length(which(snp==0))
  id.0.5<-length(which(snp==0.5))
  maf.1<-id.1/m
  maf.0.5<-id.0.5/m
  maf.0<-id.0/m
  ma1<-(2*id.1+id.0.5)/(2*leng)
  ma2<-(2*id.0+id.0.5)/(2*leng)
  maf.min<-min(ma1,ma2)
  return(list(maf.1,maf.0,maf.0.5,maf.min))
}


pve.fun<-function(beta,maf){
  
  pve<-(maf$p1-maf$p1^2+0.25*maf$p3-0.25*maf$p3^2-maf$p1*maf$p3)*beta^2
  
  return(pve)
}

if(flagREMLE==1){
  
  remle1<-emma.REMLE(Y, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
  
  
  
}
if(flagREMLE==0){
  remle1<-emma.MLE(Y, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
  
}

remle1.deltahat.g<-remle1$delta
remle1.B1<-emma.maineffects.B(Z=NULL,K,remle1.deltahat.g)
C2<-remle1.B1$mC

rm(remle1.B1)
gc()

if(flagREMLE==1){
  ys=Y;xs=mydata[,]; K=1; Z=C2; X0=W; ngrids=100;llim=-10; ulim=10; esp=1e-10
  
  stopifnot(K == 1)
  ys <- Z%*%ys
  xs <- Z%*%xs
  X0 <- Z%*%X0
  
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  n <- nrow(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)
  MLE0<-emma.REMLE0.c(ys,X0)
  
  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  vgs <- vector(length=t)
  ves <- vector(length=t)
  deltas <- vector(length=t)
  bhats<-vector(length=t)
  var.bhats.ratio <- vector(length=t)
  d <- vector(length=t)
  stats <- vector(length=t)
  ps <- vector(length=t)
  
  cl.cores <- detectCores()
  if ((cl.cores<=2)||(is.null(CLO)==FALSE)){
    cl.cores<-1
  }else if(cl.cores>2){
    if(cl.cores>10){
      cl.cores<-10
    }else {  
      cl.cores <- detectCores()-1
    }
  }   
  
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl)
  
  ff=foreach(i=1:t,.combine = 'rbind')%dopar%{
    
    vids <- !is.na(xs[,i])
    xv <- xs[vids,i]
    yv <- ys[vids]
    x0v<-X0[vids,]
    MLE1 <- emma.REMLE.c (yv, x0v, K=1, xv, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
    
    if(length(MLE1$vg)!=0){
      ML1s[i]<-MLE1$REML
      ML0s[i]<-MLE0$REML
      vgs[i]<-MLE1$vg
      ves[i]<-MLE1$ve
      deltas[i] <- MLE1$delta
      nv<-length(MLE1$etas)
      Lam<-diag(c(1/(MLE1$delta*MLE1$lambda+1),rep(1,nv-1)))
      Lam1 <- diag(c(1/(MLE1$delta*MLE1$lambda+1)^2,rep(1,nv-1)))
      temp <- crossprod(xv,MLE1$U_R)
      bhats[i] <- MLE1$delta*temp%*%Lam%*%MLE1$etas
      var.bhats.ratio[i] <- MLE1$delta^2*temp%*%Lam%*%t(temp)%*%temp%*%Lam%*%t(temp)+MLE1$delta*temp%*%Lam1%*%t(temp)
      d[i] <- (1-var.bhats.ratio[i])*((1-var.bhats.ratio[i])>=0) 
      stats[i]<- 2*(MLE1$REML-MLE0$REML)
      ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2
    }else{
      ps[i]<-1 
    }
    ff<-c(ps[i],bhats[i],deltas[i],d[i],ML1s[i],ML0s[i],stats[i],vgs[i],ves[i])
  }
  stopCluster(cl)
  row.names(ff)<-NULL
  
  REML.LRT.c2<-  list(ps=ff[,1],bhats=ff[,2],deltas=ff[,3],d=ff[,4],ML1s=ff[,5],ML0s=ff[,6],stats=ff[,7],vbs=ff[,8],ves=ff[,9])
}
if(flagREMLE==0){
  ys=Y;xs=mydata[,];K=1; Z=C2; X0=W; ngrids=100;llim=-10;ulim=10;esp=1e-10
  stopifnot(K == 1)
  ys <- Z%*%ys 
  xs <- Z%*%xs
  X0 <- Z%*%X0
  
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  n <- nrow(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)
  MLE0<-emma.MLE0.c(ys,X0)
  
  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  vgs <- vector(length=t)
  ves <- vector(length=t)
  deltas<-vector(length=t)
  bhats<-vector(length=t)
  var.bhats.ratio <- vector(length=t)
  d <- vector(length=t)
  
  stats <- vector(length=t)
  ps <- vector(length=t)
  
  cl.cores <- detectCores()
  if((cl.cores<=2)||(is.null(CLO)==FALSE)){
    cl.cores<-1
  }else if(cl.cores>2){
    if(cl.cores>10){
      cl.cores<-10
    }else {  
      cl.cores <- detectCores()-1
    }
  }   
  
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl)
  
  ff=foreach(i=1:t,.combine = 'rbind')%dopar%{
    vids <- !is.na(xs[,i])
    xv <- xs[vids,i]
    yv <- ys[vids]
    x0v<-X0[vids,]
    MLE1 <- emma.MLE.c (yv, x0v, K=1, xv, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
    
    if(length(MLE1$vg)!=0){
      ML1s[i]<-MLE1$ML
      ML0s[i]<-MLE0$ML
      vgs[i]<-MLE1$vg
      ves[i]<-MLE1$ve
      deltas[i]<-MLE1$delta
      nv<-length(MLE1$etas)
      Lam<-diag(c(1/(MLE1$delta*MLE1$lambda+1),rep(1,nv-1)))
      Lam1 <- diag(c(1/(MLE1$delta*MLE1$lambda+1)^2,rep(1,nv-1)))
      temp <- crossprod(xv,MLE1$U_R)
      bhats[i] <- MLE1$delta*temp%*%Lam%*%MLE1$etas
      var.bhats.ratio[i] <- MLE1$delta^2*temp%*%Lam%*%t(temp)%*%temp%*%Lam%*%t(temp)+MLE1$delta*temp%*%Lam1%*%t(temp)
      d[i] <- (1-var.bhats.ratio[i])*((1-var.bhats.ratio[i])>=0)
      stats[i]<- 2*(MLE1$ML-MLE0$ML)
      ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2
    }else{
      ps[i]<-1
    }
    ff<-c(ps[i],bhats[i],deltas[i],d[i],ML1s[i],ML0s[i],stats[i],vgs[i],ves[i])
  }
  stopCluster(cl)
  row.names(ff)<-NULL
  REML.LRT.c2<-list(ps=ff[,1],bhats=ff[,2],deltas=ff[,3],d=ff[,4],ML1s=ff[,5],ML0s=ff[,6],stats=ff[,7],vbs=ff[,8],ves=ff[,9])
}

rm(Z,xs)
gc()

REML.LRT.c2.new<-data.frame(REML.LRT.c2)
mafall<-apply(snp1[,],1,maf.fun)

rm(snp1)
gc()

mafall1<-unlist(mafall)
mafall2<-matrix(mafall1,nrow = 4)
mafall3<-t(mafall2)
mafall4<-data.frame(mafall3)
names(mafall4)<-c("p1","p2","p3","maf")
MAF<-mafall4$maf

var.bb<-apply((C2%*%mydata[,]),2,var)*REML.LRT.c2.new$bhats^2
pve.allr2<-var.bb/apply(cbind(matrix(var(C2%*%Y),nrow = dim(REML.LRT.c2.new)[1]),(var.bb+REML.LRT.c2.new$ves)),1,max)


rm(C2,mydata)
gc()

parms<-data.frame(chr.locus=xnames,REML.LRT.c2.new,MAF,pve.allr2)
names(parms)<-NULL
parms<-as.matrix(parms)
ress<-parms[,1:3]
ress1<-ress[ress[,3]!=1,]
resp<-as.matrix(ress1[,3])
pmin<-min(resp[resp!=0])
locsub<-which(resp==0)

if(length(locsub)!=0){
  subvalue<-10^(1.1*log10(pmin))
  ress1[locsub,3]<-subvalue
  ress1<-ress1
}else{
  ress1<-ress1
}
rowsnp <- dim(ress1)[1]
snpname <- numeric()
snpname<-as.matrix(paste("rs",c(1:rowsnp),sep = ""))
snpname<-snpname
pe<-as.numeric(parms[,6])
pee<-sum(pe)
newp<-0.05/pee
mannewp<-newp
manstandchoice<-1
ps<-as.matrix(parms[,3])
alpha<-1
obs.x<-sort(ps)
newobs.x<-obs.x[obs.x<alpha]
n<-length(newobs.x)
es<-(1:n)/(n+1)
x<--log10(es)
y<--log10(newobs.x*2)
tt<-qchisq(newobs.x*2,df=1,lower.tail=F)
GIF.median<-median(tt)/0.4549
GIF.mean<-mean(tt)
GIF.regress<-crossprod(x,y)/crossprod(x,x)
vif<-GIF.mean
VIF<-matrix("",dim(parms)[1],1)
VIF[1]<-vif
parmeter<-cbind(parms[,1:5],parms[,10:13])
parmeter[,3]<--log10(parmeter[,3])
parmeter[which(abs(parmeter)>1e-4)]<-round(parmeter[which(abs(parmeter)>1e-4)],4)
parmeter[which( abs(parmeter)<1e-4)]<-as.numeric(sprintf("%.4e", parmeter[which( abs(parmeter)<1e-4)]))
parmeter[which(abs(parmeter[,5])>1e-4),5]<-round(parmeter[which(abs(parmeter[,5])>1e-4),5],4)
parmeter[which(abs(parmeter[,5])<1e-4),5]<-as.numeric(sprintf("%.4e",parmeter[which(abs(parmeter[,5])<1e-4),5]))

if(inputform==1){
  meadd<-matrix(" ",nrow(parms),1)
  if(newp>1e-4){
    newp<-round(newp,4)
  }else{
    newp<-sprintf("%.4e",newp)
  }
  meadd[1,1]<-newp
  parmsShow<-cbind(genRaw[-1,1],parmeter,genRaw[-1,4],meadd)
  colnames(parmsShow)<-c("Marker","Chromosome","Marker position (bp)","-log10(P)","SNP effect","QTN-to-residual variance ratio","QTN variance","Residual variance","MAF","r2 (%)","Genotype for code 1","Critical P-value")
  
}
if(inputform==2){
  outATCG<-matrix(outATCG,,1)
  meadd<-matrix(" ",nrow(parms),1)
  if(newp>1e-4){
    newp<-round(newp,4)
  }else{
    newp<-sprintf("%.4e",newp)
  }
  meadd[1,1]<-newp
  parmsShow<-cbind(genRaw[-1,1],parmeter,outATCG,meadd)
  colnames(parmsShow)<-c("Marker","Chromosome","Marker position (bp)","-log10(P)","SNP effect","QTN-to-residual variance ratio","QTN variance","Residual variance","MAF","r2 (%)","Genotype for code 1","Critical P-value")
  
}
if(inputform==3){
  outATCG<-matrix(outATCG,,1)
  outATCG<-unlist(strsplit(outATCG,""))
  outATCG<-matrix(outATCG[c(TRUE,FALSE)],,1)
  
  meadd<-matrix(" ",nrow(parms),1)
  if(newp>1e-4){
    newp<-round(newp,4)
  }else{
    newp<-sprintf("%.4e",newp)
  }
  meadd[1,1]<-newp
  parmsShow<-cbind(genRaw[-1,1],parmeter,outATCG,meadd)
  colnames(parmsShow)<-c("Marker","Chromosome","Marker position (bp)","-log10(P)","SNP effect","QTN-to-residual variance ratio","QTN variance","Residual variance","MAF","r2 (%)","Genotype for code 1","Critical P-value")
}

bpnumber <- numeric()
chrnum <- unique(ress1[,1])
for(i in 1:length(chrnum))
{
  bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(ress1[,1]==chrnum[i])))))
}
parms <- data.frame(ress1,snpname,bpnumber)
colnames(parms)<-c("Chromosome","Position","P-value","SNPname","BPnumber")
parms<-parms[,-2]

mannewp <- as.matrix(mannewp)
rowbl<-matrix("",(nrow(parms)-1),1)
mannepr<-rbind(mannewp,rowbl)
colnames(mannepr)<-"Manhattan p-value"

parms<-cbind(as.matrix(parms),mannepr)
parms<-as.data.frame(parms,stringsAsFactors=FALSE)
parms[,c(1,2,4)]<-sapply(parms[,c(1,2,4)],as.numeric)
ress1<-as.data.frame(ress1[,-(1:2)])

colnames(ress1)<-"p-value"

Xemma<-data.frame(chr.locus=xnames,REML.LRT.c2.new)
vid<-which(as.numeric(Xemma$ps)<=svpal)

if(length(vid)!=0){
  if(length(vid)==1){
    snp.emma.opt<-matrix(gen[vid,],1,)
    xname.emma.opt<-matrix(snp.emma.opt[,1:2],1,)
    mafall4.opt<-mafall4[vid,]
    snp4<-matrix(snp.emma.opt[,3:dim(snp.emma.opt)[2]],1,)
    xdata<-t(snp4)
    xdata<-matrix(xdata,,1)
  }else{
    snp.emma.opt<-as.matrix(gen[vid,])
    xname.emma.opt<-snp.emma.opt[,1:2]
    mafall4.opt<-mafall4[vid,]
    snp4<-snp.emma.opt[,3:dim(snp.emma.opt)[2]]
    xdata<-t(snp4)
  }
  
  xdata<-t(snp4)
  ydata<-Y
  u1<-ebayes_EM(x=W,z=xdata,y=ydata)
  emma.lod<-likelihood(xxn=W,xxx=xdata,yn=ydata,bbo=u1$u)
  idslod<-which(emma.lod>=svmlod)
  
  if(length(idslod)!=0){
    maf.snp.4<-mafall4.opt[idslod,]
    if(length(idslod)==1){
      chrlocus<-matrix(xname.emma.opt[idslod,],1,)
    }else{
      chrlocus<-as.matrix(xname.emma.opt[idslod,])
    }
    
    
    pve.all.1<-pve.fun(u1$u[idslod],maf.snp.4)
    pve.all<-pve.all.1/as.vector(max(var(Y),(sum(pve.all.1)+u1$sigma2)))*100
    
    qtneffect<-matrix(u1$u[idslod],,1)
    lodscore<-matrix(emma.lod[idslod],,1)
    log10P <- as.matrix(-log10(1-pchisq(lodscore*4.605,1)))
    maff<-matrix(maf.snp.4$maf,,1)
    r2<-matrix(pve.all,,1)
    wanbefore<-cbind(qtneffect,lodscore,log10P,r2,maff)
    
    wanbefore[which(abs(wanbefore)>1e-4)]<-round(wanbefore[which(abs(wanbefore)>1e-4)],4)
    wanbefore[which(abs(wanbefore)<1e-4)]<-as.numeric(sprintf("%.4e", wanbefore[which(abs(wanbefore)<1e-4)]))
    wanbefore <- matrix(wanbefore,,5)
    
    wan<-cbind(chrlocus,wanbefore)
    phenotype.var<-var(Y)
    sigma2<-u1$sigma2
    pee<-matrix("",dim(wan)[1],1)
    vess<-matrix("",dim(wan)[1],1)
    pee[1]<-round(phenotype.var,4)
    vess[1]<-round(sigma2,4)
 
    if(inputform==1){
      genRaw<-as.data.frame(genRaw)
      genraw<-genRaw[-1,1:4]
      
      wan_len<-dim(wan)[1]
      marker<-character()
      snp<-character()
      
      for(i in 1:wan_len){
        chr_pos<-which(genraw[,2]==wan[i,1])
        new_matrix<-genraw[chr_pos,]
        posi_pos<-which(new_matrix[,3]==wan[i,2])
        mark<-matrix(new_matrix[posi_pos,1],1,)
        marker<-rbind(marker,mark)
        sn<-matrix(new_matrix[posi_pos,4],1,)
        snp<-rbind(snp,sn)
      }
    }
    if(inputform==2){
      
      genRaw<-as.data.frame(genRaw)
      genraw<-genRaw[-1,1:4]
      
      wan_len<-dim(wan)[1]
      marker<-character()
      snp<-character()
      for(i in 1:wan_len){
        chr_pos<-which(genraw[,2]==wan[i,1])
        new_matrix<-genraw[chr_pos,]
        posi_pos<-which(new_matrix[,3]==wan[i,2])
        mark<-matrix(new_matrix[posi_pos,1],1,)
        marker<-rbind(marker,mark)
        sn<-matrix(new_matrix[posi_pos,4],1,)
        snp<-rbind(snp,sn)
      }
      
    }
    if(inputform==3){
      genRaw<-as.data.frame(genRaw)
      genraw<-genRaw[-1,c(1,3,4,12)]
      
      wan_len<-dim(wan)[1]
      marker<-character()
      snp<-character()
      for(i in 1:wan_len){
        chr_pos<-which(genraw[,2]==wan[i,1])
        new_matrix<-genraw[chr_pos,]
        posi_pos<-which(new_matrix[,3]==wan[i,2])
        mark<-matrix(new_matrix[posi_pos,1],1,)
        marker<-rbind(marker,mark)
        sn<-matrix(new_matrix[posi_pos,4],1,)
        snp<-rbind(snp,sn)
      }
    }
    
    wan<-cbind(marker,wan,snp,vess,pee)
    colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","-log10(P)","r2 (%)","MAF","Genotype  for code 1","Var_Error","Var_phen(total)")
    wan<-as.data.frame(wan)
    
   }  
  }  
  parmsShow<-parmsShow[,-c(6,7,8,9,10,12)]
  parmsShow<-parmsShow[,c(1,2,3,5,4,6)]

  output<-list(result1=parmsShow,result2=wan,Manhattan=parms,QQ=ress1)

  return(output)
  }
 }
}
