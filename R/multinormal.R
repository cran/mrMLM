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
  
  iter<-0;err<-1000;iter_max<-200;err_max<-1e-8
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

Plot<-function(plotresult=NULL,color1=NULL,color2=NULL,p_stand=NULL,method=NULL,type=NULL){
  
  Manhattan<-function(plotresult,color1,color2){
    parms<-as.data.frame(plotresult)
    mannewp<-as.numeric(parms[1,5])
    svgwline<-round(-log10(mannewp),4)  
    standline<-svgwline
    manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P-value",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = standline)
  }
  
  
  QQplot1<-function(plotresult,p_stand,color1,color2){

    p_value<-as.matrix(plotresult)
    pvalue<-matrix(p_value,,1)
    observed<-sort(pvalue[,1])
    observed<-observed/2
    observed<-observed[which(observed!=0)]
    newobserved<-observed[which(observed<(p_stand/2))]
    lobs<--(log10(newobserved))
    expected<-c(1:length(newobserved))
    lexp<--(log10(expected/(length(pvalue)+1)))
    plot(lexp,lobs,xlim=c(0,max(lexp)),ylim=c(0,max(lobs)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color2)
    abline(0,1,col=color1)
  }
  
  QQplot2<-function(plotresult,color1,color2){

    ress1<-as.data.frame(plotresult)
    pvalue<-as.matrix(ress1)
    ps<-pvalue[,1]
    obs.x<-sort(ps)
    newobs.x<-obs.x[obs.x<1]
    n<-length(newobs.x)
    es<-(1:n)/(n+1)
    x<--log10(es)
    y<--log10(newobs.x)
    y<-y-0.3
    plot(x,y,xlim=c(0.3,max(x)),ylim=c(0.3,max(y)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color2)
    abline(0,1,col=color1)
    
  }
  
  
  LOD<-function(fileplot=NULL,color1,method=NULL){
    
    data<-as.matrix(plotresult)
    data<-as.data.frame(data,stringsAsFactors = F)
    gen<-data[,1:2]
    resulty<-data[,3:5]
    resultkq<-as.matrix(resulty)
    resultk<-which(resultkq=="",arr.ind = TRUE)
    resultq<-resulty[1:(resultk[1]-1),] 
    
    if(nrow(resultq)>1){
      result<-resultq
    }else{
      result<-t(as.matrix(resultq))
    }
    
    galaxyy<-as.data.frame(result)
    galaxyy<-sapply(galaxyy,as.numeric)
    chr_pos <- gen[,1:2]
    chr_pos<-sapply(chr_pos,as.numeric)
    
    chr_num <- length(unique(chr_pos[,1]))
    chr <- matrix(0,chr_num,1)
    pos <- matrix(0,chr_num,1)
    for(i in 1:chr_num)
    {
      
      temp <- numeric()
      temp <- length(which(chr_pos[,1]==i))
      if(i==1)
      {
        pos[i] <- temp
        chr[i] <- chr_pos[pos[i],2]
      }else{
        pos[i] <- pos[i-1] + temp
        chr[i] <- chr_pos[pos[i],2]
      }
    }
    
    pos_acc <- matrix(0,chr_num,1)
    for(i in 1:chr_num)
    {
      if(i==1){
        pos_acc[i] <- chr[i]
      }else{
        pos_acc[i] <- pos_acc[i-1] + chr[i]
      }
    }
    
    newres_pos <- galaxyy[,2]
    res_sumpos <- pos_acc[galaxyy[which(galaxyy[,1]>1),1]-1] + galaxyy[which(galaxyy[,1]>1),2] 
    newres_pos[which(galaxyy[,1]>1)] <- res_sumpos 
    pospic<-c(newres_pos)
    lodpic<-c(galaxyy[,3])  
    mm<-round(max(pospic)/4000)
    mm<-as.numeric(format(mm,digits = 1,scientific = TRUE))
    pospicx<-pospic/mm
    if(pospicx[1]<20){
      pospicx[1]<-pospicx[1]+20
    }
    pos_acc1<-pos_acc/mm
    resdf1 <- data.frame(pospicx,lodpic)
    
    pp <- ggplot(data=resdf1, aes(x=pospicx, y=lodpic)) +
      geom_bar(stat="identity", width=0.5, fill="white", linetype="solid",color=color1)
    
    pp <- pp + geom_vline(xintercept=c(0,pos_acc1),linetype="dashed",alpha=0.2)
    pp <- pp  + scale_x_continuous(expand=c(0,0),limits=c(0,(pos_acc1[dim(pos_acc1)[1]]+100))) +
      scale_y_continuous(expand=c(0,0))
    pp <- pp + xlab(paste("Genome position (",mm,"bp)",sep = "")) + ylab("LOD score") + ggtitle("") + theme_classic()
    pp <- pp + theme(axis.title.y = element_text( vjust = 2,hjust=0.5,size = 14),
                     axis.title.x = element_text(vjust = -0.5,hjust=0.5,size = 14))
    
    pp <- pp + theme(panel.background = element_rect(fill = "white"))
    pp <- pp + theme(text=element_text(family="mono"))
    pp <- pp + theme(axis.line.y = element_line(colour = "black", linetype = "solid"),
                     axis.line.x = element_line(colour = "black", linetype = "solid"))
    print(pp)
    
  }  
  
  if(type=="Manhattan"){
    
    Manhattan(plotresult,color1,color2)
    
  }else if(type=="qq"){
    
    if(method=="FASTmrEMMA"){
      QQplot2(plotresult,color1,color2)  
    }else{
      QQplot1(plotresult,p_stand,color1,color2) 
    }
    
  }else if(type=="LOD"){
    
    LOD(plotresult,color1)
    
  }  
}


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
    dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }
  
  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    snps <- snps[rowSums(is.na(snps))==0,]
  }
  
  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1
  
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }
  }
  return(K)
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

#restricted likelihood
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