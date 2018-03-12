
mrMLMFun<-function(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svrad,svmlod,Genformat,CLO){
  inputform<-Genformat
  
  if(is.null(kk)){
    envgen <- gen
    if(exists("envgen")==FALSE)
    {
      warning("Please input correct genotypic dataset !")
    }else{
      envgen<-envgen[,3:(ncol(envgen))]
      envgen<-t(envgen)
      m<-ncol(envgen)
      n<-nrow(envgen)
      kk1<-matrix(0,n,n)
      for(k in 1:m){
        z<-as.matrix(envgen[,k])
        kk1<-kk1+z%*%t(z)
      }
      cc<-mean(diag(kk1))
      kk1<-kk1/cc
      kk<-as.matrix(kk1)
    }
  } 
  
  if(is.null(psmatrix)){
    flagps<-1
  }else{
    flagps<-0
  }
  
  if(is.null(svpal)==TRUE||is.null(svrad)==TRUE||is.null(svmlod)==TRUE){
    warning("Please set parameters!")
  }
  
  if((svpal<0)||(svpal>1))
  {
    warning("Please input critical P-value between 0 and 1!")
  }
  if(svrad<0)
  {
    warning("Please input search radius of candidate gene: > 0 !")
  }
  if(svmlod<0)
  {
    warning("Please input critical LOD score: > 0 !")
  }
  if(exists("gen")==FALSE)
  {
    warning("Please input correct genotypic dataset !")
  }
  if(exists("phe")==FALSE)
  {
    warning("Please input correct phenotypic dataset !")
  }
  if(exists("kk")==FALSE)
  {
    warning("Please input correct kinship (K) dataset !")
  }
  if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(ncol(gen)!=(nrow(phe)+2)))
  {
    warning("Sample size in genotypic dataset doesn't equal to the sample size in phenotypic dataset !")
  }
  
  if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2)))&&(svpal>=0)&&(svpal<=1)&&(svrad>0)&&(svmlod>=0))
  {
    parmsShow<-NULL
    wan<-NULL
    parms<-NULL
    parms.pchange<-NULL
    mannewp<-NULL
    
    mixed<-function(x,y,kk){
      
      loglike<-function(theta){
        lambda<-exp(theta)
        logdt<-sum(log(lambda*delta+1))
        h<-1/(lambda*delta+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,q,1)
        xx<-matrix(0,q,q)
        for(i in 1:q){
          yx[i]<-sum(yu*h*xu[,i])
          for(j in 1:q){
            xx[i,j]<-sum(xu[,i]*h*xu[,j])
          }
        }
        loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
        return(-loglike)
      }
      
      fixed<-function(lambda){
        h<-1/(lambda*delta+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,q,1)
        xx<-matrix(0,q,q)
        for(i in 1:q){
          yx[i]<-sum(yu*h*xu[,i])
          for(j in 1:q){
            xx[i,j]<-sum(xu[,i]*h*xu[,j])
          }
        }
        beta<-solve(xx,yx)
        sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
        sigma2<-drop(sigma2)
        var<-diag(solve(xx)*sigma2)
        stderr<-sqrt(var)
        return(c(beta,stderr,sigma2))
      }
      
      qq<-eigen(kk)
      delta<-qq[[1]]
      uu<-qq[[2]]
      q<-ncol(x)
      n<-ncol(kk)
      vp<-var(y)
      yu<-t(uu)%*%y
      xu<-t(uu)%*%x
      theta<-0
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
      lambda<-exp(parm$par)
      conv<-parm$convergence
      fn1<-parm$value
      fn0<-loglike(-Inf)
      lrt<-2*(fn0-fn1)
      hess<-parm$hessian
      parmfix<-fixed(lambda)
      beta<-parmfix[1:q]
      stderr<-parmfix[(q+1):(2*q)]
      sigma2<-parmfix[2*q+1]
      lod<-lrt/4.61
      p_value<-1-pchisq(lrt,1)
      sigma2g<-lambda*sigma2
      goodness<-(vp-sigma2)/vp
      par<-data.frame(lrt,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
      return(par)
    }
    
    loglike<-function(theta){
      xi<-exp(theta)
      tmp0<-zz*xi+1
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      logdt2<-log(det(tmp0))
      loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
      return(-loglike)
    }
    
    fixed<-function(xi){
      tmp0<-zz*xi+diag(1)
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      zHy<-zy-zz%*%tmp%*%zy
      zHx<-zx-zx%*%tmp%*%zz
      zHz<-zz-zz%*%tmp%*%zz
      beta<-solve(xHx,yHx)
      tmp2<-solve(xHx)
      sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
      gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
      var<-abs((xi*diag(1)-xi*zHz*xi)*as.numeric(sigma2))
      stderr<-sqrt(diag(var))
      result<-list(gamma,stderr,beta,sigma2)
      return(result)
    }
    
    name<-gen[,1:2]
    gen<-gen[,3:(ncol(gen))]
    gen<-t(gen)
    n<-nrow(gen)
    m<-ncol(gen)
    if((flagps==1)||(exists("psmatrix")==FALSE))
    {
      x<-matrix(1,n,1)
    }else if(flagps==0)
    {
      x<-cbind(matrix(1,n,1),psmatrix)
    }
    ll<-numeric()
    s<-ncol(x)
    kk<-as.matrix(kk)
    qq<-eigen(kk)
    delta<-qq[[1]]
    uu<-qq[[2]]
    xu<-t(uu)%*%x
    
    for(ii in 1:1)
    {
      yy<-phe[,1]
      y<-as.matrix(yy)
      parm<-mixed(x=x,y=y,kk=kk)
      lambda<-parm$lambda[1]
      h<-1/(delta*lambda+1)
      yu<-t(uu)%*%y
      xx<-matrix(0,s,s)
      for(i in 1:s){
        for(j in 1:s){
          xx[i,j]<-sum(xu[,i]*h*xu[,j])
        }
      }
      yy<-sum(yu*h*yu)
      yx<-matrix(0,s,1)
      for(i in 1:s){
        yx[i]<-sum(yu*h*xu[,i])
      }
      
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
      
      if((flagps==1)||(is.null("psmatrix")))
      {
        ff=foreach(k=1:m, .multicombine=TRUE, .combine = 'rbind')%dopar%
        {
          browser()
          z<-as.matrix(gen[,k])
          zu<-t(uu)%*%z
          zy<-as.matrix(sum(yu*h*zu))
          zz<-as.matrix(sum(zu*h*zu))
          zx<-matrix(0,s,1)
          for(i in 1:s){
            zx[i]<-sum(xu[,i]*h*zu)
          }
          theta<-c(0)
          par<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10)
          xi<-exp(par$par)
          conv<-par$convergence
          fn1<-par$value
          hess<-par$hessian
          parmfix<-fixed(xi)
          gamma<-parmfix[[1]]
          stderr<-parmfix[[2]]
          
          beta<-parmfix[[3]]
          
          sigma2<-parmfix[[4]]
          lambda<-xi
          sigma2g<-lambda*sigma2
          fn0<-loglike(-Inf)
          lrt<-2*(fn0-fn1)
          p_lrt<-1-pchisq(lrt,1)
          wald<-(gamma/stderr)^2
          p_wald<-1-pchisq(wald,1)
          parm0<-c(ii,name[k,1],name[k,2],beta,sigma2,sigma2g,gamma,stderr,wald,p_wald)
          
        }
        stopCluster(cl)
        
        ll<-rbind(ll,ff)
      }else if(flagps==0){
        ff=foreach(k=1:m, .multicombine=TRUE, .combine = 'rbind')%dopar%
        {
          browser()
          z<-as.matrix(gen[,k])
          zu<-t(uu)%*%z
          zy<-as.matrix(sum(yu*h*zu))
          zz<-as.matrix(sum(zu*h*zu))
          zx<-matrix(0,s,1)
          for(i in 1:s){
            zx[i]<-sum(xu[,i]*h*zu)
          }
          theta<-c(0)
          par<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10)
          xi<-exp(par$par)
          conv<-par$convergence
          fn1<-par$value
          hess<-par$hessian
          parmfix<-fixed(xi)
          gamma<-parmfix[[1]]
          stderr<-parmfix[[2]]
          
          beta<-parmfix[[3]][1]
          
          sigma2<-parmfix[[4]]
          lambda<-xi
          sigma2g<-lambda*sigma2
          fn0<-loglike(-Inf)
          lrt<-2*(fn0-fn1)
          p_lrt<-1-pchisq(lrt,1)
          wald<-(gamma/stderr)^2
          p_wald<-1-pchisq(wald,1)
          parm0<-c(ii,name[k,1],name[k,2],beta,sigma2,sigma2g,gamma,stderr,wald,p_wald)
          
        }
        stopCluster(cl)
        
        ll<-rbind(ll,ff)
      }
    }
    
    parms<-ll
    parms<-matrix(parms,,10)
    
    gen<-t(gen)
    chr_pos<-parms[,2:3]
    pfit<-which(parms[,10]<=(svpal))
    pfit<-as.matrix(pfit)
    pfitrow<-nrow(pfit)
    no_p<-cbind((1:(nrow(parms))),parms[,10])
    no_porder<-order(no_p[,2])
    no_p<-no_p[no_porder,]
    choose_orderp<-no_p[1:pfitrow,]
    orderno<-no_p[1:pfitrow,1]
    orderno<-as.matrix(orderno)
    sigma2g_SNPerr<-cbind(parms[,6],parms[,8])
    correct_each<-matrix(1,(nrow(sigma2g_SNPerr)),1)-sigma2g_SNPerr[,2]*sigma2g_SNPerr[,2]/sigma2g_SNPerr[,1]
    k0<-which(correct_each<0)
    k0<-as.matrix(k0)
    if(nrow(k0)>0){
      correct_each[k0,1]<-matrix(0,(nrow(k0)),1)
    }
    correct_sum<-sum(correct_each)
    newp<-0.05/correct_sum
    mannewp<-newp
    manstandchoice<-1
    no_porder<-which(no_p[,2]<=newp)
    no_porder<-as.matrix(no_porder)
    no_porderrow<-nrow(no_porder)
    gg<-orderno
    for (ii in 1:(nrow(orderno)-1)){
      for (jj in (ii+1):(nrow(orderno))){
        ci<- chr_pos[orderno[ii],1]
        cj<- chr_pos[orderno[jj],1]
        if (ci==cj){
          ye<-abs(chr_pos[orderno[ii],2]-chr_pos[orderno[jj],2])
          if (ye<=((svrad)*1000)){
            gg[jj,1]<-0
          }
        }
      }
    }
    
    parms.pchange<-parms
    parmsp<-as.matrix(parms.pchange[,10])
    locsub<-which(parmsp==0)
    if(length(locsub)!=0){
      pmin<-min(parmsp[parmsp!=0])
      subvalue<-10^(1.1*log10(pmin))
      parms.pchange[locsub,10]<-subvalue
    }else{
      parms.pchange<-parms
    }
    
    if(inputform==1){
      #output result1 using mrMLM numeric format
      parmsShow<-parms[,-1]
      meadd<-matrix(1,nrow(parms),1)
      meadd[which(parms[,10]<newp),1]<-sprintf("%.4e",newp)
      meadd[which(parms[,10]>=newp),1]<-"  "
      tempparms<-parms[,4:10]
      tempparms[,7]<--log10(tempparms[,7])
      tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
      tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
      parmsShow<-cbind(genRaw[-1,1],parms[,2:3],tempparms,genRaw[-1,4],meadd)
      colnames(parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","-log10(P)","Genotype for code 1","Significance")
    }
    if(inputform==2){
      #output result1 using mrMLM character format
      parmsShow<-parms[,-1]
      outATCG<-matrix(outATCG,,1)
      meadd<-matrix(1,nrow(parms),1)
      meadd[which(parms[,10]<newp),1]<-sprintf("%.4e",newp)
      meadd[which(parms[,10]>=newp),1]<-"  "
      tempparms<-parms[,4:10]
      tempparms[,7]<--log10(tempparms[,7])
      tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
      tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
      parmsShow<-cbind(genRaw[-1,1],parms[,2:3],tempparms,outATCG,meadd)
      colnames(parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","-log10(P)","Genotype  for code 1","Significance")
    }
    if(inputform==3){
      #output result1 using TASSEL format
      parmsShow<-parms[,-1]
      outATCG<-matrix(outATCG,,1)
      outATCG<-unlist(strsplit(outATCG,""))
      outATCG<-matrix(outATCG[c(TRUE,FALSE)],,1)
      meadd<-matrix(1,nrow(parms),1)
      meadd[which(parms[,10]<newp),1]<-sprintf("%.4e",newp)
      meadd[which(parms[,10]>=newp),1]<-"  "
      tempparms<-parms[,4:10]
      tempparms[,7]<--log10(tempparms[,7])
      tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
      tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
      parmsShow<-cbind(genRaw[-1,1],parms[,2:3],tempparms,outATCG,meadd)
      colnames(parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","-log10(P)","Genotype  for code 1","Significance")
    }
    rowsnp <- dim(parms)[1]
    snpname <- numeric()
    snpname <- as.matrix(paste("rs",c(1:rowsnp),sep=""))
    
    bpnumber <- numeric()
    chrnum <- unique(parms[,2])
    for(i in 1:length(chrnum))
    {
      bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(parms[,2]==chrnum[i])))))
    }
    parms <- data.frame(parms.pchange,snpname,bpnumber)
    colnames(parms)<-c("Trait","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P-value","SNPname","BPnumber")
    parms<-parms[,-c(1,3,4,5,6,7,8,9)]
    
    mannewp <- as.matrix(mannewp)
    rowbl<-matrix("",(nrow(parms)-1),1)
    mannepr<-rbind(mannewp,rowbl)
    
    colnames(mannepr)<-"Manhattan p-value"
    
    parms<-cbind(as.matrix(parms),mannepr)
    parms<-as.data.frame(parms,stringsAsFactors=FALSE)
    parms[,c(1,2,4)]<-sapply(parms[,c(1,2,4)],as.numeric)
    
    parms.pchange<-as.data.frame(parms.pchange[,-(1:9)])
    colnames(parms.pchange)<-"p-value"
    
    gg<-as.matrix(gg)
    misfit<-numeric()
    kk<- numeric()
    kk0<- numeric()
    l0<- numeric()
    bong<-no_porderrow
    if (bong>0){
      g0<-gg[1:no_porderrow,1]
      g0<-as.matrix(g0)
      kk0<-no_porderrow
      no_porderrow<-which(g0>0)
      no_porderrow<-as.matrix(no_porderrow)
      g0<-g0[no_porderrow,1]
      g0<-as.matrix(g0)
      xxx0<-gen[g0,]
      if(dim(g0)[1]==1){
        xxx0<-as.matrix(xxx0)
      }
      if(dim(g0)[1]>1)
      {
        xxx0<-as.matrix(xxx0)
        xxx0<-t(xxx0)
      }
      phe<-as.matrix(phe)
      if((flagps==1)||(exists("psmatrix")==FALSE))
      {
        par<-likelihood(matrix(1,(nrow(xxx0)),1),xxx0,phe,bbo=NULL)
        lod<-par
      }else if(flagps==0)
      {
        temp<-cbind(matrix(1,(nrow(xxx0)),1),psmatrix)
        par<-likelihood(temp,xxx0,phe,bbo=NULL)
        lod<-par
      }
      kk<-which(lod>=1.5)
      kk<-as.matrix(kk)
      kk1<-which(lod<1.5)
      kk1<-as.matrix(kk1)
      if ((nrow(kk1))>0){
        misfit<-g0[kk1,1]
        misfit<-as.matrix(misfit)
      }
      if ((nrow(kk))>0){
        g0<-as.matrix(g0)
        g0<-g0[kk,1]
        xx0<-xxx0[,kk]
        lo<-lod[kk,1]
      }
      if ((nrow(kk))==0){kk<-0}
    }
    if (bong==0){
      kk0<-0
      kk<-0
    }
    nleft<-as.matrix(gg[(kk0+1):(nrow(gg)),1])
    if ((length(misfit))>0){gg<-rbind(nleft,misfit)}
    if ((length(misfit))==0){gg<-nleft}
    a1<-which(gg>0)
    a1<-as.matrix(a1)
    a2<-gg[a1,1]
    a2<-as.matrix(a2)
    xx<-t(gen[a2,])
    xx<-as.matrix(xx)
    if((flagps==1)||(exists("psmatrix")==FALSE))
    {
      if (length(kk)>1){xin<-cbind(matrix(1,(nrow(xx)),1),xx0)}
      if (length(kk)==1){
        if(kk==0){
          xin<- matrix(1,(nrow(xx)),1)
        }
        if(kk>0){
          xin<-cbind(matrix(1,(nrow(xx)),1),xx0)
        }
      }
    }else if(flagps==0)
    {
      temp<-cbind(matrix(1,(nrow(xx)),1),psmatrix)
      if (length(kk)>1){xin<-cbind(temp,xx0)}
      if (length(kk)==1){
        if(kk==0){
          xin<-temp
        }
        if(kk>0){
          xin<-cbind(temp,xx0)
        }
      }
    }
    xin<-as.matrix(xin)
    par1<-ebayes_EM(xin,xx,phe)
    par<-par1$wang
    
    w2<-which(par[,1]<=0.01)
    
    if(length(w2)!=0){
      w2<-as.matrix(w2)
      ww<- numeric()
      if ((nrow(w2))>0){
        orderno<-a2[w2,1]
        orderno<-as.matrix(orderno)
        x3<-cbind(xin,xx[,w2])
        x3<-as.matrix(x3)
        lodfix<-matrix(x3[,1],nrow(x3),)
        lodrand<-matrix(x3[,2:(ncol(x3))],nrow(x3),)
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          lod<-likelihood(lodfix,lodrand,phe,bbo=NULL)
        }else if(flagps==0)
        {
          temp<-cbind(psmatrix,lodfix)
          lod<-likelihood(temp,lodrand,phe,bbo=NULL)
        }
        w3<-which(lod[,1]>=(svmlod))
        w3<-as.matrix(w3)
        if ((kk[1])>0){
          g0<-as.matrix(g0)
          orderno<-rbind(g0,orderno)
          orderno<-as.matrix(orderno)
        }
        if ((w3[1])>0){
          if((flagps==1)||(exists("psmatrix")==FALSE))
          {
            lo<-lod[w3,1]
            ww<-orderno[w3,]
          }else if(flagps==0)
          {
            lo<-lod[w3,1]
            no_loc<-w3-ncol(psmatrix)
            ww<-orderno[no_loc,]
          }
        }
        if ((nrow(w3))==0){ww<-0}
      }
      if ((nrow(w2))==0){
        g0<-as.matrix(g0)
        lo<-as.matrix(lo)
        yang<-which(lo>=(svmlod))
        yang<-as.matrix(yang)
        if ((nrow(yang))>0){
          ww<-g0[yang,1]
          lo<-lo[yang,1]
        }
        if ((nrow(yang))==0){ww<-0}
      }
      ww<-as.matrix(ww)
      needww<-ww
      if (length(ww)>=1){
        
        if (length(ww)>1){
          if((flagps==1)||(exists("psmatrix")==FALSE))
          {
            ex<-cbind(matrix(1,(nrow(xx)),1),t(gen[ww,]))
          }else if(flagps==0)
          {
            ex<-cbind(cbind(matrix(1,(nrow(xx)),1),psmatrix),t(gen[ww,]))
          }
          
        }else{
          if((flagps==1)||(exists("psmatrix")==FALSE))
          {
            ex<-cbind(matrix(1,(nrow(xx)),1),as.matrix(gen[ww,]))
          }else if(flagps==0)
          {
            ex<-cbind(cbind(matrix(1,(nrow(xx)),1),psmatrix),as.matrix(gen[ww,]))
          }  
        }
        ex<-as.matrix(ex)
        cui<-det(t(ex)%*%ex)
        p1<-rep(1,ncol(ex))
        p2<-diag(p1)
        if (cui<1e-6){bbbb<-solve(t(ex)%*%ex+p2*0.01)%*%t(ex)%*%phe}
        if (cui>=1e-6){ bbbb<-solve(t(ex)%*%ex)%*%t(ex)%*%phe }
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          eeff<-bbbb[2:(nrow(bbbb)),1]
        }else if(flagps==0)
        {
          eeff<-bbbb[(2+ncol(psmatrix)):(nrow(bbbb)),1]
        }
        
        eeff<-as.matrix(eeff)
        er<-as.numeric()
        her<-as.numeric()
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          excol<-ncol(ex)
          for(i in 1:(excol-1))
          {
            em<-ex[,(1+i)]
            as1<-length(which(em==1))/nrow(ex)
            as2<-1-as1
            er<-rbind(er,(1-(as1-as2)*(as1-as2))*eeff[i]*eeff[i])
          }
          v0<-(1/(nrow(ex)-1))*(t(phe-ex%*%bbbb)%*%(phe-ex%*%bbbb))
          
          if(var(phe)>=(sum(er)+v0)){
            her<-(er/as.vector(var(phe)))*100 
          }else{
            her<-(er/as.numeric(sum(er)+v0))*100
          }
          
        }else if(flagps==0)
        {
          excol<-ncol(ex)
          for(i in 1:(excol-1-ncol(psmatrix)))
          {
            em<-ex[,(1+ncol(psmatrix)+i)]
            as1<-length(which(em==1))/nrow(ex)
            as2<-1-as1
            er<-rbind(er,(1-(as1-as2)*(as1-as2))*eeff[i]*eeff[i])
          }
          v0<-(1/(nrow(ex)-1))*(t(phe-ex%*%bbbb)%*%(phe-ex%*%bbbb))
          
          if(var(phe)>=(sum(er)+v0)){
            her<-(er/as.vector(var(phe)))*100 
          }else{
            her<-(er/as.numeric(sum(er)+v0))*100
          }
        }
        
        vee<-round(v0,4)
        pee<-round(var(y),4)
        
        if(nrow(her)>1){
          vees<-matrix("",nrow = nrow(her),1)
          pees<-matrix("",nrow = nrow(her),1)
          pees[1,1]<-pee
          vees[1,1]<-vee
        }else{
          pees<-as.matrix(pee)
          vees<-as.matrix(vee)
        }
        
        X<-t(gen)
        XX1<-X
        x<-XX1[3:nrow(XX1),]
        X1<-as.matrix(x)
        xxxx<-as.matrix(X1[,ww])
        
        xxmaf<-t(xxxx)
 
        maf.fun<-function(snp){
          leng<-length(snp)
          snp1<-length(which(snp==1))
          snp11<-length(which(snp==-1))
          snp0<-length(which(snp==0))
          ma1<-(2*snp1+snp0)/(2*leng)
          ma2<-(2*snp11+snp0)/(2*leng)
          maf<-min(ma1,ma2)
          return(maf)
        }
        
        maf<-apply(xxmaf,1,maf.fun)
        maf<-as.matrix(round(maf,4))
        
        eeff[which(abs(eeff)>=1e-4)] <- round(eeff[which(abs(eeff)>=1e-4)],4)
        eeff[which(abs(eeff)<1e-4)] <- as.numeric(sprintf("%.4e",eeff[which(abs(eeff)<1e-4)]))
        lo[which(abs(lo)>=1e-4)] <- round(lo[which(abs(lo)>=1e-4)],4)
        lo[which(abs(lo)<1e-4)] <- as.numeric(sprintf("%.4e",lo[which(abs(lo)<1e-4)]))
        her[which(abs(her)>=1e-4)] <- round(her[which(abs(her)>=1e-4)],4)
        her[which(abs(her)<1e-4)] <- as.numeric(sprintf("%.4e",her[which(abs(her)<1e-4)]))
        log10P <- as.matrix(-log10(1-pchisq(lo*4.605,1)))
        
        if (length(ww)>1){
          wan<-data.frame(parmsShow[needww,1],chr_pos[ww,],eeff,lo,log10P,her,maf,parmsShow[needww,11],vees,pees)
        }else{
          
          wan<-data.frame(parmsShow[needww,1],t(as.matrix(chr_pos[ww,])),eeff,lo,log10P,her,maf,parmsShow[needww,11],vees,pees)  
        }
        colnames(wan)<-c("RS#","Chromosome","Marker Position (bp)","QTN effect","LOD score","-log10(P)","r2 (%)","MAF","Genotype  for code 1","Var_Error","Var_phen (total)")
      }
    }
    parmsShow<-parmsShow[,-c(4,5,6,8,9,12)]
    
    output<-list(result1=parmsShow,result2=wan,Manhattan=parms,QQ=parms.pchange)
    return(output) 
  }
}