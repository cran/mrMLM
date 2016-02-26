mrMLM<-function(){
  mrenv <- new.env()
  gnewtable<-function (items, multiple = FALSE, chosencol = 1, icon.FUN = NULL, 
                       filter.column = NULL, filter.labels = NULL, filter.FUN = NULL, 
                       handler = NULL, action = NULL, container = NULL, ..., toolkit = guiToolkit()) 
  {
    if (!missing(items)) {
      if (is.vector(items)) 
        items <- data.frame(.= items, stringsAsFactors = FALSE)
      if (is.matrix(items)) 
        items <- data.frame(items, stringsAsFactors = FALSE)
    }
    widget <- .gtable(toolkit, items = items, multiple = multiple, 
                      chosencol = chosencol, icon.FUN = icon.FUN, filter.column = filter.column, 
                      filter.labels = filter.labels, filter.FUN = filter.FUN, 
                      handler = handler, action = action, container = container, 
                      ...)
    obj <- new("gTable", widget = widget, toolkit = toolkit)
    return(obj)
  }
  
  window<-gwindow(title="Multi-locus Random-SNP-effect Mixed Linear Model (mrMLM)",visible=TRUE,width=1240,height=730,expand=TRUE)
  plotwin<-gwindow("Manhattan Plot",visible=FALSE,width=600,height=360)
  gpw<-ggroup(container=plotwin)
  ggpw<-ggraphics(container=gpw)
  plotwin1<-gwindow("Q-Q Plot",visible=FALSE,width=600,height=360)
  gpw1<-ggroup(container=plotwin1)
  ggpw1<-ggraphics(container=gpw1)
  choicekk<-gwindow("Choice kinship",visible=FALSE,width=250,height=150)
  gkk<-ggroup(container=choicekk,expand=FALSE)
  includeps<-gwindow("Include population structure?",visible=FALSE,width=300,height=150)
  gps<-ggroup(container=includeps,expand=FALSE)
  
  lyt<-glayout(container=window,spacing=13)
  
  genotype<-gbutton("Genotype",container=lyt)
  phenotype<-gbutton("Phenotype",container=lyt)
  kinship<-gbutton("Kinship",container=lyt)
  population<-gbutton("Population Structure",container=lyt)
  manhattan<-gbutton("Manhattan Plot",container=lyt)
  qqplot<-gbutton("QQ Plot",container=lyt)
  
  savefile<-gbutton(" Save ",container=lyt)
  run<-gbutton("Run",container=lyt)
  exit<-gbutton("Exit",container=lyt)
  gwline<-glabel("Critical value for -logP",container=lyt)
  gwedit<-gedit("3",width=20,coerce.with=as.numeric,container=lyt)
  svgwline<-svalue(gwedit)
  gwstandp<-glabel("Critical P-value for QQ plot",container=lyt)
  gwedit1<-gedit("0.992",width=20,coerce.with=as.numeric,container=lyt)
  svgwstandp<-svalue(gwedit1)
  
  lyt[1,1]<-genotype
  lyt[2,1]<-phenotype
  lyt[3,1]<-kinship
  lyt[4,1]<-population
  lyt[7,1]<-run
  lyt[8,1]<-savefile
  lyt[11,1]<-gwline
  lyt[12,1]<-gwedit
  lyt[13,1]<-manhattan
  lyt[14,1]<-gwstandp
  lyt[15,1]<-gwedit1
  lyt[16,1]<-qqplot
  lyt[19,1]<-exit
  
  nb1<-gnotebook(tab.pos=3,closebuttons=TRUE,dontCloseThese=TRUE,container=lyt,expand=TRUE)
  size(nb1)<-c(680,540)
  tb<-gnewtable("     
                1. mrMLM is a R software package for genome-wide association studies based on a multi-locus random-SNP-effect mixed linear model.
                
                2. Please cite: Wang Shi-Bo, Feng Jian-Ying, Ren Wen-Long, Huang Bo, Zhou Ling, Wen Yang-Jun, Zhang Jin, Jim M. Dunwell, Xu Shizhong (*), Zhang Yuan-Ming (*). 2016. 
                Improving power and accuracy of genome-wide association studies via a multi-locus mixed linear model methodology.Scientific Reports 6: 19444. 
                
                3. The software package is developed by Wen-Long Ren, Shi-Bo Wang & Yuan-Ming Zhang.
                
                
                Version 1.0, Realeased January 2016",multiple=TRUE,container=nb1,expand=TRUE,label="About the software")
  
  
  lyt[1:20,2,expand=TRUE]<-nb1
  
  staprogress<-gtkButton()
  lyt[21,2,expand=TRUE]<-staprogress
  
  
  addHandlerClicked(genotype,handler=function(h,...){
    mrenv$flagps<-1
    input1<-gfile(text="Select a file...",type="open",
                  filter=list("All files"=list(patterns=c("*")),
                              "CSV files"=list(patterns=c("*.csv"))))
    
    if(is.na(input1))
    {
      gmessage("Please input correct genotype data !","Warning",icon="warning")
      return
    }else{
      mrenv$gen<-as.matrix(read.csv(input1,header=F))
      tbdfe1<-gdfedit(mrenv$gen,container=nb1,expand=TRUE,label="Genotype") 
    }
  })
  
  
  addHandlerClicked(kinship,handler=function(h,...){
    if(isExtant(choicekk)==FALSE)
    {
      choicekk<-gwindow("Choice kinship",visible=FALSE,width=250,height=150)
      gkk<-ggroup(container=choicekk,expand=FALSE)
    }
    lytkk<-glayout(container=gkk,spacing=13)
    mrenv$okkk<-gbutton("     OK    ",container=lytkk)
    cancelkk<-gbutton(" Cancel ",container=lytkk)
    mrenv$radiokk<-gradio(c("Input kinship file directly","Compute by this program"),selected=1,horizontal=FALSE,container=lytkk)
    lytkk[2:3,2:5]<-mrenv$radiokk
    lytkk[5,3]<-mrenv$okkk
    lytkk[5,5]<-cancelkk
    visible(choicekk)<-TRUE
    addHandlerClicked(mrenv$okkk,handler=function(h,...){
      if(svalue(mrenv$radiokk)=="Input kinship file directly"){
        input2<-gfile(text="Select a file...",type="open",
                      filter=list("All files"=list(patterns=c("*")),
                                  "CSV files"=list(patterns=c("*.csv"))))
        if(is.na(input2))
        {
          gmessage("Please input correct kinship data !","Warning",icon="warning")
          return
        }else{
          mrenv$kk<-as.matrix(read.csv(input2,header=F)) 
          tbdfe2<-gdfedit(mrenv$kk,container=nb1,expand=TRUE,label="Kinship")
          dispose(choicekk)
        }
      }else{
        envgen <- mrenv$gen
        if(exists("envgen")==FALSE)
        {
          gmessage("Please input correct genotype data !","Warning",icon="warning")
          return
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
          mrenv$kk<-kk1
          tbdfe2<-gdfedit(mrenv$kk,container=nb1,expand=TRUE,label="Kinship")
          dispose(choicekk)
        }     
      }
    })
    addHandlerClicked(cancelkk,handler=function(h,...){
      dispose(choicekk)
    })
  })
  
  addHandlerClicked(phenotype,handler=function(h,...){
    input3<-gfile(text="Select a file...",type="open",
                  filter=list("All files"=list(patterns=c("*")),
                              "CSV files"=list(patterns=c("*.csv"))))
    if(is.na(input3))
    {
      gmessage("Please input correct phenotype data !","Warning",icon="warning")
      return
    }else{
      mrenv$phe<-as.matrix(read.csv(input3,header=F)) 
      tbdfe3<-gdfedit(mrenv$phe,container=nb1,expand=TRUE,label="Phenotype")
    }
  })
  
  addHandlerClicked(population,handler=function(h,...){
    if(isExtant(includeps)==FALSE)
    {
      includeps<-gwindow("Include population structure?",visible=FALSE,width=300,height=150)
      gps<-ggroup(container=includeps,expand=FALSE)
    }
    lytps<-glayout(container=gps,spacing=13)
    okps<-gbutton("     OK    ",container=lytps)
    cancelps<-gbutton(" Cancel ",container=lytps)
    radiops<-gradio(c("Do not need population structure","Input population structure file directly"),selected=1,horizontal=FALSE,container=lytps)
    lytps[2:3,2:5]<-radiops
    lytps[5,3]<-okps
    lytps[5,5]<-cancelps
    visible(includeps)<-TRUE
    addHandlerClicked(okps,handler=function(h,...){
      if(svalue(radiops)=="Input population structure file directly"){
        mrenv$flagps<-0
        input4<-gfile(text="Select a file...",type="open",
                      filter=list("All files"=list(patterns=c("*")),
                                  "CSV files"=list(patterns=c("*.csv"))))
        if(is.na(input4))
        {
          gmessage("Please input correct population data !","Warning",icon="warning")
          return
        }else{
          mrenv$psmatrix<-as.matrix(read.csv(input4,header=F)) 
          tbdfe4<-gdfedit(mrenv$psmatrix,container=nb1,expand=TRUE,label="Population Structure")
          dispose(includeps)
        }
      }else{
        mrenv$flagps<-1
        enabled(population)<-FALSE
        dispose(includeps)
      }
    })
    addHandlerClicked(cancelps,handler=function(h,...){
      dispose(includeps)
    })
  })
  
  addHandlerClicked(exit,handler=function(h,...){
    gconfirm("Yes or no?",handler=function(h,...){dispose(window)})
  })
  
  addHandlerClicked(run,handler=function(h,...){
    gen<-mrenv$gen
    phe<-mrenv$phe
    kk<-mrenv$kk
    flagps<-mrenv$flagps
    psmatrix<-mrenv$psmatrix
    if(exists("gen")==FALSE)
    {
      gmessage("Please input correct genotype data !","Warning",icon="warning")
      return
    }
    if(exists("phe")==FALSE)
    {
      gmessage("Please input correct phenotype data !","Warning",icon="warning")
      return
    }
    if(exists("kk")==FALSE)
    {
      gmessage("Please input correct kinship data !","Warning",icon="warning")
      return
    }
    if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(ncol(gen)!=(nrow(phe)+2)))
    {
      gmessage("Sample size between genotype and phenotype is inconsistent!","Error",icon="error")
      return
    }
    
    if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2))))
    {
      progress_bar <- gtkProgressBar ( )
      staprogress$add(progress_bar)
      progress_bar$setText ( "Please be patient ..." )
      progress_bar$setFraction(2/100)
      
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
      {yy<-phe[,1]
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
       
       qq<-numeric()
       for(k in 1:m){
         progress_bar$setFraction((2+(90/m)*k)/100)
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
         if((flagps==1)||(exists("psmatrix")==FALSE))
         {
           beta<-parmfix[[3]]
         }else if(flagps==0)
         {
           beta<-parmfix[[3]][1]
         }
         sigma2<-parmfix[[4]]
         lambda<-xi
         sigma2g<-lambda*sigma2
         fn0<-loglike(-Inf)
         lrt<-2*(fn0-fn1)
         p_lrt<-1-pchisq(lrt,1)
         wald<-(gamma/stderr)^2
         p_wald<-1-pchisq(wald,1)
         parm0<-c(ii,name[k,1],name[k,2],beta,sigma2,sigma2g,gamma,stderr,wald,p_wald)
         qq<-rbind(qq,parm0)
       }
       ll<-rbind(ll,qq)
      }
      mrenv$parms<-ll
      mrenv$parms<-matrix(mrenv$parms,,10)
      
      multinormal<-function(y,mean,sigma)
      {
        pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
        return (pdf_value)
      }
      #LOD value test
      likelihood<-function(xxn,xxx,yn)
      {
        nq<-ncol(xxx)
        ns<-nrow(yn)
        at1<-0
        ww1<-1:ncol(xxx)
        ww1<-as.matrix(ww1)
        at1<-dim(ww1)[1]
        lod<-matrix(rep(0,nq),nq,1)
        if(at1>0.5)
          ad<-cbind(xxn,xxx[,ww1])
        else
          ad<-xxn
        #if(abs(det(crossprod(ad,ad)))<1e-6)
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
            ij<-which(sub!=sub[i+1])
            ad1<-ad[,ij]
            #if(abs(det(crossprod(ad1,ad1)))<1e-6)
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
      
      #2010 EM_Bayes
      ebayes_EM<-function(x,z,y)
      {
        n<-nrow(z);k<-ncol(z)
        if(abs(min(eigen(crossprod(x,x))$values))<1e-6)
          b<-solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y)
        else
          b<-solve(crossprod(x,x))%*%crossprod(x,y)
        v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
        u<-matrix(rep(0,k),k,1)
        v<-matrix(rep(0,k),k,1)
        s<-matrix(rep(0,k),k,1)
        for(i in 1:k)
        {
          zz<-z[,i]
          s[i]<-((crossprod(zz,zz))^(-1))*v0
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
        iter<-0;err<-1000;iter_max<-100;err_max<-1e-8
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
          }else
          {
            if(abs(min(eigen(xtv%*%x)$values))<1e-6){
              b<-solve((xtv%*%x)+diag(ncol(x))*1e-8)%*%(xtv%*%y)
            }
            else{
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
          wang[i]<- p
        }
        return (wang)
      }
      
      gen <- t(gen)
      cccc<-mrenv$parms[,2:3]
      h0<-which(mrenv$parms[,10]<=0.01)
      h0<-as.matrix(h0)
      hh<-nrow(h0)
      aa<-cbind((1:(nrow(mrenv$parms))),mrenv$parms[,10])
      aa0<-order(aa[,2])
      aa<-aa[aa0,]
      cc<-aa[1:hh,]
      name<-aa[1:hh,1]
      name<-as.matrix(name)
      w0<-cbind(mrenv$parms[,6],mrenv$parms[,8])
      ww0<-matrix(1,(nrow(w0)),1)-w0[,2]*w0[,2]/w0[,1]
      k0<-which(ww0<0)
      k0<-as.matrix(k0)
      if (nrow(k0)>0){ww0[k0,1]<-matrix(0,(nrow(k0)),1)}
      nn<-sum(ww0)
      pp<-0.05/nn
      aa0<-which(aa[,2]<=pp)
      aa0<-as.matrix(aa0)
      a0<-nrow(aa0)
      gg<-name
      for (ii in 1:(nrow(name)-1)){
        for (jj in (ii+1):(nrow(name))){
          ci<- cccc[name[ii],1]
          cj<- cccc[name[jj],1]
          if (ci==cj){
            ye<-abs(cccc[name[ii],2]-cccc[name[jj],2])
            if (ye<=20000){gg[jj,1]=0}
          }
        }
      }
      progress_bar$setFraction(95/100)
      gg<-as.matrix(gg)
      misfit<-numeric()
      kk<- numeric
      kk0<- numeric
      l0<- numeric
      bong<-a0
      if (bong>0){
        g0<-gg[1:a0,1]
        g0<-as.matrix(g0)
        kk0<-a0
        a0<-which(g0>0)
        a0<-as.matrix(a0)
        g0<-g0[a0,1]
        g0<-as.matrix(g0)
        xxx0<-t(gen[g0,])
        xxx0<-as.matrix(xxx0)
        phe<-as.matrix(phe)
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          par<-likelihood(matrix(1,(nrow(xxx0)),1),xxx0,phe)
          lod<-par
        }else if(flagps==0)
        {
          temp<-cbind(matrix(1,(nrow(xxx0)),1),psmatrix)
          par<-likelihood(temp,xxx0,phe)
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
      ne<-as.matrix(gg[(kk0+1):(nrow(gg)),1])
      if ((length(misfit))>0){gg<-rbind(ne,misfit)}
      if ((length(misfit))==0){gg<-ne}
      a1<-which(gg>0)
      a1<-as.matrix(a1)
      a2<-gg[a1,1]
      a2<-as.matrix(a2)
      xx<-t(gen[a2,])
      xx<-as.matrix(xx)
      if((flagps==1)||(exists("psmatrix")==FALSE))
      {
        if (length(kk)>1){xin<-cbind(matrix(1,(nrow(xx)),1),xx0)}
        if (length(kk)==1){xin<- matrix(1,(nrow(xx)),1)}
      }else if(flagps==0)
      {
        temp<-cbind(matrix(1,(nrow(xx)),1),psmatrix)
        if (length(kk)>1){xin<-cbind(temp,xx0)}
        if (length(kk)==1){xin<-temp}
      }
      xin<-as.matrix(xin)
      par<-ebayes_EM(xin,xx,phe)
      w2<-which(par[,1]<=0.01)
      w2<-as.matrix(w2)
      ww<- numeric
      if ((nrow(w2))>0){
        name<-a2[w2,1]
        name<-as.matrix(name)
        x3<-cbind(xin,xx[,w2])
        x3<-as.matrix(x3)
        lodfix<-matrix(x3[,1],nrow(x3),)
        lodrand<-matrix(x3[,2:(ncol(x3))],nrow(x3),)
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          lod<-likelihood(lodfix,lodrand,phe)
        }else if(flagps==0)
        {
          temp<-cbind(psmatrix,lodfix)
          lod<-likelihood(temp,lodrand,phe)
        }
        w3<-which(lod[,1]>=3)
        w3<-as.matrix(w3)
        if ((kk[1])>0){
          g0<-as.matrix(g0)
          name<-rbind(g0,name)
          name<-as.matrix(name)
          if ((w3[1])>0){
            lo<-lod[w3,1]
            tpww<-which(w3<dim(name)[1])
            w3<-as.matrix(w3[tpww])
            ww<-name[w3,]
          }
          if ((nrow(w3))==0){ww<-0}
        }
        if (length(kk)==1){ww<-0}
      }
      if ((nrow(w2))==0){
        g0<-as.matrix(g0)
        lo<-as.matrix(lo)
        yang<-which(lo>=3)
        yang<-as.matrix(yang)
        if ((nrow(yang))>0){
          ww<-g0[yang,1]
          lo<-lo[yang,1]
        }
        if ((nrow(yang))==0){ww<-0}
      }
      ww<-as.matrix(ww)
      if (length(ww)>1){
        ex<-cbind(matrix(1,(nrow(xx)),1),t(gen[ww,]))
        ex<-as.matrix(ex)
        cui<-det(t(ex)%*%ex)
        p1<-rep(1,ncol(ex))
        p2<-diag(p1)
        if (cui<1e-6){bbbb<-solve(t(ex)%*%ex+p2*0.01)%*%t(ex)%*%phe}
        if (cui>=1e-6){ bbbb<-solve(t(ex)%*%ex)%*%t(ex)%*%phe }
        eeff<-bbbb[2:(nrow(bbbb)),1]
        mrenv$wan<-cbind(cccc[ww,],eeff,lo)
        colnames(mrenv$wan)<-c("Chromosome","Position","Effect","Lod")
      }
      wan<-mrenv$wan
      if(exists("wan")==FALSE||is.null(wan)==TRUE)
      {
        gmessage("There is no result meets the requirements !","Info",icon="info")
      }else{
        tbdfe4<-gdfedit(wan,container=nb1,expand=TRUE,label="Result")
      }
      progress_bar$setFraction(100/100)
      progress_bar$setText("All done.")
    }
    return
  })
  
  addHandlerClicked(savefile,handler=function(h,...){
    wan<-mrenv$wan
    if(exists("wan")==FALSE||is.null(wan)==TRUE)
    {
      gmessage("There is no result meets the requirements !","Info",icon="info")
    }else{
      output<-gfile(text="Save a file...",type="save",
                    filter=list("All files"=list(patterns=c("*")),
                                "CSV files"=list(patterns=c("*.csv"))))
      write.table(wan,output,sep = ",",row.names=FALSE,col.names = FALSE) 
    }
  })
  
  addHandlerClicked(manhattan,handler=function(h,...){
    if((exists("svgwline")==FALSE)||(svalue(gwedit)<=0))
    {
      gmessage("Please input correct genomewideline value!","Warning",icon="warning")
      return
    }else{
      svgwline<-svalue(gwedit)
      mrenv$standline<-svgwline
      if(isExtant(plotwin)==FALSE)
      {
        plotwin<-gwindow("Manhattan Plot",visible=FALSE,width=600,height=360)
        gpw<-ggroup(container=plotwin)
        ggpw<-ggraphics(container=gpw)
      }
      addHandlerChanged(ggpw, handler=function(h,...) {
        parms<-as.data.frame(mrenv$parms)
        plotman<-manhattan(parms,chr = "V2",bp ="V3",p ="V10",snp="v5",suggestiveline=FALSE,genomewideline = mrenv$standline)
      })
      visible(plotwin)<-TRUE
    }
  })
  
  addHandlerClicked(qqplot,handler=function(h,...){
    if((exists("svgwstandp")==FALSE)||(svalue(gwedit1)<=0))
    {
      gmessage("Please input correct standard P-value!","Warning",icon="warning")
      return
    }else{
      svgwstandp<-svalue(gwedit1)
      mrenv$standp<-svgwstandp
      if(isExtant(plotwin1)==FALSE)
      {
        plotwin1<-gwindow("Q-Q Plot",visible=FALSE,width=600,height=360)
        gpw1<-ggroup(container=plotwin1)
        ggpw1<-ggraphics(container=gpw1)
      }
      addHandlerChanged(ggpw1, handler=function(h,...) {
        pvalue<-matrix(mrenv$parms[,10],,1)
        observed<-sort(pvalue[,1])
        newobserved<-observed[which(observed<mrenv$standp)]
        lobs<--(log10(newobserved))
        expected<-c(1:length(newobserved))
        lexp<--(log10(expected/(length(expected)+1)))
        plot(lexp,lobs,xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'))
        abline(0,1,col="red")
      })
      visible(plotwin1)<-TRUE
    }
  })
}