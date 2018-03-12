ISIS<-function(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,CLO){
  
  pvalue<-svpal
  svlod<-svmlod
  inputform<-Genformat
  
  if(is.null(psmatrix)){
    flagps<-1
  }else{ 
    flagps<-0
  }
  
  y<-phe
  ps<-psmatrix
  
  
  if(is.null(svpal)==TRUE||is.null(svmlod)==TRUE){
    warning("Please set parameters!")
  }
  
  if((svpal<0)||(svpal>1))
  {
    warning("Please input critical P-value between 0 and 1!")
  }
  if(svmlod<0)
  {
    warning("Please input critical LOD score: > 0 !")
  }
  if(is.null(gen)==TRUE)
  {
    warning("Please input correct genotypic dataset !","Warning",icon="warning")
    
  }
  if(is.null(y)==TRUE)
  {
    warning("Please input correct phenotypic dataset !","Warning",icon="warning")
    
  }
  
  if((is.null(gen)==FALSE)&&(is.null(y)==FALSE)&&(ncol(gen)!=(nrow(y)+2)))
  {
    warning("Sample size in genotypic dataset doesn't equal to the sample size in phenotypic dataset!","Error",icon="error")
    
  }
  
  if((is.null(gen)==FALSE)&&(is.null(y)==FALSE)&&((ncol(gen)==(nrow(y)+2)))&&(svpal>=0)&&(svpal<=1)&&(svmlod>=0))
  {
    
    wan<-NULL
    result<-NULL
    X<-t(gen)
    
    set.seed(1)
    X1<-X[3:nrow(X),]
    sig<-seq(1:ncol(X1))
    x<-data.frame(X1)
    y<-as.matrix(y)
    le<-length(sig)
    xnew<-x[sig]
    pval<-pvalue
    y1<-matrix(nrow=nrow(y),ncol=ncol(y))
    
    if (is.null(ps)==FALSE)
    {
      ps1<-cbind(matrix(1,nrow=nrow(y)),ps)
      vhat<-solve(crossprod(ps1,ps1))%*%crossprod(ps1,y)
      vhat1<-vhat[-1]
      y1<-y-ps%*%vhat1
    }else{
      y1<-y 
    }
    y1<-as.matrix(y1)
    
    xxx<-xnew
    mat<-vector()
    matcor<-vector()
    
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
    
    
    mm=foreach(i=1:le, .multicombine=TRUE, .combine = 'rbind')%dopar%
    {
      if (var(xxx[,i])>0){
        mat[i]<-cor.test(xxx[,i],y1)$p.value
        matcor[i]<-abs(cor(xxx[,i],y1))
        m<-c(mat[i],matcor[i])
      }else{
        mat[i]<-1
        matcor[i]<-0
        m<-c(mat[i],matcor[i])
      }
    }
    rownames(mm)<-NULL
    mat<-mm[,1];matcor<-mm[,2]
    
    stopCluster(cl)
    
    
    if(length(which(mat<pval))<=nrow(y)){
      ee<-as.vector(which(mat<pval))
    }else{
      n1<-nrow(y1)-1
      ee<-as.vector(which(rank(matcor)>=(le-n1), arr.ind=T))
    }
    
    xxxnew<-X1[,sig[ee]]
    yyy<-y1
    
    cvfit1 <- ncvreg::cv.ncvreg(scale(xxxnew), yyy, family="gaussian",penalty="SCAD",gamma=3.7,warn=FALSE)
    
    
    fit1 <- cvfit1$fit
    obj11 <- (as.vector(fit1$beta[,cvfit1$min]))
    obj1<-obj11[-1]
    if (length(which(abs(obj1)!=0))!=0)
    {
      sig1a<-which(abs(obj1)!=0)
      sig1b<-sig[-(sig[ee][sig1a])]
      yyy1<-y1-(X1[,sig[ee][sig1a]]%*%as.matrix(obj1[sig1a]))
      xxx1<-X1[,sig1b]
      mat1<-vector()
      for (i in 1:length(sig1b))
      {
        if (var(xxx1[,i])>0){
          
          mat1[i]<-abs(cor(xxx1[,i],yyy1))
        }else{
          mat1[i]<-0
          
        }
      }
      
      n2<-nrow(yyy1)-1
      ee1<-as.vector(which(rank(mat1)>=(ncol(xxx1)-n2), arr.ind=T))
      
      xxxnew1<-X1[,sig1b[ee1]]
      cvfit2 <- ncvreg::cv.ncvreg(scale(xxxnew1), yyy1, family="gaussian",penalty="SCAD",gamma=3.7,warn=FALSE)
      
      
      fit2 <- cvfit2$fit
      obj22 <- (as.vector(fit2$beta[,cvfit2$min]))
      obj2<-obj22[-1]
      sig1c<-sig1b[ee1][which(abs(obj2)!=0)]
      sigg<-sort(c(sig[ee][sig1a],sig1c))
      
    }else{
      
      sigg<-sig[ee]
    }
    le1<-length(sigg)
    
    if(le1!=0){
      
      ###########if result just have one column##############modified 2017.3.22###############
      if(le1==1){
        xxxnew11<-matrix(X1[,sigg],,1)
      }else{
        xxxnew11<-X1[,sigg]
      }
	  
      z<-matrix()
      
      if (is.null(ps)==TRUE)
      {
        z<-matrix(1,nrow(X1),1)
      }else{
        z<-cbind(matrix(1,nrow(X1),1),ps) 
      }
      
      u1<-ebayes_EM(z,xxxnew11,y)
      obj3<-u1$u 
      result1<-matrix(0,ncol(X1)*1,ncol=1,nrow=ncol(X1))
      for (i in 1: le1)
      {
        result1[(sigg)[i],1]=obj3[i]
      }
      Res<- t(as.matrix((rowSums(result1)/ncol(result1))))
      Res1<-as.vector(Res)	
      le2<-length(which(abs(Res1)>1e-5))
      
      if(le2!=0){
        
        sig1<-which(abs(Res1)>1e-5)
        bbo<-matrix(0,le2,1)
        for (i in 1:le2){
          bbo[i,]=Res1[sig1[i]]
        }
        
        her1<-vector(length=le2)
        for (i in 1:le2){
          p1<-length(as.vector(which(X1[,sig1[i]]==1)))/length(X1[,sig1[i]])
          p2<-1-p1
          her1[i]=((p1+p2)-(p1-p2)^2)*(Res1[sig1[i]])^2
        }
        
        if(var(y)>=sum(her1)+u1$sigma2){
          her<-(her1/as.vector(var(y)))*100  
          
        }else{
          her<-(her1/(sum(her1)+u1$sigma2))*100 
        } 
        
        if(length(sig1)!= 0){
          
          if(length(sig1)==1){
            xxxx<-as.matrix(X1[,sig1])
            
          }else{
            xxxx<-X1[,sig1]
          }
          
          yn<-as.matrix(y)
          xxn<-z
          
          lod<-likelihood(xxn,xxxx,yn,bbo)
          slod<-cbind(sig1,lod,her)
          
          if(length(which(slod[,2]>=svlod))>=1){
            
            if(length(which(slod[,2]>=svlod))==1){
              sslod<-t(as.matrix(slod[which(slod[,2]>=svlod),]))
              sig1<-slod[which(slod[,2]>=svlod),1]
            }else if(length(which(slod[,2]>=svlod))>1){
              sslod<-slod[which(slod[,2]>=svlod),]
              sig1<-sslod[,1]
            }
            xxxx<-as.matrix(X1[,sig1])
            lod<-sslod[,2]  
            her<-sslod[,3]
            
            ii<-as.vector(sig1)
            qqq<-matrix(0,nrow=length(ii),ncol=6)
            qqq[,1]=as.matrix(ii)
            for (j in 1:length(ii)){
              qqq[j,2]=X[1,ii[j]]
              qqq[j,3]=X[2,ii[j]]
              qqq[j,4]=result1[ii[j],]
              
              qqq[j,5]=lod[j]
              qqq[j,6]=her[j]
            }
            id<-which(qqq[,5]==0)
            
            if(length(id)!=dim(qqq)[1]){
              
              if(length(id)!=0){
                qqq1<-qqq[-id,]
              }else{
                qqq1<-qqq
              }
              #######revised 2017 3.4##############################
              if(length(sig1)==1){
                xxmaf<-t(xxxx)
                xxmaf<-matrix(xxmaf,1,)
                result<-matrix(qqq1[,-1],1,)
              }else{
                xxmaf<-t(xxxx)
                result<-as.matrix(qqq1[,-1])
              }
              
              leng.maf<-dim(xxmaf)[2]
             
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
              vee<-round(u1$sigma2,4)
              pee<-round(var(y),4)
              
              if(nrow(qqq1)>1){
                
                vees<-matrix("",nrow = nrow(result),1)
                pees<-matrix("",nrow = nrow(result),1)
                pees[1,1]<-pee
                vees[1,1]<-vee
                result<-as.matrix(qqq1[,-1])
                result<-result
                temp<-as.matrix(result[,3:5])
                temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
                temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
                wan<-cbind(result[,1:2],temp)
                
              }else{
                pees<-as.matrix(pee)
                vees<-as.matrix(vee) 
                result<-t(as.matrix(qqq1[,-1]))
                result<-result
                temp<-t(as.matrix(result[,3:5]))
                temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
                temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
                wan<-cbind(t(as.matrix(result[,1:2])),temp)
                
              }
              
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
              
              wan<-cbind(marker,wan,maf,snp,vees,pees)
              tempwan <- wan
              lodscore1 <- as.numeric(tempwan[,5])
              log10P <- as.matrix(round(-log10(1-pchisq(lodscore1*4.605,1)),4))
              
              if(nrow(tempwan)>1){
                tempwan1 <- cbind(tempwan[,1:5],log10P,tempwan[,6:10])
              }else{
                tempwan1 <- cbind(t(as.matrix(tempwan[,1:5])),log10P,t(as.matrix(tempwan[,6:10])))  
              }
              
              wan <- tempwan1
              
              colnames(wan)<-c("RS#","Chromosome","Marker Position (bp)","QTN effect","LOD score","-log10(P)","r2 (%)","MAF","Genotype  for code 1","Var_Error","Var_phen (total)")
            }
          }
        }
      }  
    }
    
    if(nrow(result)>1){
      r1<-as.matrix(result[,c(1,2,4)]) 
    }else{
      r1<-t(as.matrix(result[,c(1,2,4)]))  
    }
    r2<-as.matrix(gen[,1:2])
    
    rowbl<-nrow(r2)-nrow(r1)
    bl<-matrix("",rowbl,3)
    r12<-rbind(r1,bl)
    result<-cbind(r2,r12)
    
    colnames(result)<-c("Chromosome","Marker Position (bp)","Chromosome(detected)","Marker Position (bp)(detected)","LOD score(detected)")
    
    output<-list(result=wan,plot=result)
    
    return(output)
    
  }
}