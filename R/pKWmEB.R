pKWmEB<-function(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,CLO){
  
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
  
  
  if((flagps==1)||(exists("psmatrix")==FALSE))
  {
    phe<-phe
  }else if(flagps==0)
  {
    phe<-phe
    fixps <- cbind(matrix(1,nrow(phe),1),psmatrix)
    cui<-det(t(fixps)%*%fixps)
    p1<-rep(1,ncol(fixps))
    p2<-diag(p1)
    if (cui<1e-6){bbps<-solve(t(fixps)%*%fixps+p2*0.01)%*%t(fixps)%*%phe}
    if (cui>=1e-6){ bbps<-solve(t(fixps)%*%fixps)%*%t(fixps)%*%phe }
    bbps <- bbps[2:(nrow(bbps)),1]
    phe <- as.matrix(phe) - as.matrix(psmatrix)%*%as.matrix(bbps)
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
  
  if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2)))&&(svpal>=0)&&(svpal<=1)&&(svmlod>=0))
  {
    
    parmsShow<-NULL
    wan<-NULL
    parms.pchange<-NULL
    parmsm<-NULL
    
    K.data <- kk
    Y.data <- phe
    rawgen <- gen
    rawphe <- Y.data
    
    gene.data <- rawgen[,3:dim(rawgen)[2]]
    nsample <- dim(gene.data)[2]
    fix <- matrix(1,nsample,1)
    sam <- nsample
    
    Y.data <- matrix(Y.data,nsample,1)
    
    n<-dim(Y.data)[1]
    W.orig<-matrix(1,n,1)
    
    W <- W.orig
    K <- K.data
    YY <- Y.data
    
    p_value <- svpal
    
    ffpptotal <- numeric()
    gglartotal <- numeric()
    pvaluetotal <- numeric()
    for(ii in 1:1){
      remle2<-emma.REMLE(YY[,ii], W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
      
      remle1.B1<-emma.maineffects.B(Z=NULL,K,remle2$delta)
      C2<-remle1.B1$mC
      
      Y_c <- C2%*%YY[,ii]
      W_c <- C2%*%W
      G_c <- C2%*%t(gene.data)
      
      GGG <- t(G_c)
      
      allrowmean <- rowMeans(GGG)
      nnG <- nrow(GGG)
      
      for(jj in 1:nnG)
      {
        GGG[jj,which(GGG[jj,]>=allrowmean[jj])] <- 1
        GGG[jj,which(GGG[jj,]<allrowmean[jj])] <- -1
      }
      
      gentran <- GGG
      phetran <- Y_c
      nn <- dim(gentran)[1]
      bb<-numeric()
      cc <- numeric()
      ff <- numeric()
      
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
      
      bb=foreach(i=1:nn, .combine = 'rbind')%dopar%
        
      {
        requireNamespace("coin")
        requireNamespace("lars")
        newphe <- cbind(matrix(c(1:sam),,1),phetran)
        ph <- unique(newphe[,2])
        newph <- newphe[match(ph,newphe[,2],0L),]
        newy <- newph[,2]
        sob <- newph[,1]
        temp <- as.matrix(gentran[i,sob])
        xy <- cbind(temp,newy)
        b <- unique(xy[,1])
        
        temp <- factor(temp)
        snp <- data.frame(newy,temp)
        
        kw <- kruskal_test(newy~temp, data = snp,distribution = "asymptotic")
        kw <- pvalue(kw)
        aa <- kw[1] 
        
        
      }
      stopCluster(cl)
      kk <- matrix(seq(1:nn),nn,1)
      bb <- matrix(bb,nn,1)
      cc <- cbind(ii,kk,bb)
      pvaluetotal <- cc[,2:3]
      ff <- cc[which(cc[,3] < p_value),]
      ffpptotal <- ff
      pvaluetotal <- pvaluetotal
      
      ############lars###########################
      gg <- numeric()
      nchoice <- ff[,2]   
      genchoice <- gene.data[nchoice,]
      newpheno <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
      
      aall <- lars(t(genchoice),newpheno,type="lar",use.Gram=FALSE)
      bb2 <- aall$beta[nrow(aall$beta),]
      var <- unlist(aall[[8]])
      
      tempnn <- dim(ff)[1]
      if(tempnn<=150)
      {
        if(tempnn>=nsample)
        {
          tempnn <- nsample - 1
        }else if(tempnn <nsample)
        {
          tempnn <- dim(ff)[1]
        }
        var1 <- var[1:tempnn]
        bb2 <- bb2[abs(var1)]
        gg <- as.matrix(nchoice[abs(var1)])
        ############Empirical Bayes##################
        ggbayes <- numeric()
        optloci <- gg
        optgen <- gene.data[optloci,]
        newphebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
        bbeff <- ebayes_EM(fix,t(optgen),newphebayes)
        lod <- likelihood(fix,t(optgen),newphebayes,bbeff$u)
        optlod <- which(lod>svmlod)
        if(length(optlod)>0){
          locich <- optloci[optlod]
          ggbayes <- cbind(ii,locich,matrix(rawgen[locich,1:2],,2),bbeff$u[optlod],lod[optlod],bbeff$sigma2)
        }
        gglartotal <- ggbayes
        
      }else if((tempnn > 150)&&(nsample > 150))
      {
        if(tempnn>=nsample)
        {
          tempnn <- nsample - 1
        }else if(tempnn <nsample)
        {
          tempnn <- dim(ff)[1]
        }
        var1 <- var[1:tempnn]
        bb2 <- bb2[abs(var1)]
        gg <- as.matrix(nchoice[abs(var1)])
        
        aic <- numeric()
        hhbayes50 <- numeric()
        hhbayes100 <- numeric()
        hhbayes150 <- numeric()
        
        ggbayes <- numeric()
        ggbayes50 <- numeric()
        ggbayes100 <- numeric()
        ggbayes150 <- numeric()
        
        optloci <- gg
        optloci50 <- as.matrix(optloci[1:50])
        optloci100 <- as.matrix(optloci[1:100])
        optloci150 <- as.matrix(optloci[1:150])
        
        ##################choose 50 number variable from lars######################
        optgen50 <- gene.data[optloci50,]
        phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
        bbeff50 <- ebayes_EM(fix,t(optgen50),phebayes)
        lod50 <- likelihood(fix,t(optgen50),phebayes,bbeff50$u)
        
        optlod50 <- which(lod50>svmlod)
        if(length(optlod50)>0){
          locich50 <- optloci50[optlod50]
          ggbayes50 <- cbind(ii,locich50,matrix(rawgen[locich50,1:2],,2),bbeff50$u[optlod50],lod50[optlod50],bbeff50$sigma2)
          hhbayes50 <- rbind(hhbayes50,ggbayes50)  
        }
        
        ##################choose 100 number variable from lars#####################
        optgen100 <- gene.data[optloci100,]
        phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
        bbeff100 <- ebayes_EM(fix,t(optgen100),phebayes)
        lod100 <- likelihood(fix,t(optgen100),phebayes,bbeff100$u)
        
        optlod100 <- which(lod100>svmlod)
        if(length(optlod100)>0){
          locich100 <- optloci100[optlod100]
          ggbayes100 <- cbind(ii,locich100,matrix(rawgen[locich100,1:2],,2),bbeff100$u[optlod100],lod100[optlod100],bbeff100$sigma2)
          hhbayes100 <- rbind(hhbayes100,ggbayes100)  
        }
        
        ##################choose 150 number variable from lars#####################
        optgen150 <- gene.data[optloci150,]
        phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
        bbeff150 <- ebayes_EM(fix,t(optgen150),phebayes)
        lod150 <- likelihood(fix,t(optgen150),phebayes,bbeff150$u)
        
        optlod150 <- which(lod150>svmlod)
        if(length(optlod150)>0){
          locich150 <- optloci150[optlod150]
          ggbayes150 <- cbind(ii,locich150,matrix(rawgen[locich150,1:2],,2),bbeff150$u[optlod150],lod150[optlod150],bbeff150$sigma2)
          hhbayes150 <- rbind(hhbayes150,ggbayes150)  
        }
        
        ####################################AIC#####################################
        if(length(optlod50)==0)
        {
          lmres1 <- lm(phebayes~fix)
          aic1 <- AIC(lmres1)
        }
        if(length(optlod100)==0)
        {
          lmres2 <- lm(phebayes~fix)
          aic2 <- AIC(lmres2)
        }
        if(length(optlod150)==0)
        {
          lmres3 <- lm(phebayes~fix)
          aic3 <- AIC(lmres3)
        }
        
        if(length(optlod50)==1)
        {
          xx1 <- as.matrix(gene.data[ggbayes50[,2],])
          lmres1 <- lm(phebayes~xx1)
          aic1 <- AIC(lmres1)
        }
        if(length(optlod100)==1)
        {
          xx2 <- as.matrix(gene.data[ggbayes100[,2],])
          lmres2 <- lm(phebayes~xx2)
          aic2 <- AIC(lmres2)
        }
        if(length(optlod150)==1)
        {
          xx3 <- as.matrix(gene.data[ggbayes150[,2],])
          lmres3 <- lm(phebayes~xx3)
          aic3 <- AIC(lmres3)
        }
        
        if(length(optlod50)>1)
        {
          xx1 <- t(gene.data[ggbayes50[,2],])
          lmres1 <- lm(phebayes~xx1)
          aic1 <- AIC(lmres1)
        }
        if(length(optlod100)>1)
        {
          xx2 <- t(gene.data[ggbayes100[,2],])
          lmres2 <- lm(phebayes~xx2)
          aic2 <- AIC(lmres2)
        }
        if(length(optlod150)>1)
        {
          xx3 <- t(gene.data[ggbayes150[,2],])
          lmres3 <- lm(phebayes~xx3)
          aic3 <- AIC(lmres3)
        }
        
        aic <- rbind(aic,matrix(c(ii,aic1,aic2,aic3),1,4))
		
        ############################################################################
        if(aic1==min(aic1,aic2,aic3))
        {
          ggbayes <- ggbayes50
        }else if(aic2==min(aic1,aic2,aic3)){
          ggbayes <- ggbayes100
        }else if(aic3==min(aic1,aic2,aic3)){
          ggbayes <- ggbayes150
        }
        gglartotal <- ggbayes
      }
    }
    
    gglartotal <- gglartotal
    
    
    
    if(inputform==1){
      #output result1 using mrMLM numeric format
      parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
      tempparms<-parmsShow
      tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
      tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
      kong<-matrix("",nrow(tempparms),1)
      parmsShow<-data.frame(genRaw[-1,1],gen[,1:2],kong,tempparms,genRaw[-1,4])
      colnames(parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","SNP effect","-log10(P)","Genotype for code 1")
      
    }
    if(inputform==2){
      #output result1 using mrMLM character format
      parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
      outATCG<-matrix(outATCG,,1)
      tempparms<-parmsShow
      tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
      tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
      kong<-matrix("",nrow(tempparms),1)
      parmsShow<-data.frame(genRaw[-1,1],gen[,1:2],kong,tempparms,outATCG)
      colnames(parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","SNP effect","-log10(P)","Genotype  for code 1")
      
    }
    if(inputform==3){
      #output result1 using TASSEL format
      parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
      outATCG<-matrix(outATCG,,1)
      outATCG<-unlist(strsplit(outATCG,""))
      outATCG<-matrix(outATCG[c(TRUE,FALSE)],,1)
      tempparms<-parmsShow
      tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
      tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
      kong<-matrix("",nrow(tempparms),1)
      parmsShow<-data.frame(genRaw[-1,1],gen[,1:2],kong,tempparms,outATCG)
      colnames(parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","SNP effect","-log10(P)","Genotype  for code 1")
    }
    
    parms.pchange<-as.matrix(pvaluetotal[,2])
    parmsp<-parms.pchange
    locsub<-which(parmsp==0)
    if(length(locsub)!=0){
      pmin<-min(parmsp[parmsp!=0])
      subvalue<-10^(1.1*log10(pmin))
      parms.pchange[locsub]<-subvalue
    }else{
      parms.pchange<-parmsp
    }
    colnames(parms.pchange)<-"p-value"
    
    bpnumber <- numeric()
    chrnum <- unique(gen[,1])
    for(i in 1:length(chrnum))
    {
      bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(gen[,1]==chrnum[i])))))
    }
    rowsnp <- dim(gen)[1]
    snpname <- numeric()
    snpname <- as.matrix(paste("rs",c(1:rowsnp),sep=""))
    parmstemp <- parmsShow[,-4]
    parmstemp[,4] <-parms.pchange
    parms <- data.frame(parmstemp,snpname,bpnumber)
    colnames(parms)<-c("RS#","Chromosome","Marker Position (bp)","P-value","Genotype  for code 1","SNPname","BPnumber")
    parms<-parms[,-c(1,3,5)]
    
    mannewp <- as.matrix(0.05/as.numeric(nrow(gen)))
    rowbl<-matrix("",(nrow(parms)-1),1)
    mannepr<-rbind(mannewp,rowbl)
    
    colnames(mannepr)<-"Manhattan p-value"
    
    parmsm<-cbind(as.matrix(parms),mannepr)
    parmsm<-as.data.frame(parmsm,stringsAsFactors=FALSE)
    parmsm[,c(1,2,4)]<-sapply(parmsm[,c(1,2,4)],as.numeric)
    
    finalres <- gglartotal
    if(length(finalres[,2])>1){
      
      if((flagps==1)||(exists("psmatrix")==FALSE))
      {
        ex<-cbind(fix,t(gene.data[finalres[,2],]))
      }else if(flagps==0)
      {
        ex<-cbind(cbind(fix,psmatrix),t(gene.data[finalres[,2],]))
      }
      
    }else{
      
      if((flagps==1)||(exists("psmatrix")==FALSE))
      {
        ex<-cbind(fix,as.matrix(gene.data[finalres[,2],]))
      }else if(flagps==0)
      {
        ex<-cbind(cbind(fix,psmatrix),as.matrix(gene.data[finalres[,2],]))
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
      
      if(var(phe)>=sum(er)+v0){
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
      
      if(var(phe)>=sum(er)+v0){
        her<-(er/as.vector(var(phe)))*100 
      }else{
        
        her<-(er/as.numeric(sum(er)+v0))*100 
      } 
    }
    
    X<-t(gene.data)
    XX1<-X
    x<-XX1[3:nrow(XX1),]
    X1<-as.matrix(x)
    xxxx<-as.matrix(X1[,finalres[,2]])
    
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
    
    eeff <- finalres[,5]
    lo <- finalres[,6]
    eeff[which(abs(eeff)>=1e-4)] <- round(eeff[which(abs(eeff)>=1e-4)],4)
    eeff[which(abs(eeff)<1e-4)] <- as.numeric(sprintf("%.4e",eeff[which(abs(eeff)<1e-4)]))
    lo[which(abs(lo)>=1e-4)] <- round(lo[which(abs(lo)>=1e-4)],4)
    lo[which(abs(lo)<1e-4)] <- as.numeric(sprintf("%.4e",lo[which(abs(lo)<1e-4)]))
    her[which(abs(her)>=1e-4)] <- round(her[which(abs(her)>=1e-4)],4)
    her[which(abs(her)<1e-4)] <- as.numeric(sprintf("%.4e",her[which(abs(her)<1e-4)]))
    needrs <- genRaw[-1,1]
    needrs <- as.matrix(needrs[finalres[,2]])
    needgenofor <- as.character()
    if(inputform==1)
    {
      needgenofor <- genRaw[-1,4]
      needgenofor <- as.matrix(needgenofor[finalres[,2]])
    }
    if(inputform==2)
    {
      needgenofor <- outATCG
      needgenofor <- as.matrix(needgenofor[finalres[,2]])
    }
    if(inputform==3)
    {
      needgenofor <- outATCG
      needgenofor <- as.matrix(needgenofor[finalres[,2]])
    }
    
    phevartotal<-var(Y.data)
    if(finalres[1,7]>=1e-4){finalres[1,7]<-round(finalres[1,7],4)}
    if(finalres[1,7]<1e-4){finalres[1,7]<-as.numeric(sprintf("%.4e",finalres[1,7]))}
    if(phevartotal>=1e-4){phevartotal<-round(phevartotal,4)}
    if(phevartotal<1e-4){phevartotal<-as.numeric(sprintf("%.4e",phevartotal))}
    tempvar <- dim(as.matrix(lo))[1]
    if(tempvar==1)
    {
      wan<-data.frame(needrs,t(as.matrix(gen[finalres[,2],1:2])),as.matrix(eeff),as.matrix(lo),her,maf,needgenofor,as.matrix(finalres[,7]),phevartotal) 
    }else if(tempvar>1)
    {
      wan<-data.frame(needrs,gen[finalres[,2],1:2],eeff,lo,her,maf,needgenofor,rbind(finalres[1,7],as.matrix(rep("",(tempvar-1)))),rbind(phevartotal,as.matrix(rep("",(tempvar-1))))) 
    }
    
    tempwan <- wan
    lodscore1 <- as.numeric(tempwan[,5])
    log10P <- as.matrix(round(-log10(1-pchisq(lodscore1*4.605,1)),4))
    tempwan1 <- cbind(tempwan[,1:5],log10P,tempwan[,6:10])
    wan <- tempwan1
    
    colnames(wan)<-c("RS#","Chromosome","Marker Position (bp)","QTN effect","LOD score","-log10(P)","r2 (%)","MAF","Genotype  for code 1","Var_Error","Var_Phen(total)")
    wan<-as.data.frame(wan)
    
    output<-list(result1=parmsShow,result2=wan,Manhattan=parmsm,QQ=parms.pchange)
    
    return(output)
    
  }
}