mrMLM<-function(fileGen=NULL,filePhe=NULL,fileKin=NULL,filePS=NULL,Genformat=NULL,method=NULL,
                Likelihood=NULL,trait=NULL,SearchRadius=NULL,CriLOD=NULL,SelectVariable=NULL,
                Bootstrap=NULL,DrawPlot=NULL,Plotformat=NULL,Resolution=NULL,dir=NULL){
  
  if(DrawPlot==TRUE){
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
  }
  
  svrad<-SearchRadius;svmlod<-CriLOD;lars1<-SelectVariable
  
  if(Genformat=="Num"){Genformat<-1}else if(Genformat=="Cha"){Genformat<-2}else if(Genformat=="Hmp"){Genformat<-3}
  
  Plotformat1<-paste("*.",Plotformat,sep="");Plotformat2<-paste("*.",Plotformat,sep="")
  
  readraw<-ReadData(fileGen,filePhe,fileKin,filePS,Genformat)
  PheName<-readraw$phename
  CLO<-readraw$CLO
  
  print("Running in progress, please be patient...")
  
  for (i in trait){
    
    InputData<-inputData(readraw,Genformat,method,i)
    
    reMR<-NULL;reFMR<-NULL;reFME<-NULL;rePLA<-NULL;rePKW<-NULL;reISIS<-NULL
    re1MR<-NULL;re1FMR<-NULL;re1FME<-NULL;re1PLA<-NULL;re1PKW<-NULL;re1ISIS<-NULL
    remanMR<-NULL;reqqMR<-NULL;remanFMR<-NULL;reqqFMR<-NULL;remanFME<-NULL;reqqFME<-NULL;replPLA<-NULL;remanPKW<-NULL;reqqPKW<-NULL; replISIS<-NULL
    
    TRY1<-try({
      
      if("mrMLM"%in%method){
        outMR<-mrMLMFun(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,0.01,svrad,svmlod,Genformat,CLO)
        if(is.null(outMR$result2)==FALSE){
          me<-matrix("mrMLM",nrow(outMR$result2),1)
          tr<-matrix(i,nrow(outMR$result2),1)
          trna<-matrix(PheName[i,],nrow(outMR$result2),1)
          colnames(me)<-"Method"
          colnames(tr)<-"Trait ID"
          colnames(trna)<-"Trait name"
          reMR<-cbind(tr,trna,me,as.matrix(outMR$result2))
        }
        
        me1<-matrix("mrMLM",nrow(outMR$result1),1)
        tr1<-matrix(i,nrow(outMR$result1),1)
        tr1na<-matrix(PheName[i,],nrow(outMR$result1),1)
        
        colnames(me1)<-"Method"
        colnames(tr1)<-"Trait ID"
        colnames(tr1na)<-"Trait name"
        re1MR<-cbind(tr1,tr1na,me1,as.matrix(outMR$result1))
        
        remanMR<-outMR$Manhattan
        reqqMR<-outMR$QQ
        
        if(DrawPlot==TRUE){
          
          if(Resolution=="Low"){
            manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
          }else if(Resolution=="High"){
            manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-300 
          }
          
          if(Plotformat1=="*.png"){
            png(paste(dir,"/",i,"_mrMLM_Manhattan.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat1=="*.tiff"){
            tiff(paste(dir,"/",i,"_mrMLM_Manhattan.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat1=="*.jpeg"){
            jpeg(paste(dir,"/",i,"_mrMLM_Manhattan.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat1=="*.pdf"){
            pdf(paste(dir,"/",i,"_mrMLM_Manhattan.pdf",sep=""),width=10)
          }
          
          Plot(plotresult=remanMR,color1="red",color2="blue",0.95,method="mrMLM",type="Manhattan")
          dev.off()
          
          if(Plotformat2=="*.png"){
            png(paste(dir,"/",i,"_mrMLM_qq.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat2=="*.tiff"){
            tiff(paste(dir,"/",i,"_mrMLM_qq.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat2=="*.jpeg"){
            jpeg(paste(dir,"/",i,"_mrMLM_qq.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat2=="*.pdf"){
            pdf(paste(dir,"/",i,"_mrMLM_qq.pdf",sep=""),width=10)
          }
          
          Plot(plotresult=reqqMR,color1="red",color2="blue",0.95,method="mrMLM",type="qq")
          dev.off()
        }
      }
      
      
    },silent=FALSE)  
    
    
    if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){  
      TRY2<-try({
        
        if("FASTmrMLM"%in%method){
          outFMR<-FASTmrMLM(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,0.01,svrad,svmlod,Genformat,CLO)   
          if(is.null(outFMR$result2)==FALSE){
            me<-matrix("FASTmrMLM",nrow(outFMR$result2),1)
            tr<-matrix(i,nrow(outFMR$result2),1)
            trna<-matrix(PheName[i,],nrow(outFMR$result2),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            reFMR<-cbind(tr,trna,me,as.matrix(outFMR$result2))
          }
          me1<-matrix("FASTmrMLM",nrow(outFMR$result1),1)
          tr1<-matrix(i,nrow(outFMR$result1),1)
          tr1na<-matrix(PheName[i,],nrow(outFMR$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1FMR<-cbind(tr1,tr1na,me1,as.matrix(outFMR$result1))
          
          remanFMR<-outFMR$Manhattan
          reqqFMR<- outFMR$QQ
          
          if(DrawPlot==TRUE){
            
            if(Resolution=="Low"){
              manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
            }else if(Resolution=="High"){
              manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-300 
            }
            if(Plotformat1=="*.png"){
              png(paste(dir,"/",i,"_FASTmrMLM_Manhattan.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(dir,"/",i,"_FASTmrMLM_Manhattan.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_FASTmrMLM_Manhattan.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(dir,"/",i,"_FASTmrMLM_Manhattan.pdf",sep=""),width=10)
            }
            
            Plot(plotresult=remanFMR,color1="red",color2="blue",0.95,method="FASTmrMLM",type="Manhattan")
            dev.off()
            
            if(Plotformat2=="*.png"){
              png(paste(dir,"/",i,"_FASTmrMLM_qq.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.tiff"){
              tiff(paste(dir,"/",i,"_FASTmrMLM_qq.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_FASTmrMLM_qq.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.pdf"){
              pdf(paste(dir,"/",i,"_FASTmrMLM_qq.pdf",sep=""),width=10)
            }
            
            Plot(plotresult=reqqFMR,color1="red",color2="blue",0.95,method="FASTmrMLM",type="qq")
            dev.off()
          }
        }
        
      },silent=FALSE)  
    }
    
    
    if ('try-error' %in% class(TRY2)|| !('try-error' %in% class(TRY2))){
      
      TRY3<-try({
        
        if("FASTmrEMMA"%in%method){
          outFME<-FASTmrEMMA(InputData$doFME$gen,InputData$doFME$phe,InputData$doFME$outATCG,InputData$doFME$genRaw,InputData$doFME$kk,InputData$doFME$psmatrix,0.005,svmlod,Genformat,Likelihood,CLO)   
          
          if(is.null(outFME$result2)==FALSE){
            me<-matrix("FASTmrEMMA",nrow(outFME$result2),1)
            tr<-matrix(i,nrow(outFME$result2),1)
            trna<-matrix(PheName[i,],nrow(outFME$result2),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            reFME<-cbind(tr,trna,me,as.matrix(outFME$result2))
          }
          
          me1<-matrix("FASTmrEMMA",nrow(outFME$result1),1)
          tr1<-matrix(i,nrow(outFME$result1),1)
          tr1na<-matrix(PheName[i,],nrow(outFME$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1FME<-cbind(tr1,tr1na,me1,as.matrix(outFME$result1))
          
          remanFME<-outFME$Manhattan
          reqqFME<-outFME$QQ
          
          if(DrawPlot==TRUE){
            if(Resolution=="Low"){
              manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
            }else if(Resolution=="High"){
              manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-300 
            }
            if(Plotformat1=="*.png"){
              png(paste(dir,"/",i,"_FASTmrEMMA_Manhattan.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(dir,"/",i,"_FASTmrEMMA_Manhattan.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_FASTmrEMMA_Manhattan.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(dir,"/",i,"_FASTmrEMMA_Manhattan.pdf",sep=""),width=10)
            }
            
            Plot(plotresult=remanFME,color1="red",color2="blue",0.95,method="FASTmrEMMA",type="Manhattan")
            dev.off()
            
            if(Plotformat2=="*.png"){
              png(paste(dir,"/",i,"_FASTmrEMMA_qq.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.tiff"){
              tiff(paste(dir,"/",i,"_FASTmrEMMA_qq.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_FASTmrEMMA_qq.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.pdf"){
              pdf(paste(dir,"/",i,"_FASTmrEMMA_qq.pdf",sep=""),width=10)
            }
            
            Plot(plotresult=reqqFME,color1="red",color2="blue",0.95,method="FASTmrEMMA",type="qq")
            dev.off()
          }
        }
      },silent=FALSE)  
      
    } 
    
    
    if ('try-error' %in% class(TRY3)|| !('try-error' %in% class(TRY3))){ 
      
      TRY4<-try({
        
        if("pLARmEB"%in%method){
          outPLA<-pLARmEB(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,CriLOD,lars1,Genformat,Bootstrap,CLO)   
          if(is.null(outPLA$result)==FALSE){
            me<-matrix("pLARmEB",nrow(outPLA$result),1)
            tr<-matrix(i,nrow(outPLA$result),1)
            trna<-matrix(PheName[i,],nrow(outPLA$result),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            rePLA<-cbind(tr,trna,me,as.matrix(outPLA$result))
          }
          replPLA<-outPLA$plot
          
          if(DrawPlot==TRUE){
            if(Resolution=="Low"){
              manwidth<-960;manhei<-240;manwordre<-12;manfigurere<-72
            }else if(Resolution=="High"){
              manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-600 
            }
            if(Plotformat1=="*.png"){
              png(paste(dir,"/",i,"_pLARmEB_LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(dir,"/",i,"_pLARmEB_LOD.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_pLARmEB_LOD.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(dir,"/",i,"_pLARmEB_LOD.pdf",sep=""),width=12)
            }
            
            Plot(plotresult=replPLA,color1="red",color2="blue",0.95,method="pLARmEB",type="LOD")
            dev.off()
          }
        }
      },silent=FALSE)  
      
    } 
    
    
    if ('try-error' %in% class(TRY4)|| !('try-error' %in% class(TRY4))){  
      
      TRY5<-try({
        
        if("pKWmEB"%in%method){
          outPKW<-pKWmEB(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,0.05,svmlod,Genformat,CLO)  
          
          if(is.null(outPKW$result2)==FALSE){
            me<-matrix("pKWmEB",nrow(outPKW$result2),1)
            tr<-matrix(i,nrow(outPKW$result2),1)
            trna<-matrix(PheName[i,],nrow(outPKW$result2),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            rePKW<-cbind(tr,trna,me,as.matrix(outPKW$result2))
          }
          
          me1<-matrix("pKWmEB",nrow(outPKW$result1),1)
          tr1<-matrix(i,nrow(outPKW$result1),1)
          tr1na<-matrix(PheName[i,],nrow(outPKW$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1PKW<-cbind(tr1,tr1na,me1,as.matrix(outPKW$result1))
          
          remanPKW<-outPKW$Manhattan
          reqqPKW<-outPKW$QQ
          
          if(DrawPlot==TRUE){
            if(Resolution=="Low"){
              manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
            }else if(Resolution=="High"){
              manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-300 
            }
            if(Plotformat1=="*.png"){
              png(paste(dir,"/",i,"_pKWmEB_Manhattan.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(dir,"/",i,"_pKWmEB_Manhattan.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_pKWmEB_Manhattan.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(dir,"/",i,"_pKWmEB_Manhattan.pdf",sep=""),width=10)
            }
            
            Plot(plotresult=remanPKW,color1="red",color2="blue",0.95,method="pKWmEB",type="Manhattan")
            dev.off()
            
            if(Plotformat2=="*.png"){
              png(paste(dir,"/",i,"_pKWmEB_qq.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.tiff"){
              tiff(paste(dir,"/",i,"_pKWmEB_qq.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_pKWmEB_qq.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.pdf"){
              pdf(paste(dir,"/",i,"_pKWmEB_qq.pdf",sep=""),width=10)
            }
            
            Plot(plotresult=reqqPKW,color1="red",color2="blue",0.95,method="pKWmEB",type="qq")
            dev.off()
          }
        }
      },silent=FALSE)  
    } 
    
    
    if ('try-error' %in% class(TRY5)|| !('try-error' %in% class(TRY5))){   
      
      TRY6<-try({ 
        
        if("ISIS EM-BLASSO"%in%method){
          outISIS<-ISIS(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,0.01,svmlod,Genformat,CLO) 
          if(is.null(outISIS$result)==FALSE){
            me<-matrix("ISIS EM-BLASSO",nrow(outISIS$result),1)
            tr<-matrix(i,nrow(outISIS$result),1)
            trna<-matrix(PheName[i,],nrow(outISIS$result),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            reISIS<-cbind(tr,trna,me,as.matrix(outISIS$result))
          }
          replISIS<-outISIS$plot
          
          if(DrawPlot==TRUE){
            if(Resolution=="Low"){
              manwidth<-960;manhei<-240;manwordre<-12;manfigurere<-72
            }else if(Resolution=="High"){
              manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-600 
            }
            if(Plotformat1=="*.png"){
              png(paste(dir,"/",i,"_ISIS EM-BLASSO_LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(dir,"/",i,"_ISIS EM-BLASSO_LOD.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_ISIS EM-BLASSO_LOD.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(dir,"/",i,"_ISIS EM-BLASSO_LOD.pdf",sep=""),width=12)
            }
            
            Plot(plotresult=replISIS,color1="red",color2="blue",0.95,method="ISIS EM-BLASSO",type="LOD")
            dev.off()
          }
        }
        
      },silent=FALSE)  
    }
    
    
    if ('try-error' %in% class(TRY6)|| !('try-error' %in% class(TRY6))){    
      
      TRY7<-try({ 
        
        output1<-list(re1MR,re1FMR,re1FME,re1PKW)
        output1<-do.call(rbind,output1)
        output<-list(reMR,reFMR,reFME,rePLA,rePKW,reISIS)
        output<-do.call(rbind,output)
        write.table(output,paste(dir,"/",i,"_Final result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
        write.table(output1,paste(dir,"/",i,"_intermediate result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
      },silent=FALSE)  
      
    }
  }
  
}


inputData<-function(readraw,Genformat=NULL,method=NULL,trait=NULL){
  
  doMR<-NULL;doFME<-NULL
  
  if("mrMLM"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2)
  }
  
  if("FASTmrMLM"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2)
  }  
  
  if("FASTmrEMMA"%in%method){
    doFME<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=1)  
  }
  
  if("pLARmEB"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2) 
  }
  if("pKWmEB"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2) 
  }  
  
  if("ISIS EM-BLASSO"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2) 
  }
  
  output<-list(doMR=doMR,doFME=doFME) 
  return(output)
  
}


DoData<-function(genRaw=NULL,Genformat=NULL,pheRaw1q=NULL,kkRaw=NULL,psmatrixRaw=NULL,trait=NULL,type=NULL){
  inputform<-Genformat
  pheRaw1qq<-as.matrix(pheRaw1q[,2:ncol(pheRaw1q)])
  pheRaw1<-cbind(pheRaw1q[,1],pheRaw1qq[,trait])
  pheRaw2<-pheRaw1[-1,]
  pheRaw3<-as.data.frame(pheRaw2,stringsAsFactors=FALSE)
  pheRaw4<-as.matrix(pheRaw3[is.na(pheRaw3[,2])==F,])
  pheRawthem<-matrix(c(pheRaw1[1,1]," "),1,)
  pheRaw<-rbind(pheRawthem,pheRaw4)
  row.names(pheRaw)<-NULL
  pheRaw<-as.matrix(pheRaw)
  
  if(type==1&&inputform==1){
    genRawz<-genRaw[-1,-c(1:4)]
    genRawz2<-gsub("0","0.5",genRawz)
    genRawz3<-gsub("-1","0",genRawz2)
    genRawz4<-cbind(genRaw[-1,c(1:4)],genRawz3)
    genRaw<-rbind(genRaw[1,],genRawz4)
  }else{
    genRaw<-genRaw
  }
  if(inputform==1){
    nameGen <- as.matrix(genRaw[1,],1,)
    namePhe <- as.matrix(pheRaw[,1],,1)
    sameName <- intersect(nameGen,namePhe)
    ##########To find the location of the same name
    locGen <- match(sameName,nameGen)
    locPhe <- match(sameName,namePhe)
    ##########Produce new genotype matrix and phenotype matrix
    hapName <- matrix(c("rs#","chrom","pos","genotype for code 1"),1,)
    hapHave <- intersect(nameGen,hapName)
    locHap <- match(hapHave,nameGen)
    newGenloc <- c(locHap,locGen)
    newPheloc <- locPhe
    newGen <- as.matrix(genRaw[-1,newGenloc])
    newPhe <- as.matrix(pheRaw[newPheloc,])
    nnhap <- length(hapHave)
    rownewGen <- dim(newGen)[1]
    colnewGen <- dim(newGen)[2]
    rownewPhe <- dim(newPhe)[1]
    ###########To show on the table ----newGen
    newGen <-rbind(genRaw[1,newGenloc],newGen)
    ###########To be computed ----gen
    locChr <- as.numeric(which(newGen[1,]=="chrom"))
    locPos <- as.numeric(which(newGen[1,]=="pos"))
    needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
    needGen <- newGen[,needloc]
    
    genq<-as.matrix(needGen[-1,])
    gen<-big.matrix(nrow(genq),ncol(genq),type="double",shared = FALSE)
    gen[,]<-genq[,]
    rm(newGen,needGen,genq)
    gc()
    ###########To show on the table ----newPhe
    pheRaw[1,2]<-"  "
    newPhe<-rbind(pheRaw[1,],newPhe)
    ###########To be computed ----phe
    phe<-as.matrix(newPhe[-1,-1])
    phe<-matrix(as.numeric(phe),nrow=nrow(phe))
    outATCG<-NULL
  }else if(inputform==2){
    ##########To find the same individual ID between genotype and phenotype
    nameGen <- as.matrix(genRaw[1,],1,)
    namePhe <- as.matrix(pheRaw[,1],,1)
    sameName <- intersect(nameGen,namePhe)
    ##########To find the location of the same name
    locGen <- match(sameName,nameGen)
    locPhe <- match(sameName,namePhe)
    ##########Produce new genotype matrix and phenotype matrix
    hapName <- matrix(c("rs#","chrom","pos"),1,)
    hapHave <- intersect(nameGen,hapName)
    locHap <- match(hapHave,nameGen)
    newGenloc <- c(locHap,locGen)
    newPheloc <- locPhe
    newGen <- as.matrix(genRaw[-1,newGenloc])
    newPhe <- as.matrix(pheRaw[newPheloc,])
    ##########Transfer ATCG to numeric
    nnhap <- length(hapHave)
    rownewGen <- dim(newGen)[1]
    colnewGen <- dim(newGen)[2]
    rownewPhe <- dim(newPhe)[1]
    computeGen <- newGen[,(nnhap+1):colnewGen]
    colComGen <- ncol(computeGen)
    referSam <- as.vector(computeGen[,1])
    ATCGloc <- c(which(computeGen[,1]=="A"),which(computeGen[,1]=="T"),which(computeGen[,1]=="C"),which(computeGen[,1]=="G"))
    NNRRloc <- setdiff(c(1:rownewGen),ATCGloc)
    for(i in 2:colComGen)
    {
      if(length(NNRRloc)>0){
        referSam[NNRRloc] <- as.vector(computeGen[NNRRloc,i])
        ATCGlocLoop <- c(which(computeGen[NNRRloc,i]=="A"),which(computeGen[NNRRloc,i]=="T"),which(computeGen[NNRRloc,i]=="C"),which(computeGen[NNRRloc,i]=="G"))
        NNRRloc <- setdiff(NNRRloc,NNRRloc[ATCGlocLoop])
      }else{
        break
      }
    }
    for(i in 1:rownewGen)
    {
      tempSel1 <- as.vector(c(which(computeGen[i,]=="A"),which(computeGen[i,]=="T"),which(computeGen[i,]=="C"),which(computeGen[i,]=="G")))
      tempSel2 <- as.vector(c(which(computeGen[i,]==referSam[i])))
      notRef <- setdiff(tempSel1,tempSel2)
      notATCG <- setdiff(c(1:colComGen),tempSel1)
      computeGen[i,tempSel2] <- as.numeric(1)
      
      if(type==1){
        computeGen[i,notRef] <- as.numeric(0)
        computeGen[i,notATCG] <- as.numeric(0.5)
      }else{
        computeGen[i,notRef] <- as.numeric(-1)
        computeGen[i,notATCG] <- as.numeric(0)
      }
    }
    outATCG<-as.matrix(referSam)
    ###########To show on the table ----newGen
    newGen <- cbind(newGen[,1:nnhap],computeGen)
    newGen <-rbind(genRaw[1,newGenloc],newGen)
    rm(computeGen)
    gc()
    ###########To be computed ----gen
    locChr <- as.numeric(which(newGen[1,]=="chrom"))
    locPos <- as.numeric(which(newGen[1,]=="pos"))
    needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
    needGen<-newGen[,needloc]
    
    genq<-as.matrix(needGen[-1,])
    gen<-big.matrix(nrow(genq),ncol(genq),type="double",shared = FALSE)
    gen[,]<-genq[,]
    rm(newGen,needGen,genq)
    gc()
    ###########To show on the table ----newPhe
    pheRaw[1,2]<-"  "
    newPhe<-rbind(pheRaw[1,],newPhe)
    ###########To be computed ----phe
    phe<-as.matrix(newPhe[-1,-1])
    phe<-matrix(as.numeric(phe),nrow=nrow(phe))
  }else if(inputform==3){
    ##########To find the same individual ID between genotype and phenotype
    nameGen<-as.matrix(genRaw[1,],1,)
    namePhe<-as.matrix(pheRaw[,1],,1)
    sameName<-intersect(nameGen,namePhe)
    ##########To find the location of the same name
    locGen<-match(sameName,nameGen)
    locPhe<-match(sameName,namePhe)
    ##########Produce new genotype matrix and phenotype matrix
    hapName<-matrix(c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panel","QCcode"),1,)
    hapHave<-intersect(nameGen,hapName)
    locHap<-match(hapHave,nameGen)
    newGenloc<-c(locHap,locGen)
    newPheloc<-locPhe
    newGen<-as.matrix(genRaw[-1,newGenloc])
    newPhe<-as.matrix(pheRaw[newPheloc,])
    ##########Transfer ATCG to numeric
    nnhap<-length(hapHave)
    rownewGen<-dim(newGen)[1]
    colnewGen<-dim(newGen)[2]
    rownewPhe<-dim(newPhe)[1]
    computeGen<-newGen[,(nnhap+1):colnewGen]
    colComGen<-ncol(computeGen)
    referSam<-as.vector(computeGen[,1])
    ATCGloc<-c(which(computeGen[,1]=="AA"),which(computeGen[,1]=="TT"),which(computeGen[,1]=="CC"),which(computeGen[,1]=="GG"))
    NNRRloc<-setdiff(c(1:rownewGen),ATCGloc)
    for(i in 2:colComGen)
    {
      if(length(NNRRloc)>0){
        referSam[NNRRloc]<-as.vector(computeGen[NNRRloc,i])
        ATCGlocLoop<-c(which(computeGen[NNRRloc,i]=="AA"),which(computeGen[NNRRloc,i]=="TT"),which(computeGen[NNRRloc,i]=="CC"),which(computeGen[NNRRloc,i]=="GG"))
        NNRRloc<-setdiff(NNRRloc,NNRRloc[ATCGlocLoop])
      }else{
        break
      }
    }
    for(i in 1:rownewGen)
    {
      tempSel1<-as.vector(c(which(computeGen[i,]=="AA"),which(computeGen[i,]=="TT"),which(computeGen[i,]=="CC"),which(computeGen[i,]=="GG")))
      tempSel2<-as.vector(c(which(computeGen[i,]==referSam[i])))
      notRef<-setdiff(tempSel1,tempSel2)
      notATCG<-setdiff(c(1:colComGen),tempSel1)
      computeGen[i,tempSel2]<-as.numeric(1)
      if(type==1){
        computeGen[i,notRef]<-as.numeric(0)
        computeGen[i,notATCG]<-as.numeric(0.5)
      }else{
        computeGen[i,notRef]<-as.numeric(-1)
        computeGen[i,notATCG]<-as.numeric(0)
      }
    }
    outATCG<-as.matrix(referSam)
    ###########To show on the table ----newGen
    newGen<-cbind(newGen[,1:nnhap],computeGen)
    newGen<-rbind(genRaw[1,newGenloc],newGen)
    
    rm(computeGen)
    gc()
    ###########To be computed ----gen
    locChr<-as.numeric(which(newGen[1,]=="chrom"))
    locPos<-as.numeric(which(newGen[1,]=="pos"))
    needloc<-c(locChr,locPos,(nnhap+1):colnewGen)
    needGen<-newGen[,needloc]
    
    genq<-as.matrix(needGen[-1,])
    gen<-big.matrix(nrow(genq),ncol(genq),type="double",shared = FALSE)
    gen[,]<-genq[,]
    rm(newGen,needGen,genq)
    gc()
    ###########To show on the table ----newPhe
    pheRaw[1,2]<-"  "
    newPhe<-rbind(pheRaw[1,],newPhe)
    ###########To be computed ----phe
    phe<-as.matrix(newPhe[-1,-1])
    phe<-matrix(as.numeric(phe),nrow=nrow(phe))
    
  }
  
  if(is.null(kkRaw)){
    kk<-NULL
  }else{
    kkPre<-as.matrix(kkRaw[-1,-1])
    nameKin<-as.matrix(kkRaw[-1,1])
    sameGenKin<-intersect(sameName,nameKin)
    locKin<-match(sameGenKin,nameKin)
    kk<-kkPre[locKin,locKin]
    kk<-matrix(as.numeric(kk),nrow=nrow(kk))
  }
  
  
  if(is.null(psmatrixRaw)){
    psmatrix<-NULL
  }else{
    nnpprow<-dim(psmatrixRaw)[1]
    nnppcol<-dim(psmatrixRaw)[2]
    psmatrixRaw[1,2:nnppcol]<-"  "
    psmatrixPre<-psmatrixRaw[3:nnpprow,]
    namePop<-as.matrix(psmatrixPre[,1])
    sameGenPop<-intersect(sameName,namePop)
    locPop<-match(sameGenPop,namePop)
    ##revised
    filtername<-as.vector(psmatrixRaw[2,2:nnppcol])
    selectpsmatrix<-matrix(as.numeric(psmatrixPre[locPop,-1]),nrow = length(locPop))
    psum<-apply(selectpsmatrix,1,sum)
    psum<-round(psum)
    sumps<-sum(psum)
    m<-dim(selectpsmatrix)[1]
    if(sumps>=m){
      combovalue<-"Q1"
      coldelet<-unlist(str_extract_all(combovalue,"[0-9]+"))
      coldelet<-as.numeric(coldelet)
      psmatrix<-as.matrix(selectpsmatrix[,-coldelet])
      psmatrixRaw<-as.matrix(psmatrixRaw[,-(coldelet+1)])
    }else{
      psmatrix<-selectpsmatrix
    }
  }
  genRaw<-genRaw[,1:12]
  doresult<-list(gen=gen,phe=phe,outATCG=outATCG,genRaw=genRaw,kk=kk,psmatrix=psmatrix)
  return(doresult)
}


ReadData<-function(fileGen=NULL,filePhe=NULL,fileKin=NULL,filePS=NULL,Genformat=NULL){
  kkRaw<-NULL
  psmatrixRaw<-NULL
  inputform<-Genformat
  CLO<-NULL
  if(!is.null(fileGen)){
    if(is.character(fileGen)==TRUE){
      genRaw<-fread(fileGen,header = FALSE,stringsAsFactors=T)
    }else{
      genRaw<-fileGen
      CLO<-1
    }
    genRaw<-as.matrix(genRaw)
  }
  wnameGen <- as.matrix(genRaw[1,],1,)
  if(inputform==1){
    
    titlenameGen<-wnameGen[1:4,]
    hapName<-c("rs#","chrom","pos","genotype for code 1")
    
    if(all(titlenameGen==hapName)==FALSE){
      warning("please check the individual's name in genotypic file")
    }
  }
  
  if(inputform==2){
    
    titlenameGen<-wnameGen[1:3,]
    hapName<-c("rs#","chrom","pos")
    
    if(all(titlenameGen==hapName)==FALSE){
      warning("please check the individual's name in genotypic file")
    }
  }
  
  if(inputform==3){
    
    titlenameGen<-wnameGen[1:11,]
    hapName<- c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panel","QCcode")
    
    if(all(titlenameGen==hapName)==FALSE){
      warning("please check the individual's name in genotypic file")
    }
  }
  
  
  
  if(!is.null(filePhe)){
    
    if(is.character(filePhe)==TRUE){
      pheRaw1q<-fread(filePhe,header=F, stringsAsFactors=T) 
    }else{
      pheRaw1q<-filePhe
      CLO<-1
    }
    pheRaw1q<-as.matrix(pheRaw1q)
  }
  
  
  wnamePhe <- as.matrix(pheRaw1q[,1],,1)
  wsameName <- intersect(wnameGen,wnamePhe)
  wlocGen <- match(wsameName,wnameGen)
  if(is.null(wlocGen)){
    warning("please check the individual's name (ID) in genotypic and phenotypic files")
  }
  
  if(!is.null(fileKin)){
    kkRaw<-fread(fileKin,header = FALSE,stringsAsFactors=T)
    kkRaw<-as.matrix(kkRaw)
    nnkk<-dim(kkRaw)[1]
    kkRaw[1,2:nnkk]<-"  "
  }
  
  if(!is.null(filePS)){
    psmatrixRaw<-fread(filePS,header = FALSE,stringsAsFactors=T)
    psmatrixRaw<-as.matrix(psmatrixRaw)
  }
  
  
  phename<-as.matrix(pheRaw1q[1,2:ncol(pheRaw1q)])
  
  output<-list(genRaw=genRaw,pheRaw1q=pheRaw1q,kkRaw=kkRaw,psmatrixRaw=psmatrixRaw,phename=phename,CLO=CLO)
  return(output)
}

