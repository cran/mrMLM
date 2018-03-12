mrMLM<-function(fileGen=NULL,filePhe=NULL,fileKin=NULL,filePS=NULL,Genformat=NULL,method=NULL,
             Likelihood=NULL,trait=NULL,SearchRadius=NULL,CriLOD=NULL,SelectVariable=NULL,
             Bootstrap=NULL,DrawPlot=NULL,Plotformat=NULL,Resolution=NULL){
  
  svrad<-SearchRadius;svmlod<-CriLOD;lars1<-SelectVariable
  
  if(Genformat=="Num"){Genformat<-1}else if(Genformat=="Cha"){Genformat<-2}else if(Genformat=="Hmp"){Genformat<-3}
 
  Plotformat1<-paste("*.",Plotformat,sep="");Plotformat2<-paste("*.",Plotformat,sep="")
  
  readraw<-ReadData(fileGen,filePhe,fileKin,filePS,Genformat)
  PheName<-readraw$phename
  CLO<-readraw$CLO
  
  print("Running in progress, please be patient...")
  
  for (i in trait){
    
    InputData<-InputData(readraw,Genformat,method,i)
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
            png(paste(i,"_mrMLM_Manhattan.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat1=="*.tiff"){
            tiff(paste(i,"_mrMLM_Manhattan.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat1=="*.jpeg"){
            jpeg(paste(i,"_mrMLM_Manhattan.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat1=="*.pdf"){
            pdf(paste(i,"_mrMLM_Manhattan.pdf"),width=10)
          }
          
          Plot(plotresult=remanMR,color1="red",color2="blue",0.95,method="mrMLM",type="Manhattan")
          dev.off()
          
          if(Plotformat2=="*.png"){
            png(paste(i,"_mrMLM_qq.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat2=="*.tiff"){
            tiff(paste(i,"_mrMLM_qq.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat2=="*.jpeg"){
            jpeg(paste(i,"_mrMLM_qq.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
          }else if(Plotformat2=="*.pdf"){
            pdf(paste(i,"_mrMLM_qq.pdf"),width=10)
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
              png(paste(i,"_FASTmrMLM_Manhattan.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(i,"_FASTmrMLM_Manhattan.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(i,"_FASTmrMLM_Manhattan.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(i,"_FASTmrMLM_Manhattan.pdf"),width=10)
            }
            
            Plot(plotresult=remanFMR,color1="red",color2="blue",0.95,method="FASTmrMLM",type="Manhattan")
            dev.off()
            
            if(Plotformat2=="*.png"){
              png(paste(i,"_FASTmrMLM_qq.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.tiff"){
              tiff(paste(i,"_FASTmrMLM_qq.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.jpeg"){
              jpeg(paste(i,"_FASTmrMLM_qq.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.pdf"){
              pdf(paste(i,"_FASTmrMLM_qq.pdf"),width=10)
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
              png(paste(i,"_FASTmrEMMA_Manhattan.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(i,"_FASTmrEMMA_Manhattan.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(i,"_FASTmrEMMA_Manhattan.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(i,"_FASTmrEMMA_Manhattan.pdf"),width=10)
            }
            
            Plot(plotresult=remanFME,color1="red",color2="blue",0.95,method="FASTmrEMMA",type="Manhattan")
            dev.off()
            
            if(Plotformat2=="*.png"){
              png(paste(i,"_FASTmrEMMA_qq.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.tiff"){
              tiff(paste(i,"_FASTmrEMMA_qq.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.jpeg"){
              jpeg(paste(i,"_FASTmrEMMA_qq.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.pdf"){
              pdf(paste(i,"_FASTmrEMMA_qq.pdf"),width=10)
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
              manwidth<-1000;manhei<-6000;manwordre<-30;manfigurere<-300 
            }
            if(Plotformat1=="*.png"){
              png(paste(i,"_pLARmEB_LOD.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(i,"_pLARmEB_LOD.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(i,"_pLARmEB_LOD.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(i,"_pLARmEB_LOD.pdf"),width=12)
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
              png(paste(i,"_pKWmEB_Manhattan.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(i,"_pKWmEB_Manhattan.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(i,"_pKWmEB_Manhattan.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(i,"_pKWmEB_Manhattan.pdf"),width=10)
            }
            
            Plot(plotresult=remanPKW,color1="red",color2="blue",0.95,method="pKWmEB",type="Manhattan")
            dev.off()
            
            if(Plotformat2=="*.png"){
              png(paste(i,"_pKWmEB_qq.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.tiff"){
              tiff(paste(i,"_pKWmEB_qq.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.jpeg"){
              jpeg(paste(i,"_pKWmEB_qq.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat2=="*.pdf"){
              pdf(paste(i,"_pKWmEB_qq.pdf"),width=10)
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
              manwidth<-1000;manhei<-6000;manwordre<-30;manfigurere<-300 
            }
            if(Plotformat1=="*.png"){
              png(paste(i,"_ISIS EM-BLASSO_LOD.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.tiff"){
              tiff(paste(i,"_ISIS EM-BLASSO_LOD.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(i,"_ISIS EM-BLASSO_LOD.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
            }else if(Plotformat1=="*.pdf"){
              pdf(paste(i,"_ISIS EM-BLASSO_LOD.pdf"),width=12)
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
     
        write.table(output,paste(i,"_Final result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
        
        write.table(output1,paste(i,"_intermediate result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
        
      },silent=FALSE)  
      
    }
  }
  
}


InputData<-function(readraw,Genformat=NULL,method=NULL,trait=NULL){

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
    genRaw[which(genRaw==0)]<-0.5
    genRaw[which(genRaw==-1)]<-0 
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
    gen<-as.matrix(needGen[-1,])
    gen<-matrix(as.numeric(gen),nrow=nrow(gen))
    ###########To show on the table ----newPhe
    pheRaw[1,2]<-"  "
    newPhe<-rbind(pheRaw[1,],newPhe)
    ###########To be computed ----phe
    phe<-as.matrix(newPhe[-1,-1])
    phe<-matrix(as.numeric(phe),nrow=nrow(phe))
    shownewGen<-newGen[-1,]
    colnames(shownewGen)<-newGen[1,]
    shownewGen<-as.data.frame(shownewGen)
    shownewPhe<-newPhe[-1,]
    colnames(shownewPhe)<-c(newPhe[1,1],"   ")
    shownewPhe<-as.data.frame(shownewPhe)
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
    ###########To be computed ----gen
    locChr <- as.numeric(which(newGen[1,]=="chrom"))
    locPos <- as.numeric(which(newGen[1,]=="pos"))
    needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
    needGen<-newGen[,needloc]
    gen<-as.matrix(needGen[-1,])
    gen<-matrix(as.numeric(gen),nrow=nrow(gen))
    ###########To show on the table ----newPhe
    pheRaw[1,2]<-"  "
    newPhe<-rbind(pheRaw[1,],newPhe)
    ###########To be computed ----phe
    phe<-as.matrix(newPhe[-1,-1])
    phe<-matrix(as.numeric(phe),nrow=nrow(phe))
    shownewGen<-newGen[-1,]
    colnames(shownewGen)<-newGen[1,]
    shownewGen<-as.data.frame(shownewGen)
    shownewPhe<-newPhe[-1,]
    colnames(shownewPhe)<-c(newPhe[1,1],"   ")
    shownewPhe<-as.data.frame(shownewPhe)
    
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
    ###########To be computed ----gen
    locChr<-as.numeric(which(newGen[1,]=="chrom"))
    locPos<-as.numeric(which(newGen[1,]=="pos"))
    needloc<-c(locChr,locPos,(nnhap+1):colnewGen)
    needGen<-newGen[,needloc]
    gen<-as.matrix(needGen[-1,])
    gen<-matrix(as.numeric(gen),nrow=nrow(gen))
    ###########To show on the table ----newPhe
    pheRaw[1,2]<-"  "
    newPhe<-rbind(pheRaw[1,],newPhe)
    ###########To be computed ----phe
    phe<-as.matrix(newPhe[-1,-1])
    phe<-matrix(as.numeric(phe),nrow=nrow(phe))
    shownewGen<-newGen[-1,]
    colnames(shownewGen)<-newGen[1,]
    shownewGen<-as.data.frame(shownewGen)
    shownewPhe<-newPhe[-1,]
    colnames(shownewPhe)<-c(newPhe[1,1],"   ")
    shownewPhe<-as.data.frame(shownewPhe)
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

