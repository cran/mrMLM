mrMLM<-function()
{
  inte_env <- new.env()
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
  inte_window<-gwindow(title="Methodologies for linkage and association analysis",visible=TRUE,width=420,height=300,expand=TRUE)
  inte_lyt<-glayout(container=inte_window,spacing=13)
  mrMLMbutton<-gbutton("Multi-locus Random-SNP-effect Mixed Linear Model (mrMLM) for GWAS",container=inte_lyt)
  GCIMbutton<-gbutton("Genome-wide Composite Interval Mapping (GCIM) for Linkage Analysis",container=inte_lyt)
  inte_lyt[5,1,expand=TRUE]<-mrMLMbutton
  inte_lyt[6,1,expand=TRUE]<-glabel("Cited: Wang et al. 2016 Sci Rep 6:19444")
  inte_lyt[13,1,expand=TRUE]<-GCIMbutton
  inte_lyt[14,1,expand=TRUE]<-glabel("Cited: Wang et al. 2016 Sci Rep 6:29951")
  
  mrMLMSub<-function(){
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
    
    window<-gwindow(title="Multilocus Random-SNP-effect Mixed Linear Model (mrMLM)",visible=TRUE,width=1260,height=730,expand=TRUE)
    importwin<-gwindow("Input Dataset",visible=FALSE,width=250,height=420)
    gimpwin<-ggroup(container=importwin,expand=FALSE)
    
    plotwin<-gwindow("Manhattan Plot",visible=FALSE,width=600,height=360)
    gpw<-ggroup(container=plotwin,horizontal = FALSE)
    ggpw<-ggraphics(container=gpw)
    plotwin1<-gwindow("Q-Q Plot",visible=FALSE,width=600,height=360)
    gpw1<-ggroup(container=plotwin1,horizontal=FALSE)
    ggpw1<-ggraphics(container=gpw1)
    choicekk<-gwindow("Choose Kinship option",visible=FALSE,width=320,height=150)
    gkk<-ggroup(container=choicekk,expand=FALSE)
    includeps<-gwindow("Include population structure?",visible=FALSE,width=320,height=150)
    gps<-ggroup(container=includeps,expand=FALSE)
    parsetwin<-gwindow("Parameter Setting",visible=FALSE,width=260,height=280)
    gpar<-ggroup(container=parsetwin,expand=FALSE)
    choicesave<-gwindow("Save as ...",visible=FALSE,width=250,height=150)
    gcsave<-ggroup(container=choicesave,expand=FALSE)
    
    lyt<-glayout(container=window,spacing=13)
    
    importdata<-gbutton("Input Dataset",container=lyt)
    parset<-gbutton("Parameter Setting",container=lyt)
    manhattan<-gbutton("Manhattan Plot",container=lyt)
    qqplot<-gbutton("QQ Plot",container=lyt)
    
    helpfile<-gbutton("User Manual",container=lyt)
    savefile<-gbutton(" Save ",container=lyt)
    run<-gbutton("Run",container=lyt)
    exit<-gbutton("Exit",container=lyt)
    gwline<-glabel("Critical value for Manhattan Plot",container=lyt)
    gwedit<-gedit("3",width=20,coerce.with=as.numeric,container=lyt)
    svgwline<-svalue(gwedit)
    gwstandp<-glabel("Critical P-value for QQ plot",container=lyt)
    gwedit1<-gedit("0.992",width=20,coerce.with=as.numeric,container=lyt)
    svgwstandp<-svalue(gwedit1)
    
    lyt[1,1]<-importdata
    lyt[4,1]<-parset
    lyt[5,1]<-run
    lyt[6,1]<-savefile
    lyt[9,1]<-gwline
    lyt[10,1]<-gwedit
    lyt[11,1]<-manhattan
    lyt[12,1]<-gwstandp
    lyt[13,1]<-gwedit1
    lyt[14,1]<-qqplot
    lyt[17,1]<-helpfile
    lyt[18,1]<-exit
    
    nb1<-gnotebook(tab.pos=3,closebuttons=TRUE,dontCloseThese=TRUE,container=lyt,expand=TRUE)
    size(nb1)<-c(680,540)
    tb<-gnewtable("     
                  1. mrMLM is a R software package for genome-wide association studies based on a multi-locus random-SNP-effect mixed linear model.
                  
                  2. Please cite: Wang Shi-Bo, Feng Jian-Ying, Ren Wen-Long, Huang Bo, Zhou Ling, Wen Yang-Jun, Zhang Jin, Jim M. Dunwell, Xu Shizhong (*), Zhang Yuan-Ming (*).
                  2016. Improving power and accuracy of genome-wide association studies via a multi-locus mixed linear model methodology.Scientific Reports 6: 19444. 
                  
                  3. The software package is developed by Wen-Long Ren, Shi-Bo Wang, Bo Huang & Yuan-Ming Zhang.
                  
                  
                  Version 1.3, Realeased July 2016",multiple=TRUE,container=nb1,expand=TRUE,label="About the software")
    
    
    lyt[1:20,2,expand=TRUE]<-nb1
    
    staprogress<-gtkButton()
    lyt[21,2,expand=TRUE]<-staprogress
    
    
    addHandlerClicked(importdata,handler=function(h,...){
      if(isExtant(importwin)==FALSE)
      {
        importwin<-gwindow("Input Dataset",visible=FALSE,width=250,height=420)
        gimpwin<-ggroup(container=importwin,expand=FALSE)
      }
      lytimp<-glayout(container=gimpwin,spacing=13)
      impchoose<-glabel("1. Choose dataset format",container=lytimp)
      impfile1<-glabel("2. Input Genotypic and Phenotypic files",container=lytimp)
      impprepare<-glabel("3. Sort & Transform for dataset",container=lytimp)
      impfile2<-glabel("4. Input Kinship and Population-structure files",container=lytimp)
      radioimp<-gradio(c("mrMLM numeric format","mrMLM character format","Hapmap (TASSEL) format"),selected=3,horizontal=FALSE,container=lytimp)
      genotype<-gbutton("Genotype",container=lytimp)
      phenotype<-gbutton("Phenotype",container=lytimp)
      kinship<-gbutton("Kinship",container=lytimp)
      population<-gbutton("Population Structure",container=lytimp)
      preimp<-gbutton("Do",container=lytimp)
      lytimp[1,2:5]<-impchoose
      lytimp[2:4,2:5]<-radioimp
      lytimp[5,2:5]<-impfile1
      lytimp[6,2:4]<-genotype
      lytimp[7,2:4]<-phenotype
      lytimp[8,2:5]<-impprepare
      lytimp[9,2:4]<-preimp
      lytimp[10,2:5]<-impfile2
      lytimp[11,2:4]<-kinship
      lytimp[12,2:4]<-population
      visible(importwin)<-TRUE
      
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
          mrenv$genRaw<-as.matrix(read.csv(input1,header=FALSE))
          showgenRaw<-mrenv$genRaw[-1,]
          colnames(showgenRaw)<-mrenv$genRaw[1,]
          showgenRaw<-as.data.frame(showgenRaw)
          tbdfe1<-gdfedit(showgenRaw,container=nb1,expand=TRUE,label="Raw_Genotype") 
        }
      })
      
      addHandlerClicked(phenotype,handler=function(h,...){
        input2<-gfile(text="Select a file...",type="open",
                      filter=list("All files"=list(patterns=c("*")),
                                  "CSV files"=list(patterns=c("*.csv"))))
        if(is.na(input2))
        {
          gmessage("Please input correct phenotype data !","Warning",icon="warning")
          return
        }else{
          mrenv$pheRaw<-as.matrix(read.csv(input2,header=FALSE)) 
          showpheRaw<-mrenv$pheRaw[-1,]
          colnames(showpheRaw)<-c(mrenv$pheRaw[1,1],"   ")
          showpheRaw<-as.data.frame(showpheRaw)
          tbdfe2<-gdfedit(showpheRaw,container=nb1,expand=TRUE,label="Raw_Phenotype")
        }
      })
      
      addHandlerClicked(preimp,handler=function(h,...){
        if(svalue(radioimp)=="mrMLM numeric format"){
          mrenv$inputform<-1
          nameGen <- as.matrix(mrenv$genRaw[1,],1,)
          namePhe <- as.matrix(mrenv$pheRaw[,1],,1)
          mrenv$sameName <- intersect(nameGen,namePhe)
          ##########To find the location of the same name 
          locGen <- match(mrenv$sameName,nameGen)
          locPhe <- match(mrenv$sameName,namePhe)
          ##########Produce new genotype matrix and phenotype matrix
          hapName <- matrix(c("rs#","chrom","pos","genotype for code 1"),1,)
          hapHave <- intersect(nameGen,hapName)
          locHap <- match(hapHave,nameGen)
          newGenloc <- c(locHap,locGen)
          newPheloc <- locPhe
          newGen <- as.matrix(mrenv$genRaw[-1,newGenloc])
          newPhe <- as.matrix(mrenv$pheRaw[newPheloc,])
          nnhap <- length(hapHave)
          rownewGen <- dim(newGen)[1]
          colnewGen <- dim(newGen)[2]
          rownewPhe <- dim(newPhe)[1]
          ###########To show on the table ----newGen
          mrenv$newGen <-rbind(mrenv$genRaw[1,newGenloc],newGen)
          ###########To be computed ----gen
          locChr <- as.numeric(which(mrenv$newGen[1,]=="chrom"))
          locPos <- as.numeric(which(mrenv$newGen[1,]=="pos"))
          mrenv$needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
          mrenv$needGen <- mrenv$newGen[,mrenv$needloc]
          mrenv$gen<-as.matrix(mrenv$needGen[-1,])
          mrenv$gen<-matrix(as.numeric(mrenv$gen),nrow=nrow(mrenv$gen))
          ###########To show on the table ----newPhe
          mrenv$pheRaw[1,2]<-"  "
          mrenv$newPhe<-rbind(mrenv$pheRaw[1,],newPhe)
          ###########To be computed ----phe
          mrenv$phe<-as.matrix(mrenv$newPhe[-1,-1])
          mrenv$phe<-matrix(as.numeric(mrenv$phe),nrow=nrow(mrenv$phe))
          shownewGen<-mrenv$newGen[-1,]
          colnames(shownewGen)<-mrenv$newGen[1,]
          shownewGen<-as.data.frame(shownewGen)
          shownewPhe<-mrenv$newPhe[-1,]
          colnames(shownewPhe)<-c(mrenv$newPhe[1,1],"   ")
          shownewPhe<-as.data.frame(shownewPhe)
          tbdfe3<-gdfedit(shownewGen,container=nb1,expand=TRUE,label="Genotype")
          tbdfe4<-gdfedit(shownewPhe,container=nb1,expand=TRUE,label="Phenotype")
        }else if(svalue(radioimp)=="mrMLM character format"){
          mrenv$inputform<-2
          ##########To find the same name between genotype and phenotype
          nameGen <- as.matrix(mrenv$genRaw[1,],1,)
          namePhe <- as.matrix(mrenv$pheRaw[,1],,1)
          mrenv$sameName <- intersect(nameGen,namePhe)
          ##########To find the location of the same name 
          locGen <- match(mrenv$sameName,nameGen)
          locPhe <- match(mrenv$sameName,namePhe)
          ##########Produce new genotype matrix and phenotype matrix
          hapName <- matrix(c("rs#","chrom","pos"),1,)
          hapHave <- intersect(nameGen,hapName)
          locHap <- match(hapHave,nameGen)
          newGenloc <- c(locHap,locGen)
          newPheloc <- locPhe
          newGen <- as.matrix(mrenv$genRaw[-1,newGenloc])
          newPhe <- as.matrix(mrenv$pheRaw[newPheloc,])
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
            computeGen[i,notRef] <- as.numeric(-1)
            computeGen[i,notATCG] <- as.numeric(0)
          }
          mrenv$outATCG<-referSam
          ###########To show on the table ----newGen
          newGen <- cbind(newGen[,1:nnhap],computeGen)
          mrenv$newGen <-rbind(mrenv$genRaw[1,newGenloc],newGen)
          ###########To be computed ----gen
          locChr <- as.numeric(which(mrenv$newGen[1,]=="chrom"))
          locPos <- as.numeric(which(mrenv$newGen[1,]=="pos"))
          mrenv$needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
          mrenv$needGen<-mrenv$newGen[,mrenv$needloc]
          mrenv$gen<-as.matrix(mrenv$needGen[-1,])
          mrenv$gen<-matrix(as.numeric(mrenv$gen),nrow=nrow(mrenv$gen))
          ###########To show on the table ----newPhe
          mrenv$pheRaw[1,2]<-"  "
          mrenv$newPhe<-rbind(mrenv$pheRaw[1,],newPhe)
          ###########To be computed ----phe
          mrenv$phe<-as.matrix(mrenv$newPhe[-1,-1])
          mrenv$phe<-matrix(as.numeric(mrenv$phe),nrow=nrow(mrenv$phe))
          shownewGen<-mrenv$newGen[-1,]
          colnames(shownewGen)<-mrenv$newGen[1,]
          shownewGen<-as.data.frame(shownewGen)
          shownewPhe<-mrenv$newPhe[-1,]
          colnames(shownewPhe)<-c(mrenv$newPhe[1,1],"   ")
          shownewPhe<-as.data.frame(shownewPhe)
          tbdfe3<-gdfedit(shownewGen,container=nb1,expand=TRUE,label="Genotype")
          tbdfe4<-gdfedit(shownewPhe,container=nb1,expand=TRUE,label="Phenotype")
        }else if(svalue(radioimp)=="Hapmap (TASSEL) format"){
          mrenv$inputform<-3
          ##########To find the same name between genotype and phenotype
          nameGen<-as.matrix(mrenv$genRaw[1,],1,)
          namePhe<-as.matrix(mrenv$pheRaw[,1],,1)
          mrenv$sameName<-intersect(nameGen,namePhe)
          ##########To find the location of the same name 
          locGen<-match(mrenv$sameName,nameGen)
          locPhe<-match(mrenv$sameName,namePhe)
          ##########Produce new genotype matrix and phenotype matrix
          hapName<-matrix(c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panel","QCcode"),1,)
          hapHave<-intersect(nameGen,hapName)
          locHap<-match(hapHave,nameGen)
          newGenloc<-c(locHap,locGen)
          newPheloc<-locPhe
          newGen<-as.matrix(mrenv$genRaw[-1,newGenloc])
          newPhe<-as.matrix(mrenv$pheRaw[newPheloc,])   
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
            computeGen[i,notRef]<-as.numeric(-1)
            computeGen[i,notATCG]<-as.numeric(0)
          }
          mrenv$outATCG<-referSam
          ###########To show on the table ----mrenv$newGen
          newGen<-cbind(newGen[,1:nnhap],computeGen)
          mrenv$newGen<-rbind(mrenv$genRaw[1,newGenloc],newGen)
          ###########To be computed ----mrenv$gen
          locChr<-as.numeric(which(mrenv$newGen[1,]=="chrom"))
          locPos<-as.numeric(which(mrenv$newGen[1,]=="pos"))
          mrenv$needloc<-c(locChr,locPos,(nnhap+1):colnewGen)
          mrenv$needGen<-mrenv$newGen[,mrenv$needloc]
          mrenv$gen<-as.matrix(mrenv$needGen[-1,])
          mrenv$gen<-matrix(as.numeric(mrenv$gen),nrow=nrow(mrenv$gen))
          ###########To show on the table ----mrenv$newPhe
          mrenv$pheRaw[1,2]<-"  "
          mrenv$newPhe<-rbind(mrenv$pheRaw[1,],newPhe)
          ###########To be computed ----mrenv$phe
          mrenv$phe<-as.matrix(mrenv$newPhe[-1,-1])
          mrenv$phe<-matrix(as.numeric(mrenv$phe),nrow=nrow(mrenv$phe))
          shownewGen<-mrenv$newGen[-1,]
          colnames(shownewGen)<-mrenv$newGen[1,]
          shownewGen<-as.data.frame(shownewGen)
          shownewPhe<-mrenv$newPhe[-1,]
          colnames(shownewPhe)<-c(mrenv$newPhe[1,1],"   ")
          shownewPhe<-as.data.frame(shownewPhe)
          tbdfe3<-gdfedit(shownewGen,container=nb1,expand=TRUE,label="Genotype")
          tbdfe4<-gdfedit(shownewPhe,container=nb1,expand=TRUE,label="Phenotype")
        }
      })
      
      addHandlerClicked(kinship,handler=function(h,...){
        if(isExtant(choicekk)==FALSE)
        {
          choicekk<-gwindow("Choose Kinship option",visible=FALSE,width=320,height=150)
          gkk<-ggroup(container=choicekk,expand=FALSE)
        }
        lytkk<-glayout(container=gkk,spacing=13)
        mrenv$okkk<-gbutton("     OK    ",container=lytkk)
        mrenv$cancelkk<-gbutton(" Cancel ",container=lytkk)
        mrenv$radiokk<-gradio(c("Input the Kinship matrix file","Calculate the Kinship matrix by this software"),selected=1,horizontal=FALSE,container=lytkk)
        lytkk[2:3,2:5]<-mrenv$radiokk
        lytkk[5,2]<-mrenv$okkk
        lytkk[5,5]<-mrenv$cancelkk
        visible(choicekk)<-TRUE
        addHandlerClicked(mrenv$okkk,handler=function(h,...){
          if(svalue(mrenv$radiokk)=="Input the Kinship matrix file"){
            input3<-gfile(text="Select a file...",type="open",
                          filter=list("All files"=list(patterns=c("*")),
                                      "CSV files"=list(patterns=c("*.csv"))))
            if(is.na(input3))
            {
              gmessage("Please input correct kinship data !","Warning",icon="warning")
              return
            }else{
              mrenv$kkRaw<-read.csv(input3,header=FALSE)
              nnkk<-dim(mrenv$kkRaw)[1]
              mrenv$kkRaw[1,2:nnkk]<-"  "
              tbdfe5<-gdfedit(mrenv$kkRaw,container=nb1,expand=TRUE,label="Kinship")
              kkPre<-as.matrix(mrenv$kkRaw[-1,-1])
              nameKin<-as.matrix(mrenv$kkRaw[-1,1])
              sameGenKin<-intersect(mrenv$sameName,nameKin)
              locKin<-match(sameGenKin,nameKin)
              mrenv$kk<-kkPre[locKin,locKin]
              mrenv$kk<-matrix(as.numeric(mrenv$kk),nrow=nrow(mrenv$kk))
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
              mrenv$kk<-as.matrix(kk1)
              rowsize<-dim(mrenv$kk)[1]
              aa<-as.character()
              for(i in 1:(rowsize+1))
              {
                a<-paste("V",i,sep="")
                aa<-c(aa,a)
              }
              mrenv$kkShow<-cbind(matrix(mrenv$newPhe[-1,1],,1),round(mrenv$kk,5))
              tempFirst<-rep("  ",rowsize)
              tempFirst<-c(as.character(rowsize),tempFirst)
              mrenv$kkShow<-as.matrix(rbind(tempFirst,mrenv$kkShow))
              colnames(mrenv$kkShow)<-aa
              rownames(mrenv$kkShow)<-c(1:(rowsize+1))
              tbdfe5<-gdfedit(mrenv$kkShow,container=nb1,expand=TRUE,label="Kinship")
              dispose(choicekk)
            }     
          }
        })
        addHandlerClicked(mrenv$cancelkk,handler=function(h,...){
          dispose(choicekk)
        })
      })
      
      
      addHandlerClicked(population,handler=function(h,...){
        if(isExtant(includeps)==FALSE)
        {
          includeps<-gwindow("Include population structure?",visible=FALSE,width=320,height=150)
          gps<-ggroup(container=includeps,expand=FALSE)
        }
        lytps<-glayout(container=gps,spacing=13)
        okps<-gbutton("     OK    ",container=lytps)
        cancelps<-gbutton(" Cancel ",container=lytps)
        radiops<-gradio(c("Include","No"),selected=2,horizontal=FALSE,container=lytps)
        lytps[2:3,2:5]<-radiops
        lytps[5,2]<-okps
        lytps[5,5]<-cancelps
        visible(includeps)<-TRUE
        addHandlerClicked(okps,handler=function(h,...){
          if(svalue(radiops)=="Include"){
            mrenv$flagps<-0
            input4<-gfile(text="Select a file...",type="open",
                          filter=list("All files"=list(patterns=c("*")),
                                      "CSV files"=list(patterns=c("*.csv"))))
            if(is.na(input4))
            {
              gmessage("Please input correct population data !","Warning",icon="warning")
              return
            }else{
              mrenv$psmatrixRaw<-as.matrix(read.csv(input4,header=FALSE))
              nnpprow<-dim(mrenv$psmatrixRaw)[1]
              nnppcol<-dim(mrenv$psmatrixRaw)[2]
              mrenv$psmatrixRaw[1,2:nnppcol]<-"  "
              psmatrixPre<-mrenv$psmatrixRaw[3:nnpprow,]
              namePop<-as.matrix(psmatrixPre[,1])
              sameGenPop<-intersect(mrenv$sameName,namePop)
              locPop<-match(sameGenPop,namePop)
              mrenv$psmatrix<-psmatrixPre[locPop,-1:-2]
              mrenv$psmatrix<-matrix(as.numeric(mrenv$psmatrix),nrow=nrow(mrenv$psmatrix))
              tbdfe6<-gdfedit(mrenv$psmatrixRaw,container=nb1,expand=TRUE,label="Population Structure")
              dispose(includeps)
              dispose(importwin)
            }
          }else{
            mrenv$flagps<-1
            enabled(population)<-FALSE
            dispose(includeps)
            dispose(importwin)
          }
        })
        addHandlerClicked(cancelps,handler=function(h,...){
          dispose(includeps)
        })
      })
      
    })
    
    addHandlerClicked(parset,handler=function(h,...){
      if(isExtant(parsetwin)==FALSE)
      {
        parsetwin<-gwindow("Parameter Setting",visible=FALSE,width=260,height=280)
        gpar<-ggroup(container=parsetwin,expand=FALSE)
      }
      lytpar<-glayout(container=gpar,spacing=13)
      mrenv$pvallabel<-glabel("1. Critical P-value in rMLM:",container=lytpar)
      mrenv$pvaledit<-gedit("0.01",width=20,coerce.with=as.numeric,container=lytpar)
      mrenv$radlabel<-glabel("2. Search radius of candidate gene (kb):",container=lytpar)
      mrenv$radedit<-gedit("20",width=20,coerce.with=as.numeric,container=lytpar)
      mrenv$mlodlabel<-glabel("3. Critical LOD score in mrMLM:",container=lytpar)
      mrenv$mlodedit<-gedit("3",width=20,coerce.with=as.numeric,container=lytpar)
      mrenv$okpar<-gbutton("     OK    ",container=lytpar)
      mrenv$cancelpar<-gbutton(" Cancel ",container=lytpar)
      lytpar[1,1:5]<-mrenv$pvallabel
      lytpar[2,1:5]<-mrenv$pvaledit
      lytpar[3,1:5]<-mrenv$radlabel
      lytpar[4,1:5]<-mrenv$radedit
      lytpar[5,1:5]<-mrenv$mlodlabel
      lytpar[6,1:5]<-mrenv$mlodedit
      lytpar[7,1]<-mrenv$okpar
      lytpar[7,4]<-mrenv$cancelpar
      visible(parsetwin)<-TRUE
      addHandlerClicked(mrenv$okpar,handler=function(h,...){
        mrenv$svpal<-svalue(mrenv$pvaledit)
        mrenv$svrad<-svalue(mrenv$radedit)
        mrenv$svmlod<-svalue(mrenv$mlodedit)
        if((mrenv$svpal<0)||(mrenv$svpal>1))
        {
          gmessage("Please input critical P-value more than 0 and less than 1!","Warning",icon="warning")
          return
        }
        if(mrenv$svrad<0)
        {
          gmessage("Please input search radius of candidate gene more than 0!","Warning",icon="warning")
          return
        }
        if(mrenv$svmlod<0)
        {
          gmessage("Please input critical LOD score more than 0!","Warning",icon="warning")
          return
        }
        if((mrenv$svpal>0)&&(mrenv$svpal<1)&&(mrenv$svrad>=0)&&(mrenv$svmlod>=0))
        {
          dispose(parsetwin)
        }
      })
      
      addHandlerClicked(mrenv$cancelpar,handler=function(h,...){
        mrenv$svpal<-svalue(mrenv$pvaledit)
        mrenv$svrad<-svalue(mrenv$radedit)
        mrenv$svmlod<-svalue(mrenv$mlodedit)
        dispose(parsetwin)
      })
    })
    
    addHandlerClicked(helpfile,handler=function(h,...){
      RShowDoc("Instruction.pdf",package="mrMLM")  
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
              #ij<-which(sub!=sub[i+1])
              ij<-which(sub!=sub[i+ncol(xxn)])
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
            wang[i]<-p
          }
          return (wang)
        }
        
        gen<-t(gen) 
        chr_pos<-mrenv$parms[,2:3]
        pfit<-which(mrenv$parms[,10]<=(mrenv$svpal))
        pfit<-as.matrix(pfit)
        pfitrow<-nrow(pfit)
        no_p<-cbind((1:(nrow(mrenv$parms))),mrenv$parms[,10])
        no_porder<-order(no_p[,2])
        no_p<-no_p[no_porder,]
        choose_orderp<-no_p[1:pfitrow,]
        orderno<-no_p[1:pfitrow,1]
        orderno<-as.matrix(orderno)
        sigma2g_SNPerr<-cbind(mrenv$parms[,6],mrenv$parms[,8])
        correct_each<-matrix(1,(nrow(sigma2g_SNPerr)),1)-sigma2g_SNPerr[,2]*sigma2g_SNPerr[,2]/sigma2g_SNPerr[,1]
        k0<-which(correct_each<0)
        k0<-as.matrix(k0)
        if(nrow(k0)>0){
          correct_each[k0,1]<-matrix(0,(nrow(k0)),1)
        }
        correct_sum<-sum(correct_each)
        newp<-0.05/correct_sum
        mrenv$mannewp<-newp
        mrenv$manstandchoice<-1
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
              if (ye<=((mrenv$svrad)*1000)){
                gg[jj,1]<-0
              }
            }
          }
        }
        progress_bar$setFraction(95/100)
        if(mrenv$inputform==1){
          #output result1 using mrMLM numeric format
          mrenv$parmsShow<-mrenv$parms[,-1]
          meadd<-matrix(1,nrow(mrenv$parms),1)
          meadd[which(mrenv$parms[,10]<newp),1]<-round(newp,10)
          meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
          mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],mrenv$parms[,2:3],mrenv$parms[,4:10],mrenv$genRaw[-1,4],meadd)
          colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype for code 1","Significance")
          tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
        }
        if(mrenv$inputform==2){
          #output result1 using mrMLM character format
          mrenv$parmsShow<-mrenv$parms[,-1]
          mrenv$outATCG<-matrix(mrenv$outATCG,,1)
          meadd<-matrix(1,nrow(mrenv$parms),1)
          meadd[which(mrenv$parms[,10]<newp),1]<-round(newp,10)
          meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
          mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],mrenv$parms[,2:3],mrenv$parms[,4:10],mrenv$outATCG,meadd)
          colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype  for code 1","Significance")
          tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
        }
        if(mrenv$inputform==3){
          #output result1 using TASSEL format
          mrenv$parmsShow<-mrenv$parms[,-1]
          mrenv$outATCG<-matrix(mrenv$outATCG,,1)
          mrenv$outATCG<-unlist(strsplit(mrenv$outATCG,""))
          mrenv$outATCG<-matrix(mrenv$outATCG[c(TRUE,FALSE)],,1)
          meadd<-matrix(1,nrow(mrenv$parms),1)
          meadd[which(mrenv$parms[,10]<newp),1]<-round(newp,10)
          meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
          mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],mrenv$parms[,2:3],mrenv$parms[,4:10],mrenv$outATCG,meadd)
          colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype  for code 1","Significance")
          tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
        }
        
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
        par<-ebayes_EM(xin,xx,phe)
        w2<-which(par[,1]<=0.01)
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
            lod<-likelihood(lodfix,lodrand,phe)
          }else if(flagps==0)
          {
            temp<-cbind(psmatrix,lodfix)
            lod<-likelihood(temp,lodrand,phe)
          }
          w3<-which(lod[,1]>=(mrenv$svmlod))
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
          #if (length(kk)==1){ww<-0}
        }
        if ((nrow(w2))==0){
          g0<-as.matrix(g0)
          lo<-as.matrix(lo)
          yang<-which(lo>=(mrenv$svmlod))
          yang<-as.matrix(yang)
          if ((nrow(yang))>0){
            ww<-g0[yang,1]
            lo<-lo[yang,1]
          }
          if ((nrow(yang))==0){ww<-0}
        }
        ww<-as.matrix(ww)
        mrenv$needww<-ww
        if (length(ww)>1){
          if((flagps==1)||(exists("psmatrix")==FALSE))
          {
            ex<-cbind(matrix(1,(nrow(xx)),1),t(gen[ww,]))
          }else if(flagps==0)
          {
            ex<-cbind(cbind(matrix(1,(nrow(xx)),1),psmatrix),t(gen[ww,]))
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
            her<-er/as.numeric(sum(er)+v0)
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
            her<-er/as.numeric(sum(er)+v0)
          }
          
          mrenv$wan<-data.frame(mrenv$parmsShow[mrenv$needww,1],chr_pos[ww,],eeff,lo,her,mrenv$parmsShow[mrenv$needww,11])
          colnames(mrenv$wan)<-c("RS#","Chromosome","Position","QTN effect","LOD score","R2","Genotype  for code 1")
          
        }
        wan<-mrenv$wan
        if(exists("wan")==FALSE||is.null(wan)==TRUE)
        {
          gmessage("There is no result meets the requirements in the second step!","Info",icon="info")
        }else{
          tbdfe8<-gdfedit(wan,container=nb1,expand=TRUE,label="Result2")
        }
        progress_bar$setFraction(100/100)
        progress_bar$setText("All done.")
      }
      return
    })
    
    addHandlerClicked(savefile,handler=function(h,...){
      if(isExtant(choicesave)==FALSE)
      {
        choicesave<-gwindow("Save as ...",visible=FALSE,width=250,height=150)
        gcsave<-ggroup(container=choicesave,expand=FALSE)
      }
      lytsavere<-glayout(container=gcsave,spacing=13)
      mrenv$oksa<-gbutton("     OK    ",container=lytsavere)
      mrenv$cancelsa<-gbutton(" Cancel ",container=lytsavere)
      mrenv$radiosa<-gradio(c("Result1","Result2"),selected=1,horizontal=FALSE,container=lytsavere)
      lytsavere[2:3,2:5]<-mrenv$radiosa
      lytsavere[5,3]<-mrenv$oksa
      lytsavere[5,5]<-mrenv$cancelsa
      visible(choicesave)<-TRUE
      
      addHandlerClicked(mrenv$oksa,handler=function(h,...){
        if(svalue(mrenv$radiosa)=="Result1"){
          parms<-mrenv$parmsShow
          if(exists("parms")==FALSE||is.null(parms)==TRUE)
          {
            gmessage("There is something wrong in the first step!","Info",icon="info")
            return
          }else{
            output<-gfile(text="Save a file...",type="save",
                          filter=list("All files"=list(patterns=c("*")),
                                      "CSV files"=list(patterns=c("*.csv"))))
            write.table(parms,output,sep = ",",row.names=FALSE,col.names = TRUE) 
          }
        }else{
          wan<-mrenv$wan
          if(exists("wan")==FALSE||is.null(wan)==TRUE)
          {
            gmessage("There is no result meets the requirements in the second step!","Info",icon="info")
            return
          }else{
            output<-gfile(text="Save a file...",type="save",
                          filter=list("All files"=list(patterns=c("*")),
                                      "CSV files"=list(patterns=c("*.csv"))))
            write.table(wan,output,sep = ",",row.names=FALSE,col.names = TRUE) 
          }
        }
      })
      
      addHandlerClicked(mrenv$cancelsa,handler=function(h,...){
        dispose(choicesave)
      })
    })  
    
    addHandlerClicked(manhattan,handler=function(h,...){
      if((exists("svgwline")==FALSE)||(svalue(gwedit)<=0))
      {
        gmessage("Please input correct genomewideline value!","Warning",icon="warning")
        return
      }else{
        svgwline<-svalue(gwedit)
        if(mrenv$manstandchoice==1)
        {
          mrenv$standline<--log10(mrenv$mannewp)
          mrenv$manstandchoice<-mrenv$manstandchoice+1
        }else{
          mrenv$standline<-svgwline
        }
        if(isExtant(plotwin)==FALSE)
        {
          plotwin<-gwindow("Manhattan Plot",visible=FALSE,width=600,height=360)
          gpw<-ggroup(container=plotwin,horizontal = FALSE)
          ggpw<-ggraphics(container=gpw)
        }
        
        plotlyt <- glayout(container = gpw)
        plotlabel1 <- glabel("Width (px):", container = plotlyt)
        plotedit1 <- gedit("960", container = plotlyt)
        widvalue <- as.numeric(svalue(plotedit1))
        plotlabel2 <- glabel("Height (px):", container = plotlyt)
        plotedit2 <- gedit("600", container = plotlyt)
        heightvalue <- as.numeric(svalue(plotedit2))
        plotlabel3 <- glabel("Point size (1/72 inch, ppi):", container = plotlyt)
        plotedit3 <- gedit("12", container = plotlyt)
        pointsizevalue <- as.numeric(svalue(plotedit3))
        plotbt <- gbutton(" Save ", container = plotlyt)
        plotmancl <- gbutton("  Cancel  ", container = plotlyt)
        combo_box1 <- gcombobox(selected=1,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotlyt)
        combo_box2 <- gcombobox(selected=3,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotlyt)
        plotlabel4 <- glabel("Chromosome color (odd):", container = plotlyt)
        plotlabel5 <- glabel("Chromosome color (even):", container = plotlyt)
        plotlyt[2, 1] <- plotlabel1
        plotlyt[2, 2] <- plotedit1
        plotlyt[3, 1] <- plotlabel2
        plotlyt[3, 2] <- plotedit2
        plotlyt[4, 1] <- plotlabel3
        plotlyt[4, 2] <- plotedit3
        plotlyt[2, 4] <- plotlabel4
        plotlyt[2, 5] <- combo_box1
        plotlyt[3, 4] <- plotlabel5
        plotlyt[3, 5] <- combo_box2
        plotlyt[4, 4] <- plotbt
        plotlyt[4, 5] <- plotmancl
        
        visible(plotwin)<-TRUE
        
        rowsnp <- dim(mrenv$parms)[1]
        snpname <- numeric()
        for(i in 1:rowsnp)
        {
          snpname <- rbind(snpname,paste("rs",i,sep=""))
        }
        
        bpnumber <- numeric()
        chrnum <- unique(mrenv$parms[,2])
        
        for(i in 1:length(chrnum))
        {
          bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(mrenv$parms[,2]==chrnum[i])))))
        }
        
        addHandlerChanged(ggpw, handler=function(h,...) {
          color1 <- svalue(combo_box1)
          color2 <- svalue(combo_box2)
          
          parms <- data.frame(mrenv$parms,snpname,bpnumber)
          colnames(parms)<-c("Trait","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","SNPname","BPnumber")
          manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P_wald",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = mrenv$standline)
        })
        addhandlerclicked(plotbt, handler = function(h, ...) {
          widvalue <- as.numeric(svalue(plotedit1))
          heightvalue <- as.numeric(svalue(plotedit2))
          pointsizevalue <- as.numeric(svalue(plotedit3))
          color1 <- svalue(combo_box1)
          color2 <- svalue(combo_box2)
          
          output <- gfile(text = "Save a file...", type = "save", 
                          filter = list(`All files` = list(patterns = c("*")), 
                                        `TIFF files` = list(patterns = c("*.tiff")),
                                        `PNG files` = list(patterns = c("*.png")),
                                        `JPEG files` = list(patterns = c("*.jpeg"))))
          if((length(grep(".png",output))==1)||(length(grep(".PNG",output))==1)){
            png(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue)
          }else if((length(grep(".tiff",output))==1)||(length(grep(".TIFF",output))==1)){
            tiff(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue)
          }else if((length(grep(".jpeg",output)==1))||(length(grep(".JPEG",output))==1)){
            jpeg(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue)  
          }else{
            gmessage("Please input correct image format !")
          }
          parms <- data.frame(mrenv$parms,snpname,bpnumber)
          colnames(parms)<-c("Trait","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","SNPname","BPnumber")
          manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P_wald",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = mrenv$standline)
          dev.off()
          
        })
        addHandlerClicked(plotmancl,handler=function(h,...){
          dispose(plotwin)
        })
      }
    })
    
    qqplotfun <- function(p_value,p_stand){
      pvalue<-matrix(p_value,,1)
      observed<-sort(pvalue[,1])
      observed<-observed/2
      observed<-observed[which(observed!=0)]
      newobserved<-observed[which(observed<(p_stand/2))]
      lobs<--(log10(newobserved))
      expected<-c(1:length(newobserved))
      lexp<--(log10(expected/(length(pvalue)+1)))
      plot(lexp,lobs,xlim=c(0,max(lexp)),ylim=c(0,max(lobs)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'))
      abline(0,1,col="red")
    }
    
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
          gpw1<-ggroup(container=plotwin1,horizontal=FALSE)
          ggpw1<-ggraphics(container=gpw1)
        }
        plotqqlyt <- glayout(container = gpw1)
        plotqqlabel1 <- glabel("Width (px):", container = plotqqlyt)
        plotqqedit1 <- gedit("960", container = plotqqlyt)
        widqqvalue <- as.numeric(svalue(plotqqedit1))
        plotqqlabel2 <- glabel("Height (px):", container = plotqqlyt)
        plotqqedit2 <- gedit("600", container = plotqqlyt)
        heightqqvalue <- as.numeric(svalue(plotqqedit2))
        plotqqlabel3 <- glabel("Point size (1/72 inch, ppi):", container = plotqqlyt)
        plotqqedit3 <- gedit("12", container = plotqqlyt)
        pointsizeqqvalue <- as.numeric(svalue(plotqqedit3))
        plotqqbt <- gbutton(" Save ", container = plotqqlyt)
        plotqqcl <- gbutton("  Cancel  ", container = plotqqlyt)
        
        plotqqlyt[2, 1] <- plotqqlabel1
        plotqqlyt[2, 2] <- plotqqedit1
        plotqqlyt[3, 1] <- plotqqlabel2
        plotqqlyt[3, 2] <- plotqqedit2
        plotqqlyt[4, 1] <- plotqqlabel3
        plotqqlyt[4, 2] <- plotqqedit3
        plotqqlyt[3, 3] <- plotqqbt
        plotqqlyt[4, 3] <- plotqqcl 
        visible(plotwin1)<-TRUE
        
        addHandlerChanged(ggpw1, handler=function(h,...) {
          qqplotfun(mrenv$parms[,10],mrenv$standp)
        })
        
        addhandlerclicked(plotqqbt, handler = function(h, ...) {
          widqqvalue <- as.numeric(svalue(plotqqedit1))
          heightqqvalue <- as.numeric(svalue(plotqqedit2))
          pointsizeqqvalue <- as.numeric(svalue(plotqqedit3))
          
          output <- gfile(text = "Save a file...", type = "save", 
                          filter = list(`All files` = list(patterns = c("*")), 
                                        `TIFF files` = list(patterns = c("*.tiff")),
                                        `PNG files` = list(patterns = c("*.png")),
                                        `JPEG files` = list(patterns = c("*.jpeg"))))
          if((length(grep(".png",output))==1)||(length(grep(".PNG",output))==1)){
            png(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue)
          }else if((length(grep(".tiff",output))==1)||(length(grep(".TIFF",output))==1)){
            tiff(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue)
          }else if((length(grep(".jpeg",output)==1))||(length(grep(".JPEG",output))==1)){
            jpeg(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue)  
          }else{
            gmessage("Please input correct image format !")
          }
          qqplotfun(mrenv$parms[,10],mrenv$standp)
          dev.off()
          
        })
        
        addHandlerClicked(plotqqcl,handler=function(h,...){
          dispose(plotwin1)
        })
        
      }
    })
  }
  
  addHandlerClicked(mrMLMbutton,handler=function(h,...){
    
    mrMLMSub()
    
  })
  
  
  
  GCIM<-function(){
    
    envim<-new.env()
    Sys.setenv(LANGUAGE="en")
    
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
    
    
    
    window<-gwindow("Genome-wide Composite Interval Mapping (GCIM)",visible=TRUE,width=1240,height=650,expand=TRUE)
    wgroup<-ggroup(container = window)
    
    
    icimw<-gwindow("QTLIciMaping Format",visible = FALSE,width =80,height = 140,expand=FALSE)
    icimg<-ggroup(container = icimw,expand=FALSE)
    cartw<-gwindow("WinQTLCart  Format",visible = FALSE,width = 80,height = 120,expand=FALSE)
    cartg<-ggroup(container = cartw)
    gcimw<-gwindow("GCIM Format",visible = FALSE,width = 85,height = 160,expand=FALSE)
    gcimg<-ggroup(container = gcimw)
    
    covw<-gwindow("Include Covariate ?",visible = FALSE,width = 250,height = 130,expand=FALSE)
    covg<-ggroup(container = covw)
    
    plotwin<-gwindow("Plot of LOD Score against Genome Postion",visible=FALSE,width=960,height=220,expand=TRUE)
    gpw<-ggroup(container=plotwin,horizontal = FALSE,spacing = 10)
    
    aboutw<-gwindow("About GCIM",visible=FALSE,width=600,height=360,expand=FALSE)
    aboutg<-ggroup(container=aboutw)
    
    
    lyt<-glayout(container =wgroup,spacing =8)
    frame<-gframe(container = lyt)
    framlyt<-glayout(container = frame,spacing = 10)
    
    framlabel_trait<-glabel("   Trait  Selection:",container = framlyt)
    framedit_trait<-gedit(" 1 ",coerce.with=as.numeric,container = framlyt)
    envim$ii<-svalue(framedit_trait)
    framlabel_inter<-glabel("Walk Speed for Genome-wide Scanning (cM):",container = framlyt)
    framedit_inter<-gedit("1",coerce.with=as.numeric,container = framlyt)
    cl<-svalue(framedit_inter)
    framlabel_lod<-glabel("Critical LOD score:",container = framlyt)
    framedit_lod<-gedit("2.5",coerce.with=as.numeric,container = framlyt)
    sLOD<-svalue(framedit_lod)
    covbt<-gbutton(" Covariate ",container =framlyt)
    runbt<-gbutton("          Run             ",container =framlyt)
    clearbt<-gbutton("  Clear       ",container = framlyt)
    exitbt<-gbutton("  Exit      ",container = framlyt)
    framlyt[1,2]<-framlabel_lod
    framlyt[1,3:4]<-framedit_lod
    framlyt[1,5]<-framlabel_inter
    framlyt[1,6:7]<-framedit_inter
    framlyt[2,2]<-framlabel_trait
    framlyt[2,3:4]<-framedit_trait
    framlyt[2,5]<-covbt
    framlyt[2,6:7]<-runbt
    
    framlyt[1, 20:28] <-clearbt 
    framlyt[2, 20:28] <-exitbt
    
    
    
    
    nb1<-gnotebook(tab.pos=3,closebuttons=TRUE,dontCloseThese=TRUE,container=lyt,expand=TRUE)
    
    
    size(nb1)<-c(1200,650)
    
    progress_bar <- gtkProgressBar()
    lyt[1:2,1:20,expand=TRUE]<-frame
    
    lyt[3:20,1:20,expand=TRUE]<-nb1
    lyt[21,1:20,expand=TRUE]<-progress_bar
    
    tb<-gnewtable("     
                  1. GCIM is a R software package for Linkage analysis based on Genome-wide Composite Interval Mapping.
                  
                  2. Please cite: Wang Shi-Bo, Wen Yang-Jun, Ren Wen-Long, Ni Yuan-Li, Zhang Jin, Feng Jian-Ying, Zhang Yuan-Ming.  
                  
                  Mapping small-effect and linked quantitative trait loci for complex traits in backcross or DH populations via
                  
                  a multi-locus GWAS methodology. Scientific Reports 2016, 6: 29951.
                  
                  3. The software package is developed by Yuan-Li Ni,Wen-Long Ren, Shi-Bo Wang & Yuan-Ming Zhang.
                  
                  
                  Version 1.3, Realeased July 2016",multiple=TRUE,container=nb1,expand=TRUE,label="About the software")
    
    
    font(tb)<-c(size="x-large")
    
    addhandlerclicked(exitbt,handler = function(h,...){
      gconfirm("Yes or no ?",handler=function(h,...){dispose(window)})
      
    })
    
    
    fuicim<-function(h,...){
      if(isExtant(icimw)==FALSE){
        icimw<-gwindow("QTLIciMaping Format",visible = FALSE,width =80,height = 140)
        icimg<-ggroup(container = icimw,expand=FALSE)
      }
      iclyt<-glayout(container = icimg,spacing = 8)
      iclabel1<-glabel("1. Please input QTLIciMaping format data",container = iclyt)
      icimbt<-gbutton("QTLIciMaping_Format",container = iclyt)
      
      iclabel<-glabel("2. Please Select Mapping Population Type",container = iclyt)
      icradio<-gradio(c("BC1 ( F1 X P1 )","BC2 ( F1 X P2 )","DH","RIL","Chromosome Segment Substitution Line (CSSL)"),container =iclyt,selected = 1,horizontal = FALSE )
      icok<-gbutton(" DO ",container = iclyt)
      iccancel<-gbutton(" Cancel",container = iclyt )
      
      iclabelmd<-glabel("3. Please Select QTL-effect Model",container = iclyt)
      icradiomd<-gradio(c("Random Model","Fixed Model"),container = iclyt,selected = 1,horizontal = FALSE)
      iclyt[1,1:6]<-iclabel1
      iclyt[2,1:3]<- icimbt
      iclyt[4,1:6]<-iclabel
      iclyt[5:9,1:6]<-icradio
      iclyt[11,1:5]<-iclabelmd
      iclyt[12:13,1:3]<-icradiomd
      iclyt[14,1:3]<-icok
      iclyt[15,1:3]<-iccancel
      visible(icimw)<-TRUE
      addhandlerclicked(icimbt,handler = function(h,...){
        input1<-gfile(text = "select a directory...",type = "open",
                      filter = list("All files" = list(patterns = c("*")),
                                    "Excel files"=list(patterns=c("*.xlsx"))))
        if(is.na(input1)){
          gmessage("Please input QTLIciMaping format data!",title = "warning",icon="warning")
          return
        }else{
          envim$geo<-read.xlsx(input1,sheet = "Genotype",colNames = FALSE)
          envim$pos<-read.xlsx(input1,sheet = "LinkageMap",colNames = FALSE)
          envim$pho<-read.xlsx(input1,sheet = "Phenotype",colNames = FALSE)
          
          gdf2<-gdfedit(as.data.frame(envim$geo),expand=TRUE,container=nb1,label="Raw_Genotype")
          gdf3<-gdfedit(as.data.frame(envim$pho),expand=TRUE,container=nb1,label="Raw_phenotype") 
          gdf4<-gdfedit(as.data.frame(envim$pos),expand=TRUE,container=nb1,label="Raw_Linkage Map")
        }
      })
      
      addhandlerclicked(icok,handler = function(h,...){
        
        if((is.null(envim$geo)==TRUE)&&(is.null(envim$pho)==TRUE)&&(is.null(envim$pos)==TRUE)){
          gmessage("Please input QTLIciMaping format data!",title = "warning",icon="warning")
          return
        }else{
          envim$geo<-as.matrix(envim$geo)
          envim$pos<-as.matrix(envim$pos)
          envim$pho<-as.matrix(envim$pho)
          geo<-envim$geo
          pos<-envim$pos
          pho<-envim$pho
          
          if(svalue(icradio)=="BC1 ( F1 X P1 )"&&svalue(icradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1 
            envim$flagRIL<-0
            
            gen_0<-geo[,-1]
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("1","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)

            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            dispose(icimw)
            
          }else if(svalue(icradio)=="BC1 ( F1 X P1 )"&&svalue(icradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            
            
            gen_0<-geo[,-1]
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("1","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            dispose(icimw)
            
          }else if(svalue(icradio)=="BC2 ( F1 X P2 )"&&svalue(icradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1 
            envim$flagRIL<-0
            
            gen_0<-geo[,-1]
            gen_1<-gsub("-1","99",gen_0)
            gen_11<-gsub("0","-1",gen_1)
            # gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            
            
            dispose(icimw)
            
          } else if(svalue(icradio)=="BC2 ( F1 X P2 )"&&svalue(icradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            
            
            gen_0<-geo[,-1]
            gen_1<-gsub("-1","99",gen_0)
            gen_11<-gsub("0","-1",gen_1)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            
            dispose(icimw)
            
          }else if(svalue(icradio)=="DH"&&svalue(icradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-0
            
            gen_0<-geo[,-1]
            
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("0","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            
            dispose(icimw) 
            
          } else if(svalue(icradio)=="DH"&&svalue(icradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            
            
            gen_0<-geo[,-1]
            
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("0","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            
            dispose(icimw)
            
            
          }else if(svalue(icradio)=="RIL"&&svalue(icradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-1
            
            
            gen_0<-geo[,-1]
            
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("0","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            
            dispose(icimw)
            
          }else if(svalue(icradio)=="RIL"&&svalue(icradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-1
            
            gen_0<-geo[,-1]
            
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("0","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            
            dispose(icimw)
          }else if(svalue(icradio)=="Chromosome Segment Substitution Line (CSSL)"&&svalue(icradiomd)=="Random Model"){
            envim$flagCSSL<-1
            envim$flag<-1
            
            gen_0<-geo[,-1]
            
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("0","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            
            dispose(icimw) 
            
            
          }else{
            envim$flagCSSL<-1
            envim$flag<-0
            
            
            gen_0<-geo[,-1]
            
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("0","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            
            pos_0<-pos[,-1]
            gen_con<-cbind(pos_0,gen_11)
            gen<-matrix(as.numeric(gen_con),nrow(gen_con),ncol(gen_con))
            
            phett<-t(pho)
            phe_m<-as.matrix(phett[-1,])
            
            phe_sub<-gsub(-100,NA,phe_m)
            phe_00<-matrix(as.numeric(phe_sub),nrow(phe_m),ncol(phe_m))
            
            seq_indiv<-seq(1,nrow(phe_00))
            seq_indiv1<-c("genotype",seq_indiv)
            seq_indiv1<-matrix(seq_indiv1,1,)
            geo1<-cbind(geo[,1],gen_11)
            
            seq_indiv2<-c("phenotype",seq_indiv)
            seq_indiv2<-matrix(seq_indiv2,,1)
            
            phename<-matrix(phett[1,],1,)
            phe<-rbind(phename,phe_00)
            pheRaw<-cbind(seq_indiv2,phe)
            genRaw<-rbind(seq_indiv1,geo1)
            colname_mapRaw1<-c("marker","chr","pos")
            colname_mapRaw1<-matrix(colname_mapRaw1,1,)
            mapRaw1<-rbind(colname_mapRaw1,pos)
            envim$pheRaw<-pheRaw
            envim$genRaw<-genRaw
            envim$mapRaw1<-mapRaw1
            
            pheRaw_show<-pheRaw[-1,]
            genRaw_show<-genRaw[-1,]
            mapRaw_show<-mapRaw1[-1,]
            colnames(pheRaw_show)<-pheRaw[1,]
            colnames(genRaw_show)<-genRaw[1,]
            colnames(mapRaw_show)<-mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(genRaw_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(pheRaw_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(mapRaw_show),expand=TRUE,container=nb1,label="Linkage Map")
            
            dispose(icimw) 
            
            
          }
        }
      })
      
      addhandlerclicked(iccancel,handler=function(h,...){
        dispose(icimw)
      })
    }
    
    
    fuwinqtl<-function(h,...){
      if(isExtant(cartw)==FALSE){
        cartw<-gwindow("WinQTLCart  Format",visible = FALSE,width = 80,height = 120)
        cartg<-ggroup(container = cartw)
      }
      cartlyt<-glayout(container = cartg,spacing = 8)
      cartlabel1<-glabel("1. Please input WinQTLCart format file",container = cartlyt )
      cartbt<-gbutton("WinQTLCart_Format",container = cartlyt)
      
      cartlabel<-glabel("2. Please Select Mapping Population Type",container = cartlyt)
      cartradio<-gradio(c("BC1 ( F1 X P1 )","BC2 ( F1 X P2 )","DH","RIL","Chromosome Segment Substitution Line (CSSL)"),container =cartlyt,selected = 1,horizontal = FALSE )
      cartok<-gbutton(" DO ",container = cartlyt)
      cartcancel<-gbutton(" Cancel",container = cartlyt )
      cartlabelmd<-glabel("3. Please Select QTL-effect Model",container = cartlyt)
      cartradiomd<-gradio(c("Random Model","Fixed Model"),container =cartlyt,selected = 1,horizontal = FALSE )
      
      cartlyt[1,2:5]<-cartlabel1
      cartlyt[2,2:3]<-cartbt
      cartlyt[3,2:5]<-cartlabel
      cartlyt[4:7,2:6]<-cartradio
      cartlyt[8,2:4]<-cartlabelmd
      cartlyt[9:10,2:4]<-cartradiomd
      cartlyt[11,2:3]<-cartok
      cartlyt[12,2:3]<-cartcancel
      visible(cartw)<-TRUE
      addhandlerclicked(cartbt,handler = function(h,...){
        input2<-gfile(text = "select a directory...",type = "open",
                      filter = list("All files" = list(patterns = c("*")),
                                    "Excel files"=list(patterns=c("*.mcd"))))
        if(is.na(input2)){
          gmessage("Please input WinQTLCart data!",title = "warning",icon="warning")
          return
        }else{
          envim$y_jun00<-scan(input2,what = "",sep = "\n")
          envim$y_jun3<-scan(input2,what = "",sep = "")
          tb<-gnewtable(envim$y_jun00,multiple=TRUE,container=nb1,expand=TRUE,label="WinQTLCart")
        }
      })
      
      addhandlerclicked(cartok,handler = function(h,...){
        
        if((is.null(envim$y_jun00)==TRUE)&&(is.null(envim$y_jun3)==TRUE)){
          gmessage("Please input WinQTLCart data!",title = "warning",icon="warning")
          return
        }else{
          
          y_jun00<-envim$y_jun00
          y_jun3<-envim$y_jun3
          if(svalue(cartradio)=="BC1 ( F1 X P1 )"&&svalue(cartradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-0
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp1<-gsub("1","-1",marker_snp)
            marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp2,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
            
          }else if(svalue(cartradio)=="BC1 ( F1 X P1 )"&&svalue(cartradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp1<-gsub("1","-1",marker_snp)
            marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp2,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
            
          }else if(svalue(cartradio)=="BC2 ( F1 X P2 )"&&svalue(cartradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-0
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp1<-gsub("0","-1",marker_snp)
            #marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp1,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            envim$genRaw<-as.data.frame(envim$genRaw)
            envim$pheRaw<-as.data.frame(envim$pheRaw)
            envim$mapRaw1<-as.data.frame(envim$mapRaw1)
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
            
          }else if(svalue(cartradio)=="BC2 ( F1 X P2 )"&&svalue(cartradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            
            
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp1<-gsub("0","-1",marker_snp)
            #marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp1,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
          }else if(svalue(cartradio)=="DH"&&svalue(cartradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-0
            
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp1<-gsub("0","-1",marker_snp)
            marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp2,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
            
            
          }else if(svalue(cartradio)=="DH"&&svalue(cartradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp1<-gsub("0","-1",marker_snp)
            marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp2,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
            
            
          }else if(svalue(cartradio)=="RIL"&&svalue(cartradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-1
            
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp10<-gsub("1","99",marker_snp)
            marker_snp1<-gsub("0","-1",marker_snp10)
            marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp2,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg1<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
            
            
          }else if(svalue(cartradio)=="RIL"&&svalue(cartradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-1
            
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp10<-gsub("1","99",marker_snp)
            marker_snp1<-gsub("0","-1",marker_snp10)
            marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp2,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
          }else if(svalue(cartradio)=="Chromosome Segment Substitution Line (CSSL)"&&svalue(cartradiomd)=="Random Model"){
            envim$flagCSSL<-1
            envim$flag<-1
            
            
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp1<-gsub("0","-1",marker_snp)
            marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp2,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            dispose(cartw)
            
          }else{
            envim$flagCSSL<-1
            envim$flag<-0
            
            
            
            start_dex<-grep("-start",y_jun3,fixed = TRUE)
            stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
            
            chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
            
            chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
            chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
            chrdata_dexlen<-length(chrdata_dex)
            
            chr_num<-numeric()
            chrname_num<-numeric()
            chr_numfirst<-numeric()
            markername0<-numeric()
            chr_pos<-numeric()
            chrRaw_name<-numeric()
            chr_Rawnumfirst<-numeric()
            for(i in 1:(chrdata_dexlen-1)){
              chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
              marker_name<-numeric()
              marker_pos<-numeric()
              
              for(j in 0:(trunc(length(chr_each)/2)-1) ){
                marker_name0<-chr_each[2*j+1]
                marker_name<-cbind(marker_name,marker_name0)
                marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
                marker_pos<-cbind(marker_pos,marker_pos0)
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                  marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                  marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                  
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                  marker_posm<-100*marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                  marker_posm<-marker_pos
                  
                  markerlen<-length(marker_posm)
                  marker_pos1<-marker_posm[1:(markerlen-1)]
                  marker_pos2<-c(0,marker_pos1)
                  marker_possum<-cumsum(marker_pos2)
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                  marker_possum<-100*marker_pos
                }
                if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                  marker_possum<-marker_pos
                }
              }
              markername0<-cbind(markername0,marker_name)
              markername<-matrix(markername0,,1)
              marker_possum0<-matrix(marker_possum,,1)
              chr_pos<-rbind(chr_pos,marker_possum0)
              chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
              chr_data<-na.omit(chr_a)
              chr_datalen<-length(chr_data)
              chr_num<-rbind(chr_num,chr_datalen)
              
              chrRawname<-chrdata[(chrdata_dex[i]+1)]
              chrname<-str_extract_all(chrRawname,"[0-9]+")
              chrRawname<-matrix(chrRawname,,1)
              chrRaw_name<-rbind(chrRaw_name,chrRawname)
              
              chrname_num<-rbind(chrname_num,chrname)
              chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
              chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
              chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
              chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
            }
            
            chr_leng<-length(chr_pos)
            
            chr_numtwo<-cbind(chr_numfirst,chr_pos)
            
            marker_dex<-grep("markers",y_jun3,fixed = TRUE)
            marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
            marker_snp1<-gsub("0","-1",marker_snp)
            marker_snp2<-gsub("2","1",marker_snp1)
            marker_snpnum<-gsub("*","99",marker_snp2,fixed = TRUE)
            
            snpa<-suppressWarnings(as.numeric(marker_snpnum))
            snpdata<-na.omit(snpa)
            indi_num<-length(snpdata)/chr_leng
            snp_data<-numeric()
            
            if_indi<-y_jun3[marker_dex[1]-1]
            if(if_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-t(snp_data)
            }else{
              for(i in 0:(chr_leng-1)){
                snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
                snp_data<-rbind(snp_data,snp_eve)
              }
              snp_data1<-snp_data
            }
            
            trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
            
            for(i in 1:length(trait_total)){
              if(trait_total[i]=="."){
                trait_total[i]<-"0"
              }
            }
            trait_dex<-grep("traits",trait_total)
            traita<-suppressWarnings(as.numeric(trait_total))
            traitdata<-na.omit(traita)
            trait_num<-length(traitdata)/indi_num
            trait_data<-numeric()
            
            iftrait_indi<-trait_total[trait_dex[1]-1]
            if(iftrait_indi=="individuals"){
              
              for(i in 0:(indi_num-1)){
                trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
                for(j in 1:length( trait_bbb)){
                  if(trait_bbb[j]==0){
                    trait_bbb[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_bbb,1,)
                trait_data<-rbind(trait_data,trait_eve)
              }
              
            }else{
              
              for(i in 0:(trait_num-1)){
                trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
                for(j in 1:length(trait_aaa)){
                  if(trait_aaa[j]==0){
                    trait_aaa[j]<-NA
                  }
                }
                trait_eve<-matrix(trait_aaa,,1)
                trait_data<-cbind(trait_data,trait_eve)
              }
              
            }
            
            if(length(start_dex)==3){yygg1<-0}
            if(length(start_dex)==4){
              if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
                cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
                cov_dex<-grep("otraits",cov_total)
                cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
                bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
                
                otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
                
                
                if(y_jun3[start_dex[4]+1]=="otraits"){
                  
                  covnumonly<-numeric()
                  for( i in 0:(otrait_indi-1)){
                    cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                    covnumonly<-rbind(covnumonly,cov_each)
                  }
                  covnum<-covnumonly[,-1] 
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                  
                }
                if(y_jun3[start_dex[4]+1]=="individuals"){
                  covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                  covnum<-numeric()
                  otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                  for(m in 1:otrait_indi){
                    
                    coveach<-numeric()
                    for(n in 0:(indi_num-1)){
                      cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                      
                      coveach<-cbind(coveach,cov_each)
                    }
                    covnum<-rbind(covnum,coveach)
                  }
                  covnumonly<-cbind(otrait_name,covnum)
                  yygg<-numeric()
                  for(i in 1:nrow(covnum)){
                    
                    otrait_ind<-unique(covnum[i,])
                    cov_col<-length(otrait_ind)-1
                    
                    col_each<-numeric()
                    for(j in 1:length(covnum[i,])){
                      
                      if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                        cov_0<-matrix(-1,1,cov_col)
                        
                      }else{
                        cov_0<-matrix(0,1,cov_col)
                        covnum_loc<-which(otrait_ind[]==covnum[i,j])
                        cov_0[1,covnum_loc]<-1
                      }
                      col_each<-rbind(col_each,cov_0)
                      
                    }
                    yygg<-cbind(yygg,col_each)
                    
                  } 
                }
                
              }
            }
            
            seq_indi<-seq(1,nrow(trait_data))
            seq_indi1<-c("genotype",seq_indi)
            seq_indi1<-matrix(seq_indi1,1,)
            snp1<-cbind(markername,snp_data1)
            seq_indi2<-c("phenotype",seq_indi)
            seq_indi2<-matrix(seq_indi2,,1)
            num_trait<-ncol(trait_data)
            seq_trait<-seq(1,num_trait)
            seq_trait<-matrix(seq_trait,1,)
            trait_data00<-rbind(seq_trait,trait_data)
            
            colnames_mapname<-c("marker","chr","pos")
            colnames_mapname<-matrix(colnames_mapname,1,)
            mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
            envim$mapRaw1<-rbind(colnames_mapname,mapRaw1)
            
            envim$genRaw<-rbind(seq_indi1,snp1)
            envim$pheRaw<-cbind(seq_indi2,trait_data00)
            
            gen_show<-envim$genRaw[-1,]
            phe_show<-envim$pheRaw[-1,]
            map_show<-envim$mapRaw1[-1,]
            colnames(gen_show)<-envim$genRaw[1,]
            colnames(phe_show)<-envim$pheRaw[1,]
            colnames(map_show)<-envim$mapRaw1[1,]
            
            gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
            gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
            gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
            if(exists("yygg")==TRUE){
              envim$yygg<-yygg
              envim$cov_en<-t(covnumonly)
              enabled(covbt)<-FALSE
              gdf6<-gdfedit(as.data.frame(envim$cov_en),expand=TRUE,container=nb1,label="Covariate")
            }
            
            dispose(cartw)
            
          }
          
        }
      })
      addhandlerclicked(cartcancel,handler = function(h,...){
        dispose(cartw)
      })
    }
    
    
    fugcim<-function(h,...){
      if(isExtant(gcimw)==FALSE){
        gcimw<-gwindow("GCIM Format",visible = FALSE,width = 85,height = 160)
        gcimg<-ggroup(container = gcimw)
      }
      
      gclyt<-glayout(container = gcimg,spacing = 8)
      gclabel1<-glabel("1. Please input Genotypic, Phenotypic and Linkage Map files",container = gclyt)
      gcbtgen<-gbutton("Genotype",container = gclyt)
      gcbtphe<-gbutton("Phenotype",container = gclyt)
      gcbtmap<-gbutton("Linkage Map",container = gclyt)
      gclabel<-glabel("2. Please Select Mapping Population Type",container = gclyt)
      gcradio<-gradio(c("BC1 ( F1 X P1 )","BC2 ( F1 X P2 )","DH","RIL","Chromosome Segment Substitution Line (CSSL)"),container = gclyt,selected = 1,horizontal = FALSE)
      frame<-gframe("Marker Genotype Table",container = gclyt)
      framlyt<-glayout(container =frame )
      framlabel_AA<-glabel("AA",container = framlyt)
      framedit_AA<-gedit("1",container = framlyt,width = 10,coerce.with=as.character)
      framlabel_Aa<-glabel("Aa",container = framlyt)
      framedit_Aa<-gedit("-1",container = framlyt,width = 10,coerce.with=as.character)
      framlabel_aa<-glabel("aa",container = framlyt)
      framedit_aa<-gedit("-1",container = framlyt,width = 10,coerce.with=as.character)
      framlabel_<-glabel("Missing",container = framlyt)
      framedit_<-gedit("99",container = framlyt,width = 10,coerce.with=as.character)
      framlyt[1,1]<-framlabel_AA
      framlyt[1,4:6]<-framedit_AA
      framlyt[2,1]<-framlabel_Aa
      framlyt[2,4:6]<-framedit_Aa
      framlyt[3,1]<-framlabel_aa
      framlyt[3,4:6]<-framedit_aa
      framlyt[4,1]<-framlabel_
      framlyt[4,4:6]<-framedit_
      value_AA<-svalue(framedit_AA) 
      value_Aa<-svalue(framedit_Aa) 
      value_aa<-svalue(framedit_aa) 
      value_<-svalue(framedit_) 
      
      
      gclabelmd<-glabel("3. Please Select QTL-effect Model",container = gclyt)
      gcradiomd<-gradio(c("Random Model","Fixed Model"),container = gclyt)
      gcdobt<-gbutton("DO",container = gclyt)
      gccancel<-gbutton(" Cancel ",container = gclyt)
      gclyt[1,1:6]<-gclabel1
      gclyt[2,1:3]<-gcbtgen
      gclyt[3,1:3]<-gcbtphe
      gclyt[4,1:3]<-gcbtmap
      gclyt[5,1:6]<-gclabel
      gclyt[6:9,1:6]<-gcradio
      gclyt[10:13,1:6]<-frame
      
      gclyt[14,1:5]<-gclabelmd
      gclyt[15:16,1:4]<-gcradiomd
      gclyt[17,1:3]<-gcdobt
      gclyt[18,1:3]<-gccancel
      visible(gcimw)<-TRUE
      
      if(svalue(gcradio)=="BC1 ( F1 X P1 )"){
        enabled(framedit_aa) <- FALSE
      }
      addhandlerchanged(gcradio,handler=function(h,...){
        if(svalue(gcradio)=="BC1 ( F1 X P1 )"){
          enabled(framedit_Aa) <- TRUE
          enabled(framedit_aa) <- FALSE
          enabled(framedit_AA) <- TRUE
          svalue(framedit_Aa)<--1
          envim$value_Aa<-"-1"
        }
        else if(svalue(gcradio)=="BC2 ( F1 X P2 )"){
          enabled(framedit_Aa) <- TRUE
          enabled(framedit_aa) <- TRUE
          enabled(framedit_AA) <- FALSE
          svalue(framedit_Aa)<-1
          envim$value_Aa<-"1"
        }
        else if(svalue(gcradio)=="DH"){
          enabled(framedit_Aa) <- FALSE
          enabled(framedit_aa) <- TRUE
          enabled(framedit_AA) <- TRUE
        }
        else  if(svalue(gcradio)=="Chromosome Segment Substitution Line (CSSL)"){
          enabled(framedit_Aa) <- FALSE
          enabled(framedit_aa) <- TRUE
          enabled(framedit_AA) <- TRUE
        }else{
          enabled(framedit_Aa) <- FALSE
          enabled(framedit_aa) <- TRUE
          enabled(framedit_AA) <- TRUE
        }
        
      })
      addhandlerclicked(gcbtgen,handler = function(h,...){
        input3<-gfile(text="Select a directory...",type="open",
                      filter=list("All files"=list(patterns=c("*")),
                                  "CSV files"=list(patterns=c("*.csv"))))
        if(is.na(input3)){
          gmessage("Please input phenotype data!",title = "warning",icon="warning")
          return
        }else{
          envim$genRaw1<-read.csv(input3,header = F)
          
          genRaw1_1<-as.matrix(envim$genRaw1)
          genR_show<-genRaw1_1[-1,]
          colnames(genR_show)<-genRaw1_1[1,]
          
          gdf3<-gdfedit(as.data.frame(genR_show),expand=TRUE,container=nb1,label="Raw_Genotype")
          
        }
      })
      addhandlerclicked(gcbtphe,handler = function(h,...){
        input4<-gfile(text="Select a directory...",type="open",
                      filter=list("All files"=list(patterns=c("*")),
                                  "CSV files"=list(patterns=c("*.csv"))))
        if(is.na(input4)){
          gmessage("Please input phenotype data!",title = "warning",icon="warning")
          return
        }else{
          envim$pheRaw<-read.csv(input4,header = F)
          
          
          pheRaw1_1<-as.matrix(envim$pheRaw)
          pheR_show<-pheRaw1_1[-1,]
          colnames(pheR_show)<-pheRaw1_1[1,]
          
          gdf4<-gdfedit(as.data.frame(pheR_show),expand=TRUE,container=nb1,label="Raw_Phenotype")
          
        }
      })
      addhandlerclicked(gcbtmap,handler = function(h,...){
        input5<-gfile(text="Select a directory...",type="open",
                      filter=list("All files"=list(patterns=c("*")),
                                  "CSV files"=list(patterns=c("*.csv"))))
        if(is.na(input5)){
          gmessage("Please input linkage map data!",title = "warning",icon="warning")
          return
        }else{
          envim$mapRaw11<-read.csv(input5,header = F)
          
          mapRaw1_1<-as.matrix(envim$mapRaw11)
          mapR_show<-mapRaw1_1[-1,]
          colnames(mapR_show)<-mapRaw1_1[1,]
          gdf7<-gdfedit(as.data.frame(mapR_show),expand=TRUE,container=nb1,label="Raw_Linkage Map")
        }
      })
      addhandlerclicked(gcdobt,handler = function(h,...){
        enabled(covbt)<-TRUE
        
        genRaw1<-envim$genRaw1
        pheRaw<-envim$pheRaw
        mapRaw1<-envim$mapRaw11
        
        if(is.null(genRaw1)==TRUE){
          gmessage("Please input genotype data!",title = "warning",icon="warning")
          return
        }
        if(is.null(pheRaw)==TRUE){
          gmessage("Please input phenotype data!",title = "warning",icon="warning")
          return
        }
        if(is.null(mapRaw1)==TRUE){
          gmessage("Please input linkage map data!",title = "warning",icon="warning")
          return
        }
        if((is.null(genRaw1)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)){
          genRaw1<-as.matrix(genRaw1)
          pheRaw<-as.matrix(pheRaw)
          mapRaw1<-as.matrix(mapRaw1)
          
          if((svalue(gcradio)=="BC1 ( F1 X P1 )")&&(svalue(gcradiomd)=="Random Model")){
            
            enabled(framedit_aa)<-FALSE
            
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-0
            envim$value_AA<-svalue(framedit_AA) 
            envim$value_Aa<-svalue(framedit_Aa)
            envim$value_aa<-svalue(framedit_aa)
            envim$value_<-svalue(framedit_) 
            
            if(envim$value_AA=="1"&&envim$value_Aa=="-1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }else{
              
              
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              genRaw_Aa<-gsub(envim$value_Aa,"-1",genRaw_)
              genRaw_AA<-gsub(envim$value_AA,"1",genRaw_Aa)
              
              envim$genRaw<-genRaw_AA
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              
              
              dispose(gcimw) 
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
            
            
            
          }else if(svalue(gcradio)=="BC1 ( F1 X P1 )"&&svalue(gcradiomd)=="Fixed Model"){
            enabled(framedit_aa)<-FALSE
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            envim$value_AA<-svalue(framedit_AA) 
            envim$value_Aa<-svalue(framedit_Aa)
            envim$value_aa<-svalue(framedit_aa)
            envim$value_<-svalue(framedit_) 
            
            
            if(envim$value_AA=="1"&&envim$value_Aa=="-1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }else{
              
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              genRaw_Aa<-gsub(envim$value_Aa,"-1",genRaw_)
              genRaw_AA<-gsub(envim$value_AA,"1",genRaw_Aa)
              
              envim$genRaw<-genRaw_AA
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
            
            
          }else if((svalue(gcradio)=="BC2 ( F1 X P2 )")&&(svalue(gcradiomd)=="Random Model")){
            
            enabled(framedit_AA)<-FALSE
            
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-0
            envim$value_AA<-svalue(framedit_AA) 
            envim$value_Aa<-svalue(framedit_Aa)
            envim$value_aa<-svalue(framedit_aa)
            envim$value_<-svalue(framedit_) 
            
            if(envim$value_Aa=="1"&&envim$value_aa=="-1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
              
            }else{
              
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              genRaw_aa<-gsub(envim$value_aa,"-1",genRaw_)
              genRaw_Aa<-gsub(envim$value_Aa,"1",genRaw_aa)
              
              envim$genRaw<-genRaw_Aa
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
            
            
            
          }else if((svalue(gcradio)=="BC2 ( F1 X P2 )")&&(svalue(gcradiomd)=="Fixed Model")){
            
            enabled(framedit_AA)<-FALSE
            
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            envim$value_AA<-svalue(framedit_AA) 
            envim$value_Aa<-svalue(framedit_Aa)
            envim$value_aa<-svalue(framedit_aa)
            envim$value_<-svalue(framedit_) 
            
            if(envim$value_Aa=="1"&&envim$value_aa=="-1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
              
            }else{
              
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              genRaw_aa<-gsub(envim$value_aa,"-1",genRaw_)
              genRaw_Aa<-gsub(envim$value_Aa,"1",genRaw_aa)
              
              envim$genRaw<-genRaw_Aa
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
            
            
          }else if(svalue(gcradio)=="DH"&&svalue(gcradiomd)=="Random Model"){
            enabled(framedit_Aa)<-FALSE
            
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-0
            
            envim$value_AA<-svalue(framedit_AA) 
            
            envim$value_aa<-svalue(framedit_aa) 
            envim$value_<-svalue(framedit_) 
            
            if(envim$value_AA=="1"&&envim$value_aa=="-1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }else{
              
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              genRaw_aa<-gsub(envim$value_aa,"-1",genRaw_)
              genRaw_AA<-gsub(envim$value_AA,"1",genRaw_aa)
              
              envim$genRaw<-genRaw_AA
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              dispose(icimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
          }else if(svalue(gcradio)=="DH"&&svalue(gcradiomd)=="Fixed Model"){
            enabled(framedit_Aa)<-FALSE
            
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-0
            envim$value_AA<-svalue(framedit_AA) 
            
            envim$value_aa<-svalue(framedit_aa) 
            envim$value_<-svalue(framedit_)
            
            if(envim$value_AA=="1"&&envim$value_aa=="-1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }else{
              
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              genRaw_aa<-gsub(envim$value_aa,"-1",genRaw_)
              genRaw_AA<-gsub(envim$value_AA,"1",genRaw_aa)
              
              envim$genRaw<-genRaw_AA
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              dispose(icimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
            
            
          }else if(svalue(gcradio)=="RIL"&&svalue(gcradiomd)=="Random Model"){
            envim$flagCSSL<-0
            envim$flag<-1
            envim$flagRIL<-1
            envim$value_AA<-svalue(framedit_AA) 
            envim$value_Aa<-svalue(framedit_Aa) 
            envim$value_aa<-svalue(framedit_aa) 
            envim$value_<-svalue(framedit_) 
            
            if(envim$value_aa=="-1"&&envim$value_AA=="1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }else{
              
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              
              genRaw_aa<-gsub(envim$value_aa,"-1",genRaw_)
              genRaw_AA<-gsub(envim$value_AA,"1",genRaw_aa)
              
              envim$genRaw<-genRaw_AA
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              dispose(icimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
          }else if(svalue(gcradio)=="RIL"&&svalue(gcradiomd)=="Fixed Model"){
            envim$flagCSSL<-0
            envim$flag<-0
            envim$flagRIL<-1
            
            envim$value_AA<-svalue(framedit_AA) 
            envim$value_Aa<-svalue(framedit_Aa) 
            envim$value_aa<-svalue(framedit_aa) 
            envim$value_<-svalue(framedit_) 
            
            if(envim$value_aa=="-1"&&envim$value_AA=="1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }else{
              
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              
              genRaw_aa<-gsub(envim$value_aa,"-1",genRaw_)
              genRaw_AA<-gsub(envim$value_AA,"1",genRaw_aa)
              
              envim$genRaw<-genRaw_AA
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              dispose(icimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
             }
            }else if(svalue(gcradio)=="Chromosome Segment Substitution Line (CSSL)"&&svalue(gcradiomd)=="Random Model"){
            enabled(framedit_Aa)<-FALSE
            envim$flagCSSL<-1
            envim$flag<-1
            
            envim$value_AA<-svalue(framedit_AA) 
            
            envim$value_aa<-svalue(framedit_aa) 
            envim$value_<-svalue(framedit_)
            
            if(envim$value_aa=="-1"&&envim$value_AA=="1"&&envim$value_=="99"){
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              enabled(framedit_inter)<-FALSE
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }else{
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              genRaw_aa<-gsub(envim$value_aa,"-1",genRaw_)
              genRaw_AA<-gsub(envim$value_AA,"1",genRaw_aa)
              
              envim$genRaw<-genRaw_AA
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              enabled(framedit_inter)<-FALSE
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
          }else{
            enabled(framedit_Aa)<-FALSE
            envim$flagCSSL<-1
            envim$flag<-0
            envim$value_AA<-svalue(framedit_AA) 
            
            envim$value_aa<-svalue(framedit_aa) 
            envim$value_<-svalue(framedit_)
            
            if(envim$value_aa=="-1"&&envim$value_AA=="1"&&envim$value_=="99"){
              enabled(framedit_inter)<-FALSE
              envim$genRaw<-genRaw1
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }else{
              genRaw_<-gsub(envim$value_,"99",genRaw1)
              genRaw_aa<-gsub(envim$value_aa,"-1",genRaw_)
              genRaw_AA<-gsub(envim$value_AA,"1",genRaw_aa)
              
              envim$genRaw<-genRaw_AA
              envim$pheRaw<-pheRaw
              envim$mapRaw1<-mapRaw1
              
              enabled(framedit_inter)<-FALSE
              
              gen_show<-envim$genRaw[-1,]
              phe_show<-envim$pheRaw[-1,]
              map_show<-envim$mapRaw1[-1,]
              colnames(gen_show)<-envim$genRaw[1,]
              colnames(phe_show)<-envim$pheRaw[1,]
              colnames(map_show)<-envim$mapRaw1[1,]
              
              gdf3<-gdfedit(as.data.frame(gen_show),expand=TRUE,container=nb1,label="Genotype")
              gdf4<-gdfedit(as.data.frame(phe_show),expand=TRUE,container=nb1,label="Phenotype")
              gdf5<-gdfedit(as.data.frame(map_show),expand=TRUE,container=nb1,label="Linkage Map")
              dispose(gcimw)
              infowin<-gwindow("Info",width =200,height = 150 ,visible=FALSE)
              infolyt<-glayout(container = infowin)
              infobt<-gbutton("Continue",container =infolyt )
              infolyt[5:6,4:6]<-infobt
              addhandlerclicked(infobt,handler = function(h,...)dispose(infowin))
              visible(infowin)<-TRUE
            }
            
          }
        }
      })
      addhandlerclicked(gccancel,handler = function(h,...){
        dispose(gcimw)
      })
    }
    
    addhandlerclicked(covbt,function(h,...){
      if(isExtant(covw)==FALSE){
        covw<-gwindow("Include Covariate ?",visible = FALSE,width = 250,height = 130)
        covg<-ggroup(container = covw)
      }
      lytcov<-glayout(container=covg,spacing=13)
      okcov<-gbutton("     OK    ",container=lytcov)
      cancelcov<-gbutton(" Cancel ",container=lytcov)
      radiocov<-gradio(c("Included","No"),selected=1,horizontal=FALSE,container=lytcov)
      lytcov[2:4,2:5]<-radiocov
      lytcov[5,2]<-okcov
      lytcov[5,7]<-cancelcov
      visible(covw)<-TRUE
      addhandlerclicked(okcov,handler = function(h,...){
        if(svalue(radiocov)=="Included"){
          
          envim$flagcov<-1
          
          input4<-gfile(text = "select a directory...",type = "open",
                        filter = list("All files" = list(patterns = c("*")),
                                      "Excel files"=list(patterns=c("*.csv"))))
          if(is.na(input4)){
            gmessage("Please input covariates data!",title = "warning",icon="warning")
          }else{
            envim$cov_en<-as.matrix(read.csv(input4,header = FALSE))
            
            cov_enShow<-envim$cov_en[-1,]
            colnames(cov_enShow)<-envim$cov_en[1,]
            
            gdf4<-gdfedit(as.data.frame(cov_enShow),expand=TRUE,container = nb1,label="Covariate")
            
            cov_en1<-envim$cov_en[-1,2:ncol(envim$cov_en)]
            covnum<-t(cov_en1)
            yygg1<-numeric()
            for(i in 1:nrow(covnum)){
              
              otrait_ind<-unique(covnum[i,])
              cov_col<-length(otrait_ind)-1
              
              col_each<-numeric()
              for(j in 1:length(covnum[i,])){
                
                if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                  cov_0<-matrix(-1,1,cov_col)
                  
                }else{
                  cov_0<-matrix(0,1,cov_col)
                  covnum_loc<-which(otrait_ind[]==covnum[i,j])
                  cov_0[1,covnum_loc]<-1
                }
                col_each<-rbind(col_each,cov_0)
                
              }
              yygg1<-cbind(yygg1,col_each)
              
            } 
            
            envim$yygg1<-yygg1
            dispose(covw)
          }
          
        }else{
          envim$flagcov<-0
          
          dispose(covw)
        }
      })
      addhandlerclicked(cancelcov,function(h,...){
        dispose(covw)
        
      })
    })
    
    addhandlerclicked(runbt,handler = function(h,...){
      cl<-svalue(framedit_inter)
      sLOD<-svalue(framedit_lod)
      envim$ii<-svalue(framedit_trait)
      pheRaw<-envim$pheRaw
      mapRaw1<-envim$mapRaw1
      genRaw<-envim$genRaw
      
      flagRIL<-envim$flagRIL
      flag<-envim$flag
      
      yygg<-envim$yygg
      if(is.null(genRaw)==TRUE){
        gmessage("Please input correct genotype data!","Warning",icon="warning")
        return
      }
      if(is.null(pheRaw)==TRUE){
        gmessage("Please input correct phenotype data!","Warning",icon="warning")
        return
      }
      if(is.null(mapRaw1)==TRUE){
        gmessage("Please input correct linkage map data!","Warning",icon="warning")
        return
      }
      if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(envim$ii>(ncol(pheRaw)-1))){
        gmessage("It is more than Trait numbers ","Warning",icon="warning")
        return
      }
      if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(envim$ii<=(ncol(pheRaw)-1))&&(cl<0)){
        gmessage("Please input Walk Speed more than 0!","Warning",icon="warning")
        return
      }
      if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(envim$ii<=(ncol(pheRaw)-1))&&(cl>0)&&(sLOD<0)){
        gmessage("Please input critical LOD score more than 0!","Warning",icon="warning")
        return
      }
      if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(envim$ii<=(ncol(pheRaw)-1))&&(cl>0)&&(sLOD>0)){
        
        progress_bar$setText ( "Please be patient ..." )
        progress_bar$setFraction(2/100)
        

        pheRaw<-as.matrix(pheRaw)
        mapRaw1<-as.matrix(mapRaw1)
        genRaw<-as.matrix(genRaw)
        
        if(is.null(envim$yygg1)==FALSE){
          envim$cov_en<-as.matrix(envim$cov_en)
          envim$yygg1<-as.matrix(envim$yygg1)
          covname<-envim$cov_en[2:nrow(envim$cov_en),1]
          yygg1<-cbind(covname,envim$yygg1)
          
          mapRaw10<-mapRaw1[-1,]
          
          chr_name<-unique(mapRaw10[,2])
          chr_secon<-mapRaw10[,2]
          
          mm<-numeric()
          map_chr<-numeric()
          for(i in 1:length(chr_name)){
            chr_i<-which(chr_secon[]==chr_name[i])
            len<-matrix(length(chr_i),,1)
            mm<-rbind(mm,len)
            chr_name[i]<-i
          }
          
          map_chr<-numeric()
          for(i in 1:length(chr_name)){
            chr_name<-as.numeric(chr_name)
            
            chr_pos1<-matrix(chr_name[i],mm[i],1)
            map_chr<-rbind(map_chr,chr_pos1)
            
          }
          
          map_marker<-matrix(mapRaw10[,1],,1)
          map_pos<-matrix(mapRaw10[,3],,1)
          mapRaw<-cbind(map_marker,map_chr,map_pos)
          
          nameMap<-matrix(mapRaw[,1],,1)
          nameGenrow<-matrix(genRaw[1,],1,)
          nameGencol<-matrix(genRaw[,1],,1)
          namePhe<-as.matrix(pheRaw[,1],,1)
          nameCov<-matrix(yygg1[,1],,1)
          if(nameGenrow[2]==" 1"){
            phee<-pheRaw[-1,-1]
            phe<-matrix(as.numeric(phee),nrow(phee),ncol(phee))
            mapname<-mapRaw
            genn<-genRaw[-1,-1]
            genn<-matrix(as.numeric(genn),nrow(genn),ncol(genn))
            mapnametwo<-matrix(as.numeric(mapname[,2:3]),nrow(mapname),2)
            mx<-cbind(mapnametwo,genn)
            chrRaw_name<-unique(mapRaw10[,2])
            chr_name<-chr_name
            yygg<-yygg1
            envim$mapname<-mapname
            envim$chrRaw_name<-chrRaw_name
            envim$chr_name<-chr_name
          }else{
            sameName_MG<-intersect(nameMap,nameGencol)
            sameName_PG<-intersect(namePhe,nameGenrow)
            locPhe<-match(sameName_PG,namePhe)
            locMap<-match(sameName_MG,nameMap)
            locGen_PG<-match(sameName_PG,nameGenrow)
            locGen_MG<-match(sameName_MG,nameGencol)
            
            locCov<-match(sameName_PG,nameCov)
            newyygg<-as.matrix(yygg1[locCov,])
            yyggChar<-newyygg[,-1]
            yygg<-matrix(as.numeric(yyggChar),nrow(yyggChar),ncol(yyggChar))
            
            newPhe<-as.matrix(pheRaw[locPhe,])
            newMap<-as.matrix(mapRaw[locMap,])
            
            newGenrow<-as.matrix(genRaw[,locGen_PG])
            newGen<-as.matrix(newGenrow[locGen_MG,])
            
            gen_two<-newMap[,2:3]
            genChar<-cbind(gen_two,newGen)
            if(ncol(newPhe)==2)
            {
              pheChar<-as.matrix(newPhe[,2:ncol(newPhe)])  
            }else if(ncol(newPhe)>2){
              pheChar<-newPhe[,2:ncol(newPhe)]
            }
            
            chrRaw_name<-unique(mapRaw10[,2])
            chr_name<-chr_name
            mx<-matrix(as.numeric(genChar),nrow(genChar),ncol(genChar))
            phe<-matrix(as.numeric(pheChar),nrow(pheChar),ncol(pheChar))
            mapname<-newMap
            
            envim$mapname<-mapname
            
            envim$chr_name<-chr_name
          }
        } else{
          mapRaw10<-mapRaw1[-1,]
          
          chr_name<-unique(mapRaw10[,2])
          chr_secon<-mapRaw10[,2]
          
          mm<-numeric()
          map_chr<-numeric()
          for(i in 1:length(chr_name)){
            chr_i<-which(chr_secon[]==chr_name[i])
            len<-matrix(length(chr_i),,1)
            mm<-rbind(mm,len)
            chr_name[i]<-i
          }
          
          map_chr<-numeric()
          for(i in 1:length(chr_name)){
            chr_name<-as.numeric(chr_name)
            
            chr_pos1<-matrix(chr_name[i],mm[i],1)
            map_chr<-rbind(map_chr,chr_pos1)
            
          }
          
          map_marker<-matrix(mapRaw10[,1],,1)
          map_pos<-matrix(mapRaw10[,3],,1)
          mapRaw<-cbind(map_marker,map_chr,map_pos)
          
          nameMap<-matrix(mapRaw[,1],,1)
          nameGenrow<-matrix(genRaw[1,],1,)
          nameGencol<-matrix(genRaw[,1],,1)
          namePhe<-as.matrix(pheRaw[,1],,1)
          
          if(nameGenrow[2]==" 1"){
            phee<-pheRaw[-1,-1]
            phe<-matrix(as.numeric(phee),nrow(phee),ncol(phee))
            mapname<-mapRaw
            genn<-genRaw[-1,-1]
            genn<-matrix(as.numeric(genn),nrow(genn),ncol(genn))
            mapnametwo<-matrix(as.numeric(mapname[,2:3]),nrow(mapname),2)
            mx<-cbind(mapnametwo,genn)
            chrRaw_name<-unique(mapRaw10[,2])
            chr_name<-chr_name
            
            envim$mapname<-mapname
            envim$chrRaw_name<-chrRaw_name
            envim$chr_name<-chr_name
          }else{
            sameName_MG<-intersect(nameMap,nameGencol)
            sameName_PG<-intersect(namePhe,nameGenrow)
            locPhe<-match(sameName_PG,namePhe)
            locMap<-match(sameName_MG,nameMap)
            locGen_PG<-match(sameName_PG,nameGenrow)
            locGen_MG<-match(sameName_MG,nameGencol)
            
            newPhe<-as.matrix(pheRaw[locPhe,])
            newMap<-as.matrix(mapRaw[locMap,])
            
            newGenrow<-as.matrix(genRaw[,locGen_PG])
            newGen<-as.matrix(newGenrow[locGen_MG,])
            
            gen_two<-newMap[,2:3]
            genChar<-cbind(gen_two,newGen)
            
            if(ncol(newPhe)==2)
            {
              pheChar<-as.matrix(newPhe[,2:ncol(newPhe)])  
            }else if(ncol(newPhe)>2){
              pheChar<-newPhe[,2:ncol(newPhe)]
            }
            
            chrRaw_name<-unique(mapRaw10[,2])
            chr_name<-chr_name
            mx<-matrix(as.numeric(genChar),nrow(genChar),ncol(genChar))
            phe<-matrix(as.numeric(pheChar),nrow(pheChar),ncol(pheChar))
            mapname<-newMap
            
            envim$mapname<-mapname
            
            envim$chr_name<-chr_name
          }
          
        }
        gentwo<-mx[,1:2]
        t_gen<-mx[,3:ncol(mx)]
        sumphe<-apply(phe, 1,sum)
        deletRow<-which(is.na(sumphe)==TRUE)
        
        if(length(deletRow)>0){
          
          deletRow<-which(is.na(sumphe)==TRUE)
          phe<-as.matrix(phe[-deletRow,])
          t_gen1<-t_gen[,-deletRow]
          yygg<-yygg[-deletRow,]
          phe<-phe
          mx<-cbind(gentwo,t_gen1)
          
        }else{
          mx<-mx
          phe<-phe
          yygg<-yygg
        }
        
        envim$mx<-mx
        
        
        mx<-as.matrix(mx)
        phe<-as.matrix(phe[,envim$ii])
        
        map<-mx[,1:2]
        geno<-t(mx[,3:(ncol(mx))])
        n_sam<-nrow(geno)
        
        
        if(is.null(yygg)==FALSE){
          fx<-cbind(matrix(1,n_sam,1),yygg)
        }
        if(is.null(yygg)==TRUE){
          fx<-matrix(1,n_sam,1)
        }
        
        
        gg1<-1
        gg2<--1
        gg0<-99
        
        if((envim$flagCSSL)==0)
        {
          ##BC DH population
          r_tr<-matrix(0,2,2)
          tmatr<-function(r,flagRIL)
          {if (flagRIL==0)
          {r_tr[1,1:2]<-cbind((1-r),r)
          r_tr[2,1:2]<-cbind(r,(1-r))
          }
            if (flagRIL==1)
            {r_tr[1,1:2]<-cbind(1/(1+2*r),2*r/(1+2*r))
            r_tr[2,1:2]<-cbind(2*r/(1+2*r),1/(1+2*r))
            }
            return(r_tr)
          }
          
          
          con_pro<-function(j_1,abc,chr,uu1,uu2,geno,dd2,r01,r02,r0t,flagRIL,cl,gg1,gg2,gg0)
          {az1<-matrix(1,1,2)
          nni<-0
          for (itt in uu1:(uu2-1))
          {nni<-nni+1
          dd1<-diag(1,2,2)
          if (geno[j_1,itt]==gg1){dd1[-1,-1]<-0}
          if (geno[j_1,itt]==gg2){dd1[-2,-2]<-0}
          if (geno[j_1,itt]==gg0){dd1<-dd1*0.5}
          if ((abs(itt-abc))<1e-6)
          {az1<-az1%*%dd1
          r_tr<-tmatr(r01,flagRIL)
          az1<-az1%*%r_tr%*%dd2
          r_tr<-tmatr(r02,flagRIL)
          az1<-az1%*%r_tr
          }
          if ((abs(itt-abc))>1e-6)
          {az1<-az1%*%dd1
          r_tr<-tmatr(r0t[itt-uu1+1],flagRIL)
          az1<-az1%*%r_tr
          }
          }
          dd1<-diag(1,2,2)
          if (geno[j_1,uu2]==gg1){dd1[-1,-1]<-0}
          if (geno[j_1,uu2]==gg2){dd1[-2,-2]<-0}
          if (geno[j_1,uu2]==gg0){dd1<-dd1*0.5}
          az1<-az1%*%dd1
          sw<-az1%*%matrix(1,(ncol(az1)),1)
          return(sw)
          }
          
          
          mapinsert<-function(map,cl)
          {k<-0
          k1<-0
          mp<-numeric()
          for (ichr in 1:nrow(as.matrix(unique(map[,1]))))
          { 
            
            q1<-as.matrix(which(map[,1]==ichr))
            for (i in 2:(nrow(q1)))
            {rr<-map[q1[i],2]-map[q1[i-1],2]
            ll<-floor(rr/cl)
            q2<-rr-ll*cl
            if (q2>0){ll<-ll+1} 
            ss<-rr/ll
            k<-k+1
            for (j in 1:ll)
            {k1<-k1+1
            q3<-cbind((map[q1[i-1],2]+(j-1)*ss),ichr,(i-1),k,k1)
            mp<-rbind(mp,q3)
            }
            }
            k1<-k1+1
            q4<-cbind(map[q1[nrow(q1)],2],ichr,(nrow(q1)-1),k,k1)
            mp<-rbind(mp,q4)
          }
          return(mp)
          }
          
          
          markerinsert<-function(mp,geno,map,cl,gg1,gg2,gg0)
          {markerpart<-NULL
          for (ire in 1:nrow(mp))
          {   

            if (mp[ire,2]==1)
            {wp1<-1
            wp2<-nrow(as.matrix(which(map[,1]==1)))
            }
            if (mp[ire,2]>1)
            {wp1<-nrow(as.matrix(which(map[,1]<=(mp[ire,2]-1))))+1
            wp2<-nrow(as.matrix(which(map[,1]<=mp[ire,2])))
            }
            ichr<-mp[ire,2]
            q5<-as.matrix(which(map[,1]==ichr))
            r0t<-matrix(0,1,q5-1)
            for (ii in 1:(nrow(q5)-1))
            {r0t[ii]<-(1-exp(-2*(map[q5[ii+1],2]-map[q5[ii],2])/100))/2
            }
            r01<-(1-exp(-2*(mp[ire,1]-map[q5[mp[ire,3]],2])/100))/2
            r02<-(1-exp(-2*(-mp[ire,1]+map[q5[mp[ire,3]+1],2])/100))/2
            abc<-mp[ire,4]+mp[ire,2]-1
            ac4<-matrix(0,1,2)
            cn00<-matrix(0,1,2)
            dd2<-diag(1,2,2)
            dd2[-1,-1]<-0
            ac4[1]<-con_pro(1,abc,mp[ire,2],wp1,wp2,geno,dd2,r01,r02,r0t,flagRIL,cl,gg1,gg2,gg0)
            dd2<-diag(1,2,2)
            dd2[-2,-2]<-0
            ac4[2]<-con_pro(1,abc,mp[ire,2],wp1,wp2,geno,dd2,r01,r02,r0t,flagRIL,cl,gg1,gg2,gg0)
            cn00[1,1:2]<-ac4/sum(ac4)
            markerpart<-cbind(markerpart,cn00[1,1]-cn00[1,2])
          }
          return(markerpart)
          }
          
          
          mp<-mapinsert(map,cl)
          
          nq<-nrow(mp)
          mapp<-map
          mpp<-mp
          genoo<-geno
          
          markerall<-matrix(0,nrow(genoo),nrow(mpp))
          for (i in 1:nrow(as.matrix(unique(mapp[,1]))))
          {
            
            progress_bar$setFraction((2+((30/nrow(as.matrix(unique(mapp[,1]))))*i))/100)
            
            location<-as.matrix(which(mapp[,1]==i))
            pmap<-mapp[location,]
            pgeno<-genoo[,location]
            for (j in 1:(ncol(pgeno)-1))
            {ppmap<-pmap[j:(j+1),]
            ppgeno<-pgeno[,j:(j+1)]
            ppmap[,1]<-matrix(1,nrow(ppmap),1)
            map<-ppmap
            mp<-mapinsert(map,cl)
            ploc<-as.matrix(which(mpp[,2]==i & mpp[,1]>=map[1,2] & mpp[,1]<=map[2,2]))
            
            if (nrow(mp)>2)
            {for (ii in 1:n_sam)
            {geno<-t(as.matrix(ppgeno[ii,]))
            if ((geno[1]!=gg0) && (geno[2]!=gg0))
            {markerall[ii,ploc]<-markerinsert(mp,geno,map,cl,gg1,gg2,gg0)
            }
            }
            }
            
            if (nrow(mp)==2)
            {w1<-which(ppgeno==gg1)
            w2<-which(ppgeno==gg2)
            ppgeno2<-ppgeno
            ppgeno2[w1]<-1
            ppgeno2[w2]<--1
            markerall[,ploc]<-ppgeno2
            }
            
            
            }
          }
          
          map<-mapp
          mp<-mpp
          geno<-genoo
          
          hm0<-matrix(0,nrow(geno),ncol(geno))
          for (i in 1:nrow(as.matrix(unique(map[,1]))))
          {
            progress_bar$setFraction((32+((10/nrow(as.matrix(unique(map[,1]))))*i))/100)
            hm<-as.matrix(which(map[,1]==i))
            for (j in 1:(nrow(hm)-1))
            {for (ii in 1:n_sam)
            {hmm<-as.matrix(cbind(geno[ii,hm[j]],geno[ii,hm[j+1]]))
            if (nrow(as.matrix(which(hmm==gg0)))==1)
            {if (as.matrix(which(hmm==gg0))==1)
            {hm0[ii,hm[j]]<-2
            hm0[ii,hm[j+1]]<-1
            }
              if (as.matrix(which(hmm==gg0))==2)
              {hm0[ii,hm[j]]<-1
              hm0[ii,hm[j+1]]<-2
              }
            }
            if (nrow(as.matrix(which(hmm==gg0)))==2)
            {hm0[ii,hm[j]]<-2
            hm0[ii,hm[j+1]]<-2
            }
            
            }
            }
          }
          
          
          for (i in 1:n_sam)
          {
            progress_bar$setFraction((42+(20/n_sam)*i)/100)
            for (j in 1:nrow(as.matrix(unique(mapp[,1]))))
            {am<-as.matrix(which(mapp[,1]==j))
            pos0<-mapp[am,]
            hmm0<-t(as.matrix(hm0[i,am]))
            if (nrow(as.matrix(which(hmm0==1)))>0)
            {for (ii in 1:(nrow(as.matrix(which(hmm0==1)))+1))
            {amm<-as.matrix(which(hmm0==1))
            if (ii==1)
            {if (nrow(as.matrix(which(hmm0[1,1:amm[1]]==2)))>0)
            {geno<-t(as.matrix(genoo[i,am[1]:am[amm[1]]]))
            map<-pos0[1:amm[1],]
            loc<-as.matrix(which(mpp[,2]==map[1,1] & mpp[,1]>=map[1,2] & mpp[,1]<=map[nrow(map),2]))
            map[,1]<-map[,1]-map[1,1]+1
            mp<-mapinsert(map,cl)
            markerall[i,loc]<-markerinsert(mp,geno,map,cl,gg1,gg2,gg0)
            }
            }
            if (ii>1 & ii<(nrow(as.matrix(which(hm0[i,am]==1)))+1))
            {if (nrow(as.matrix(which(hmm0[1,amm[ii-1]:amm[ii]]==2)))>0)
            {geno<-t(as.matrix(genoo[i,am[amm[ii-1]]:am[amm[ii]]]))
            map<-pos0[amm[ii-1]:amm[ii],]
            loc<-as.matrix(which(mpp[,2]==map[1,1] & mpp[,1]>=map[1,2] & mpp[,1]<=map[nrow(map),2]))
            map[,1]<-map[,1]-map[1,1]+1
            mp<-mapinsert(map,cl)
            markerall[i,loc]<-markerinsert(mp,geno,map,cl,gg1,gg2,gg0)
            }
            }
            if (ii==(nrow(as.matrix(which(hm0[i,am]==1)))+1))
            {if (nrow(as.matrix(which(hmm0[1,amm[nrow(amm)]:ncol(hmm0)]==2)))>0)
            {geno<-t(as.matrix(genoo[i,am[amm[nrow(amm)]]:am[nrow(am)]]))
            map<-pos0[amm[nrow(amm)]:nrow(pos0),]
            loc<-as.matrix(which(mpp[,2]==map[1,1] & mpp[,1]>=map[1,2] & mpp[,1]<=map[nrow(map),2]))
            map[,1]<-map[,1]-map[1,1]+1
            mp<-mapinsert(map,cl)
            markerall[i,loc]<-markerinsert(mp,geno,map,cl,gg1,gg2,gg0)
            }
            }
            }
            }
            }
          }
          map<-mapp
          mp<-mpp
          geno<-genoo
          gen<-cbind(mpp[,2],mpp[,1],t(markerall))
        }
        
        if(envim$flagCSSL==1)
        {
          enabled(framedit_inter)<-FALSE
          pl<-which(geno==gg0)
          geno[pl]<-0
          gen<-cbind(map,t(geno))
        }
        
        
        ori<-NULL
        for (j in 1:(nrow(map)))
        {
          progress_bar$setFraction((62+(10/nrow(map))*j)/100)
          ta<-as.matrix(which(gen[,1]==map[j,1]))
          tb<-gen[ta,2]
          cori<-ta[as.matrix(which(tb==map[j,2])),1]
          ori<-rbind(ori,cori)
        }
        corie<-matrix(0,nrow(gen),1)
        corie[ori,1]<-matrix(1,nrow(ori),1)
        gen<-cbind(gen[,1:2],corie,gen[,3:(ncol(gen))])
        
        wg<-gen
        
        iw<-as.matrix(which(gen[,3]==1))
        gk<-t(gen[iw,4:(ncol(gen))])
        m<-ncol(gk)
        n<-nrow(gk)
        kk<-matrix(0,n,n)
        for(k in 1:m){
          z<-as.matrix(gk[,k])
          kk<-kk+z%*%t(z)
        }
        cc<-mean(diag(kk))
        kk<-kk/cc
        gen<-cbind(gen[,1:2],gen[,4:(ncol(gen))])
        
        
        ##fixed model
        fix<-function(x,gen,y,kk){
          
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
            if(abs(min(eigen(xx)$values))<1e-6)
              loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx+diag(ncol(xx))*0.01)%*%yx)-0.5*log(det(xx+diag(ncol(xx))*0.01))
            else
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
            if(abs(min(eigen(xx)$values))<1e-6)
              beta<-solve(xx+diag(ncol(xx))*0.01,yx)
            else
              beta<-solve(xx,yx)
            if(abs(min(eigen(xx)$values))<1e-6)       
              sigma2<-(yy-t(yx)%*%solve(xx+diag(ncol(xx))*0.01)%*%yx)/(n-q)
            else
              sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
            sigma2<-drop(sigma2)
            if(abs(min(eigen(xx)$values))<1e-6)
              vertue<-solve(xx+diag(ncol(xx))*0.01)
            else
              vertue<-solve(xx)
            var<-diag(vertue*sigma2)
            stderr<-sqrt(var)
            return(c(beta,stderr,sigma2))
          }
          qq<-eigen(as.matrix(kk))
          delta<-qq[[1]]
          uu<-qq[[2]]    
          qx<-ncol(x)
          n<-length(y)
          yu<-t(uu)%*%y
          
          tempx<-x
          parmm<-NULL
          for (i in 1:nrow(gen)) 
          {x<-tempx
          z<-gen[i,3:(ncol(gen))]
          qz<-ncol(z)  
          x<-cbind(x,z)
          q<-ncol(x)
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
          poly.lod<-lrt/4.61
          poly.p<-1-pchisq(lrt,1)
          sigma2g<-lambda*sigma2
          g<-beta[-c(1:qx)]
          g.err<-stderr[-c(1:qx)]
          b<-beta[c(1:qx)]
          b.err<-stderr[c(1:qx)]
          wald<-g^2/g.err^2
          p<-1-pchisq(wald,1)
          parmm<-rbind(parmm,c(b[1],sigma2,lambda,sigma2g,poly.lod,poly.p,g,g.err,wald,p))
          }
          return(parmm)
        }
        
        
        ## random model
        random<-function(fx,gen,phe,kk)
        {
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
          x<-fx
          
          s<-ncol(x)
          kk<-as.matrix(kk)
          qq<-eigen(kk)
          delta<-qq[[1]]
          uu<-qq[[2]]
          xu<-t(uu)%*%x
          y<-as.matrix(phe)
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
            parm0<-c(k,name[k,1],name[k,2],beta[1],sigma2,sigma2g,gamma,stderr,wald,p_wald)
            qq<-rbind(qq,parm0)
          }
          parms<-qq
          return(parms)
        }
        
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
              ij<-which(sub!=sub[(i+ncol(xxn))])
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
            wang[i]<-p
          }
          return (wang)
        }
        
        if (envim$flag==1)
        {code<-random(fx=fx,gen=gen,phe=phe,kk=kk)
        }
        if (envim$flag==0)
        {code<-fix(x=fx,gen=gen,y=phe,kk=kk)
        }
        
        
        x0<-t(gen[,3:ncol(gen)])
        y<-phe
        bb<-code
        bb<-as.matrix(bb)
        aa<-numeric()
        for (i in 1:nrow(as.matrix(unique(gen[,1])))){
          progress_bar$setFraction((72+(10/nrow(as.matrix(unique(gen[,1]))))*i)/100)
          mc<-which(gen[,1]==i)
          mc<-as.matrix(mc)
          for (j in 1:(nrow(mc)-2))
          {if (bb[mc[j+1],10]<bb[mc[j],10] & bb[mc[j+1],10]<bb[mc[j+2],10]){aa<-rbind(aa,mc[j+1])}
          }
          if (bb[mc[1],10]<bb[mc[2],10]){aa<-rbind(aa,mc[1])}
          if (bb[mc[nrow(mc)],10]<bb[mc[nrow(mc)-1],10]){aa<-rbind(aa,mc[nrow(mc)])}
        }
        
        
        xx<-x0[,aa]
        par<-ebayes_EM(fx,xx,y)
        cc<-which(par[,1]<=0.01)
        name<-as.matrix(aa[cc,1])
        xxx<-as.matrix(x0[,name])
        y<-as.matrix(y)
        lod<-likelihood(fx,xxx,y)
        dd<-as.matrix(which(lod[,1]>=sLOD))
        na<-name[dd]
        wow<-cbind(fx,xxx[,dd])
        bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
        ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
        galaxy<-cbind(gen[na,1],gen[na,2],ef,lod[dd,])
        
        c1<-as.matrix(galaxy)
        
        xx1<-numeric()
        for (i in 1:nrow(c1)){
          ng1<-as.matrix(which(gen[,1]==c1[i,1]))
          ng2<-gen[ng1,]
          ng3<-as.matrix(which(ng2[,2]==c1[i,2]))
          xx1<-rbind(xx1,ng2[ng3,])
        }
        
        xx1<-as.matrix(xx1)
        x1<-xx1[,3:(ncol(xx1))]
        x1<-t(x1)
        
        wow<-cbind(fx,x1)
        bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
        ef<-as.matrix(bbbb[(ncol(fx)+1):nrow(bbbb),1])
        y<-y-x1%*%ef

        
        if (envim$flag==1)
        {code<-random(fx=fx,gen=gen,phe=y,kk=kk)
        }
        if (envim$flag==0)
        {code<-fix(x=fx,gen=gen,y=y,kk=kk)
        }
        
        x0<-t(gen[,3:(ncol(gen))])
        bb<-code
        bb<-as.matrix(bb)
        aa<-numeric()
        for (i in 1:nrow(as.matrix(unique(gen[,1])))){
          progress_bar$setFraction((82+(8/nrow(as.matrix(unique(gen[,1]))))*i)/100)
          mc<-which(gen[,1]==i)
          mc<-as.matrix(mc)
          for (j in 1:(nrow(mc)-2))
          {if (bb[mc[j+1],10]<bb[mc[j],10] & bb[mc[j+1],10]<bb[mc[j+2],10]){aa<-rbind(aa,mc[j+1])}
          }
          if (bb[mc[1],10]<bb[mc[2],10]){aa<-rbind(aa,mc[1])}
          if (bb[mc[nrow(mc)],10]<bb[mc[nrow(mc)-1],10]){aa<-rbind(aa,mc[nrow(mc)])}
        }
        
        mi<-code[aa,2:3]
        
        style<-numeric()
        for (i in 1:nrow(mi))
        {
          progress_bar$setFraction((90+(3/nrow(mi))*i)/100)
          for (j in 1:nrow(xx1))
          {if (mi[i,1]==xx1[j,1] & mi[i,2]==xx1[j,2])
          {style<-rbind(style,aa[i])
          }
          }
        }
        aa<-as.matrix(setdiff(aa,style))
        xx<-x0[,aa]
        par<-ebayes_EM(fx,xx,y)
        cc<-as.matrix(which(par[,1]<=0.01))
        if (nrow(cc)>0)
        {
          name<-as.matrix(aa[cc,1])
          xxx<-as.matrix(x0[,name])
          y<-as.matrix(y)
          lod<-likelihood(fx,xxx,y)
          dd<-as.matrix(which(lod[,1]>=sLOD))
          if (nrow(dd)>0)
          {na<-as.matrix(name[dd])
          wow<-cbind(fx,xxx[,dd])
          bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
          ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
          galaxy<-cbind(gen[na,1],gen[na,2],ef,lod[dd,])
          
          pp<-as.matrix(phe[,1])
          x2<-as.matrix(cbind(x1,xxx[,dd]))
          na2<-rbind(xx1[,1:2],galaxy[,1:2])
          lod1<-likelihood(fx,x2,pp)
          dd1<-which(lod1[,1]>=sLOD)
          
          woww<-cbind(fx,x2[,dd1])
          bbbb<-solve(t(woww)%*%woww)%*%t(woww)%*%pp
          eff<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
          galaxyy<-cbind(na2[dd1,],eff,lod1[dd1,])
          }
          if (nrow(dd)==0)
          {galaxyy<-galaxy
          woww<-wow
          }
        }
        
        if (nrow(cc)==0)
        {galaxyy<-galaxy
        woww<-wow
        }
        progress_bar$setFraction(95/100)
        pp<-as.matrix(phe[,1])
        
        va<-galaxyy[,3]*galaxyy[,3]
        ve<-(1/(n_sam-1))*t(pp-woww%*%bbbb)%*%(pp-woww%*%bbbb)
        vp<-(1/(n_sam-1))*t(pp-mean(pp))%*%(pp-mean(pp))
        vy<-(sum(va)+ve)
        
        if (vy>=vp){
          heredity<-va/vy
          pv<-vy}
        if (vy<vp){
          heredity<-va/vp
          pv<-vp}
        
        va<-matrix(va,,1)
        heredity<-100*heredity
        heredity<-matrix(heredity,,1)
        
        if(is.null(mapname)==FALSE){
          map<-as.numeric(mapname[,2:3])
          map<-matrix(map,nrow(mapname),2)
          galaxytwo<-galaxyy[,1:2]
          left_marker<-numeric()
          right_marker<-numeric()
          
          
          for( i in 1:nrow(galaxyy)){
            
            allchr<-as.vector(map[which(map[,1]==galaxytwo[i,1]),2])
            chr_loc<-which(map[,1]==galaxytwo[i,1])
            allmarker<-mapname[chr_loc,1]
            chose_left<-(map[which(map[,1]==galaxytwo[i,1]),2]<galaxytwo[i,2])
            max_left<-max(allchr[chose_left[]==TRUE])
            chr_loclen<-length(chr_loc)
            if(max_left==-Inf){
              
              leftmarker<-matrix(mapname[chr_loc[1],1],,1)
            }else{
              leftloc<-which(allchr[]==max_left)
              
              leftmarker<-matrix(allmarker[leftloc],,1)
            }
            
            chose_right<-(map[which(map[,1]==galaxytwo[i,1]),2]>galaxytwo[i,2])
            min_right<-min(allchr[chose_right[]==TRUE])
            if(min_right==Inf){
              rightmarker<-matrix(mapname[chr_loc[chr_loclen],1],,1)
            }else{
              rightloc<-which(allchr[]==min_right)
              
              rightmarker<-matrix(allmarker[rightloc],,1)
            }
            
            
            left_marker<-rbind(left_marker,leftmarker)
            right_marker<-rbind(right_marker,rightmarker)
          }
          
        }else{
          left_marker<-matrix("------",nrow(galaxyy),1)
          right_marker<-matrix("------",nrow(galaxyy),1) 
        }
        
        
        
        
        if((is.null(chrRaw_name)==FALSE)&&(is.null(chr_name)==FALSE)){
          
          chr_name<-chr_name
          chrRaw_name<-chrRaw_name
          galaxyysec<-galaxyy[,1]
          galaxyylast<-galaxyy[,2:ncol(galaxyy)]
          
          
          chrName<-numeric()
          for( i in 1:length(galaxyysec)){
            chrLoc<-which(chr_name[]==galaxyysec[i])
            chrName0<-matrix(chrRaw_name[chrLoc],,1)
            chrName<-rbind(chrName,chrName0)
          }
          
          galaxyy_A<-cbind(chrName,galaxyylast)
        }else{
          galaxyy_A<-galaxyy
        }
        
        galaxyy_A<-as.matrix(galaxyy_A)
        vee<-matrix("",nrow(galaxyy_A),1)
        
        vee[1,1]<-ve
        
        vee<-matrix(vee,,1)
        vpp<-matrix("",nrow(galaxyy_A),1)
        
        vpp[1,1]<-pv
        vpp<-matrix(vpp,,1)
        traitid<-matrix(envim$ii,nrow(galaxyy_A),1)
        
        envim$galaxyy<-galaxyy
        envim$result<-cbind(traitid,galaxyy_A,left_marker,right_marker,va,heredity,vee,vpp)
        colnames(envim$result)<-c("TraitID","Chr","Position (cM)","Additive Effect","LOD","Left_Marker","Right_Marker","Var_Genet_(i)","r2 (%)","Var_Error",
                                  "Var_Phen (total)")
        
        
        if(is.null(envim$result)==TRUE){
          gmessage("There is no result meets the requirements in the last step!","Info",icon="info")
          return
        }else{
          tbdfe8<-gdfedit(envim$result,container=nb1,expand=TRUE,label="Result")
        }
        progress_bar$setFraction(100/100)
        progress_bar$setText("All done.")
      }
      return
    })  
    
    addHandlerClicked(clearbt,handler=function(h,...){
      progress_bar$setFraction(0/100)
      progress_bar$setText("")
      mm<-length(nb1)
      for(i in 1:(mm-1)){
        if(is.null(envim$genRaw)==FALSE){
          genRaw <- NULL
          genRaw1 <- NULL
          flag <- NULL
          flagRIL <- NULL
          flagCSSL <- NULL
          rm(genRaw,envir = as.environment(envim))
          rm(genRaw1,envir = as.environment(envim))
          rm(flag,envir = as.environment(envim))
          rm(flagRIL,envir = as.environment(envim))
          rm(flagCSSL,envir = as.environment(envim))
        }else if(is.null(envim$pheRaw)==FALSE){
          pheRaw <- NULL
          rm(pheRaw,envir = as.environment(envim))
          
        }else if(is.null(envim$mapRaw1)==FALSE){
          mapRaw1 <- NULL
          mapRaw11 <- NULL
          rm(mapRaw1,envir = as.environment(envim))
          rm(mapRaw11,envir = as.environment(envim))
        }else if(is.null(envim$cov_en)==FALSE){
          yygg <- NULL
          cov_en <- NULL
          yygg1 <- NULL
          flagcov <- NULL
          rm(yygg,envir = as.environment(envim))
          rm(cov_en,envir = as.environment(envim))
          rm(yygg1,envir = as.environment(envim))
          rm(flagcov,envir = as.environment(envim))
        }
        
        else if(is.null(envim$result)==FALSE){
          result <- NULL
          chr_name <- NULL
          galaxyy <- NULL
          chrRaw_name <- NULL
          mapname <- NULL
          mx <- NULL
          ii <- NULL
          rm(result,envir = as.environment(envim))
          rm(chr_name,envir = as.environment(envim))
          rm(galaxyy,envir = as.environment(envim))
          rm(chrRaw_name,envir = as.environment(envim))
          rm(mapname,envir = as.environment(envim))
          rm(mx,envir = as.environment(envim))
          rm(ii,envir = as.environment(envim))
        } else if(is.null(envim$y_jun3)==FALSE){
          y_jun00 <- NULL
          flag <- NULL
          flagRIL <- NULL
          flagCSSL <- NULL
          y_jun3 <- NULL
          rm(y_jun00,envir = as.environment(envim))
          rm(flag,envir = as.environment(envim))
          rm(flagRIL,envir = as.environment(envim))
          rm(flagCSSL,envir = as.environment(envim))
          rm(y_jun3,envir = as.environment(envim))
        } else if(is.null(envim$geo)==FALSE){
          geo <- NULL
          rm(geo,envir = as.environment(envim))
        }else if(is.null(envim$pos)==FALSE){
          flag <- NULL
          flagRIL <- NULL
          flagCSSL <- NULL
          pos <- NULL
          rm(flag,envir = as.environment(envim))
          rm(flagRIL,envir = as.environment(envim))
          rm(flagCSSL,envir = as.environment(envim))
          rm(pos,envir = as.environment(envim))
        }else if(is.null(envim$pho)==FALSE){
          pho <- NULL
          rm(pho,envir = as.environment(envim))
        }
        dispose(nb1)
      }
      
    })
    
    
    
    fusave<-function(h,...){
      
      result<-envim$result
      if(is.null(result)==TRUE)
      {
        gmessage("There is no result meets the requirements!","Info",icon="info")
        return
      }else{
        output<-gfile(text="Save a file...",type="save",
                      filter=list("All files"=list(patterns=c("*")),
                                  "CSV files"=list(patterns=c("*.csv"))))
        write.table(result,output,sep = ",",row.names=FALSE,col.names = TRUE) 
      }
    }
    
    
    
    fuplot<-function(h,...){
      if(is.null(envim$mx)==TRUE){
        gmessage("Please input correct genotype data !","Warning",icon="warning")
        return
      }
      if(is.null(envim$galaxyy)==TRUE){
        gmessage("There is no result meets the requirements!","Info",icon="info")
        return
      }
      
      if((is.null(envim$mx)==FALSE)&&(is.null(envim$galaxyy)==FALSE)&&(is.null(envim$chr_name)==FALSE)){
        
        plotwin<-gwindow("Plot of LOD Score against Genome Postion",visible=FALSE,width=960,height=220)
        gpw<-ggroup(container=plotwin,horizontal = FALSE,spacing = 10)
        ggpw<-ggraphics(container=gpw)
        addSpring(gpw)
        plotlyt<-glayout(container =gpw )
        plotfram<-gframe(container =  plotlyt)
        plotlabel<-glabel("Width (mm):",container = plotlyt)
        plotedit<-gedit("960",container = plotlyt)
        widvalue<-svalue(plotedit)
        plotlabel1<-glabel("Height (mm):",container = plotlyt)
        plotedit1<-gedit("240",container = plotlyt)
        
        plotlabel2<-glabel("Precision (dpi):",container = plotlyt)
        plotedit2<-gedit("300",container = plotlyt)
        prervalue<-svalue(plotedit2)
        plotbt<-gbutton(" Save ",container =plotlyt )
        plotlyt[1,2]<-plotlabel
        plotlyt[1,3]<-plotedit
        plotlyt[1,5]<-plotlabel1
        plotlyt[1,6]<-plotedit1
        plotlyt[1,8]<-plotlabel2
        plotlyt[1,9]<-plotedit2
        plotlyt[2,3:4]<-plotbt
        visible(plotwin) <- TRUE
        addHandlerChanged(ggpw, handler=function(h,...) {
          chr_pos <- envim$mx[,1:2]
          chr_num <- length(envim$chr_name)
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
          
          newres_pos <- envim$galaxyy[,2]
          res_sumpos <- pos_acc[envim$galaxyy[which(envim$galaxyy[,1]>1),1]-1] + envim$galaxyy[which(envim$galaxyy[,1]>1),2] 
          newres_pos[which(envim$galaxyy[,1]>1)] <- res_sumpos 
          pospic<-c(newres_pos)
          lodpic<-c(envim$galaxyy[,4])
          resdf <- data.frame(pospic,lodpic)
          envim$pp <- ggplot(data=resdf, aes(x=pospic, y=lodpic)) +
            geom_bar(stat="identity", width=0.5, fill="white", linetype="solid",color="black")
          envim$pp <- envim$pp + geom_vline(xintercept=c(0,pos_acc),linetype="dashed",alpha=0.2)
          envim$pp <- envim$pp  + scale_x_continuous(expand=c(0,0),limits=c(0,pos_acc[dim(pos_acc)[1]])) +
            scale_y_continuous(expand=c(0,0))
          envim$pp <- envim$pp + xlab("Genome position (cM)") + ylab("LOD score") + ggtitle("") + theme_classic()
          envim$pp <- envim$pp + theme(axis.title.y = element_text( vjust = 2,hjust=0.5),
                                       axis.title.x = element_text(vjust = -2,hjust=0.5)) 
          envim$pp <- envim$pp + coord_fixed(ratio=16)
          envim$pp <- envim$pp + theme(panel.background = element_rect(fill = "white"))
          envim$pp <- envim$pp + theme(text=element_text(family="mono"))
          envim$pp <- envim$pp + theme(axis.line.y = element_line(colour = "black", linetype = "solid"),
                                       axis.line.x = element_line(colour = "black", linetype = "solid"))
          print(envim$pp)
        })
        
        addhandlerclicked(plotbt,handler = function(h,...){
          widvalue<-as.numeric(svalue(plotedit))
          heightvalue<-as.numeric(svalue(plotedit1))
          prervalue<-as.numeric(svalue(plotedit2))
          if(is.null(envim$pp)==TRUE)
          {
            gmessage("There is no result meets the requirements!","Info",icon="info")
            return
          }else{
            output<-gfile(text="Save a file...",type="save",
                          filter=list("All files"=list(patterns=c("*")),
                                      "TIF files"=list(patterns=c("*.tif"))))
            plot(envim$pp)
            ggsave(output,width=widvalue,height=heightvalue,units="mm",limitsize = FALSE,dpi = prervalue)
          }
        })
      }
    }
    
    
    fuquit<-function(h,...){
      gconfirm("Yes or no ?",handler=function(h,...){dispose(window)})
    }
    
    fuhelp<-function(h,...){
      RShowDoc("Instruction1.pdf",package="mrMLM")
      
    }
    
    aICIM<- gaction(label="QTLIciMaping",handler=fuicim)
    aWINQTL<- gaction(label="WinQTLCart",handler=fuwinqtl)
    aGCIM<-gaction(label = "GCIM",handler = fugcim)
    aSave <- gaction(label="Save",handler=fusave,icon = "save")
    aQuit<- gaction(label="Quit", handler=fuquit,icon = "quit")
    ahelp<-gaction(label = "User__Manual",handler = fuhelp)
    
    aqqPlot<-gaction(label = "Plot",handler = fuplot,key.accel="Ctrl+P",icon = "plot")
    
    
    ml<-list(File=list(
      
      save=aSave,
      sep=list(separator=TRUE),
      quit=aQuit),
      
      Data_Format=list(
        GCIM_Format=aGCIM,
        sep=list(separator=TRUE),
        icim=aICIM,
        sep=list(separator=TRUE),
        winqtl=aWINQTL
      ),
      
      Graphic=list(
        qqPlot=aqqPlot),
      Help=list(
        help=ahelp 
      )
    )
    gmenu(ml,container=window)
    
  }

  addHandlerClicked(GCIMbutton,handler=function(h,...){
    
    GCIM()
    
  })
  
}






