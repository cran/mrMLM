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
inte_window<-gwindow(title="Methodologies for multi-locus GWAS and multi-locus linkage analysis",visible=TRUE,width=578,height=300,expand=TRUE)
inte_lyt<-glayout(container=inte_window,spacing=13)
mrMLMbutton<-gbutton("multi-locus random-SNP-effect Mixed Linear Model (mrMLM)",container=inte_lyt)
FASTbutton<-gbutton("FAST multi-locus random-SNP-effect EMMA (FASTmrEMMA)",container=inte_lyt)
GCIMbutton<-gbutton("Genome-wide Composite Interval Mapping (GCIM)",container=inte_lyt)
gwasmed<-glabel("Multi-locus GWAS methodologies:")
linkmed<-glabel("Multi-locus linkage analysis methodologies:")
inte_lyt[1,1,expand=TRUE]<-gwasmed
inte_lyt[2,1,expand=TRUE]<-mrMLMbutton
inte_lyt[3,1,expand=TRUE]<-glabel("Cited: Wang et al. 2016. Sci Rep 6:19444")
inte_lyt[4,1,expand=TRUE]<-FASTbutton
inte_lyt[5,1,expand=TRUE]<-glabel("Cited: Wen et al. 2016. Briefings in Bioinformatics, DOI: 10.1093/bib/bbw145")
inte_lyt[8,1,expand=TRUE]<-linkmed
inte_lyt[9,1,expand=TRUE]<-GCIMbutton
inte_lyt[10,1,expand=TRUE]<-glabel("Cited: Wang et al. 2016. Sci Rep 6:29951")
inte_lyt[13,1,expand=TRUE]<-glabel("Version 2.0, December 2016")
font(gwasmed)<-c(size="x-large")
font(linkmed)<-c(size="x-large")

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
  
  window<-gwindow(title="multi-locus random-SNP-effect Mixed Linear Model (mrMLM)",visible=TRUE,width=1260,height=730,expand=TRUE)
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
  includeps<-gwindow("Include population structure (Q matrix)?",visible=FALSE,width=400,height=150)
  gps<-ggroup(container=includeps,expand=FALSE)
  parsetwin<-gwindow("Parameter Settings",visible=FALSE,width=260,height=280)
  gpar<-ggroup(container=parsetwin,expand=FALSE)
  choicesave<-gwindow("Save as ...",visible=FALSE,width=250,height=150)
  gcsave<-ggroup(container=choicesave,expand=FALSE)
  
  lyt<-glayout(container=window,spacing=20)
  
  importdata<-gbutton("Input Datasets",container=lyt)
  parset<-gbutton("Parameter Settings",container=lyt)
  manhattan<-gbutton("Manhattan Plot",container=lyt)
  qqplot<-gbutton("QQ Plot",container=lyt)
  
  helpfile<-gbutton("User Manual",container=lyt)
  savefile<-gbutton(" Save ",container=lyt)
  run<-gbutton("Run",container=lyt)
  exit<-gbutton("Exit",container=lyt)
  
  lyt[1,1]<-importdata
  lyt[4,1]<-parset
  lyt[5,1]<-run  
  lyt[6,1]<-savefile
  lyt[9,1]<-manhattan
  lyt[10,1]<-qqplot
  lyt[13,1]<-helpfile
  lyt[14,1]<-exit
  
  nb1<-gnotebook(tab.pos=3,closebuttons=TRUE,dontCloseThese=TRUE,container=lyt,expand=TRUE)
  size(nb1)<-c(680,540)
  tb<-gnewtable("     
                1.mrMLM is a R software package for multi-locus genome-wide association studies based on a multi-locus 
                  random-SNP-effect mixed linear model.
                
                2.Please cite: Wang Shi-Bo, Feng Jian-Ying, Ren Wen-Long, Huang Bo, Zhou Ling, Wen Yang-Jun, Zhang Jin, 
                  Jim M. Dunwell, Xu Shizhong (*), Zhang Yuan-Ming (*).2016. Improving power and accuracy of genome-wide 
                  association studies via a multi-locus mixed linear model methodology. Scientific Reports 6:19444. 
                
                3.Please cite: Zhang Yuan-Ming et al. Mapping quantitative trait loci using naturally occurring genetic
                  variance among commercial inbred line of maize (Zea mays L.). Genetics 2005, 169:2267-2275.
                
                4.The software package is developed by Wen-Long Ren, Shi-Bo Wang, Bo Huang & Yuan-Ming Zhang.
                  ",multiple=TRUE,container=nb1,expand=TRUE,label="About the software")
  
  font(tb)<-c(size="x-large")
  lyt[1:14,2,expand=TRUE]<-nb1
  
  staprogress<-gtkButton()
  lyt[15,2,expand=TRUE]<-staprogress
  
  
  addHandlerClicked(importdata,handler=function(h,...){
    if(isExtant(importwin)==FALSE)
    {
      importwin<-gwindow("Input Dataset",visible=FALSE,width=250,height=420)
      gimpwin<-ggroup(container=importwin,expand=FALSE)
    }
    lytimp<-glayout(container=gimpwin,spacing=13)
    impchoose<-glabel("1. Choose dataset formats",container=lytimp)
    impfile1<-glabel("2. Input Genotypic and Phenotypic files",container=lytimp)
    impprepare<-glabel("3. Sort & Transform for datasets",container=lytimp)
    impfile2<-glabel("4. Input Kinship (K matrix) and Population-structure (Q matrix) files",container=lytimp)
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
        gmessage("Please input correct genotype dataset !","Warning",icon="warning")
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
        gmessage("Please input correct phenotype dataset !","Warning",icon="warning")
      }else{
        mrenv$pheRaw1<-as.matrix(read.csv(input2,header=FALSE)) 
        mrenv$pheRaw2<-mrenv$pheRaw1[-1,]
        mrenv$pheRaw3<-as.data.frame(mrenv$pheRaw2,stringsAsFactors=FALSE)
        mrenv$pheRaw4<-as.matrix(mrenv$pheRaw3[is.na(mrenv$pheRaw3[,2])==F,])
        mrenv$pheRawthem<-matrix(c(mrenv$pheRaw1[1,1]," "),1,)
        mrenv$pheRaw<-rbind(mrenv$pheRawthem,mrenv$pheRaw4)
        row.names(mrenv$pheRaw)<-NULL
        mrenv$pheRaw<-as.matrix(mrenv$pheRaw)
        showpheRaw<-mrenv$pheRaw1[-1,]
        colnames(showpheRaw)<-c(mrenv$pheRaw1[1,1],"   ")
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
      mrenv$radiokk<-gradio(c("Input Kinship matrix file","Calculate Kinship matrix by this software"),selected=1,horizontal=FALSE,container=lytkk)
      notekk<-glabel("Note: Please select the 1st option if no. of markers is more than 50,000.",container=lytkk)
      lytkk[2:3,2:5]<-mrenv$radiokk
      lytkk[5,2]<-mrenv$okkk
      lytkk[5,5]<-mrenv$cancelkk
      lytkk[6,1:5]<-notekk
      visible(choicekk)<-TRUE
      addHandlerClicked(mrenv$okkk,handler=function(h,...){
        if(svalue(mrenv$radiokk)=="Input Kinship matrix file"){
          input3<-gfile(text="Select a file...",type="open",
                        filter=list("All files"=list(patterns=c("*")),
                                    "CSV files"=list(patterns=c("*.csv"))))
          if(is.na(input3))
          {
            gmessage("Please input correct kinship dataset !","Warning",icon="warning")
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
            gmessage("Please input correct genotype dataset !","Warning",icon="warning")
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
        includeps<-gwindow("Include population structure (Q matrix)?",visible=FALSE,width=400,height=150)
        gps<-ggroup(container=includeps,expand=FALSE)
      }
      lytps<-glayout(container=gps,spacing=13)
      okps<-gbutton("     OK    ",container=lytps)
      cancelps<-gbutton(" Cancel ",container=lytps)
      radiops<-gradio(c("Not included in the model","Included"),selected=1,horizontal=FALSE,container=lytps)
      lytps[2:3,2:5]<-radiops
      lytps[5,2]<-okps
      lytps[5,5]<-cancelps
      visible(includeps)<-TRUE
      addHandlerClicked(okps,handler=function(h,...){
        if(svalue(radiops)=="Included"){
          mrenv$flagps<-0
          input4<-gfile(text="Select a file...",type="open",
                        filter=list("All files"=list(patterns=c("*")),
                                    "CSV files"=list(patterns=c("*.csv"))))
          if(is.na(input4))
          {
            gmessage("Please input correct population dataset !","Warning",icon="warning")
          }else{
            mrenv$psmatrixRaw<-as.matrix(read.csv(input4,header=FALSE))
            nnpprow<-dim(mrenv$psmatrixRaw)[1]
            nnppcol<-dim(mrenv$psmatrixRaw)[2]
            mrenv$psmatrixRaw[1,2:nnppcol]<-"  "
            psmatrixPre<-mrenv$psmatrixRaw[3:nnpprow,]
            namePop<-as.matrix(psmatrixPre[,1])
            sameGenPop<-intersect(mrenv$sameName,namePop)
            locPop<-match(sameGenPop,namePop)
            ##revised
            filtername<-as.vector(mrenv$psmatrixRaw[2,2:nnppcol])
            selectpsmatrix<-matrix(as.numeric(psmatrixPre[locPop,-1]),nrow = nrow(psmatrixPre))
            psum<-apply(selectpsmatrix,1,sum)
            psum<-round(psum)
            sumps<-sum(psum)
            m<-dim(selectpsmatrix)[1]
            if(sumps==m){
              filterps<-gwindow("Filter",visible=FALSE,width=300,height=120)
              filergp<-ggroup(container=filterps,expand=FALSE)
              filterlyt<-glayout(container = filergp)
              filterlabel<-glabel("Please choose the column that should be deleted in population structure Q matrix",container = filterlyt)
              filtercombo<-gcombobox(filtername,editable = TRUE,container =filterlyt)
              filterlabel1<-glabel("Filter:",container = filterlyt)
              filterok<-gbutton(" OK ",container = filterlyt)
              filtercancel<-gbutton("Cancel",container = filterlyt)
              filterlyt[2,2:5]<-filterlabel
              filterlyt[4,2]<-filterlabel1
              filterlyt[4,3]<-filtercombo
              filterlyt[5,2]<-filterok
              filterlyt[5,5]<-filtercancel
              visible(filterps)<-TRUE
              addhandlerclicked(filterok,handler = function(h,...){
                combovalue<-svalue(filtercombo)
                coldelet<-unlist(str_extract_all(combovalue,"[0-9]+"))
                coldelet<-as.numeric(coldelet)
                mrenv$psmatrix<-as.matrix(selectpsmatrix[,-coldelet])
                mrenv$psmatrixRaw<-as.matrix(mrenv$psmatrixRaw[,-(coldelet+1)])
                tbdfe6<-gdfedit(mrenv$psmatrixRaw,container=nb1,expand=TRUE,label="Population Structure")
                dispose(filterps)
                dispose(includeps)
                dispose(importwin)
              })
              addhandlerclicked(filtercancel,handler = function(h,...){
                dispose(filterps)
              })
              
            }else{
              mrenv$psmatrix<-selectpsmatrix
              tbdfe6<-gdfedit(mrenv$psmatrixRaw,container=nb1,expand=TRUE,label="Population Structure")
              
              dispose(includeps)
              dispose(importwin)
            }
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
      parsetwin<-gwindow("Parameter Settings",visible=FALSE,width=260,height=280)
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
        gmessage("Please input critical P-value between 0 and 1!","Warning",icon="warning")
      }
      if(mrenv$svrad<0)
      {
        gmessage("Please input search radius of candidated gene: > 0 !","Warning",icon="warning")
      }
      if(mrenv$svmlod<0)
      {
        gmessage("Please input critical LOD score: > 0 !","Warning",icon="warning")
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
      gmessage("Please input correct genotype dataset !","Warning",icon="warning")
    }
    if(exists("phe")==FALSE)
    {
      gmessage("Please input correct phenotype dataset !","Warning",icon="warning")
    }
    if(exists("kk")==FALSE)
    {
      gmessage("Please input correct kinship dataset !","Warning",icon="warning")
    }
    if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(ncol(gen)!=(nrow(phe)+2)))
    {
      gmessage("Sample sizes between genotypic and phenotypic datasets do not equal !","Error",icon="error")
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
      
      
      mrenv$parms.pchange<-mrenv$parms
      parmsp<-as.matrix(mrenv$parms.pchange[,10])
      locsub<-which(parmsp==0)
      if(length(locsub)!=0){
        pmin<-min(parmsp[parmsp!=0])
        subvalue<-10^(1.1*log10(pmin))
        mrenv$parms.pchange[locsub,10]<-subvalue
      }else{
        mrenv$parms.pchange<-mrenv$parms
      }
      
      progress_bar$setFraction(95/100)
      if(mrenv$inputform==1){
        #output result1 using mrMLM numeric format
        mrenv$parmsShow<-mrenv$parms[,-1]
        meadd<-matrix(1,nrow(mrenv$parms),1)
        meadd[which(mrenv$parms[,10]<newp),1]<-sprintf("%.4e",newp)
        meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
        tempparms<-mrenv$parms[,4:10]
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],mrenv$parms[,2:3],tempparms,mrenv$genRaw[-1,4],meadd)
        colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype for code 1","Significance")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      if(mrenv$inputform==2){
        #output result1 using mrMLM character format
        mrenv$parmsShow<-mrenv$parms[,-1]
        mrenv$outATCG<-matrix(mrenv$outATCG,,1)
        meadd<-matrix(1,nrow(mrenv$parms),1)
        meadd[which(mrenv$parms[,10]<newp),1]<-sprintf("%.4e",newp)
        meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
        tempparms<-mrenv$parms[,4:10]
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],mrenv$parms[,2:3],tempparms,mrenv$outATCG,meadd)
        colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype  for code 1","Significance")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      if(mrenv$inputform==3){
        #output result1 using TASSEL format
        mrenv$parmsShow<-mrenv$parms[,-1]
        mrenv$outATCG<-matrix(mrenv$outATCG,,1)
        mrenv$outATCG<-unlist(strsplit(mrenv$outATCG,""))
        mrenv$outATCG<-matrix(mrenv$outATCG[c(TRUE,FALSE)],,1)
        meadd<-matrix(1,nrow(mrenv$parms),1)
        meadd[which(mrenv$parms[,10]<newp),1]<-sprintf("%.4e",newp)
        meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
        tempparms<-mrenv$parms[,4:10]
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],mrenv$parms[,2:3],tempparms,mrenv$outATCG,meadd)
        colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Marker Position (bp)","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype  for code 1","Significance")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      
      rowsnp <- dim(mrenv$parms)[1]
      mrenv$snpname <- numeric()
      mrenv$snpname <- as.matrix(paste("rs",c(1:rowsnp),sep=""))
      
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
          her<-(er/as.numeric(sum(er)+v0))*100
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
          her<-(er/as.numeric(sum(er)+v0))*100
        }
        eeff[which(abs(eeff)>=1e-4)] <- round(eeff[which(abs(eeff)>=1e-4)],4)
        eeff[which(abs(eeff)<1e-4)] <- as.numeric(sprintf("%.4e",eeff[which(abs(eeff)<1e-4)]))
        lo[which(abs(lo)>=1e-4)] <- round(lo[which(abs(lo)>=1e-4)],4)
        lo[which(abs(lo)<1e-4)] <- as.numeric(sprintf("%.4e",lo[which(abs(lo)<1e-4)]))
        her[which(abs(her)>=1e-4)] <- round(her[which(abs(her)>=1e-4)],4)
        her[which(abs(her)<1e-4)] <- as.numeric(sprintf("%.4e",her[which(abs(her)<1e-4)]))
        mrenv$wan<-data.frame(mrenv$parmsShow[mrenv$needww,1],chr_pos[ww,],eeff,lo,her,mrenv$parmsShow[mrenv$needww,11])
        colnames(mrenv$wan)<-c("RS#","Chromosome","Marker Position (bp)","QTN effect","LOD score","r2 (%)","Genotype  for code 1")
        
      }
      wan<-mrenv$wan
      if(exists("wan")==FALSE||is.null(wan)==TRUE)
      {
        gmessage("No result meets the requirements in the second step!","Info",icon="info")
      }else{
        tbdfe8<-gdfedit(wan,container=nb1,expand=TRUE,label="Result2")
      }
      progress_bar$setFraction(100/100)
      progress_bar$setText("All done.")
    }
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
          gmessage("No result meets the requirements in the second step!","Info",icon="info")
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
    plotlabel3 <- glabel("Word resolution (1/72 inch, ppi):", container = plotlyt)
    plotedit3 <- gedit("20", container = plotlyt)
    pointsizevalue <- as.numeric(svalue(plotedit3))
    plotbt <- gbutton(" Save ", container = plotlyt)
    plotmancl <- gbutton("  Cancel  ", container = plotlyt)
    combo_box1 <- gcombobox(selected=1,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotlyt)
    combo_box2 <- gcombobox(selected=3,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotlyt)
    plotlabel4 <- glabel("Chromosome color (odd):", container = plotlyt)
    plotlabel5 <- glabel("Chromosome color (even):", container = plotlyt)
    plotlabel6 <- glabel("Figure resolution (ppi):", container = plotlyt)
    plotedit6 <- gedit("72", container = plotlyt)
    plotlabel7 <- glabel("Critical value for Manhattan Plot:",container=plotlyt)
    plotedit7 <- gwedit<-gedit(sprintf("%.6s",-log10(mrenv$mannewp)),coerce.with=as.numeric,container=plotlyt)
    plotlabel8 <- glabel("Note: The parameter settings of preview figure are for general resolution, see manual for parameter settings for high resolution figure (*.png,*.tiff,*.jpeg).", container = plotlyt)
    
    plotlyt[2, 1] <- plotlabel1
    plotlyt[2, 2] <- plotedit1
    plotlyt[2, 4] <- plotlabel2  
    plotlyt[2, 5] <- plotedit2  
    plotlyt[3, 1] <- plotlabel3
    plotlyt[3, 2] <- plotedit3
    plotlyt[3, 4] <- plotlabel6
    plotlyt[3, 5] <- plotedit6
    plotlyt[4, 1] <- plotlabel4
    plotlyt[4, 2] <- combo_box1
    plotlyt[4, 4] <- plotlabel5
    plotlyt[4, 5] <- combo_box2
    plotlyt[5, 1] <- plotlabel7
    plotlyt[5, 2] <- plotedit7
    plotlyt[5, 4] <- plotbt
    plotlyt[5, 5] <- plotmancl
    plotlyt[6, 1:5] <- plotlabel8
    
    visible(plotwin)<-TRUE
    
    bpnumber <- numeric()
    chrnum <- unique(mrenv$parms[,2])
    
    for(i in 1:length(chrnum))
    {
      bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(mrenv$parms[,2]==chrnum[i])))))
    }
    
    
    
    addHandlerChanged(ggpw, handler=function(h,...) {
      color1 <- svalue(combo_box1)
      color2 <- svalue(combo_box2)
      svgwline<-svalue(plotedit7)
      mrenv$standline<-svgwline
      
      parms <- data.frame(mrenv$parms.pchange,mrenv$snpname,bpnumber)
      colnames(parms)<-c("Trait","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","SNPname","BPnumber")
      manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P_wald",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = mrenv$standline)
    })
    addhandlerclicked(plotbt, handler = function(h, ...) {
      widvalue <- as.numeric(svalue(plotedit1))
      heightvalue <- as.numeric(svalue(plotedit2))
      pointsizevalue <- as.numeric(svalue(plotedit3))
      color1 <- svalue(combo_box1)
      color2 <- svalue(combo_box2)
      resppi <- as.numeric(svalue(plotedit6))
      svgwline<-svalue(plotedit7)
      mrenv$standline<-svgwline
      
      output <- gfile(text = "Save a file...", type = "save", 
                      filter = list(`All files` = list(patterns = c("*")), 
                                    `TIFF files` = list(patterns = c("*.tiff")),
                                    `PNG files` = list(patterns = c("*.png")),
                                    `JPEG files` = list(patterns = c("*.jpeg"))))
      if((length(grep(".png",output))==1)||(length(grep(".PNG",output))==1)){
        png(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue,res=resppi)
      }else if((length(grep(".tiff",output))==1)||(length(grep(".TIFF",output))==1)){
        tiff(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue,res=resppi)
      }else if((length(grep(".jpeg",output)==1))||(length(grep(".JPEG",output))==1)){
        jpeg(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue,res=resppi)  
      }else{
        gmessage("Please input correct image format !")
      }
      parms <- data.frame(mrenv$parms.pchange,mrenv$snpname,bpnumber)
      colnames(parms)<-c("Trait","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","SNPname","BPnumber")
      manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P_wald",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = mrenv$standline)
      dev.off()
      
    })
    addHandlerClicked(plotmancl,handler=function(h,...){
      dispose(plotwin)
    })
    
  })
  
  qqplotfun <- function(p_value,p_stand,color1,color2){
    pvalue<-matrix(p_value,,1)
    observed<-sort(pvalue[,1])
    observed<-observed/2
    observed<-observed[which(observed!=0)]
    newobserved<-observed[which(observed<(p_stand/2))]
    lobs<--(log10(newobserved))
    expected<-c(1:length(newobserved))
    lexp<--(log10(expected/(length(pvalue)+1)))
    plot(lexp,lobs,xlim=c(0,max(lexp)),ylim=c(0,max(lobs)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color1)
    abline(0,1,col=color2)
  }
  
  addHandlerClicked(qqplot,handler=function(h,...){
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
    plotqqlabel3 <- glabel("Word resolution (1/72 inch, ppi):", container = plotqqlyt)
    plotqqedit3 <- gedit("20", container = plotqqlyt)
    pointsizeqqvalue <- as.numeric(svalue(plotqqedit3))
    plotqqbt <- gbutton(" Save ", container = plotqqlyt)
    plotqqcl <- gbutton("  Cancel  ", container = plotqqlyt)
    plotqqlabel4 <- glabel("Figure resolution (ppi):", container = plotqqlyt)
    plotqqedit4 <- gedit("72", container = plotqqlyt)
    plotqqlabel5 <- glabel("Point color:", container = plotqqlyt)
    combo_box1 <- gcombobox(selected=2,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotqqlyt)
    plotqqlabel6 <- glabel("Line color:", container = plotqqlyt)
    combo_box2 <- gcombobox(selected=3,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotqqlyt)
    plotqqlabel7 <- glabel("Critical  P-value of deleting points:", container = plotqqlyt)
    plotqqedit7 <- gedit("0.95",coerce.with=as.numeric,container = plotqqlyt)
    plotqqlabel8 <- glabel("Note: The parameter settings of preview figure are for general resolution, see manual for high resolution figure (*.png,*.tiff,*.jpeg).", container = plotqqlyt)
    
    plotqqlyt[2, 1] <- plotqqlabel1
    plotqqlyt[2, 2] <- plotqqedit1
    plotqqlyt[2, 4] <- plotqqlabel2
    plotqqlyt[2, 5] <- plotqqedit2
    plotqqlyt[3, 1] <- plotqqlabel3
    plotqqlyt[3, 2] <- plotqqedit3
    plotqqlyt[3, 4] <- plotqqlabel4
    plotqqlyt[3, 5] <- plotqqedit4
    plotqqlyt[4, 1] <- plotqqlabel5
    plotqqlyt[4, 2] <- combo_box1
    plotqqlyt[4, 4] <- plotqqlabel6
    plotqqlyt[4, 5] <- combo_box2
    plotqqlyt[5, 1] <- plotqqlabel7
    plotqqlyt[5, 2] <- plotqqedit7
    plotqqlyt[5, 4] <- plotqqbt
    plotqqlyt[5, 5] <- plotqqcl 
    plotqqlyt[6, 1:5] <- plotqqlabel8
    
    visible(plotwin1)<-TRUE
    
    addHandlerChanged(ggpw1, handler=function(h,...) {
      color1 <- svalue(combo_box1)
      color2 <- svalue(combo_box2)
      svgwstandp <- svalue(plotqqedit7)
      qqplotfun(mrenv$parms.pchange[,10],svgwstandp,color1,color2)
    })
    
    addhandlerclicked(plotqqbt, handler = function(h, ...) {
      widqqvalue <- as.numeric(svalue(plotqqedit1))
      heightqqvalue <- as.numeric(svalue(plotqqedit2))
      pointsizeqqvalue <- as.numeric(svalue(plotqqedit3))
      resppi <- as.numeric(svalue(plotqqedit4))
      color1 <- svalue(combo_box1)
      color2 <- svalue(combo_box2)
      svgwstandp <- svalue(plotqqedit7)
      
      output <- gfile(text = "Save a file...", type = "save", 
                      filter = list(`All files` = list(patterns = c("*")), 
                                    `TIFF files` = list(patterns = c("*.tiff")),
                                    `PNG files` = list(patterns = c("*.png")),
                                    `JPEG files` = list(patterns = c("*.jpeg"))))
      if((length(grep(".png",output))==1)||(length(grep(".PNG",output))==1)){
        png(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue, res = resppi)
      }else if((length(grep(".tiff",output))==1)||(length(grep(".TIFF",output))==1)){
        tiff(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue, res = resppi)
      }else if((length(grep(".jpeg",output)==1))||(length(grep(".JPEG",output))==1)){
        jpeg(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue, res = resppi)  
      }else{
        gmessage("Please input correct image format !")
      }
      qqplotfun(mrenv$parms.pchange[,10],svgwstandp,color1,color2)
      dev.off()
      
    })
    
    addHandlerClicked(plotqqcl,handler=function(h,...){
      dispose(plotwin1)
    })
  })
}

addHandlerClicked(mrMLMbutton,handler=function(h,...){
  mrMLMSub()
})

GCIMSub<-function(){
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
  cartw<-gwindow("WinQTLCart Format",visible = FALSE,width = 80,height = 120,expand=FALSE)
  cartg<-ggroup(container = cartw)
  gcimw<-gwindow("GCIM Format",visible = FALSE,width = 85,height = 160,expand=FALSE)
  gcimg<-ggroup(container = gcimw)
  
  covw<-gwindow("Include Covariate?",visible = FALSE,width = 250,height = 130,expand=FALSE)
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
  
  
  size(nb1)<-c(1200,600)
  
  progress_bar <- gtkProgressBar()
  lyt[1:2,1:20,expand=TRUE]<-frame
  
  lyt[3:20,1:20,expand=TRUE]<-nb1
  lyt[21,1:20,expand=TRUE]<-progress_bar
  
  tb<-gnewtable("     
                1. GCIM is a R software package for multi-locus linkage analysis based on Genome-wide Composite Interval Mapping.
                
                2. Please cite: Wang Shi-Bo, Wen Yang-Jun, Ren Wen-Long, Ni Yuan-Li, Zhang Jin, Feng Jian-Ying, Zhang Yuan-Ming*.  
                
                   Mapping small-effect and linked quantitative trait loci for complex traits in backcross or DH populations via
                
                   a multi-locus GWAS methodology. Scientific Reports 2016, 6: 29951.
                
                3. The software package is developed by Yuan-Li Ni,Wen-Long Ren, Shi-Bo Wang & Yuan-Ming Zhang.
                
                
                   Version 2.0, Realeased Dec 2016",multiple=TRUE,container=nb1,expand=TRUE,label="About the software")
  
  
  font(tb)<-c(size="x-large")
  
  addhandlerclicked(exitbt,handler = function(h,...){
    gconfirm("Yes or No ?",handler=function(h,...){dispose(window)})
    
  })
  
  
  fuicim<-function(h,...){
    if(isExtant(icimw)==FALSE){
      icimw<-gwindow("QTLIciMaping Format",visible = FALSE,width =80,height = 140)
      icimg<-ggroup(container = icimw,expand=FALSE)
    }
    iclyt<-glayout(container = icimg,spacing = 8)
    iclabel1<-glabel("1. Please input QTLIciMaping format datasets",container = iclyt)
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
        gmessage("Please input QTLIciMaping format datasets!",title = "warning",icon="warning")
      }else{
        envim$geoo<-read.xlsx(input1,sheet = "Genotype",colNames = FALSE)
        envim$poss<-read.xlsx(input1,sheet = "LinkageMap",colNames = FALSE)
        envim$pho<-read.xlsx(input1,sheet = "Phenotype",colNames = FALSE)
        envim$parm<-read.xlsx(input1,sheet = "GeneralInfo",colNames = FALSE)
        gdf2<-gdfedit(as.data.frame(envim$geoo),expand=TRUE,container=nb1,label="Raw_Genotype")
        gdf3<-gdfedit(as.data.frame(envim$pho),expand=TRUE,container=nb1,label="Raw_phenotype") 
        gdf4<-gdfedit(as.data.frame(envim$poss),expand=TRUE,container=nb1,label="Raw_Linkage Map")
      }
    })
    
    addhandlerclicked(icok,handler = function(h,...){
      if((is.null(envim$geoo)==TRUE)&&(is.null(envim$pho)==TRUE)&&(is.null(envim$poss)==TRUE)){
        gmessage("Please input QTLIciMaping format datasets!",title = "warning",icon="warning")
      }else{
        envim$geoo<-as.matrix(envim$geoo)
        envim$poss<-envim$poss
        envim$pho<-as.matrix(envim$pho)
        geoo<-envim$geoo
        poss<-envim$poss
        pho<-envim$pho
        parm<-envim$parm
        
        if(parm[4,1]==1){
          pos.be<-numeric()
          for(i in 1:10){
            pos1<-poss[which(poss[,2]==i),]
            poss1<-pos1
            positi<-as.matrix(cumsum(poss1[,3]))
            chrr<-as.matrix(poss1[,1:2])
            poss2<-cbind(chrr,positi)
            pos.be<-rbind(pos.be,poss2)
            
          } 
        }
        if(parm[4,1]==2){
          pos.be<-poss
          
        }
        
        if(parm[5,1]==2){
          posthree<-matrix(100*pos.be[,3],,1)
          postwo<-pos.be[,1:2]
          pos<-cbind(postwo,posthree)
        }
        if(parm[5,1]==1){
          pos<-as.matrix(pos.be)
        }
        pos<-as.matrix(pos)
        geo<-geoo
        if(svalue(icradio)=="BC1 ( F1 X P1 )"&&svalue(icradiomd)=="Random Model"){
          envim$flagCSSL<-0
          envim$flag<-1 
          envim$flagRIL<-0
          gen_0<-geo[,-1]
          gen_1<-gsub("-1","99",gen_0)
          gen_2<-gsub("1","-1",gen_1)
          gen_11<-gsub("2","1",gen_2)
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
          
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
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
    }
    )
    
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
        gmessage("Please input WinQTLCart datasets!",title = "warning",icon="warning")
        return
      }else{
        envim$y_jun00<-scan(input2,what = "",sep = "\n")
        envim$y_jun3<-scan(input2,what = "",sep = "")
        tb<-gnewtable(envim$y_jun00,multiple=TRUE,container=nb1,expand=TRUE,label="WinQTLCart")
      }
    })
    
    addhandlerclicked(cartok,handler = function(h,...){
      
      if((is.null(envim$y_jun00)==TRUE)&&(is.null(envim$y_jun3)==TRUE)){
        gmessage("Please input WinQTLCart datasets!",title = "warning",icon="warning")
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
        gmessage("Please input phenotype dataset!",title = "warning",icon="warning")
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
        gmessage("Please input phenotype dataset!",title = "warning",icon="warning")
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
        gmessage("Please input linkage map dataset!",title = "warning",icon="warning")
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
        gmessage("Please input genotype dataset!",title = "warning",icon="warning")
        return
      }
      if(is.null(pheRaw)==TRUE){
        gmessage("Please input phenotype dataset!",title = "warning",icon="warning")
        return
      }
      if(is.null(mapRaw1)==TRUE){
        gmessage("Please input linkage map dataset!",title = "warning",icon="warning")
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
      covw<-gwindow("Include Covariate?",visible = FALSE,width = 250,height = 130)
      covg<-ggroup(container = covw)
    }
    lytcov<-glayout(container=covg,spacing=13)
    okcov<-gbutton("     OK    ",container=lytcov)
    cancelcov<-gbutton(" Cancel ",container=lytcov)
    radiocov<-gradio(c("Included","Not included in the model"),selected=1,horizontal=FALSE,container=lytcov)
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
          gmessage("Please input covariates dataset!",title = "warning",icon="warning")
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
      gmessage("Please input correct genotype dataset!","Warning",icon="warning")
      return
    }
    if(is.null(pheRaw)==TRUE){
      gmessage("Please input correct phenotype dataset!","Warning",icon="warning")
      return
    }
    if(is.null(mapRaw1)==TRUE){
      gmessage("Please input correct linkage map dataset!","Warning",icon="warning")
      return
    }
    if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(envim$ii>(ncol(pheRaw)-1))){
      gmessage("It is larger than the number of Traits","Warning",icon="warning")
      return
    }
    if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(envim$ii<=(ncol(pheRaw)-1))&&(cl<0)){
      gmessage("Please input Walk Speed: >0!","Warning",icon="warning")
      return
    }
    if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(envim$ii<=(ncol(pheRaw)-1))&&(cl>0)&&(sLOD<0)){
      gmessage("Please input critical LOD score: >0!","Warning",icon="warning")
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
        
        mapRaw10<-as.matrix(mapRaw1[-1,])
        
        chr_name<-unique(mapRaw10[,2])
        chr_secon<-as.matrix(mapRaw10[,2])
        
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
          phee<-as.matrix(pheRaw[-1,-1])
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
        mapRaw10<-as.matrix(mapRaw1[-1,])
        
        chr_name<-unique(mapRaw10[,2])
        chr_secon<-as.matrix(mapRaw10[,2])
        
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
          phee<-as.matrix(pheRaw[-1,-1])
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
      
      mapp1<-as.numeric(mapname[,2:3])
      mapp1<-matrix(mapp1,,2)
      chr<-length(unique(mapp1[,1]))
      ##change
      for(i in 1:chr){
        pos1<-as.matrix(mapp1[which(mapp1[,1]==i),])
        envim$delerow<-which(duplicated(pos1[,2]))
        if(length(envim$delerow)!=0){
          break
        }else{
          mapname<-mapname
          
        }
      }
      if(length(envim$delerow)!=0){
        progress_bar$setText ( "Stop Running!" )
        progress_bar$setFraction(2/100)
        gmessage("Please correct the map file to make sure whether all the marker positions are different!",title = "message",icon="error")
        
      }else{
        
        gentwo<-mx[,1:2]
        t_gen<-mx[,3:ncol(mx)]
        sumphe<-apply(phe, 1,sum)
        deletRow<-which(is.na(sumphe)==TRUE)
        
        if(length(deletRow)>0){
          
          deletRow<-which(is.na(sumphe)==TRUE)
          phe<-as.matrix(phe[-deletRow,])
          t_gen1<-t_gen[,-deletRow]
          if(is.null(yygg)==FALSE){
            yygg<-yygg[-deletRow,]
          }
          
          phe<-phe
          mx<-cbind(gentwo,t_gen1)
          
        }else{
          mx<-mx
          phe<-phe
          if(is.null(yygg)==FALSE){
            yygg<-yygg
          }
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
        tempcode<-code
        }
        if (envim$flag==0)
        {code<-fix(x=fx,gen=gen,y=phe,kk=kk)
        tempcode<-code
        tempcode[,2:3] <- gen[,1:2]
        }
        envim$res1<-tempcode
        
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
        
        if(length(dd)==1){
          na<-matrix(name[dd],1,)
          xxxm<-matrix(xxx[,dd],,1)
          wow<-cbind(fx,xxxm)
          bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
          ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
          genna1<-matrix(gen[na,1],,1)
          genna2<-matrix(gen[na,2],,1)
          genna3<-as.matrix(ef)
          genna4<-as.matrix(lod[dd,])
          galaxy<-cbind(genna1,genna2,genna3,genna4)
        }else{
          na<-as.matrix(name[dd])
          wow<-cbind(fx,xxx[,dd])
          bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
          ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
          galaxy<-cbind(gen[na,1],gen[na,2],ef,lod[dd,])
        }
        
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
        va<-round(va,4)
        heredity<-100*heredity
        heredity<-matrix(heredity,,1)
        heredity<-round(heredity,4)
        galaxyy[which(abs(galaxyy)>1e-4)]<-round(galaxyy[which(abs(galaxyy)>1e-4)],4)
        galaxyy[which(abs(galaxyy)<1e-4)]<-as.numeric(sprintf("%.4e",galaxyy[which(abs(galaxyy)<1e-4)]))
        
        
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
        vee[1,1]<-round(ve,4)
        vee<-matrix(vee,,1)
        vpp<-matrix("",nrow(galaxyy_A),1)
        vpp[1,1]<-round(pv,4)
        vpp<-matrix(vpp,,1)
        traitid<-matrix(envim$ii,nrow(galaxyy_A),1)
        
        envim$galaxyy<-galaxyy
        envim$result<-cbind(traitid,galaxyy_A,left_marker,right_marker,va,heredity,vee,vpp)
        colnames(envim$result)<-c("TraitID","Chr","Position (cM)","Additive Effect","LOD","Left_Marker","Right_Marker","Var_Genet_(i)","r2 (%)","Var_Error",
                                  "Var_Phen (total)")
        
        
        if(is.null(envim$result)==TRUE){
          gmessage("No result meets the requirements in the last step!","Info",icon="info")
          return
        }else{
          tbdfe8<-gdfedit(envim$result,container=nb1,expand=TRUE,label="Result")
        }
        progress_bar$setFraction(100/100)
        progress_bar$setText("All done.")
      }
    }
    
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
        
      }else if(is.null(envim$mapRaw11)==FALSE){
        mapRaw1 <- NULL
        mapRaw11 <- NULL
        value_AA<-NULL
        value_Aa<-NULL
        value_aa<-NULL
        value_<-NULL
        rm(value_AA,envir = as.environment(envim))
        rm(value_Aa,envir = as.environment(envim))
        rm(value_aa,envir = as.environment(envim))
        rm(value_,envir = as.environment(envim))
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
      }else if(is.null(envim$result)==FALSE){
        result <- NULL
        chr_name <- NULL
        galaxyy <- NULL
        chrRaw_name <- NULL
        mapname <- NULL
        mx <- NULL
        ii <- NULL
        delerow<-NULL
        res1<-NULL
        rm(res1,envir = as.environment(envim))
        rm(delerow,envir = as.environment(envim))
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
        mapRaw1 <- NULL
        rm(mapRaw1,envir = as.environment(envim))
        rm(y_jun00,envir = as.environment(envim))
        rm(flag,envir = as.environment(envim))
        rm(flagRIL,envir = as.environment(envim))
        rm(flagCSSL,envir = as.environment(envim))
        rm(y_jun3,envir = as.environment(envim))
      } else if(is.null(envim$geoo)==FALSE){
        geoo<- NULL
        parm<-NULL
        rm(parm,envir = as.environment(envim))
        rm(geoo,envir = as.environment(envim))
      }else if(is.null(envim$poss)==FALSE){
        flag <- NULL
        flagRIL <- NULL
        flagCSSL <- NULL
        poss <- NULL
        mapRaw1 <- NULL
        rm(mapRaw1,envir = as.environment(envim))
        rm(flag,envir = as.environment(envim))
        rm(flagRIL,envir = as.environment(envim))
        rm(flagCSSL,envir = as.environment(envim))
        rm(poss,envir = as.environment(envim))
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
      gmessage("No result meets the requirements!","Info",icon="info")
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
      gmessage("Please input correct genotype dataset!","Warning",icon="warning")
    }
    if(is.null(envim$galaxyy)==TRUE){
      gmessage("No result meets the requirements!","Info",icon="info")
    }
    
    if((is.null(envim$mx)==FALSE)&&(is.null(envim$galaxyy)==FALSE)&&(is.null(envim$chr_name)==FALSE)){
      plotwin<-gwindow("Genome-wide composite interval mapping (GCIM) figure",visible=FALSE,width=960,height=500)
      gpw<-ggroup(container=plotwin,horizontal = FALSE,spacing = 10)
      ggpw<-ggraphics(container=gpw)
      addSpring(gpw)
      plotlyt<-glayout(container =gpw )
      plotfram<-gframe(container =  plotlyt)
      
      plotlabel<-glabel("Width (px):",container = plotlyt)
      plotedit<-gedit("1500",container = plotlyt)
      plotlabel1<-glabel("Height (px):",container = plotlyt)
      plotedit1<-gedit("600",container = plotlyt)
      plotlabel2<-glabel("Word resolution (1/72 inch, ppi):",container = plotlyt)
      plotedit2<-gedit("12",container = plotlyt)
      plotbt<-gbutton(" Save ",container =plotlyt )
      combo_box1 <- gcombobox(selected=3,c("blue","black","red","yellow","green","pink","purple","gray50","brown"),editable = TRUE,container = plotlyt)
      combo_box2 <- gcombobox(selected=8,c("blue","black","red","yellow","green","pink","purple","gray50","brown"),editable = TRUE,container = plotlyt)
      plotlabel3 <- glabel("LOD line color:", container = plotlyt)
      plotlabel4 <- glabel("-log10(P) curve color:", container = plotlyt)
      plotlabel5 <- glabel("Legend and tick marks:",container = plotlyt)
      plotedit5 <- gedit("1.0",container = plotlyt)
      plotlabel6 <- glabel("LOD line size:",container = plotlyt)
      plotedit6 <- gedit("1.0",container = plotlyt)
      plotlabel7 <- glabel("Size for -log10(P) curve:",container = plotlyt)
      plotedit7 <- gedit("0.5",container = plotlyt)
      plotlabel8 <- glabel("Margin space:",container = plotlyt)
      plotedit8 <- gedit("1.5",container = plotlyt)
      plotlabel9 <- glabel("Space between tick marks and axis:",container = plotlyt)
      plotedit9 <- gedit("1.0",container = plotlyt)
      plotlabel10 <- glabel("Times for max{-log10(P)}:",container = plotlyt)
      plotedit10 <- gedit("1.5",container = plotlyt)
      plotlabel11 <- glabel("Critical LOD score:",container = plotlyt)
      plotedit11 <- gedit("2.5",container = plotlyt)
      plotlabel12 <- glabel("Note: The parameter settings of preview figure is for general resolution, see manual for parameter settings of high resolution figure (*.png,*.tiff,*.jpeg).",container = plotlyt)
      plotlabel13 <- glabel("Figure resolution (ppi):",container = plotlyt)
      plotedit13 <- gedit("72",container = plotlyt)
      
      plotlyt[1,2] <- plotlabel5
      plotlyt[1,3] <- plotedit5
      plotlyt[1,5] <- plotlabel6
      plotlyt[1,6] <- plotedit6
      plotlyt[1,8] <- plotlabel7
      plotlyt[1,9] <- plotedit7
      plotlyt[2,2] <- plotlabel8
      plotlyt[2,3] <- plotedit8
      plotlyt[2,5] <- plotlabel9
      plotlyt[2,6] <- plotedit9
      plotlyt[2,8] <- plotlabel10
      plotlyt[2,9] <- plotedit10
      plotlyt[3,2] <- plotlabel
      plotlyt[3,3] <- plotedit
      plotlyt[3,5] <- plotlabel1
      plotlyt[3,6] <- plotedit1
      plotlyt[3,8] <- plotlabel11
      plotlyt[3,9] <- plotedit11
      plotlyt[4,2] <- plotlabel2
      plotlyt[4,3] <- plotedit2
      plotlyt[4,5] <- plotlabel13
      plotlyt[4,6] <- plotedit13 
      plotlyt[5,2] <- plotlabel3
      plotlyt[5,3] <- combo_box1
      plotlyt[5,5] <- plotlabel4
      plotlyt[5,6] <- combo_box2 
      plotlyt[5,8:9] <- plotbt
      plotlyt[6,2:9] <- plotlabel12
      visible(plotwin) <- TRUE
      
      gcimFunc <- function(legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,lodthred)
      {
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
        
        firFil <- envim$res1[,2:3]
        newposadd <- as.matrix(firFil[,2])
        for(i in 1:chr_num)
        {
          temp1 <- numeric()
          temp1 <- which(firFil[,1]==i)
          if(i>1)
          {
            newposadd[temp1] <- newposadd[temp1]+pos_acc[i-1]
          }
        }
        
        newres_pos <- envim$galaxyy[,2]
        res_sumpos <- pos_acc[envim$galaxyy[which(envim$galaxyy[,1]>1),1]-1] + envim$galaxyy[which(envim$galaxyy[,1]>1),2] 
        newres_pos[which(envim$galaxyy[,1]>1)] <- res_sumpos 
        pospic<-c(newres_pos)
        lodpic<-c(envim$galaxyy[,4])
        resdf <- data.frame(pospic,lodpic)
        
        resp<-as.matrix(envim$res1[,10])
        pmin<-min(resp[resp!=0])
        locsub<-which(resp==0)
        if(length(locsub)!=0){
          subvalue<-10^(1.1*log10(pmin))
          envim$res1[locsub,10]<-subvalue
        }else{
          envim$res1<-envim$res1
        }
        
        
        negloP <- -log10(as.matrix(envim$res1[,10]))
        
        par(mar=c(2*margin_space,2*margin_space,2*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))  
        plot(pospic,lodpic,type="h",col=color1,xlab="",ylab="Logarithm of odds (LOD)",cex.axis=legend_size,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadd)),ylim=c(0,max(lodpic)))
        abline(h=lodthred)
        par(new=TRUE)
        plot(newposadd,negloP,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
        axis(side=4,cex.axis=legend_size)
        mtext(expression('-log'[10]*'(P)'),side=4,line=3*axis_space,cex=legend_size)
        abline(v=pos_acc,lty=2,col="gray")
      }
      
      addHandlerChanged(ggpw, handler=function(h,...) {
        legend_size <- as.numeric(svalue(plotedit5))
        mainline_size <- as.numeric(svalue(plotedit6))
        backline_size <- as.numeric(svalue(plotedit7))
        margin_space <- as.numeric(svalue(plotedit8))
        axis_space <- as.numeric(svalue(plotedit9))
        logPCoff <- as.numeric(svalue(plotedit10))
        lodthred <- as.numeric(svalue(plotedit11))
        color1 <- svalue(combo_box1)
        color2 <- svalue(combo_box2)
        gcimFunc(legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,lodthred)
      })
      
      addhandlerclicked(plotbt, handler = function(h, ...) {
        widqqvalue <- as.numeric(svalue(plotedit))
        heightqqvalue <- as.numeric(svalue(plotedit1))
        pointsizeqqvalue <- as.numeric(svalue(plotedit2))
        resppi <- as.numeric(svalue(plotedit13))
        
        output <- gfile(text = "Save a file...", type = "save", 
                        filter = list(`All files` = list(patterns = c("*")), 
                                      `TIFF files` = list(patterns = c("*.tiff")),
                                      `PNG files` = list(patterns = c("*.png")),
                                      `JPEG files` = list(patterns = c("*.jpeg"))))
        if((length(grep(".png",output))==1)||(length(grep(".PNG",output))==1)){
          png(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
        }else if((length(grep(".tiff",output))==1)||(length(grep(".TIFF",output))==1)){
          tiff(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
        }else if((length(grep(".jpeg",output)==1))||(length(grep(".JPEG",output))==1)){
          jpeg(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)  
        }else{
          gmessage("Please input correct image format!")
        }
        
        legend_size <- as.numeric(svalue(plotedit5))
        mainline_size <- as.numeric(svalue(plotedit6))
        backline_size <- as.numeric(svalue(plotedit7))
        margin_space <- as.numeric(svalue(plotedit8))
        axis_space <- as.numeric(svalue(plotedit9))
        logPCoff <- as.numeric(svalue(plotedit10))
        color1 <- svalue(combo_box1)
        color2 <- svalue(combo_box2)
        lodthred <- as.numeric(svalue(plotedit11))
        gcimFunc(legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,lodthred)
        dev.off()
        
      })
    }
  }
  
  
  fuquit<-function(h,...){
    gconfirm("Yes or no ?",handler=function(h,...){dispose(window)})
  }
  
  fuhelp<-function(h,...){
    RShowDoc("Instruction.pdf",package="mrMLM")
    
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
  GCIMSub()
})

FASTSub<-function(){
  mrenv <- new.env()
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
  
  window<-gwindow(title="FAST multi-locus random-SNP-effect EMMA (FASTmrEMMA)",visible=TRUE,width=1290,height=730,expand=TRUE)
  importwin<-gwindow("Input Dataset",visible=FALSE,width=250,height=420)
  gimpwin<-ggroup(container=importwin,expand=FALSE)
  plotwin<-gwindow("Manhattan Plot",visible=FALSE,width=600,height=360)
  gpw<-ggroup(container=plotwin,horizontal = FALSE)
  ggpw<-ggraphics(container=gpw)
  plotwin1<-gwindow("Q-Q Plot",visible=FALSE,width=600,height=360)
  gpw1<-ggroup(container=plotwin1,horizontal=FALSE)
  ggpw1<-ggraphics(container=gpw1)
  choicekk<-gwindow("Choose Kinship option",visible=FALSE,width=350,height=150)
  gkk<-ggroup(container=choicekk,expand=FALSE)
  includeps<-gwindow("Include population structure (Q matrix)?",visible=FALSE,width=400,height=150)
  gps<-ggroup(container=includeps,expand=FALSE)
  parsetwin<-gwindow("Parameter Settings",visible=FALSE,width=300,height=270)
  gpar<-ggroup(container=parsetwin,expand=FALSE)
  choicesave<-gwindow("Save as ...",visible=FALSE,width=250,height=150)
  gcsave<-ggroup(container=choicesave,expand=FALSE)
  lyt<-glayout(container=window,spacing=13)
  
  importdata<-gbutton("Input Datasets",container=lyt)
  parset<-gbutton("Parameter Settings",container=lyt)
  manhattan<-gbutton("Manhattan Plot",container=lyt)
  qqplot<-gbutton("QQ Plot",container=lyt)
  
  savefile<-gbutton(" Save ",container=lyt)
  run<-gbutton("Run",container=lyt)
  clear<-gbutton(" Clear ",container = lyt)
  exit<-gbutton("Exit",container=lyt)
  helpfile<-gbutton("User Manual",container=lyt)
  
  lyt[1,1]<-importdata
  lyt[4,1]<-parset
  lyt[5,1]<-run
  lyt[6,1]<-clear
  lyt[7,1]<-savefile
  lyt[11,1]<-manhattan
  lyt[12,1]<-qqplot
  lyt[16,1]<-helpfile
  lyt[17,1]<-exit
  
  nb1<-gnotebook(tab.pos=3,closebuttons=TRUE,dontCloseThese=TRUE,container=lyt,expand=TRUE)
  size(nb1)<-c(680,540)
  tb<-gnewtable("     
                1. FASTmrEMMA is a R software package for multi-locus genome-wide association studies based on a fast  
                   multi-locus random-SNP-effect EMMA.
                
                2. Please cite: Wen Yang-Jun, Zhang Hanwen, Ni Yuan-Li, Huang Bo, Zhang Jin, Feng Jian-Ying, Wang Shi-Bo, 
                   Jim M. Dunwell, Zhang Yuan-Ming*, Wu Rongling*. Methodological implementation of mixed linear models 
                   in multi-locus genome-wide association studies. Briefings in Bioinformatics, DOI: 10.1093/bib/bbw145.
                
                3. Please cite: Zhang Yuan-Ming et al. Mapping quantitative trait loci using naturally occurring genetic
                   variance among commercial inbred line of maize (Zea mays L.). Genetics 2005, 169:2267-2275.                  
                
                4. The software package is developed by Bo Huang, Yuan-Li Ni, Wen-Long Ren, Yang-Jun Wen & Yuan-Ming Zhang.
                   ",multiple=TRUE,container=nb1,expand=TRUE,label="About the software")
  
  font(tb)<-c(size="x-large")
  lyt[1:20,2,expand=TRUE]<-nb1
  progress_bar <- gtkProgressBar()
  lyt[21,2,expand=TRUE]<-progress_bar
  
  addHandlerClicked(importdata,handler=function(h,...){
    if(isExtant(importwin)==FALSE)
    {
      importwin<-gwindow("Input Dataset",visible=FALSE,width=250,height=500)
      gimpwin<-ggroup(container=importwin,expand=FALSE)
    }
    lytimp<-glayout(container=gimpwin,spacing=13)
    impchoose<-glabel("1. Choose dataset format",container=lytimp)
    impfile1<-glabel("2. Input Genotypic and Phenotypic files",container=lytimp)
    impprepare<-glabel("4. Sort & Transform for datasets",container=lytimp)
    impfile2<-glabel("5. Input Kinship (K) and Population-structure (Q matrix) files",container=lytimp)
    radioimp<-gradio(c("FASTmrEMMA numeric format","mrMLM character format","Hapmap (TASSEL) format"),selected=3,horizontal=FALSE,container=lytimp)
    genotype<-gbutton("Genotype",container=lytimp)
    phenotype<-gbutton("Phenotype",container=lytimp)
    kinship<-gbutton("Kinship",container=lytimp)
    population<-gbutton("Population Structure",container=lytimp)
    preimp<-gbutton("Do",container=lytimp)
    mllabel<-glabel("3. Please select likelihood Function",container =lytimp )
    mlradio<-gradio(c("Restricted Likelihood","Likelihood"),horizontal = FALSE,container = lytimp,selected = 1)
    
    lytimp[1,2:5]<-impchoose
    lytimp[2:4,2:5]<-radioimp
    lytimp[5,2:5]<-impfile1
    lytimp[6,2:4]<-genotype
    lytimp[7,2:4]<-phenotype
    lytimp[8,2:5]<-mllabel
    lytimp[9:10,2:4]<-mlradio
    lytimp[11,2:5]<-impprepare
    lytimp[12,2:4]<-preimp
    lytimp[13,2:5]<-impfile2
    lytimp[14,2:4]<-kinship
    lytimp[15,2:4]<-population
    visible(importwin)<-TRUE
    
    addHandlerClicked(genotype,handler=function(h,...){
      
      input1<-gfile(text="Select a file...",type="open",
                    filter=list("All files"=list(patterns=c("*")),
                                "CSV files"=list(patterns=c("*.csv"))))
      
      if(is.na(input1))
      {
        gmessage("Please input correct genotype datasets!","Warning",icon="warning")
      }else{
        mrenv$genRaw<-read.csv(input1,header=FALSE)
        mrenv$genRaw<-as.matrix(mrenv$genRaw)
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
        gmessage("Please input correct phenotype datasets!","Warning",icon="warning")
      }else{
        pheRaw1<-as.matrix(read.csv(input2,header=FALSE)) 
        pheRaw2<-pheRaw1[-1,]
        pheRaw3<-as.data.frame(pheRaw2,stringsAsFactors=FALSE)
        pheRaw4<-as.matrix(pheRaw3[is.na(pheRaw3[,2])==F,])
        pheRawthem<-matrix(c(pheRaw1[1,1]," "),1,)
        pheRaw<-rbind(pheRawthem,pheRaw4)
        row.names(pheRaw)<-NULL
        mrenv$pheRaw<-as.matrix(pheRaw)
        showpheRaw<-pheRaw1[-1,]
        colnames(showpheRaw)<-c(pheRaw1[1,1],"   ")
        showpheRaw<-as.data.frame(showpheRaw)
        tbdfe2<-gdfedit(showpheRaw,container=nb1,expand=TRUE,label="Raw_Phenotype")
        
      }
    })
    
    addHandlerClicked(preimp,handler=function(h,...){
      if((svalue(radioimp)=="FASTmrEMMA numeric format")&&(svalue(mlradio)=="Restricted Likelihood")){
        mrenv$flagREMLE<-1
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
      }else if((svalue(radioimp)=="FASTmrEMMA numeric format")&&(svalue(mlradio)=="Likelihood")){
        mrenv$inputform<-1
        mrenv$flagREMLE<-0
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
      }else if((svalue(radioimp)=="mrMLM character format")&&(svalue(mlradio)=="Restricted Likelihood")){
        mrenv$inputform<-2
        mrenv$flagREMLE<-1
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
          computeGen[i,notRef] <- as.numeric(0)
          computeGen[i,notATCG] <- as.numeric(0.5)
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
      }else if((svalue(radioimp)=="mrMLM character format")&&(svalue(mlradio)=="Likelihood")){
        mrenv$inputform<-2
        mrenv$flagREMLE<-0
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
          computeGen[i,notRef] <- as.numeric(0)
          computeGen[i,notATCG] <- as.numeric(0.5)
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
      }
      else if((svalue(radioimp)=="Hapmap (TASSEL) format")&&(svalue(mlradio)=="Restricted Likelihood")){
        mrenv$inputform<-3
        mrenv$flagREMLE<-1
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
          computeGen[i,notRef]<-as.numeric(0)
          computeGen[i,notATCG]<-as.numeric(0.5)
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
      }else if((svalue(radioimp)=="Hapmap (TASSEL) format")&&(svalue(mlradio)=="Likelihood")){
        mrenv$inputform<-3
        mrenv$flagREMLE<-0
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
          computeGen[i,notRef]<-as.numeric(0)
          computeGen[i,notATCG]<-as.numeric(0.5)
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
        choicekk<-gwindow("Choose Kinship option",visible=FALSE,width=350,height=150)
        gkk<-ggroup(container=choicekk,expand=FALSE)
      }
      lytkk<-glayout(container=gkk,spacing=13)
      mrenv$okkk<-gbutton("     OK    ",container=lytkk)
      mrenv$cancelkk<-gbutton(" Cancel ",container=lytkk)
      mrenv$radiokk<-gradio(c("Input the Kinship matrix file","Calculate the Kinship matrix by this software"),selected=1,horizontal=FALSE,container=lytkk)
      kklabel<-glabel("Note: Please select the 1st option if no. of markers is more than 50,000.",container = lytkk)
      font(kklabel)<-("bold")
      lytkk[2:3,2:5]<-mrenv$radiokk
      lytkk[5,2]<-mrenv$okkk
      lytkk[5,6]<-mrenv$cancelkk
      lytkk[6,2:7]<-kklabel
      visible(choicekk)<-TRUE
      addHandlerClicked(mrenv$okkk,handler=function(h,...){
        if(svalue(mrenv$radiokk)=="Input the Kinship matrix file"){
          input3<-gfile(text="Select a file...",type="open",
                        filter=list("All files"=list(patterns=c("*")),
                                    "CSV files"=list(patterns=c("*.csv"))))
          if(is.na(input3))
          {
            gmessage("Please input correct kinship dataset!","Warning",icon="warning")
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
          snp5<- mrenv$gen
          
          
          if(exists("snp5")==FALSE)
          {
            gmessage("Please input correct genotype dataset!","Warning",icon="warning")
          }else{
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
            snp8<-snp5[,3:dim(snp5)[2]]
            kk<-emma.kinship(snp8)
            mrenv$kk<-as.matrix(kk)
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
        includeps<-gwindow("Include population structure?",visible=FALSE,width=400,height=150)
        gps<-ggroup(container=includeps,expand=FALSE)
      }
      lytps<-glayout(container=gps,spacing=13)
      okps<-gbutton("     OK    ",container=lytps)
      cancelps<-gbutton(" Cancel ",container=lytps)
      radiops<-gradio(c("Not included in the model","Included"),selected=1,horizontal=FALSE,container=lytps)
      lytps[2:3,2:5]<-radiops
      lytps[5,2]<-okps
      lytps[5,5]<-cancelps
      visible(includeps)<-TRUE
      addHandlerClicked(okps,handler=function(h,...){
        if(svalue(radiops)=="Included"){
          mrenv$flagps<-0
          input4<-gfile(text="Select a file...",type="open",
                        filter=list("All files"=list(patterns=c("*")),
                                    "CSV files"=list(patterns=c("*.csv"))))
          if(is.na(input4))
          {
            gmessage("Please input correct population dataset!","Warning",icon="warning")
          }else{
            mrenv$psmatrixRaw<-as.matrix(read.csv(input4,header=FALSE))
            nnpprow<-dim(mrenv$psmatrixRaw)[1]
            nnppcol<-dim(mrenv$psmatrixRaw)[2]
            mrenv$psmatrixRaw[1,2:nnppcol]<-"  "
            psmatrixPre<-mrenv$psmatrixRaw[3:nnpprow,]
            namePop<-as.matrix(psmatrixPre[,1])
            sameGenPop<-intersect(mrenv$sameName,namePop)
            locPop<-match(sameGenPop,namePop)
            ##revised
            filtername<-as.vector(mrenv$psmatrixRaw[2,2:nnppcol])
            selectpsmatrix<-matrix(as.numeric(psmatrixPre[locPop,-1]),nrow = nrow(psmatrixPre))
            psum<-apply(selectpsmatrix,1,sum)
            psum<-round(psum)
            sumps<-sum(psum)
            m<-dim(selectpsmatrix)[1]
            if(sumps==m){
              filterps<-gwindow("Filter",visible=FALSE,width=300,height=120)
              filergp<-ggroup(container=filterps,expand=FALSE)
              filterlyt<-glayout(container = filergp)
              filterlabel<-glabel("Please choose the column tha should be deleted in population structure Q matrix",container = filterlyt)
              filtercombo<-gcombobox(filtername,editable = TRUE,container =filterlyt)
              filterlabel1<-glabel("Filter:",container = filterlyt)
              filterok<-gbutton("OK",container = filterlyt)
              filtercancel<-gbutton("Cancel",container = filterlyt)
              filterlyt[2,2:5]<-filterlabel
              filterlyt[4,2]<-filterlabel1
              filterlyt[4,3]<-filtercombo
              filterlyt[5,2]<-filterok
              filterlyt[5,5]<-filtercancel
              visible(filterps)<-TRUE
              addhandlerclicked(filterok,handler = function(h,...){
                combovalue<-svalue(filtercombo)
                coldelet<-unlist(str_extract_all(combovalue,"[0-9]+"))
                coldelet<-as.numeric(coldelet)
                mrenv$psmatrix<-as.matrix(selectpsmatrix[,-coldelet])
                mrenv$psmatrixRaw<-as.matrix(mrenv$psmatrixRaw[,-(coldelet+1)])
                tbdfe6<-gdfedit(mrenv$psmatrixRaw,container=nb1,expand=TRUE,label="Population Structure")
                dispose(filterps)
                dispose(includeps)
                dispose(importwin)
              })
              addhandlerclicked(filtercancel,handler = function(h,...){
                dispose(filterps)
              })
              
            }else{
              mrenv$psmatrix<-selectpsmatrix
              tbdfe6<-gdfedit(mrenv$psmatrixRaw,container=nb1,expand=TRUE,label="Population Structure")
              
              dispose(includeps)
              dispose(importwin)
            }
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
      parsetwin<-gwindow("Parameter Setting",visible=FALSE,width=300,height=270)
      gpar<-ggroup(container=parsetwin,expand=FALSE)
    }
    lytpar<-glayout(container=gpar,spacing=13)
    mrenv$pvallabel<-glabel("1. Critical P-value in the 1st step of FASTmrEMMA:",container=lytpar)
    mrenv$pvaledit<-gedit("0.005",width=20,coerce.with=as.numeric,container=lytpar)
    
    mrenv$mlodlabel<-glabel("2. Critical LOD score in the last step of FASTmrEMMA:",container=lytpar)
    mrenv$mlodedit<-gedit("3",width=20,coerce.with=as.numeric,container=lytpar)
    mrenv$okpar<-gbutton("     OK    ",container=lytpar)
    mrenv$cancelpar<-gbutton(" Cancel ",container=lytpar)
    lytpar[1,1:5]<-mrenv$pvallabel
    lytpar[2,1:5]<-mrenv$pvaledit
    
    lytpar[3,1:5]<-mrenv$mlodlabel
    lytpar[4,1:5]<-mrenv$mlodedit
    lytpar[6,1]<-mrenv$okpar
    lytpar[6,4]<-mrenv$cancelpar
    visible(parsetwin)<-TRUE
    addHandlerClicked(mrenv$okpar,handler=function(h,...){
      mrenv$svpal<-svalue(mrenv$pvaledit)
      
      mrenv$svmlod<-svalue(mrenv$mlodedit)
      if((mrenv$svpal<0)||(mrenv$svpal>1))
      {
        gmessage("Please input critical P-value between 0 and 1!","Warning",icon="warning")
      }
      
      if(mrenv$svmlod<0)
      {
        gmessage("Please input critical LOD score: >0!","Warning",icon="warning")
      }
      if((mrenv$svpal>0)&&(mrenv$svpal<1)&&(mrenv$svmlod>=0))
      {
        dispose(parsetwin)
      }
    })
    
    addHandlerClicked(mrenv$cancelpar,handler=function(h,...){
      mrenv$svpal<-svalue(mrenv$pvaledit)
      mrenv$svmlod<-svalue(mrenv$mlodedit)
      dispose(parsetwin)
    })
  })
  
  addHandlerClicked(exit,handler=function(h,...){
    gconfirm("Yes or no?",handler=function(h,...){dispose(window)})
  })
  
  addHandlerClicked(helpfile,handler=function(h,...){
    RShowDoc("Instruction.pdf",package="mrMLM")  
  })
  
  addHandlerClicked(run,handler=function(h,...){
    gen<-mrenv$gen
    phe<-mrenv$phe
    kk<-mrenv$kk
    flagps<-mrenv$flagps
    psmatrix<-mrenv$psmatrix
    if(is.null(gen)==TRUE)
    {
      gmessage("Please input correct genotype dataset!","Warning",icon="warning")
    }
    if(is.null(phe)==TRUE)
    {
      gmessage("Please input correct phenotype dataset!","Warning",icon="warning")
    }
    if(is.null(kk)==TRUE)
    {
      gmessage("Please input correct kinship dataset!","Warning",icon="warning")
    }
    if((is.null(gen)==FALSE)&&(is.null(phe)==FALSE)&&(ncol(gen)!=(nrow(phe)+2)))
    {
      gmessage("Sample size between genotypic and phenotypic datasets are not same!","Error",icon="error")
    }
    
    if((is.null(gen)==FALSE)&&(is.null(phe)==FALSE)&&(is.null(kk)==FALSE)&&((ncol(gen)==(nrow(phe)+2))))
    {
      
      progress_bar$setText ( "Please be patient ..." )
      progress_bar$setFraction(2/100)
      
      multinormal<-function(y,mean,sigma)
      {
        pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
        return (pdf_value)
      }
      #LOD value test
      likelihood<-function(xxn,xxx,yn,bbo)
      {
        nq<-ncol(xxx)
        ns<-nrow(yn)
        at1<-0
        ww1<-as.matrix(which(abs(bbo)>1e-5))
        at1<-dim(ww1)[1]
        lod<-matrix(rep(0,nq),nq,1)
        if(at1>0.5)
          ad<-cbind(xxn,xxx[,ww1])
        else
          ad<-xxn
        
        bb<-ginv(crossprod(ad,ad))%*%crossprod(ad,yn)
        
        vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns)
        ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
        
        sub<-1:ncol(ad);
        if(at1>0.5)
        {
          for(i in 1:at1)
          {
            #ij<-which(sub!=sub[i+1])
            ij<-which(sub!=sub[i+ncol(xxn)])
            ad1<-ad[,ij]
            
            bb1<-ginv(crossprod(ad1,ad1))%*%crossprod(ad1,yn) 
            vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns)
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
        b<-solve(crossprod(x,x))%*%(crossprod(x,y))
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
        iter<-0;err<-1000;iter_max<-5000;err_max<-1e-10
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
            b<-solve(xtv%*%x)%*%(xtv%*%y)
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
        #return (u)
        return (list(u=u,sigma2=sigma2))
      }
      
      #####k kinship start################
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
      
      #likelihood
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
      #restricted likelihood
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
      }####have revised 20140822
      
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
      }####have revised 20140823
      
      emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
        #t <- length(xi.1)#have no relationship with t#have revised 20140830
        delta <- exp(logdelta)
        etasq <- etas.1*etas.1
        ldelta <- delta*lambda+1
        return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(xi.1/(delta*xi.1+1))) )
      }####have revised 20140823
      
      emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
        nq <- length(etas)
        delta <-  exp(logdelta)
        return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(delta*lambda+1))))-sum(log(delta*lambda+1))) )
      }####have revised 20140823
      
      emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
        tq <- length(etas.1)
        nq <- n - t + tq#(n-t):the number of eigvalue 1;tq:the number of components of etas.1
        delta <-  exp(logdelta)
        return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
      }####have revised 20140823
      
      emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
        nq <- length(etas)
        delta <- exp(logdelta)
        etasq <- etas*etas
        ldelta <- delta*lambda+1
        return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(lambda/ldelta)) )
      }####have revised 20140823
      
      emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
        t <- t1
        tq <- length(etas.1)
        nq <- n - t + tq
        delta <- exp(logdelta)
        etasq <- etas.1*etas.1
        ldelta <- delta*lambda+1
        return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
      }####have revised 20140823
      
      emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                           esp=1e-10, eig.L = NULL, eig.R = NULL)
      {
        n <- length(y)
        t <- nrow(K)
        q <- ncol(X)
        #  stopifnot(nrow(K) == t)
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
          #    optdelta <- exp(optlogdelta)
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
          #    optdelta <- exp(optlogdelta)
        }
        maxdelta <- exp(optlogdelta[which.max(optLL)])
        #handler of grids with NaN log
        optLL=replaceNaN(optLL)  #20160728
        maxLL <- max(optLL)
        if ( is.null(Z) ) {
          maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
        }
        else {
          maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
        }
        maxvg <- maxve*maxdelta
        
        return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
        #####have revised 20140823
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
          #    optdelta <- exp(optlogdelta)
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
          ######have revised 20140823
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
          #    optdelta <- exp(optlogdelta)
        }  
        maxdelta <- exp(optlogdelta[which.max(optLL)])
        #handler of grids with NaN log
        optLL=replaceNaN(optLL)  #20160728
        maxLL <- max(optLL)
        if ( is.null(Z) ) {
          maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
        }
        else {
          maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
        }
        maxvg <- maxve*maxdelta
        ####have revised 20140823 
        return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
        ####have revised 20140823
      }
      
      #likelihood: 
      emma.eigen.L.c <- function(Z,K,complete=TRUE) {
        if ( is.null(Z) ) {
          return(emma.eigen.L.wo.Z.c(K))
        }#don't run
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
        }#don't run
        else {
          return(emma.eigen.R.w.Z.c(Z,K,X,complete))
        }
      }
      #restricted likelihood
      emma.eigen.R.wo.Z.c <- function(K, X) {
        #n <- nrow(X)
        #q <- ncol(X)
        
        #have revised 20140916
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
      }#don't run
      
      #restricted likelihood
      emma.eigen.R.w.Z.c <- function(Z, K, X, complete = TRUE) {
        if ( complete == FALSE ) {#don't run
          vids <-  colSums(Z) > 0
          Z <- Z[,vids]
          K <- K[vids,vids]
        }#don't run
        if(!is.matrix(Z)) n<-length(Z) 
        t <-1 #have revised 20140916
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
                                 complete=TRUE)[,c(1:t,(t+q+1):n)]))   # have revised 20140830
        
      }
      
      emma.delta.REML.LL.w.Z.c <- function(logdelta, lambda, etas.1, n, q, etas.2.sq ) {#t change to q#have revised 20140830
        nq<-n-q #have revised 20140830
        delta <-  exp(logdelta)
        return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
      }####have revised 20140823
      
      emma.delta.REML.dLL.w.Z.c <- function(logdelta, lambda, etas.1, n, q1, etas.2.sq ) {#t1 change to q1 #have revised 20140830
        
        q<-q1
        nq<-n-q
        delta <- exp(logdelta)
        etasq <- etas.1*etas.1
        ldelta <- delta*lambda+1
        return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
      }
      
      emma.MLE.c <- function(y, X, K=1, Z=NULL, ngrids=100, llim=-10, ulim=10,
                             esp=1e-10, eig.L = NULL, eig.R = NULL)#Z=X_c,K=1,X=W_c
      {
        if(is.matrix(y)){
          n<-nrow(y)
        }
        else{
          n <- length(y)
        }
        
        t <- 1#,K=1,t=1#have revised 20140915
        if ( is.matrix(X) ){
          q <- ncol(X)
          stopifnot(nrow(X)==n)
        }
        else{#X:n*1 vector
          q <- 1
          stopifnot(length(X)==n)
        }
        
        stopifnot(K == 1)#have revised 20140915
        
        if ( det(crossprod(X,X)) == 0 ) {
          warning("X is singular")
          return (list(ML=0,delta=0,ve=0,vg=0))
        }
        
        if ( is.null(Z) ) {#don't run
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
          #    optdelta <- exp(optlogdelta)
        }#don't run
        else {
          if ( is.null(eig.L) ) {
            eig.L <- emma.eigen.L.w.Z.c(Z,K)##have revised 20140830
          }
          if ( is.null(eig.R) ) {
            eig.R <- emma.eigen.R.w.Z.c(Z,K,X)#have revised 20140830
          }
          etas <- crossprod(eig.R$vectors,y)
          etas.1<-etas[1:t]
          etas.2<-etas[(t+1):(n-q)]#have revised 20140830
          etas.2.sq <- sum(etas.2*etas.2)
          logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
          m <- length(logdelta)
          delta <- exp(logdelta)
          Lambdas.1<-matrix(eig.R$values,t,m)
          Lambdas <- Lambdas.1 * matrix(delta,t,m,byrow=TRUE) + 1#have revised 20140830
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
          #    optdelta <- exp(optlogdelta)
        }
        
        maxdelta <- exp(optlogdelta[which.max(optLL)])
        optLL=replaceNaN(optLL)  #20160728
        maxLL <- max(optLL)
        if ( is.null(Z) ) {#don't run
          maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
        }#don't run
        else {
          maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
        }
        maxvg <- maxve*maxdelta
        return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,U_R=eig.R$vectors,etas.1=etas.1,etas=etas,lambda=eig.R$values))
      }
      
      emma.REMLE.c <- function(y, X, K=1, Z=NULL, ngrids=100, llim=-10, ulim=10,
                               esp=1e-10, eig.L = NULL, eig.R = NULL) #K=1,X=W_c,Z=X_c
      {
        if(is.matrix(y)){
          n<-nrow(y)
        }
        else{
          n <- length(y)
        }
        
        t <- 1#,K=1,t=1#have revised 20140915
        if ( is.matrix(X) ){
          q <- ncol(X)
          #stopifnot(nrow(X)==n)
        }
        else{#X:n*1 vector
          q <- 1
          #stopifnot(length(X)==n)
        }
        
        stopifnot(K == 1)#have revised 20140915
        stopifnot(nrow(X) == n)
        
        if ( det(crossprod(X,X)) == 0 ) {
          warning("X is singular")
          return (list(REML=0,delta=0,ve=0,vg=0))
        }
        
        if ( is.null(Z) ) {#don't run
          if ( is.null(eig.R) ) {
            eig.R <- emma.eigen.R.wo.Z.c(K,X)#have revised 20140830
          }
          etas <- crossprod(eig.R$vectors,y)
          logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
          m <- length(logdelta)
          delta <- exp(logdelta)
          Lambdas.1<-matrix(eig.R$values,n-q,m)
          Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
          ######have revised 20140823
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
          #    
        }#don't run
        else {
          if ( is.null(eig.R) ) {
            eig.R <- emma.eigen.R.w.Z.c(Z,K,X)#have revised 20140830
          }
          etas <- crossprod(eig.R$vectors,y)
          etas.1 <- etas[1:t]
          etas.2 <- etas[(t+1):(n-q)]#have revised 20140830
          etas.2.sq <- sum(etas.2*etas.2)
          logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
          m <- length(logdelta)
          delta <- exp(logdelta)
          Lambdas.1 <- matrix(eig.R$values,t,m) 
          Lambdas <- Lambdas.1 * matrix(delta,t,m,byrow=TRUE) + 1
          ######have revised 20140823
          Etasq <- matrix(etas.1*etas.1,t,m)#have revised 20140830
          
          dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
          optlogdelta <- vector(length=0)
          optLL <- vector(length=0)
          if ( dLL[1] < esp ) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(llim,eig.R$values,etas.1,n,q,etas.2.sq))# t change to q#have revised 20140830
          }
          if ( dLL[m-1] > 0-esp ) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(ulim,eig.R$values,etas.1,n,q,etas.2.sq))# t change to q#have revised 20140830
          }
          
          for( i in 1:(m-1) )
          {
            if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
            {
              r <- uniroot(emma.delta.REML.dLL.w.Z.c, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, q1=q, etas.2.sq = etas.2.sq )#t1=t change to q1=q#have revised 20140830
              optlogdelta <- append(optlogdelta, r$root)
              optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(r$root,eig.R$values, etas.1, n, q, etas.2.sq ))# t change to q#have revised 20140830
            }
          }
          
        }  
        maxdelta <- exp(optlogdelta[which.max(optLL)])
        optLL=replaceNaN(optLL)  #20160728
        maxLL <- max(optLL)
        if ( is.null(Z) ) {#don't run
          maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
        }# don't run
        else {
          maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
        }
        maxvg <- maxve*maxdelta
        return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,U_R=eig.R$vectors,etas.1=etas.1,etas=etas,lambda=eig.R$values))
      }
      
      emma.maineffects.B<-function(Z=NULL,K,deltahat.g,complete=TRUE){
        if( is.null(Z) ){
          return(emma.maineffects.B.Zo(K,deltahat.g))
        }
        else{
          return(emma.maineffects.B.Z(Z,K,deltahat.g,complete))
        }
      }
      
      #####
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
        # Y_c <- C%*%y
        #  W_c<-C%*%X
        M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
        etas<-crossprod(M_c,Y_c)
        #sum(etas*etas)
        
        LL <- 0.5*n*(log(n/(2*pi))-1-log(sum(etas*etas)))
        return(list(ML=LL))
      }
      
      emma.ML.LRT.c.noalpha <- function(ys, xs, K=1, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
        #K=1,Z=C:n*n,X0=W:n*c(or,1_n*1),ys:n*1
        stopifnot(K == 1)
        ys <- Z%*%ys   #Y_c=C*Y,Z=C,ys
        xs <- Z%*%xs
        X0 <- Z%*%X0
        
        ys<-as.matrix(ys)
        xs<-as.matrix(xs)
        X0<-as.matrix(X0)
        n <- nrow(ys)
        #n <- ncol(ys)
        m <- nrow(xs)
        t <- ncol(xs)
        q0 <- ncol(X0)
        #q1 <- q0 + 1
        MLE0<-emma.MLE0.c(ys,X0)
        #MLE0$ML
        
        ML1s <- vector(length=t)
        ML0s <- vector(length=t)#all of ML0s are the same,equal to MLE0<-emma.MLE0.c(ys,X0)
        vgs <- vector(length=t)
        ves <- vector(length=t)
        deltas<-vector(length=t)
        bhats<-vector(length=t)
        var.bhats.ratio <- vector(length=t)
        d <- vector(length=t)
        #have revised 20141201
        
        stats <- vector(length=t)
        ps <- vector(length=t)
        
        for (i in 1:t){
          progress_bar$setFraction((2+(88/t)*i)/100)
          vids <- !is.na(xs[,i])
          #nv <- sum(vids)
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
            d[i] <- (1-var.bhats.ratio[i])*((1-var.bhats.ratio[i])>=0) #to record me=sum(d)
            stats[i]<- 2*(MLE1$ML-MLE0$ML)
            ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2#20160619
          }else{
            ps[i]<-1
          }
        }
        return(list(ps=ps,bhats=bhats,deltas=deltas,d=d,ML1s=ML1s,ML0s=ML0s,stats=stats,vgs=vgs,ves=ves))
      } 
      
      emma.REMLE0.c <- function(Y_c,W_c){
        n <- length(Y_c)
        #stopifnot(nrow(C)==n)
        #stopifnot(ncol(C)==n)
        stopifnot(nrow(W_c)==n)
        # Y_c <- C%*%y
        #  W_c<-C%*%X
        M_c <-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
        eig <-eigen(M_c)
        t <-qr(W_c)$rank
        v <-n-t
        U_R <-eig$vector[,1:v]
        etas<-crossprod(U_R,Y_c)
        #sum(etas*etas)
        LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
        return(list(REML=LL))
      }
      emma.REML.LRT.c.noalpha <- function(ys, xs, K=1, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
        #K=1,Z=C:n*n,X0=W:n*c(or,1_n*1),ys:n*1
        stopifnot(K == 1)
        ys <- Z%*%ys   #Y_c=C*Y,Z=C,ys
        xs <- Z%*%xs
        X0 <- Z%*%X0
        
        ys<-as.matrix(ys)
        xs<-as.matrix(xs)
        X0<-as.matrix(X0)
        n <- nrow(ys)
        #n <- ncol(ys)
        m <- nrow(xs)
        t <- ncol(xs)
        q0 <- ncol(X0)
        #q1 <- q0 + 1
        MLE0<-emma.REMLE0.c(ys,X0)
        #MLE0$REML
        
        ML1s <- vector(length=t)
        ML0s <- vector(length=t)#all of ML0s are the same,equal to MLE0<-emma.MLE0.c(ys,X0)
        vgs <- vector(length=t)
        ves <- vector(length=t)
        deltas <- vector(length=t)
        bhats<-vector(length=t)
        var.bhats.ratio <- vector(length=t)
        d <- vector(length=t)
        stats <- vector(length=t)
        ps <- vector(length=t)
        
        for (i in 1:t){
          progress_bar$setFraction((2+(88/t)*i)/100)
          vids <- !is.na(xs[,i])
          #nv <- sum(vids)
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
            d[i] <- (1-var.bhats.ratio[i])*((1-var.bhats.ratio[i])>=0) #to record me=sum(d)
            stats[i]<- 2*(MLE1$REML-MLE0$REML)
            ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2
          }else{
            ps[i]<-1 
          }
        }
        return(list(ps=ps,bhats=bhats,deltas=deltas,d=d,ML1s=ML1s,ML0s=ML0s,stats=stats,vbs=vgs,ves=ves))
        
      }
      
      replaceNaN<-  function(LL) {
        index=(LL=="NaN")
        if(length(index)>0) theMin=min(LL[!index])
        if(length(index)<1) theMin="NaN"
        LL[index]=theMin
        return(LL)    
      }
      
      yraw<-matrix(phe[,1],,1)
      snp5<-gen
      xraw<-snp5
      xnames<-xraw[,1:2]
      snp1<-xraw[,3:dim(xraw)[2]]
      snp3<-matrix(snp1,nrow=dim(snp1)[1])
      mydata<-t(snp3)
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
      
      maf.fun<-function(snp){
        m<-length(snp)
        id.0<-length(which(snp==0))
        id.0.5<-length(which(snp==0.5))
        id.1<-length(which(snp==1))
        
        maf.1<-id.1/m
        maf.0.5<-id.0.5/m
        maf.0<-id.0/m
        maf.min<-min(id.1,id.0)/m
        
        return(list(maf.1,maf.0,maf.0.5,maf.min))
      }
      
      pve.fun<-function(beta,maf){
        
        pve<-(maf$p1-maf$p1^2+0.25*maf$p3-0.25*maf$p3^2-maf$p1*maf$p3)*beta^2
        
        return(pve)
      }
      
      if(mrenv$flagREMLE==1){
        
        remle1<-emma.REMLE(Y, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
        
      }
      if(mrenv$flagREMLE==0){
        remle1<-emma.MLE(Y, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
        
      }
      
      remle1.deltahat.g<-remle1$delta
      remle1.B1<-emma.maineffects.B(Z=NULL,K,remle1.deltahat.g)
      C2<-remle1.B1$mC
      
      if(mrenv$flagREMLE==1){
        REML.LRT.c2<-emma.REML.LRT.c.noalpha(ys=Y, xs=mydata, K=1, Z=C2, X0=W, ngrids=100, llim=-10, ulim=10, esp=1e-10) 
        
      }
      if(mrenv$flagREMLE==0){
        REML.LRT.c2<-emma.ML.LRT.c.noalpha(ys=Y, xs=mydata, K=1, Z=C2, X0=W, ngrids=100, llim=-10, ulim=10, esp=1e-10) 
        
      }
      REML.LRT.c2.new<-data.frame(REML.LRT.c2)
      mafall<-apply(snp1,1,maf.fun)
      mafall1<-unlist(mafall)
      mafall2<-matrix(mafall1,nrow = 4)
      mafall3<-t(mafall2)
      mafall4<-data.frame(mafall3)
      names(mafall4)<-c("p1","p2","p3","maf")
      MAF<-mafall4$maf
      pve.allr2<-(pve.fun(REML.LRT.c2$bhats,mafall4)/var(Y))*100
      mrenv$parms<-data.frame(chr.locus=xnames,REML.LRT.c2.new,MAF,pve.allr2)
      names(mrenv$parms)<-NULL
      mrenv$parms<-as.matrix(mrenv$parms)
      # colnames(mrenv$parms)<-c("Chromosome","Position","P_value","SNP effect","deltas","d","ML1s","ML0s","ML1s/ML0s","Sigma2_k","Sigma2","MAF","r2")
      ress<-mrenv$parms[,1:3]
      ress1<-ress[ress[,3]!=1,]
      resp<-as.matrix(ress1[,3])
      pmin<-min(resp[resp!=0])
      locsub<-which(resp==0)
      if(length(locsub)!=0){
        subvalue<-10^(1.1*log10(pmin))
        ress1[locsub,3]<-subvalue
        mrenv$ress1<-ress1
      }else{
        mrenv$ress1<-ress1
      }
      
      
      rowsnp <- dim(ress1)[1]
      mrenv$snpname <- numeric()
      snpname<-as.matrix(paste("rs",c(1:rowsnp),sep = ""))
      mrenv$snpname<-snpname
      pe<-as.numeric(mrenv$parms[,6])
      pee<-sum(pe)
      newp<-0.05/pee
      mrenv$mannewp<-newp
      mrenv$manstandchoice<-1
      ps<-as.matrix(mrenv$parms[,3])
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
      mrenv$vif<-GIF.mean
      VIF<-matrix("",dim(mrenv$parms)[1],1)
      VIF[1]<-mrenv$vif
      
      parmeter<-cbind(mrenv$parms[,1:5],mrenv$parms[,10:13])
      parmeter[which(abs(parmeter)>1e-4)]<-round(parmeter[which(abs(parmeter)>1e-4)],4)
      parmeter[which( abs(parmeter)<1e-4)]<-as.numeric(sprintf("%.4e", parmeter[which( abs(parmeter)<1e-4)]))
      parmeter[which(abs(parmeter[,5])>1e-4),5]<-round(parmeter[which(abs(parmeter[,5])>1e-4),5],4)
      parmeter[which(abs(parmeter[,5])<1e-4),5]<-as.numeric(sprintf("%.4e",parmeter[which(abs(parmeter[,5])<1e-4),5]))
      
      if(mrenv$inputform==1){
        #output result1 using mrMLM numeric format
        meadd<-matrix(" ",nrow(mrenv$parms),1)
        if(newp>1e-4){
          newp<-round(newp,4)
        }else{
          newp<-sprintf("%.4e",newp)
        }
        meadd[1,1]<-newp
        mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],parmeter,mrenv$genRaw[-1,4],meadd)
        colnames(mrenv$parmsShow)<-c("Marker","Chromosome","Marker position (bp)","P_value","SNP effect","QTN-to-residual variance ratio","QTN variance","Residual variance","MAF","r2 (%)","Genotype for code 1","Critical P-value")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      if(mrenv$inputform==2){
        #output result1 using mrMLM character format
        mrenv$outATCG<-matrix(mrenv$outATCG,,1)
        meadd<-matrix(" ",nrow(mrenv$parms),1)
        if(newp>1e-4){
          newp<-round(newp,4)
        }else{
          newp<-sprintf("%.4e",newp)
        }
        meadd[1,1]<-newp
        mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],parmeter,mrenv$outATCG,meadd)
        colnames(mrenv$parmsShow)<-c("Marker","Chromosome","Marker position (bp)","P_value","SNP effect","QTN-to-residual variance ratio","QTN variance","Residual variance","MAF","r2 (%)","Genotype for code 1","Critical P-value")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      if(mrenv$inputform==3){
        #output result1 using TASSEL format
        mrenv$outATCG<-matrix(mrenv$outATCG,,1)
        mrenv$outATCG<-unlist(strsplit(mrenv$outATCG,""))
        mrenv$outATCG<-matrix(mrenv$outATCG[c(TRUE,FALSE)],,1)
        
        meadd<-matrix(" ",nrow(mrenv$parms),1)
        if(newp>1e-4){
          newp<-round(newp,4)
        }else{
          newp<-sprintf("%.4e",newp)
        }
        meadd[1,1]<-newp
        
        mrenv$parmsShow<-data.frame(mrenv$genRaw[-1,1],parmeter,mrenv$outATCG,meadd)
        colnames(mrenv$parmsShow)<-c("Marker","Chromosome","Marker position (bp)","P_value","SNP effect","QTN-to-residual variance ratio","QTN variance","Residual variance","MAF","r2 (%)","Genotype for code 1","Critical P-value")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      progress_bar$setFraction(95/100)
      
      Xemma<-data.frame(chr.locus=xnames,REML.LRT.c2.new)#all
      vid<-which(as.numeric(Xemma$ps)<=mrenv$svpal)
      if(length(vid)==1){
        snp.emma.opt<-matrix(xraw[vid,],1,)
        xname.emma.opt<-matrix(snp.emma.opt[,1:2],1,)
        snp4<-matrix(snp.emma.opt[,3:dim(snp.emma.opt)[2]],1,)
        xdata<-t(snp4)
        xdata<-matrix(xdata,,1)
      }else{
        snp.emma.opt<-as.matrix(xraw[vid,])
        xname.emma.opt<-snp.emma.opt[,1:2]
        snp4<-snp.emma.opt[,3:dim(snp.emma.opt)[2]]
        xdata<-t(snp4)
      }
      
      xdata<-t(snp4)
      ydata<-Y
      u1<-ebayes_EM(x=W,z=xdata,y=ydata)
      emma.lod<-likelihood(xxn=W,xxx=xdata,yn=ydata,bbo=u1$u)
      idslod<-which(emma.lod>=mrenv$svmlod)
      
      if(length(idslod)==1){
        snp.eb.opt<-matrix(snp.emma.opt[idslod,],1,)
        snp.maf<-snp.eb.opt[,3:dim(snp.eb.opt)[2]]
        snp.maf<-matrix(snp.maf,1,)
        maf.snp<-apply(snp.maf,1,maf.fun)
        chrlocus<-matrix(xname.emma.opt[idslod,],1,)
      }else{
        snp.eb.opt<-snp.emma.opt[idslod,]
        snp.maf<-snp.eb.opt[,3:dim(snp.eb.opt)[2]]
        maf.snp<-apply(snp.maf,1,maf.fun)
        chrlocus<-as.matrix(xname.emma.opt[idslod,])
      }
      
      maf.snp.1<-unlist(maf.snp)
      maf.snp.2<-matrix(maf.snp.1,nrow = 4)
      maf.snp.3<-t(maf.snp.2)
      maf.snp.4<-data.frame(maf.snp.3)
      names(maf.snp.4)<-c("p1","p2","p3","maf")
      
      pve.all<-(pve.fun(u1$u[idslod],maf.snp.4)/var(Y))*100
      qtneffect<-matrix(u1$u[idslod],,1)
      lodscore<-matrix(emma.lod[idslod],,1)
      maff<-matrix(maf.snp.4$maf,,1)
      r2<-matrix(pve.all,,1)
      wanbefore<-cbind(qtneffect,lodscore,maff,r2)
      wanbefore[which(abs(wanbefore)>1e-4)]<-round(wanbefore[which(abs(wanbefore)>1e-4)],4)
      wanbefore[which(abs(wanbefore)<1e-4)]<-as.numeric(sprintf("%.4e", wanbefore[which(abs(wanbefore)<1e-4)]))
      wan<-cbind(chrlocus,wanbefore)
      phenotype.var<-var(Y)
      sigma2<-u1$sigma2
      pee<-matrix("",dim(wan)[1],1)
      vess<-matrix("",dim(wan)[1],1)
      pee[1]<-round(phenotype.var,4)
      vess[1]<-round(sigma2,4)
      
      wan_len<-dim(wan)[1]
      marker<-as.character()
      snp<-as.character()
      for(i in 1:wan_len){
        chr_pos<-which(mrenv$parmsShow[,2]==wan[i,1])
        new_matrix<-mrenv$parmsShow[chr_pos,]
        posi_pos<-which(new_matrix[,3]==wan[i,2])
        mark<-matrix(new_matrix[posi_pos,1],1,)
        marker<-rbind(marker,mark)
        sn<-matrix(new_matrix[posi_pos,11],1,)
        snp<-rbind(snp,sn)
      }
      progress_bar$setFraction(98/100)
      
      mrenv$wan<-cbind(marker,wan,snp,vess,pee)
      colnames(mrenv$wan)<-c("Marker","Chromosome","Marker position (bp)","QTN effect","LOD score","MAF","r2 (%)","Genotype  for code 1","Var_Error","Var_phen(total)")
      wan<-mrenv$wan
      if(exists("wan")==FALSE||is.null(wan)==TRUE)
      {
        gmessage("No result meets the requirements in the second step!","Info",icon="info")
      }else{
        tbdfe8<-gdfedit(wan,container=nb1,expand=TRUE,label="Result2")
      }
      progress_bar$setFraction(100/100)
      progress_bar$setText("All done.")
      
    }
  })
  
  addhandlerclicked(clear,handler = function(h,...){
    progress_bar$setFraction(0/100)
    progress_bar$setText("")
    m<-length(nb1)
    for(i in 1:(m-1)){
      if(is.null(mrenv$wan)==FALSE){
        wan<-NULL
        svmlod<-NULL
        svpal<-NULL
        vif<-NULL
        snpname<-NULL
        rm(snpname,envir = as.environment(mrenv))
        rm(vif,envir = as.environment(mrenv))
        rm(wan,envir = as.environment(mrenv))
        rm(svmlod,envir = as.environment(mrenv))
        rm(svpal,envir = as.environment(mrenv))
        
      }
      else if(is.null(mrenv$parmsShow)==FALSE){
        parmsShow<-NULL
        parms<-NULL
        rm(parms,envir = as.environment(mrenv))
        rm(parmsShow,envir = as.environment(mrenv))
        
      }else if(is.null(mrenv$gen)==FALSE){
        gen<-NULL
        sameName<-NULL
        needloc<-NULL
        needGen<-NULL
        newPhe<-NULL
        newGen<-NULL
        outATCG<-NULL
        inputform<-NULL
        flagREMLE<-NULL
        genRaw<-NULL
        pheRaw<-NULL
        rm(flagREMLE,envir = as.environment(mrenv))
        
        rm(gen,envir = as.environment(mrenv))
        rm(inputform,envir = as.environment(mrenv))
        rm(genRaw,envir = as.environment(mrenv))
        rm(pheRaw,envir = as.environment(mrenv))
        rm(sameName,envir = as.environment(mrenv))
        rm(needloc,envir = as.environment(mrenv))
        rm(needGen,envir = as.environment(mrenv))
        rm(newPhe,envir = as.environment(mrenv))
        rm(newGen,envir = as.environment(mrenv))
      }else if(is.null(mrenv$phe)==FALSE){
        phe<-NULL
        rm(phe,envir = as.environment(mrenv))
      }else if(is.null(mrenv$kk)==FALSE){
        kkRaw<-NULL
        kkShow<-NULL
        kk<-NULL
        rm(kk,envir = as.environment(mrenv))
        rm(kkRaw,envir = as.environment(mrenv))
        rm(kkShow,envir = as.environment(mrenv))
        
      }else if(is.null(mrenv$psmatrix)==FALSE){
        psmatrix<-NULL
        psmatrixRaw<-NULL
        flagps<-NULL
        rm(flagps,envir = as.environment(mrenv))
        rm(psmatrix,envir = as.environment(mrenv))
        rm(psmatrixRaw,envir = as.environment(mrenv))
      }
      dispose(nb1)
    }
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
          gmessage("No result meets the requirements in the second step!","Info",icon="info")
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
    plotlabel3 <- glabel("Word resolution (1/72 inch, ppi):", container = plotlyt)
    plotedit3 <- gedit("20", container = plotlyt)
    pointsizevalue <- as.numeric(svalue(plotedit3))
    plotbt <- gbutton(" Save ", container = plotlyt)
    plotmancl <- gbutton("  Cancel  ", container = plotlyt)
    combo_box1 <- gcombobox(selected=1,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotlyt)
    combo_box2 <- gcombobox(selected=3,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotlyt)
    plotlabel4 <- glabel("Chromosome color (odd):", container = plotlyt)
    plotlabel5 <- glabel("Chromosome color (even):", container = plotlyt)
    plotlabel6 <- glabel("Figure resolution (ppi):", container = plotlyt)
    plotedit6 <- gedit("72", container = plotlyt)
    plotlabel7 <- glabel("Critical value for Manhattan Plot:",container=plotlyt)
    plotedit7 <- gwedit<-gedit(sprintf("%.6s",-log10(mrenv$mannewp)),coerce.with=as.numeric,container=plotlyt)
    plotlabel8 <- glabel("Note: The parameter settings of preview figure are for general resolution, see manual for parameter settings of high resolution figure (*.png,*.tiff,*.jpeg).", container = plotlyt)
    
    svgwline<-svalue(plotedit7)
    mrenv$standline<-svgwline
    if(mrenv$manstandchoice==1)
    {
      mrenv$standline<--log10(mrenv$mannewp)
      mrenv$manstandchoice<-mrenv$manstandchoice+1
    }else{
      mrenv$standline<-svgwline
    }
    
    plotlyt[2, 1] <- plotlabel1
    plotlyt[2, 2] <- plotedit1
    plotlyt[2, 4] <- plotlabel2  
    plotlyt[2, 5] <- plotedit2  
    plotlyt[3, 1] <- plotlabel3
    plotlyt[3, 2] <- plotedit3
    plotlyt[3, 4] <- plotlabel6
    plotlyt[3, 5] <- plotedit6
    plotlyt[4, 1] <- plotlabel4
    plotlyt[4, 2] <- combo_box1
    plotlyt[4, 4] <- plotlabel5
    plotlyt[4, 5] <- combo_box2
    plotlyt[5, 1] <- plotlabel7
    plotlyt[5, 2] <- plotedit7
    plotlyt[5, 4] <- plotbt
    plotlyt[5, 5] <- plotmancl
    plotlyt[6, 1:5] <- plotlabel8
    visible(plotwin)<-TRUE
    
    ress1<-mrenv$ress1
    snpname<-mrenv$snpname
    bpnumber <- numeric()
    chrnum <- unique(ress1[,1])
    
    for(i in 1:length(chrnum))
    {
      bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(ress1[,1]==chrnum[i])))))
    }
    
    addHandlerChanged(ggpw, handler=function(h,...) {
      color1 <- svalue(combo_box1)
      color2 <- svalue(combo_box2)
      parms <- data.frame(ress1,snpname,bpnumber)
      colnames(parms)<-c("Chromosome","Position","P_wald","SNPname","BPnumber")
      manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P_wald",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = mrenv$standline)
      
    })
    addhandlerclicked(plotbt, handler = function(h, ...) {
      widvalue <- as.numeric(svalue(plotedit1))
      heightvalue <- as.numeric(svalue(plotedit2))
      pointsizevalue <- as.numeric(svalue(plotedit3))
      color1 <- svalue(combo_box1)
      color2 <- svalue(combo_box2)
      resppi <- as.numeric(svalue(plotedit6))
      output <- gfile(text = "Save a file...", type = "save", 
                      filter = list(`All files` = list(patterns = c("*")), 
                                    `TIFF files` = list(patterns = c("*.tiff")),
                                    `PNG files` = list(patterns = c("*.png")),
                                    `JPEG files` = list(patterns = c("*.jpeg"))))
      if((length(grep(".png",output))==1)||(length(grep(".PNG",output))==1)){
        png(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue,res=resppi)
      }else if((length(grep(".tiff",output))==1)||(length(grep(".TIFF",output))==1)){
        tiff(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue,res=resppi)
      }else if((length(grep(".jpeg",output)==1))||(length(grep(".JPEG",output))==1)){
        jpeg(output, width=widvalue, height=heightvalue, units= "px", pointsize = pointsizevalue,res=resppi)  
      }else{
        gmessage("Please input correct image format !")
      }
      parms <- data.frame(ress1,snpname,bpnumber)
      colnames(parms)<-c("Chromosome","Position","P_wald","SNPname","BPnumber")
      manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P_wald",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = mrenv$standline)
      dev.off()
      
    })
    addHandlerClicked(plotmancl,handler=function(h,...){
      dispose(plotwin)
    })
  })
  
  addHandlerClicked(qqplot,handler=function(h,...){
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
    plotqqlabel3 <- glabel("Word resolution (1/72 inch, ppi):", container = plotqqlyt)
    plotqqedit3 <- gedit("20", container = plotqqlyt)
    pointsizeqqvalue <- as.numeric(svalue(plotqqedit3))
    plotqqbt <- gbutton(" Save ", container = plotqqlyt)
    plotqqcl <- gbutton("  Cancel  ", container = plotqqlyt)
    plotqqlabel4 <- glabel("Figure resolution (ppi):", container = plotqqlyt)
    plotqqedit4 <- gedit("72", container = plotqqlyt)
    plotqqlabel5 <- glabel("Point color:", container = plotqqlyt)
    combo_box1 <- gcombobox(selected=2,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotqqlyt)
    plotqqlabel6 <- glabel("Line color:", container = plotqqlyt)
    combo_box2 <- gcombobox(selected=3,c("blue","black","red","yellow","green","pink","purple","gray","brown"),editable = TRUE,container = plotqqlyt)
    plotqqlabel8 <- glabel("Note: The parameter settings of preview figure are for general resolution, see manual for high resolution figure (*.png,*.tiff,*.jpeg).", container = plotqqlyt)
    
    plotqqlyt[2, 1] <- plotqqlabel1
    plotqqlyt[2, 2] <- plotqqedit1
    plotqqlyt[2, 4] <- plotqqlabel2
    plotqqlyt[2, 5] <- plotqqedit2
    plotqqlyt[3, 1] <- plotqqlabel3
    plotqqlyt[3, 2] <- plotqqedit3
    plotqqlyt[3, 4] <- plotqqlabel4
    plotqqlyt[3, 5] <- plotqqedit4
    plotqqlyt[4, 1] <- plotqqlabel5
    plotqqlyt[4, 2] <- combo_box1
    plotqqlyt[4, 4] <- plotqqlabel6
    plotqqlyt[4, 5] <- combo_box2
    plotqqlyt[5, 4] <- plotqqbt
    plotqqlyt[5, 5] <- plotqqcl 
    plotqqlyt[6, 1:5] <- plotqqlabel8
    visible(plotwin1)<-TRUE
    
    
    addHandlerChanged(ggpw1, handler=function(h,...){
      color1 <- svalue(combo_box1)
      color2 <- svalue(combo_box2)
      pvalue<-as.matrix(mrenv$ress1)
      ps<-pvalue[,3]
      obs.x<-sort(ps)
      newobs.x<-obs.x[obs.x<1]
      n<-length(newobs.x)
      es<-(1:n)/(n+1)
      x<--log10(es)
      y<--log10(newobs.x)
      y<-y-0.3
      plot(x,y,xlim=c(0.3,max(x)),ylim=c(0.3,max(y)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color1)
      abline(0,1,col=color2)
    })
    
    addhandlerclicked(plotqqbt, handler = function(h, ...) {
      widqqvalue <- as.numeric(svalue(plotqqedit1))
      heightqqvalue <- as.numeric(svalue(plotqqedit2))
      pointsizeqqvalue <- as.numeric(svalue(plotqqedit3))
      resppi <- as.numeric(svalue(plotqqedit4))
      color1 <- svalue(combo_box1)
      color2 <- svalue(combo_box2)
      
      output <- gfile(text = "Save a file...", type = "save", 
                      filter = list(`All files` = list(patterns = c("*")), 
                                    `TIFF files` = list(patterns = c("*.tiff")),
                                    `PNG files` = list(patterns = c("*.png")),
                                    `JPEG files` = list(patterns = c("*.jpeg"))))
      if((length(grep(".png",output))==1)||(length(grep(".PNG",output))==1)){
        png(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue, res = resppi)
      }else if((length(grep(".tiff",output))==1)||(length(grep(".TIFF",output))==1)){
        tiff(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue, res = resppi)
      }else if((length(grep(".jpeg",output)==1))||(length(grep(".JPEG",output))==1)){
        jpeg(output, width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue, res = resppi)  
      }else{
        gmessage("Please input correct image format !")
      }
      pvalue<-as.matrix(mrenv$ress1)
      ps<-pvalue[,3]
      obs.x<-sort(ps)
      newobs.x<-obs.x[obs.x<1]
      n<-length(newobs.x)
      es<-(1:n)/(n+1)
      x<--log10(es)
      y<--log10(newobs.x)
      y<-y-0.3
      plot(x,y,xlim=c(0.3,max(x)),ylim=c(0.3,max(y)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color1)
      abline(0,1,col=color2)
      dev.off()
      
    })
    
    addHandlerClicked(plotqqcl,handler=function(h,...){
      dispose(plotwin1)
    })
  })
}

addHandlerClicked(FASTbutton,handler=function(h,...){
  FASTSub()
})

}