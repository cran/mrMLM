mrMLM<-function(fileGen=NULL,filePhe=NULL,fileKin=NULL,filePS=NULL,PopStrType=NULL,fileCov=NULL,Genformat=NULL,method=NULL,
                Likelihood="REML",trait=NULL,SearchRadius=20,CriLOD=NULL,SelectVariable=50,
                Bootstrap=FALSE,DrawPlot=TRUE,Plotformat="tiff",dir=NULL,PC=FALSE,RAM=4){


if(DrawPlot==TRUE){

  manhattan_mrMLM<-function(data_in,data_fin,mar=c(2.9,2.8,0.7,2.8),VerLabDis=1.5,HorLabDis=1.5,
                            HorTckDis=0.2,VerTckDis=0.4,label_size=0.8,CoorLwd=5,
                            TckLen=-0.03,TckLwd=0.7,log_times=2,LOD_times=1.2,lodline){

    ###########Data process#################
    ###########intermediate result
    method<-unique(data_in[,3])
    data_method<-list(NULL)
    for(i in 1:length(method)){
      data_method[[i]]<-data_in[which(data_in[,3]==method[i]),]
    }
    logp_4method<-numeric()
    for(i in 1:length(method)){
      method_p<-data_method[[i]][,8]
      logp_4method<-cbind(logp_4method,method_p)
    }
    logp_4method<-apply(logp_4method,2,as.numeric)
    p_4method<-10^-logp_4method
    p_median<-apply(p_4method,1,median)
    locsub<-which(p_median==0)
    pmin<-min(p_median[p_median!=0])
    subvalue<-10^(1.1*log10(pmin))
    p_median[locsub]<-subvalue
    data_p<-as.matrix(p_median)
    data_num<-as.matrix(seq(1:length(p_median)))
    data_chr<-as.matrix(data_method[[1]][,5])
    data_pos<-as.matrix(data_method[[1]][,6])
    manresult<-cbind(data_chr,data_pos,data_p,data_num)
    manresult<-apply(manresult,2,as.numeric)
    colnames(manresult)<-c("Chromosome","BPnumber","P-value","SNPname")
    manresult<-as.data.frame(manresult)
    #######final result##################
    data_fin_method<-unique(data_fin[,3])
    data_fin_method_length<-1:length(unique(data_fin[,3]))
    for(r in 1:length(unique(data_fin[,3]))){
      data_fin[which(data_fin[,3]==data_fin_method[r]),3]<-r
    }
    data_fin_mark<-matrix(data_fin[,c(5,6,8,3)],,4)
    data_fin_mark<-matrix(apply(data_fin_mark,2,as.numeric),,4)
    data_fin_mark_chr<-matrix(data_fin_mark[order(data_fin_mark[,1]),],,4)
    data_fin_mark_order<-numeric()
    for(i in c(unique(data_fin_mark_chr[,1]))){
      data_fin_mark_erery_chr<-matrix(data_fin_mark_chr[which(data_fin_mark_chr[,1]==i),],,4)
      data_fin_mark_pos<-matrix(data_fin_mark_erery_chr[order(data_fin_mark_erery_chr[,2]),],,4)
      all_pos<-unique(data_fin_mark_pos[,2])
      all_pos_maxlod<-numeric()
      for(ii in 1:length(all_pos)){
        all_pos_every<-matrix(data_fin_mark_pos[which(data_fin_mark_pos[,2]==all_pos[ii]),],,4)
        lod_me<-median(all_pos_every[,3])
        all_pos_every_median<-c(all_pos_every[1,1:2],lod_me,all_pos_every[1,4])
        if(nrow(all_pos_every)>=2){
          all_pos_every_median<-c(all_pos_every[1,1:2],lod_me,max(data_fin_mark[,4])+1)
        }
        all_pos_maxlod<-rbind(all_pos_maxlod,all_pos_every_median)
      }
      data_fin_mark_order<-rbind(data_fin_mark_order,all_pos_maxlod)
    }
    snpOfInterest<-numeric()
    for(i in c(unique(data_fin_mark_order[,1]))){
      manresult_chr<-manresult[which(manresult[,1]==i),]
      data_fin_mark_order_chr<-matrix(data_fin_mark_order[which(data_fin_mark_order[,1]==i),],,4)
      mark_loc<-manresult_chr[which(manresult_chr[,2]%in%data_fin_mark_order_chr[,2]),4]
      snpOfInterest<-c(snpOfInterest,mark_loc)
    }
    bpnumber <- numeric()
    chrnum <- unique(manresult[,1])
    for(i in 1:length(chrnum))
    {
      bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(manresult[,1]==chrnum[i])))))
    }
    manresult2<-cbind(manresult[,1],bpnumber,manresult[,3:4])
    colnames(manresult2)<-c("Chromosome","BPnumber","P-value","SNPname")
    ##########prepare for data#############################
    x<-manresult2;col=c("lightgreen","lightskyblue");logp=TRUE
    chr = "Chromosome";bp ="BPnumber";p ="P-value";snp="SNPname";
    highlight<-snpOfInterest
    CHR=BP=P=index=NULL
    d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
    if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
      d$logp <- -log10(d$P)
    } else {
      d$logp <- d$P
    }
    d$pos=NA
    d$index=NA
    ind = 0
    for (i in unique(d$CHR)){
      ind = ind + 1
      d[d$CHR==i,]$index = ind
    }

    nchr = length(unique(d$CHR))
    if (nchr==1) { ## For a single chromosome
      ## Uncomment the next two linex to plot single chr results in Mb
      #options(scipen=999)
      #d$pos=d$BP/1e6
      d$pos=d$BP
      ticks=floor(length(d$pos))/2+1
      xlabel = paste('Chromosome',unique(d$CHR),'position')
      labs = ticks
    } else { ## For multiple chromosomes
      lastbase=0
      ticks=NULL
      for (i in unique(d$index)) {
        if (i==1) {
          d[d$index==i, ]$pos=d[d$index==i, ]$BP
        } else {
          lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
          d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
        }
        # Old way: assumes SNPs evenly distributed
        # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
        # New way: doesn't make that assumption
        ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
      }
      xlabel = 'Chromosomes'
      #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
      labs <- unique(d$CHR)
    }

    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)

    ########draw plot#######################

    par(mar=mar)
    def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                     xlim=c(xmin,xmax), ylim=c(0,log_times*max(d$logp)),
                     xlab=xlabel,ylab="",mgp=c(HorLabDis,0,0),cex.lab=label_size)

    dotargs <- list(NULL)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    axis(1, at=ticks, labels=labs,lwd=CoorLwd,tck=TckLen,mgp=c(2.5,HorTckDis,0.5),cex.axis=TckLwd)

    suppressWarnings(axis(2, at=seq(0,log_times*max(d$logp),ceiling(log_times*max(d$logp)/5)),lwd=CoorLwd,tck=TckLen,mgp=c(2.2,VerTckDis,0),cex.axis=TckLwd))
    mtext(expression(-log[10]('P-value')),side=2,line=VerLabDis,cex=label_size,font=1)

    # Create a vector of alternatiting colors
    col=rep(col, max(d$CHR))
    # Add points to the plot
    if (nchr==1) {
      with(d, points(pos, logp, pch=20, col=col[1]))
    } else {
      # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20))
        icol=icol+1
      }
    }
    d.highlight=d[which(d$SNP %in% highlight), ]
    highlight_LOD<-as.numeric(data_fin_mark_order[,3])
    d.highlight<-as.data.frame(cbind(d.highlight,highlight_LOD))

    ################################
    par(new=T)

    def_args <- list(xaxt='n', yaxt='n',bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                     xlim=c(xmin,xmax), ylim=c(0,LOD_times*max(highlight_LOD)),xlab="",ylab="")
    dotargs <- list(NULL)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    suppressWarnings(axis(4,mgp=c(1.4,VerTckDis,0),at=seq(0,LOD_times*max(highlight_LOD),ceiling(LOD_times*max(highlight_LOD)/5)),col="magenta",col.ticks="magenta",col.axis="magenta",lwd=CoorLwd,tck=TckLen,cex.axis=TckLwd))
    mtext("LOD score",side=4,line=VerLabDis,cex=label_size,font=1,col="magenta")
    abline(h=lodline,col="gray25",lty=2,lwd=2)
    peach_colors<-c("magenta","deepskyblue2")
    col_pos<-list(NULL)
    method_num<-sort(unique(data_fin_mark_order[,4]))

    if(max(unique(data_fin[,3]))<max(unique(data_fin_mark_order[,4]))){
      col_pos[[1]]<-which(data_fin_mark_order[,4]==max(method_num))
      col_pos[[2]]<-which(data_fin_mark_order[,4]!=max(method_num))
    }else{
      if(length(unique(data_fin[,3]))==1){
        col_pos[[1]]<-which(data_fin_mark_order[,4]==max(method_num))
      }else{
        col_pos[[1]]<-1:nrow(data_fin_mark_order)

      }
    }
    if(length(col_pos)>1&&length(col_pos[[2]])!=0){
      with(d.highlight, points(pos[col_pos[[2]]], highlight_LOD[col_pos[[2]]], col=peach_colors[2], pch=20))
      with(d.highlight, points(pos[col_pos[[2]]], highlight_LOD[col_pos[[2]]], col=peach_colors[2], pch=20,type="h",lty=2))
      with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20))
      with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20,type="h",lty=2))
    }else{
      with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20))
      with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20,type="h",lty=2))
    }
  }



  QQ_mrMLM<-function(data_in,mar=c(2.5,2.5,1,1),label_size=0.7,TckLen=-0.02,
                     CoorLwd=3,TckLwd=0.6,HorLabDis=1,HorTckDis=0.02,VerLabDis=1.1,
                     VerTckDis=0.3,P_stand=0.9){

    method<-unique(data_in[,3])
    data_method<-list(NULL)
    for(i in 1:length(method)){
      data_method[[i]]<-data_in[which(data_in[,3]==method[i]),]
    }
    logp_4method<-numeric()
    for(i in 1:length(method)){
      method_p<-data_method[[i]][,8]
      logp_4method<-cbind(logp_4method,method_p)
    }
    logp_4method<-apply(logp_4method,2,as.numeric)
    p_4method<-10^-logp_4method
    p_median<-apply(p_4method,1,median)
    locsub<-which(p_median==0)
    pmin<-min(p_median[p_median!=0])
    subvalue<-10^(1.1*log10(pmin))
    p_median[locsub]<-subvalue
    data_p<-as.matrix(p_median)
    p_value<-data_p
    pvalue<-matrix(p_value,,1)
    observed<-sort(pvalue[,1])
    observed<-observed/2
    observed<-observed[which(observed!=0)]
    newobserved<-observed[which(observed<(0.7/2))]
    lobs<--(log10(newobserved))
    expected<-c(1:length(newobserved))
    lexp<--(log10(expected/(length(pvalue)+1)))
    par(mar=mar)
    suppressWarnings(plot(lexp,lobs,xlim=c(0,max(lexp)),ylim=c(0,max(lobs)),xlab=expression('Expected -log'[10]*'(P-value)'),
                          yaxt="n",ylab="",col="blue",pch=20,cex.lab=label_size,tck=TckLen,bty="l",lwd=CoorLwd,
                          lwd.ticks=CoorLwd,cex.axis=TckLwd,mgp=c(HorLabDis,HorTckDis,0)))
    suppressWarnings(axis(2, at=seq(0,max(lobs)),lwd=CoorLwd,tck=TckLen,mgp=c(2.2,VerTckDis,0),cex.axis=TckLwd))
    mtext(expression('Observed -log'[10]*'(P-value)'),side=2,line=VerLabDis,cex=label_size,font=1)
    abline(0,1,col="red")
    box(bty="l",lwd=CoorLwd)
  }
}

if(PC==TRUE){

  #### When PC=TRUE, mrMLM v5.0 may calculate big data (millions of SNPs for thousands of individuals) on personal computer such as desktop or laptop, which have much smaller RAM than server.
  #### 2022-03-14 Wang Jing-Tian

  Genformat <- 1
  parmsShow <- NULL
  outATCG <- NULL
  svrad<-SearchRadius
  svmlod<-CriLOD
  svpal=0.01
  CLO=NULL
  lars1 <- SelectVariable
  Plotformat1<-paste("*.",Plotformat,sep="")
  Plotformat2<-paste("*.",Plotformat,sep="")


  if(RAM<4){
    BLOCK_M=5000
  }else if(RAM>=4&&RAM<=7){
    BLOCK_M=10000
  }else if(RAM>7&&RAM<=15){
    BLOCK_M=20000
  }else if(RAM>15&&RAM<=31){
    BLOCK_M=30000
  }else if(RAM>31){
    BLOCK_M=60000
  }else{
    BLOCK_M=1000
  }

  phy_match <- function(gen_fam,phy){
    phy_ID <- phy[-1,1]
    gen_ID <- t(gen_fam[,2])
    intersect_ID <- intersect(phy_ID,gen_ID)
    match_gen_ID_idex <- match(intersect_ID,gen_ID)
    match_phy_ID_idex <- match(intersect_ID,phy_ID)
    phy <- phy[c(1,(match_phy_ID_idex+1)),]

    phy_match_list <- list(phy,match_gen_ID_idex)
    return(phy_match_list)
  }
  inutpe_transform <- function(gen_bed,gen_bim,gen_fam,match_gen_ID_idex,index_left,index_right){
    gen_block <- as.matrix(cbind(gen_bim[,c(1,4)][index_left:index_right,],t(gen_bed[match_gen_ID_idex,(index_left:index_right)]-1)))
    genRaw_block <- as.matrix(rbind(t(c("rs#","chrom","pos","genotype for code 1")),gen_bim[,c(2,1,4,5)][index_left:index_right,],use.names=FALSE))
    return(list(gen_block,genRaw_block))
  }

  mrMLMFun.PC <- function(gen_bed,gen_bim,gen_fam,phy,phy_match_list,block_m=BLOCK_M,trait){
    print("Running mrMLMFun algorithm with low RAM consumption...")

    # K
    K_PC <- function(gen_bed,match_gen_ID_idex,block_m=20000){

      K_block_i <- matrix(0,length(match_gen_ID_idex),length(match_gen_ID_idex))

      #####
      block_lengthout <- round(ncol(gen_bed)/block_m,0)+1
      if(block_lengthout<3){
        block_point_left <- c(1)
        block_point_left <- ncol(gen_bed)
      }else{
        block_point_left <- round(seq(1,ncol(gen_bed),length.out=block_lengthout),0)
        block_point_right <- block_point_left[-1]-1
        block_point_right[length(block_point_right)] <- block_point_right[length(block_point_right)]+1
        block_point_left <- block_point_left[-length(block_point_left)]
      }

      block_i <- 1
      for(block_i in 1:length(block_point_left)){

        gen_block_i <- gen_bed[match_gen_ID_idex,(block_point_left[block_i] : block_point_right[block_i])]-1
        K_block_i <- K_block_i + multiplication_speed(gen_block_i,t(gen_block_i))
        rm(gen_block_i)
      }

      K <- K_block_i/mean(diag(K_block_i))

      return(K)
    }
    if(is.null(fileKin)){
      K <- K_PC(gen_bed,match_gen_ID_idex,block_m=20000)
    }else{
      fileKin <- as.matrix(fread(fileKin,header = FALSE,stringsAsFactors=T))
      fileKin[1,2:ncol(fileKin)]<-"  "
      kkPre<-as.matrix(fileKin[-1,-1])
      nameKin<-as.matrix(fileKin[-1,1])
      sameGenKin<-intersect(samename_genphy,nameKin)
      locKin<-match(sameGenKin,nameKin)
      K<-kkPre[locKin,locKin]
      K<-matrix(as.numeric(K),nrow=nrow(K))
      rm(kkPre,locKin,sameGenKin,nameKin)
    }
    gc()

    #

    #
    phe <- as.matrix(phy[,trait])
    #####
    block_lengthout <- round(ncol(gen_bed)/block_m,0)+1
    if(block_lengthout<3){
      block_point_left <- c(1)
      block_point_left <- ncol(gen_bed)
    }else{
      block_point_left <- round(seq(1,ncol(gen_bed),length.out=block_lengthout),0)
      block_point_right <- block_point_left[-1]-1
      block_point_right[length(block_point_right)] <- block_point_right[length(block_point_right)]+1
      block_point_left <- block_point_left[-length(block_point_left)]
    }


    ll_read <- matrix(0,1,10)
    genRaw_read <- matrix(0,1,4)
    block_i <- 1
    fin_block <- FALSE
    for(block_i in 1:length(block_point_left)){

      if(block_i!=length(block_point_left)){
        inpute_block <- inutpe_transform(gen_bed,gen_bim,gen_fam,phy_match_list[[2]],block_point_left[block_i],block_point_right[block_i])
        mid_result_i <- mrMLMFun_2.0(gen=inpute_block[[1]],phe=phe,genRaw=inpute_block[[2]],kk=K,fin_block=fin_block,read_ll=NULL,block_i=block_i,match_gen_ID_idex=phy_match_list[[2]])
        ll_read <- rbind(ll_read,mid_result_i$result3)
        genRaw_read <- rbind(genRaw_read,inpute_block[[2]])
        rm(inpute_block,mid_result_i)
        gc()
      }else{
        fin_block=TRUE
        ll_read <- ll_read[-1,]
        genRaw_read <- genRaw_read[-1,]
        inpute_block <- inutpe_transform(gen_bed,gen_bim,gen_fam,phy_match_list[[2]],block_point_left[block_i],block_point_right[block_i])
        total_result <- mrMLMFun_2.0(gen=inpute_block[[1]],phe=phe,genRaw=inpute_block[[2]],kk=K,fin_block=fin_block,read_ll=ll_read,read_genRaw=genRaw_read,block_i=block_i,genq_BED=gen_bed,match_gen_ID_idex=phy_match_list[[2]])

        # write.csv(total_result$result1,paste(dir,"/mid_result.csv",sep=""))
        # write.csv(total_result$result2,paste(dir,"/fin_result.csv",sep=""))

        rm(genRaw_read,ll_read)
        gc()

      }
      #print(block_i)
    }
    return(total_result)
  }
  FASTmrMLM.PC <- function(gen_bed,gen_bim,gen_fam,phy,phy_match_list,block_m=BLOCK_M,trait){
    print("Running FASTmrMLM algorithm with low RAM consumption...")
    # K
    K_FASTmrMLM_PC <- function(gen_bed,match_gen_ID_idex,block_m=20000){

      K_block_i <- matrix(0,length(match_gen_ID_idex),length(match_gen_ID_idex))
      #####
      block_lengthout <- round(ncol(gen_bed)/block_m,0)+1
      if(block_lengthout<3){
        block_point_left <- c(1)
        block_point_left <- ncol(gen_bed)
      }else{
        block_point_left <- round(seq(1,ncol(gen_bed),length.out=block_lengthout),0)
        block_point_right <- block_point_left[-1]-1
        block_point_right[length(block_point_right)] <- block_point_right[length(block_point_right)]+1
        block_point_left <- block_point_left[-length(block_point_left)]
      }

      block_i <- 1
      for(block_i in 1:length(block_point_left)){

        gen_block_i <- gen_bed[match_gen_ID_idex,(block_point_left[block_i] : block_point_right[block_i])]-1
        K_block_i <- K_block_i + multiplication_speed(gen_block_i,t(gen_block_i))
        rm(gen_block_i)
      }

      K <- K_block_i/ncol(gen_bed)

      return(K)
    }
    if(is.null(fileKin)){
      K <- K_FASTmrMLM_PC(gen_bed,match_gen_ID_idex,block_m=20000)
    }else{
      fileKin <- as.matrix(fread(fileKin,header = FALSE,stringsAsFactors=T))
      fileKin[1,2:ncol(fileKin)]<-"  "
      kkPre<-as.matrix(fileKin[-1,-1])
      nameKin<-as.matrix(fileKin[-1,1])
      sameGenKin<-intersect(samename_genphy,nameKin)
      locKin<-match(sameGenKin,nameKin)
      K<-kkPre[locKin,locKin]
      K<-matrix(as.numeric(K),nrow=nrow(K))
      rm(kkPre,locKin,sameGenKin,nameKin)
    }
    gc()

    phe <- as.matrix(phy[,trait])

    #####
    block_lengthout <- round(ncol(gen_bed)/block_m,0)+1
    if(block_lengthout<3){
      block_point_left <- c(1)
      block_point_left <- ncol(gen_bed)
    }else{
      block_point_left <- round(seq(1,ncol(gen_bed),length.out=block_lengthout),0)
      block_point_right <- block_point_left[-1]-1
      block_point_right[length(block_point_right)] <- block_point_right[length(block_point_right)]+1
      block_point_left <- block_point_left[-length(block_point_left)]
    }
    ll_read <- numeric(0)
    genRaw_read <- numeric(0)
    block_i <- 1
    fin_block <- FALSE
    for(block_i in 1:length(block_point_left)){
      if(block_i!=length(block_point_left)){

        inpute_block <- inutpe_transform(gen_bed,gen_bim,gen_fam,phy_match_list[[2]],block_point_left[block_i],block_point_right[block_i])
        mid_result_i <- FASTmrMLM_2.0(gen=inpute_block[[1]],phe=phe,genRaw=inpute_block[[2]],kk=K,fin_block=fin_block,read_ll=NULL,block_i=block_i)
        ll_read <- rbind(ll_read,mid_result_i$result3)
        genRaw_read <- rbind(genRaw_read,inpute_block[[2]])
        rm(inpute_block,mid_result_i)
        gc()
        #print(block_i)
      }else{

        fin_block=TRUE
        inpute_block <- inutpe_transform(gen_bed,gen_bim,gen_fam,phy_match_list[[2]],block_point_left[block_i],block_point_right[block_i])
        total_result <- FASTmrMLM_2.0(gen=inpute_block[[1]],phe=phe,genRaw=inpute_block[[2]],kk=K,fin_block=fin_block,read_ll=ll_read,read_genRaw=genRaw_read,block_i=block_i,genq_BED=gen_bed,match_gen_ID_idex=phy_match_list[[2]])

        # write.csv(total_result$result1,paste(dir,"/FASTmrMLM_mid_result.csv",sep=""))
        # write.csv(as.matrix(total_result$result2),paste(dir,"/FASTmrMLM_fin_result.csv",sep=""),row.names=FALSE)

        rm(genRaw_read,ll_read)
        gc()

      }
    }
    return(total_result)
  }
  FASTmrEMMA.PC <- function(gen_bed,gen_bim,gen_fam,phy,phy_match_list,block_m=BLOCK_M,trait,Likelihood=Likelihood){
    print("Running FASTmrEMMA algorithm with low RAM consumption...")
    #####
    block_lengthout <- round(ncol(gen_bed)/block_m,0)+1
    if(block_lengthout<3){
      block_point_left <- c(1)
      block_point_left <- ncol(gen_bed)
    }else{
      block_point_left <- round(seq(1,ncol(gen_bed),length.out=block_lengthout),0)
      block_point_right <- block_point_left[-1]-1
      block_point_right[length(block_point_right)] <- block_point_right[length(block_point_right)]+1
      block_point_left <- block_point_left[-length(block_point_left)]
    }

    # K
    K_FASTmrEMMA_PC <- function(gen_bed,match_gen_ID_idex){
      K_block_i <- matrix(0,length(match_gen_ID_idex),length(match_gen_ID_idex))
      block_i <- 1
      for(block_i in 1:length(block_point_left)){

        gen_block_i <- t(gen_bed[match_gen_ID_idex,(block_point_left[block_i] : block_point_right[block_i])]/2)

        ###
        flags <- matrix(as.double(rowMeans(gen_block_i,na.rm=TRUE) > 0.5),nrow(gen_block_i),ncol(gen_block_i))
        gen_block_i[!is.na(gen_block_i) & (gen_block_i == 0.5)] <- flags[!is.na(gen_block_i) & (gen_block_i == 0.5)]
        mafs <- matrix(rowMeans(gen_block_i,na.rm=TRUE),nrow(gen_block_i),ncol(gen_block_i))
        gen_block_i[is.na(gen_block_i)] <- mafs[is.na(gen_block_i)]
        rm(mafs)
        gc()



        K_block_i <- K_block_i + multiplication_speed(t(gen_block_i),gen_block_i) + multiplication_speed(t(1-gen_block_i),(1-gen_block_i))
        rm(gen_block_i)
      }

      K <- K_block_i/ncol(gen_bed)
      diag(K) <- 1
      return(K)
    }
    if(is.null(fileKin)){
      K <- K_FASTmrEMMA_PC(gen_bed,match_gen_ID_idex)
    }else{
      fileKin <- as.matrix(fread(fileKin,header = FALSE,stringsAsFactors=T))
      fileKin[1,2:ncol(fileKin)]<-"  "
      kkPre<-as.matrix(fileKin[-1,-1])
      nameKin<-as.matrix(fileKin[-1,1])
      sameGenKin<-intersect(samename_genphy,nameKin)
      locKin<-match(sameGenKin,nameKin)
      K<-kkPre[locKin,locKin]
      K<-matrix(as.numeric(K),nrow=nrow(K))
      rm(kkPre,locKin,sameGenKin,nameKin)
    }
    gc()

    phe <- as.matrix(phy[,trait])

    ll_read <- numeric(0)
    genRaw_read <- numeric(0)
    block_i <- 1
    fin_block <- FALSE
    for(block_i in 1:length(block_point_left)){
      if(block_i!=length(block_point_left)){

        inpute_block <- inutpe_transform(gen_bed,gen_bim,gen_fam,phy_match_list[[2]],block_point_left[block_i],block_point_right[block_i])
        mid_result_i <- FASTmrEMMA_2.0(gen=inpute_block[[1]],phe=phe,genRaw=inpute_block[[2]],kk=K,fin_block=fin_block,read_ll=NULL,block_i=block_i)
        ll_read <- rbind(ll_read,mid_result_i$result3)
        genRaw_read <- rbind(genRaw_read,inpute_block[[2]])
        rm(inpute_block,mid_result_i)
        gc()
        #print(block_i)
      }else{
        fin_block=TRUE
        inpute_block <- inutpe_transform(gen_bed,gen_bim,gen_fam,phy_match_list[[2]],block_point_left[block_i],block_point_right[block_i])
        total_result <- FASTmrEMMA_2.0(gen=inpute_block[[1]],phe=phe,genRaw=inpute_block[[2]],kk=K,fin_block=fin_block,read_ll=ll_read,read_genRaw=genRaw_read,block_i=block_i,genq_BED=gen_bed,match_gen_ID_idex=phy_match_list[[2]])
        rm(genRaw_read,ll_read)
        gc()
      }
    }
    return(total_result)
  }
  pKWmEB.PC <- function(gen_bed,gen_bim,gen_fam,phy,phy_match_list,block_m=BLOCK_M,trait){
    print("Running pKWmEB algorithm with low RAM consumption...")
    #####
    block_lengthout <- round(ncol(gen_bed)/block_m,0)+1
    if(block_lengthout<3){
      block_point_left <- c(1)
      block_point_left <- ncol(gen_bed)
    }else{
      block_point_left <- round(seq(1,ncol(gen_bed),length.out=block_lengthout),0)
      block_point_right <- block_point_left[-1]-1
      block_point_right[length(block_point_right)] <- block_point_right[length(block_point_right)]+1
      block_point_left <- block_point_left[-length(block_point_left)]
    }

    # K
    K_PC <- function(gen_bed,match_gen_ID_idex,block_m=20000){

      K_block_i <- matrix(0,length(match_gen_ID_idex),length(match_gen_ID_idex))

      #####
      block_lengthout <- round(ncol(gen_bed)/block_m,0)+1
      if(block_lengthout<3){
        block_point_left <- c(1)
        block_point_left <- ncol(gen_bed)
      }else{
        block_point_left <- round(seq(1,ncol(gen_bed),length.out=block_lengthout),0)
        block_point_right <- block_point_left[-1]-1
        block_point_right[length(block_point_right)] <- block_point_right[length(block_point_right)]+1
        block_point_left <- block_point_left[-length(block_point_left)]
      }

      block_i <- 1
      for(block_i in 1:length(block_point_left)){

        gen_block_i <- gen_bed[match_gen_ID_idex,(block_point_left[block_i] : block_point_right[block_i])]-1
        K_block_i <- K_block_i + multiplication_speed(gen_block_i,t(gen_block_i))
        rm(gen_block_i)
      }

      K <- K_block_i/mean(diag(K_block_i))

      return(K)
    }
    if(is.null(fileKin)){
      K <- K_PC(gen_bed,match_gen_ID_idex,block_m=20000)
    }else{
      fileKin <- as.matrix(fread(fileKin,header = FALSE,stringsAsFactors=T))
      fileKin[1,2:ncol(fileKin)]<-"  "
      kkPre<-as.matrix(fileKin[-1,-1])
      nameKin<-as.matrix(fileKin[-1,1])
      sameGenKin<-intersect(samename_genphy,nameKin)
      locKin<-match(sameGenKin,nameKin)
      K<-kkPre[locKin,locKin]
      K<-matrix(as.numeric(K),nrow=nrow(K))
      rm(kkPre,locKin,sameGenKin,nameKin)
    }
    gc()

    phe <- as.matrix(phy[,trait])

    ll_read <- numeric(0)
    genRaw_read <- numeric(0)
    block_i <- 1
    fin_block <- FALSE
    for(block_i in 1:length(block_point_left)){
      if(block_i!=length(block_point_left)){
        inpute_block <- inutpe_transform(gen_bed,gen_bim,gen_fam,phy_match_list[[2]],block_point_left[block_i],block_point_right[block_i])
        mid_result_i <- pKWmEB_2.0(gen=inpute_block[[1]],phe=phe,genRaw=inpute_block[[2]],kk=K,fin_block=fin_block,read_ll=NULL,block_i=block_i)
        ll_read <- rbind(ll_read,mid_result_i$result3)
        genRaw_read <- rbind(genRaw_read,inpute_block[[2]])
        rm(inpute_block,mid_result_i)
        gc()
      }else{
        fin_block=TRUE
        inpute_block <- inutpe_transform(gen_bed,gen_bim,gen_fam,phy_match_list[[2]],block_point_left[block_i],block_point_right[block_i])
        total_result <- pKWmEB_2.0(gen=inpute_block[[1]],phe=phe,genRaw=inpute_block[[2]],kk=K,fin_block=fin_block,read_ll=ll_read,read_genRaw=genRaw_read,block_i=block_i,match_gen_ID_idex=phy_match_list[[2]])
        # write.csv(total_result$result1,paste(dir,"/pKWmEB_mid_result.csv",sep=""))
        # write.csv(as.matrix(total_result$result2),paste(dir,"/pKWmEB_fin_result.csv",sep=""),row.names=FALSE)

        #rm(genRaw_read,ll_read)
        gc()
      }
    }
    return(total_result)
  }
  pLARmEB.PC <- function(gen_bed,gen_bim,gen_fam,phy,phy_match_list,trait){
    print("Running pLARmEB algorithm with low RAM consumption...")
    phe <- as.matrix(phy[,trait])
    total_result <- pLARmEB_2.0(phe,match_gen_ID_idex=phy_match_list[[2]],CriLOD=3)

    return(total_result)
  }

  mrMLMFun_2.0<-function(gen,phe,genRaw,kk,fin_block = FALSE,read_ll=NULL,read_genRaw=NULL,block_i,genq_BED=NULL,match_gen_ID_idex=NULL){

    inputform<-Genformat

    if(is.null(kk)){
      if(is.null(gen)==TRUE)
      {
        warning("Please input correct genotypic dataset !")
      }else{
        envgen<-t(gen[,3:ncol(gen)])
        m<-ncol(envgen)
        n<-nrow(envgen)
        kk1<-matrix(0,n,n)
        # for(k in 1:m){
        #   z<-as.matrix(envgen[,k])
        #   kk1<-kk1+z%*%t(z)
        # }
        kk1<-mrMLM::multiplication_speed(envgen,t(envgen))
        cc<-mean(diag(kk1))
        kk1<-kk1/cc
        kk<-as.matrix(kk1)
      }
      rm(envgen,kk1)
      gc()
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
      warning("Please input search radius (kb) of candidate gene: > 0 !")
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
          p<-pchisq(f,1,lower.tail = F)
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
        p_value<-pchisq(lrt,1,lower.tail = F)
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

      m<-nrow(gen)
      n<-length(match_gen_ID_idex)
      name<-gen[,1:2]
      genq<-gen[,3:ncol(gen)]
      gen<-t(genq)


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

      rm(qq)
      gc()

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
          cl.cores <- detectCores()-2
        }
      }
      if(cl.cores < 1){cl.cores <- 1}
      cl <- makeCluster(cl.cores)
      registerDoParallel(cl)

      if((flagps==1)||(is.null("psmatrix")))
      {
        k<-numeric()
        ff=foreach(k=1:m, .multicombine=TRUE, .combine = 'rbind')%dopar%
          {
            #browser()
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
            p_lrt<-pchisq(lrt,1,lower.tail = F)
            wald<-(gamma/stderr)^2
            p_wald<-pchisq(wald,1,lower.tail = F)
            parm0<-c(1,name[k,1],name[k,2],beta,sigma2,sigma2g,gamma,stderr,wald,p_wald)

          }
        stopCluster(cl)

        ll<-rbind(ll,ff)
      }else if(flagps==0){
        k<-numeric()
        ff=foreach(k=1:m, .multicombine=TRUE, .combine = 'rbind')%dopar%
          {
            #browser()
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
            p_lrt<-pchisq(lrt,1,lower.tail = F)
            wald<-(gamma/stderr)^2
            p_wald<-pchisq(wald,1,lower.tail = F)
            parm0<-c(1,name[k,1],name[k,2],beta,sigma2,sigma2g,gamma,stderr,wald,p_wald)

          }
        stopCluster(cl)

        ll<-rbind(ll,ff)
      }


      rm(uu,kk)
      gc()

      parms<-ll
      parms<-matrix(parms,,10)

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

      if(nrow(orderno)>1){
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
        colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (mrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (mrMLM)'","Genotype for code 1","Significance")
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
        colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (mrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (mrMLM)'","Genotype for code 1","Significance")
      }
      if(inputform==3){
        #output result1 using TASSEL format
        parmsShow<-parms[,-1]
        outATCG<-matrix(outATCG,,1)
        #outATCG<-unlist(strsplit(outATCG,""))
        #outATCG<-matrix(outATCG[c(TRUE,FALSE)],,1)
        meadd<-matrix(1,nrow(parms),1)
        meadd[which(parms[,10]<newp),1]<-sprintf("%.4e",newp)
        meadd[which(parms[,10]>=newp),1]<-"  "
        tempparms<-parms[,4:10]
        tempparms[,7]<--log10(tempparms[,7])
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        parmsShow<-cbind(genRaw[-1,1],parms[,2:3],tempparms,outATCG,meadd)
        colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (mrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (mrMLM)'","Genotype for code 1","Significance")
      }




      ###### ###### ###### ###### ###### ###### ###### ######
      ###### ###### ###### ###### ###### ###### ###### ######

      if(fin_block==FALSE){
        wan <- NULL
        output<-list(result1=parmsShow,result2=wan,result3=ll)
        return(output)
      }else{


        ll <- rbind(read_ll,ll)
        genRaw <- rbind(read_genRaw,genRaw[-1,])

        if(length(c(which(genRaw[,1]=="rs#")))!=1){
          genRaw <- genRaw[-c(which(genRaw[,1]=="rs#")[-1]),]
        }

        parms<-ll
        parms<-matrix(parms,,10)

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

        if(nrow(orderno)>1){
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
          colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (mrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (mrMLM)'","Genotype for code 1","Significance")
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
          colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (mrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (mrMLM)'","Genotype for code 1","Significance")
        }
        if(inputform==3){
          #output result1 using TASSEL format
          parmsShow<-parms[,-1]
          outATCG<-matrix(outATCG,,1)
          #outATCG<-unlist(strsplit(outATCG,""))
          #outATCG<-matrix(outATCG[c(TRUE,FALSE)],,1)
          meadd<-matrix(1,nrow(parms),1)
          meadd[which(parms[,10]<newp),1]<-sprintf("%.4e",newp)
          meadd[which(parms[,10]>=newp),1]<-"  "
          tempparms<-parms[,4:10]
          tempparms[,7]<--log10(tempparms[,7])
          tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
          tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
          parmsShow<-cbind(genRaw[-1,1],parms[,2:3],tempparms,outATCG,meadd)
          colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (mrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (mrMLM)'","Genotype for code 1","Significance")
        }



        ####### ######## ######

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
          xxx0<-t(genq_BED[match_gen_ID_idex,c(g0)]-1)
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
        if(nrow(a2)>1){
          xx<-genq_BED[match_gen_ID_idex,c(a2)]-1
        }else{
          xx<-genq_BED[match_gen_ID_idex,c(a2)]-1
        }
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
            #if ((nrow(w3))==0){ww<-0}change20190125

            if ((nrow(w3)!=0)&&(w3[1]>0)){
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
            #if ((nrow(yang))==0){ww<-0}change20190125
          }
          #ww<-as.matrix(ww)change20190125
          needww<-ww
          if (length(ww)>=1){
            #ww<-as.matrix(ww)chang20190125
            if (length(ww)>1){
              ww<-as.matrix(ww)#change20190125

              if((flagps==1)||(exists("psmatrix")==FALSE))
              {
                ex<-cbind(matrix(1,(nrow(xx)),1),(genq_BED[match_gen_ID_idex,c(ww)]-1))
              }else if(flagps==0)
              {
                ex<-cbind(cbind(matrix(1,(nrow(xx)),1),psmatrix),(genq_BED[match_gen_ID_idex,c(ww)]-1))
              }

            }else{
              if((flagps==1)||(exists("psmatrix")==FALSE))
              {
                ex<-cbind(matrix(1,(nrow(xx)),1),as.matrix(genq_BED[match_gen_ID_idex,c(ww)]-1))
              }else if(flagps==0)
              {
                ex<-cbind(cbind(matrix(1,(nrow(xx)),1),psmatrix),as.matrix(genq_BED[match_gen_ID_idex,c(ww)]-1))
              }
            }
            rm(genq)
            gc()

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

            #x<-gen[3:nrow(gen),]
            xxxx<-as.matrix(genq_BED[match_gen_ID_idex,ww]-1)

            #rm(x)
            gc()

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
            log10P <- as.matrix(-log10(pchisq(lo*4.605,1,lower.tail = F)))

            log10P[which(abs(log10P)>=1e-4)] <- round(log10P[which(abs(log10P)>=1e-4)],4)
            log10P[which(abs(log10P)<1e-4)] <- as.numeric(sprintf("%.4e",her[which(abs(log10P)<1e-4)]))

            if (length(ww)>1){
              wan<-data.frame(parmsShow[needww,1],chr_pos[ww,],eeff,lo,log10P,her,maf,parmsShow[needww,11])
              wan<-wan[order(wan[,2]),]
              wan<-data.frame(wan,vees,pees)
            }else{

              wan<-data.frame(parmsShow[needww,1],t(as.matrix(chr_pos[ww,])),eeff,lo,log10P,her,maf,parmsShow[needww,11],vees,pees)
            }
            colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen (total)")
          }
        }

        if(is.null(parmsShow)==FALSE){
          parmsShow<-parmsShow[,-c(4,5,6,8,9,12)]
        }

        output<-list(result1=parmsShow,result2=wan)
        return(output)



      }









      ###### ###### ###### ###### ###### ###### ###### ######
      ###### ###### ###### ###### ###### ###### ###### ######







      # rm(genRaw)
      # gc()
      #

    }
  }
  FASTmrEMMA_2.0<-function(gen,phe,genRaw,kk,fin_block = FALSE,read_ll=NULL,read_genRaw=NULL,block_i,genq_BED=NULL,match_gen_ID_idex=NULL){

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
        #stopifnot(n0+nh+n1+nNA == length(snps))
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
          #los<-intersect(which(!is.na(snps)),which(snps==0.5))
          dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
          rm(flags)
          gc()
          flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
          rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
          rm(flags,snps)
          gc()
          snps <- rbind(dsnps,rsnps)
          rm(dsnps,rsnps)
          gc()
        }
        if ( use == "all" ) {
          mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
          #losna<-which(is.na(snps))
          snps[is.na(snps)] <- mafs[is.na(snps)]
          rm(mafs)
          gc()
        }
        else if ( use == "complete.obs" ) {
          snps <- snps[rowSums(is.na(snps))==0,]
        }
        n <- ncol(snps)
        #K<-(t(snps)%*%snps+t(1-snps)%*%(1-snps))/nrow(snps)
        K<-(mrMLM::multiplication_speed(t(snps),snps)+mrMLM::multiplication_speed(t(1-snps),(1-snps)))/nrow(snps)
        diag(K) <- 1
        return(K)
      }


      if(is.null(gen)==TRUE)
      {
        warning("Please input correct genotype dataset !")
      }else{
        snp8<-gen[,3:ncol(gen)]
        kk<-emma.kinship(snp8)
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
          p<-pchisq(f,1,lower.tail = F)
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
      ################################################
      #likelihood
      FASTmrEMMA.delta.ML.LL.c<-function(logdelta,X,M,M.y,yMy,n){
        #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
        #n<-dim(M)[1]
        delta <-  exp(logdelta)
        ci<-as.numeric(crossprod(X))
        delta1<-as.numeric(t(X)%*%M%*%X)
        xMy<-as.numeric(crossprod(X,M.y))
        return(0.5*(n*((log(n/(2*pi))-log(as.numeric(yMy)-delta*(xMy)^2/(1+delta*delta1)))-1)-log(delta*ci+1)))
      }
      #dML

      FASTmrEMMA.delta.ML.dLL.c<-function(logdelta,X,M,M.y,yMy,n){
        #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
        #n<-dim(M)[1]
        delta <-  exp(logdelta)
        ci<-as.numeric(crossprod(X))
        delta1<-as.numeric(t(X)%*%M%*%X)
        xMy<-as.numeric(crossprod(X,M.y))
        return(-0.5*ci/(1+delta*ci)+0.5*n/((1+delta*delta1)*(as.numeric(yMy)*(1+delta*delta1)/(xMy^2)-delta)))
      }
      #restrict likelihood 20190902
      FASTmrEMMA.delta.REML.LL.c<-function(logdelta,X,M,M.y,yMy,v){
        #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
        #v<-n-1
        delta <-  exp(logdelta)
        #ci<-crossprod(X)
        delta1<-as.numeric(t(X)%*%M%*%X)
        xMy<-as.numeric(crossprod(X,M.y))
        return(0.5*(v*((log(v/(2*pi))-log(as.numeric(yMy)-delta*(xMy)^2/(1+delta*delta1)))-1)-log(delta*delta1+1)))
      }
      #dREML
      FASTmrEMMA.delta.REML.dLL.c<-function(logdelta,X,M,M.y,yMy,v){
        #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
        #n<-dim(M)[1]
        delta <-  exp(logdelta)
        #ci<-crossprod(X)
        delta1<-as.numeric(t(X)%*%M%*%X)
        xMy<-as.numeric(crossprod(X,M.y))
        return(-0.5*delta1/(1+delta*delta1)+0.5*v/((1+delta*delta1)*(as.numeric(yMy)*(1+delta*delta1)/(xMy^2)-delta)))
      }
      ####################
      #20190906
      FASTmrEMMA.MLE.c<-function(X,M,M.y,yMy,n, ngrids=100, llim=-10, ulim=10, esp=1e-10){

        logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
        m <- length(logdelta)
        #delta <- exp(logdelta)
        dLL<-FASTmrEMMA.delta.ML.dLL.c(logdelta,X,M,M.y,yMy,n)
        optlogdelta <- vector(length=0)
        optLL <- vector(length=0)

        if ( dLL[1] < esp ) {
          optlogdelta <- append(optlogdelta, llim)
          optLL <- append(optLL,FASTmrEMMA.delta.ML.LL.c(llim,X,M,M.y,yMy,n))
        }
        if ( dLL[m-1] > 0-esp ) {
          optlogdelta <- append(optlogdelta, ulim)
          #optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
          optLL <- append(optLL, FASTmrEMMA.delta.ML.LL.c(ulim,X,M,M.y,yMy,n))
        }
        for( i in 1:(m-1) )
        {
          if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
          {
            #r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
            #r <- uniroot(FASTmrEMMA.delta.ML.dLL.c,lower = logdelta[i],upper = logdelta[i+1],X,M,M.y,yMy,n)
            r <- uniroot(FASTmrEMMA.delta.ML.dLL.c,c(logdelta[i],logdelta[i+1]),X=X,M=M,M.y=M.y,yMy=yMy,n=n)

            optlogdelta <- append(optlogdelta, r$root)
            #optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
            optLL <- append(optLL,FASTmrEMMA.delta.ML.LL.c(r$root,X,M,M.y,yMy,n))
          }
        }
        maxdelta <- exp(optlogdelta[which.max(optLL)])
        optLL=replaceNaN(optLL)
        maxLL <- max(optLL)
        xMy<-crossprod(X,M.y)
        xMx<-crossprod(X,(M%*%X))
        maxve <-(yMy-maxdelta*(xMy)^2/(1+maxdelta*xMx))/n
        #(sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
        maxvg <- maxve*maxdelta
        #alpha<-inv()
        return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,delta1=xMx,xMy=xMy))
      }

      FASTmrEMMA.REMLE.c<-function(X,M,M.y,yMy,v, ngrids=100, llim=-10, ulim=10, esp=1e-10){

        logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
        m <- length(logdelta)
        #delta <- exp(logdelta)
        dLL<-FASTmrEMMA.delta.REML.dLL.c(logdelta,X,M,M.y,yMy,v)
        optlogdelta <- vector(length=0)
        optLL <- vector(length=0)
        if ( dLL[1] < esp ) {
          optlogdelta <- append(optlogdelta, llim)
          optLL <- append(optLL,FASTmrEMMA.delta.REML.LL.c(llim,X,M,M.y,yMy,v))
        }
        if ( dLL[m-1] > 0-esp ) {
          optlogdelta <- append(optlogdelta, ulim)
          #optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
          optLL <- append(optLL, FASTmrEMMA.delta.REML.LL.c(ulim,X,M,M.y,yMy,v))

        }
        for( i in 1:(m-1) )
        {
          if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
          {
            #r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
            r <- uniroot(FASTmrEMMA.delta.REML.dLL.c,lower = logdelta[i],upper = logdelta[i+1],X=X,M=M,M.y=M.y,yMy=yMy,v=v)

            optlogdelta <- append(optlogdelta, r$root)
            #optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
            optLL <- append(optLL,FASTmrEMMA.delta.REML.LL.c(r$root,X,M,M.y,yMy,v))
          }
        }
        maxdelta <- exp(optlogdelta[which.max(optLL)])
        optLL=replaceNaN(optLL)
        maxLL <- max(optLL)
        xMy<-crossprod(X,M.y)
        xMx<-crossprod(X,(M%*%X))
        maxve <-(yMy-maxdelta*(xMy)^2/(1+maxdelta*xMx))/v
        #(sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
        maxvg <- maxve*maxdelta
        #alpha<-inv()

        return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,delta1=xMx,xMy=xMy))
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
        return(list(ML=LL,M=M_c,n=n))

      }



      emma.REMLE0.c <- function(Y_c,W_c){#20190831

        n <- length(Y_c)
        stopifnot(nrow(W_c)==n)

        t <-qr(W_c)$rank
        v <-n-t

        M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
        etas<-crossprod(M_c,Y_c)

        LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
        return(list(REML=LL,M=M_c,v=v))

      }

      replaceNaN<-  function(LL) {
        index=(LL=="NaN")
        if(length(index)>0) theMin=min(LL[!index])
        if(length(index)<1) theMin="NaN"
        LL[index]=theMin
        return(LL)
      }


      maf.fun<-function(snp){
        leng<-length(snp)
        id.1<-length(which(snp==1))
        id.0<-length(which(snp==0))
        id.0.5<-length(which(snp==0.5))
        maf.1<-id.1/leng
        maf.0.5<-id.0.5/leng
        maf.0<-id.0/leng
        ma1<-(2*id.1+id.0.5)/(2*leng)
        ma2<-(2*id.0+id.0.5)/(2*leng)
        maf.min<-min(ma1,ma2)
        return(list(maf.1,maf.0,maf.0.5,maf.min))
      }


      pve.fun<-function(beta,maf){
        pve<-(maf$p1-maf$p1^2+0.25*maf$p3-0.25*maf$p3^2-maf$p1*maf$p3)*beta^2
        return(pve)
      }

      yraw<-matrix(phe[,1],,1)
      xnames<-gen[,1:2]
      snp1<-gen[,3:ncol(gen)]
      mydata<-t(matrix(snp1,nrow=dim(snp1)[1]))
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

      if(flagREMLE==1){
        remle1<-emma.REMLE(Y, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
      }else{
        remle1<-emma.MLE(Y, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)

      }

      remle1.deltahat.g<-remle1$delta
      remle1.B1<-emma.maineffects.B(Z=NULL,K,remle1.deltahat.g)
      C2<-remle1.B1$mC

      rm(remle1.B1)
      gc()

      if(flagREMLE==1){

        ys=Y;xs=mydata;Z=C2;X0=W;ngrids=100;llim=-10;ulim=10;esp=1e-10
        ys <- Z%*%ys
        xs <- Z%*%xs
        X0 <- Z%*%X0

        ys<-as.matrix(ys)
        xs<-as.matrix(xs)
        X0<-as.matrix(X0)

        n <- nrow(ys)
        t <- ncol(xs)
        q<- if ( is.matrix(X0) ) ncol(X0) else 1
        v<-n-q

        MLE0<-emma.REMLE0.c(ys,X0)
        ML1s <- vector(length=t)
        ML0s <- vector(length=t)
        vgs <- vector(length=t)
        ves <- vector(length=t)
        lambdas <- vector(length=t)
        bhats<-vector(length=t)
        d <- vector(length=t)

        stats <- vector(length=t)
        ps <- vector(length=t)
        M<-MLE0$M
        M.y<-M%*%ys
        yMy<-crossprod(ys,M.y)


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

        REML.LRT.c2<-foreach(i=1:t,.combine = 'rbind')%dopar%{

          #MLE1 <- emma.REMLE.c (ys, x0v, K=1, xv, qr.X0,ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)#20181112
          MLE1 <- FASTmrEMMA.REMLE.c(X=xs[,i],M,M.y,yMy,v, ngrids=100, llim=-10, ulim=10, esp=1e-10)

          if(is.na(MLE1$REML)==TRUE){
            ps[i]<-1
          }else{
            ML1s[i]<-MLE1$REML
            ML0s[i]<-MLE0$REML
            vgs[i]<-MLE1$vg
            ves[i]<-MLE1$ve
            lambdas[i] <- MLE1$delta
            ###################
            d[i] <- 1/(1+MLE1$delta*MLE1$delta1)
            #bhats[i]<-MLE1$lambda*MLE1$xMy/(1+MLE1$lambda*MLE1$delta1)
            #bhats[i]<-MLE1$delta*MLE1$xMy/d[i]
            bhats[i]<-MLE1$delta*MLE1$xMy*d[i]
            #to record me=sum(d)

            stats[i]<- 2*(MLE1$REML-MLE0$REML)
            ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2
          }

          c(ps[i],bhats[i],lambdas[i],d[i],ML1s[i],ML0s[i],stats[i],vgs[i],ves[i])
        }
        stopCluster(cl)
      }else{

        #FASTmrEMMA.ML.LRT.c <- function(ys, xs, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
        #20190910
        #Z=C,X0=W=W0,xs=x:snp,n*p
        ys=Y;xs=mydata;Z=C2;X0=W;ngrids=100;llim=-10;ulim=10;esp=1e-10
        ys <- Z%*%ys
        xs <- Z%*%xs
        X0 <- Z%*%X0
        ys<-as.matrix(ys)
        xs<-as.matrix(xs)
        X0<-as.matrix(X0)
        n <- nrow(ys)
        t <- ncol(xs)
        q<- if ( is.matrix(X0) ) ncol(X0) else 1
        v<-n-q
        MLE0<-emma.MLE0.c(ys,X0)
        ML1s <- vector(length=t)
        ML0s <- vector(length=t)
        vgs <- vector(length=t)
        ves <- vector(length=t)
        lambdas<-vector(length=t)
        bhats<-vector(length=t)
        #
        d <- vector(length=t)
        stats <- vector(length=t)
        ps <- vector(length=t)
        #n<-199
        #M<-diag(1,n)-X0%*%ginv(crossprod(X0))%*%t(X0)
        M<-MLE0$M
        M.y<-M%*%ys
        yMy<-crossprod(ys,M.y)

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

        REML.LRT.c2<-foreach(i=1:t,.combine = 'rbind')%dopar%{

          MLE1 <- FASTmrEMMA.MLE.c(X=xs[,i],M,M.y,yMy,n, ngrids=100, llim=-10, ulim=10, esp=1e-10)
          if(length(MLE1$vg)!=0){
            ML1s[i]<-MLE1$ML
            ML0s[i]<-MLE0$ML
            vgs[i]<-MLE1$vg
            ves[i]<-MLE1$ve
            lambdas[i]<-MLE1$delta
            ###################
            d[i] <- 1/(1+MLE1$delta*MLE1$delta1)
            #bhats[i]<-MLE1$lambda*MLE1$xMy/(1+MLE1$lambda*MLE1$delta1)
            #bhats[i]<-MLE1$delta*MLE1$xMy/d[i]
            bhats[i]<-MLE1$delta*MLE1$xMy*d[i]
            #to record me=sum(d)
            stats[i]<- 2*(MLE1$ML-MLE0$ML)
            ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2#20160619
          }else{
            ps[i]<-1
          }
          c(ps[i],bhats[i],lambdas[i],d[i],ML1s[i],ML0s[i],stats[i],vgs[i],ves[i])
        }
        stopCluster(cl)
      }

      rm(Z,xs)
      gc()

      REML.LRT.c2.new<-data.frame(REML.LRT.c2)

      rm(C2,mydata)
      gc()

      parms<-data.frame(chr.locus=xnames,REML.LRT.c2.new)
      names(parms)<-NULL
      parms<-as.matrix(parms)


      parmeter<-parms[,1:4]
      parmeter[,3]<--log10(parmeter[,3])
      parmeter[which(abs(parmeter)>1e-4)]<-round(parmeter[which(abs(parmeter)>1e-4)],4)
      parmeter[which( abs(parmeter)<1e-4)]<-as.numeric(sprintf("%.4e", parmeter[which( abs(parmeter)<1e-4)]))


      if(inputform==1){
        parmsShow<-cbind(genRaw[-1,1],parmeter,genRaw[-1,4])
        parmsShow<-parmsShow[,c(1,2,3,5,4,6)]
        colnames(parmsShow)<-c("Marker","Chromosome","Marker position (bp)","SNP effect (FASTmrEMMA)","'-log10(P) (FASTmrEMMA)'","Genotype for code 1")

      }


      ###### ###### ###### ###### ###### ###### ###### ######
      ###### ###### ###### ###### ###### ###### ###### ######

      if(fin_block==FALSE){
        wan <- NULL
        output<-list(result1=parmsShow,result2=wan,result3=REML.LRT.c2.new)
        return(output)
      }else{


        REML.LRT.c2.new <- rbind(read_ll,REML.LRT.c2.new)
        genRaw <- rbind(read_genRaw,genRaw[-1,])

        if(length(c(which(genRaw[,1]=="rs#")))!=1){
          genRaw <- genRaw[-c(which(genRaw[,1]=="rs#")[-1]),]
        }

        xnames <- gen_bim[,c(1,4)]
        parms<-data.frame(chr.locus=xnames,REML.LRT.c2.new)
        names(parms)<-NULL
        parms<-as.matrix(parms)


        parmeter<-parms[,1:4]
        parmeter[,3]<--log10(parmeter[,3])
        parmeter[which(abs(parmeter)>1e-4)]<-round(parmeter[which(abs(parmeter)>1e-4)],4)
        parmeter[which( abs(parmeter)<1e-4)]<-as.numeric(sprintf("%.4e", parmeter[which( abs(parmeter)<1e-4)]))


        if(inputform==1){
          parmsShow<-cbind(genRaw[-1,1],parmeter,genRaw[-1,4])
          parmsShow<-parmsShow[,c(1,2,3,5,4,6)]
          colnames(parmsShow)<-c("Marker","Chromosome","Marker position (bp)","SNP effect (FASTmrEMMA)","'-log10(P) (FASTmrEMMA)'","Genotype for code 1")

        }

        Xemma<-data.frame(chr.locus=xnames,REML.LRT.c2.new)
        vid<-which(as.numeric(Xemma[,3])<=svpal)

        if(length(vid)!=0){
          if(length(vid)==1){
            xname.emma.opt<-matrix(gen_bim[vid,c(1,4)],1,)
            xdata<-t(matrix(t(gen_bed[match_gen_ID_idex,vid]-1),1,))
            xdata<-matrix(xdata,,1)
          }else{
            xname.emma.opt<-gen_bim[vid,c(1,4)]
            xdata<-t(as.matrix(t(gen_bed[match_gen_ID_idex,vid]-1)))
          }

          ydata<-Y
          u1<-ebayes_EM(x=W,z=xdata,y=ydata)
          emma.lod<-likelihood(xxn=W,xxx=xdata,yn=ydata,bbo=u1$u)
          idslod<-which(emma.lod>=svmlod)

          if(length(idslod)!=0){
            if(length(idslod)==1){
              chrlocus<-matrix(xname.emma.opt[idslod,],1,)
            }else{
              chrlocus<-as.matrix(xname.emma.opt[idslod,])
            }

            gc()

            maf.snp.2<-matrix(t(unlist(apply((gen_bed[match_gen_ID_idex,vid][,idslod]-1),2,maf.fun))),nrow = 4)
            maf.snp.3<-t(maf.snp.2)
            maf.snp.4<-data.frame(maf.snp.3)
            names(maf.snp.4)<-c("p1","p2","p3","maf")

            pve.opt.all.1<-pve.fun(u1$u[idslod],maf.snp.4)
            pve.opt.all<-pve.opt.all.1/as.vector(max(var(Y),sum(pve.opt.all.1)+u1$sigma2))*100

            qtneffect<-matrix(u1$u[idslod],,1)
            lodscore<-matrix(emma.lod[idslod],,1)
            log10P <- as.matrix(-log10(pchisq(lodscore*4.605,1,lower.tail = F)))
            maff<-matrix(maf.snp.4$maf,,1)
            r2<-matrix(pve.opt.all,,1)
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
            colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen(total)")
            wan<-as.data.frame(wan)

          }
        }
        output<-list(result1=parmsShow,result2=wan)
        return(output)

      }
    }
  }
  FASTmrMLM_2.0<-function(gen,phe,genRaw,kk,fin_block = FALSE,read_ll=NULL,read_genRaw=NULL,block_i,genq_BED=NULL,match_gen_ID_idex=NULL){

    inputform<-Genformat
    svlod<-svmlod


    if(is.null(psmatrix)){
      flagps<-1
    }else{
      flagps<-0
    }

    if(is.null(svpal)==TRUE||is.null(svrad)==TRUE||is.null(svlod)==TRUE){
      warning("Please set parameter!")
    }

    if((svpal<0)||(svpal>1))
    {
      warning("Please input critical P-value between 0 and 1!")
    }
    if(svrad<0)
    {
      warning("Please input search radius (kb) of candidate gene: > 0 !")
    }
    if(svlod<0)
    {
      warning("Please input critical LOD score: > 0 !")
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
      warning("Sample size in genotypic dataset doesn't equal to the sample size in phenotypic dataset!")
    }

    if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2)))&&(svpal>=0)&&(svpal<=1)&&(svrad>0)&&(svmlod>=0))
    {

      parmsShow<-NULL
      wan<-NULL
      parms<-NULL
      parms.pchange<-NULL

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
          p<-pchisq(f,1,lower.tail = F)
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

      mixed1<-function(xu,yu,theta1){

        loglike<-function(theta1){
          lambda<-exp(theta1)
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
          loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))-0.5*(n-q)
          return(-loglike)
        }
        grad<-function(theta1){
          lambda<-exp(theta1)
          h<-1/(lambda*delta+1)
          d<-diag(delta,nrow(X1),nrow(X1))
          hinv<-diag(1/(lambda*delta+1),nrow(X1),nrow(X1))
          yy<-sum(yu*h*yu)
          yx<-matrix(0,q,1)
          xx<-matrix(0,q,q)
          for(i in 1:q){
            yx[i]<-sum(yu*h*xu[,i])
            for(j in 1:q){
              xx[i,j]<-sum(xu[,i]*h*xu[,j])
            }
          }
          pp=hinv- hinv%*%xu%*%solve(xx)%*%t(xu)%*%hinv
          sigma<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
          f= -0.5*{sum(diag(pp%*%d))-1/sigma*(t(yu)%*%pp%*%d%*%pp%*%yu)}
          return(c(-f))
        }
        parm<-optim(par=theta,fn=loglike,gr=grad,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
        lambda<-(parm$par)
        return(c(lambda))
      }

      lll<- function(theta){
        lambdak<-exp(theta)
        deth<-1+lambdak*g1
        tmp<-lambdak*1/deth
        yHy<-yy-zy%*%tmp%*%zy
        yHx<-yx-zx%*%tmp%*%zy
        xHx<-xx-zx%*%tmp%*%t(zx)
        logdt2<-log(deth)
        ll<- -0.5*logdt2-0.5*(n-q)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
        return(-ll)
      }
      grad2<- function(theta){
        lambdak<-exp(theta)
        deth<-1+lambdak*g1
        tmp<-lambdak*1/deth
        yHy<-yy-zy%*%tmp%*%zy
        yHx<-yx-zx%*%tmp%*%zy
        xHx<-xx-zx%*%tmp%*%t(zx)
        zHy<-zy-zz%*%tmp%*%zy
        zHx<-zx-zx%*%tmp%*%zz
        zHz<-zz-zz%*%tmp%*%zz
        sigma2<-(yHy-t(yHx)%*%solve(xHx)%*%yHx)/(n-q)
        f<- -0.5*{(zHz-t(zHx)%*%solve(xHx)%*%zHx)-(zHy-t(zHx)%*%solve(xHx)%*%yHx)^2/sigma2}
        return(c(-f))
      }


      fixed2<-function(lambdak){
        deth<-1+lambdak*g1
        tmp<-lambdak*1/deth
        yHy<-yy-zy%*%tmp%*%zy
        yHx<-yx-zx%*%tmp%*%zy
        xHx<-xx-zx%*%tmp%*%t(zx)
        zHy<-zy-zz%*%tmp%*%zy
        zHx<-zx-zx%*%tmp%*%zz
        zHz<-zz-zz%*%tmp%*%zz
        beta<-solve(xHx,yHx)
        tmp2<-solve(xHx)
        sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-q)
        gamma<-lambdak*zHy-lambdak*t(zHx)%*%tmp2%*%yHx
        var<-abs((lambdak*diag(1)-lambdak*zHz*lambdak)*as.numeric(sigma2))
        wald<-gamma^2/var
        stderr<-sqrt(diag(var))
        p_value<-pchisq(wald,1,lower.tail = F)
        result<-list(gamma,stderr,beta,sigma2,p_value,wald)
        return(result)
      }

      y<-as.matrix(phe)
      XX1<-t(gen)

      x<-XX1[3:nrow(XX1),]
      rownames(x)<-NULL
      colnames(x)<-NULL

      X1<-as.matrix(x)
      rm(x)
      gc()
      n<-nrow(X1)
      m<-ncol(X1)
      ########kinship##########
      xxx<-matrix(1,n,1)
      xxx<-matrix()
      if (is.null(psmatrix)==TRUE)
      {
        xxx<-matrix(1,n,1)
      }else{
        ps<-as.matrix(psmatrix)
        xxx<-cbind(matrix(1,n,1),ps)
      }

      qq<-eigen(kk)
      delta<-qq[[1]]
      d<-diag(delta,n,n)
      uu<-qq[[2]]
      q<-ncol(xxx)
      waving<-svrad
      xu<-t(uu)%*%xxx
      zkk<-t(uu)%*%X1
      theta1<-0
      theta<-0

      rm(kk,d,qq)
      gc()

      ll<-numeric()
      y<-as.matrix(y)
      yu<-t(uu)%*%y
      ll<-numeric()
      omeg<-mixed1(xu,yu,theta1)
      delta1<-1/sqrt(delta*exp(omeg)+1)
      d1<-diag(delta1,nrow(X1),nrow(X1))
      yc<-d1%*%yu
      yy<-sum(yc*1*yc)
      xc<-d1%*%xu
      yx<-matrix(0,q,1)
      for(i in 1:q){
        yx[i]<-sum(yc*1*xc[,i])
      }
      binv<-diag(1,nrow(X1),nrow(X1))
      xx<-matrix(0,q,q)
      for(i in 1:q){
        for(j in 1:q){
          xx[i,j]<-sum(xc[,i]*1*xc[,j])
        }
      }
      zkk1<-d1%*%zkk

      rm(d1,uu)
      gc()

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
      mat=foreach(j=1:m, .multicombine=TRUE, .combine = 'rbind')%dopar%
        {

          zc<-as.matrix(zkk1[,j])
          uu1<-as.matrix(zc)%*%t(as.matrix(zc))
          g1<-sum(diag(uu1))
          zy<-as.matrix(sum(yc*1*zc))
          zz<-as.matrix(sum(zc*1*zc))
          zx<-matrix(0,q,1)
          for(i in 1:q){
            zx[i]<-sum(xc[,i]*1*zc)
          }
          par<-tryCatch(optim(par=theta,fn=lll,hessian = TRUE,gr=grad2,method="L-BFGS-B",lower=-10,upper=10), error=function(e) optim(par=theta,fn=lll,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10))
          lambda<-exp(par$par)
          conv<-par$convergence
          fn1<-par$value
          hess<-par$hessian
          parmfix<-fixed2(lambda)
          gamma<-parmfix[[1]]
          stderr<-parmfix[[2]]
          beta<-parmfix[[3]][1,]
          sigma2<-parmfix[[4]]
          p_wald<-parmfix[[5]]
          sigma2g<-lambda*sigma2
          wald<-parmfix[[6]]
          fn0<-lll(c(-Inf))
          lrt<-2*abs(fn0-fn1)
          p_lrt<-pchisq(lrt,1,lower.tail = F)
          parm0<-c(j,beta,sigma2,sigma2g,gamma,stderr,wald,p_wald)
        }
      stopCluster(cl)

      rm(zkk,zkk1)
      gc()

      ll<-rbind(ll,mat)
      parms1<-as.matrix(ll)
      rownames(parms1)<-NULL
      newparm<-cbind(gen[,1:2],parms1[,2:8])
      parms<-newparm
      parms.pchange<-parms
      parmsp<-as.matrix(parms.pchange[,9])
      locsub<-which(parmsp==0)
      if(length(locsub)!=0){
        pmin<-min(parmsp[parmsp!=0])
        subvalue<-10^(1.1*log10(pmin))
        parms.pchange[locsub,9]<-subvalue
      }else{
        parms.pchange<-parms
      }

      if(inputform==1){
        #output result1 using mrMLM numeric format
        parmsShow<-parms
        tempparms<-parms[,3:9]
        tempparms[,7]<--log10(tempparms[,7])
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        parmsShow<-cbind(genRaw[-1,1],parms[,1:2],tempparms,genRaw[-1,4])
        colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (FASTmrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (FASTmrMLM)'","Genotype for code 1")

      }




      ###### ###### ###### ###### ###### ###### ###### ######
      ###### ###### ###### ###### ###### ###### ###### ######

      if(fin_block==FALSE){
        if(is.null(parmsShow)==FALSE){
          parmsShow<-parmsShow[,-c(4,5,6,8,9,12)]
        }
        wan <- NULL
        output<-list(result1=parmsShow,result2=wan,result3=ll)
        return(output)
      }else{



        ll <- rbind(read_ll,ll)
        genRaw <- rbind(read_genRaw,genRaw[-1,])

        if(length(c(which(genRaw[,1]=="rs#")))!=1){
          genRaw <- genRaw[-c(which(genRaw[,1]=="rs#")[-1]),]
        }



        parms1<-as.matrix(ll)
        rownames(parms1)<-NULL
        newparm<-as.matrix(cbind(as.matrix(gen_bim[,c(1,4)]),parms1[,2:8]))
        parms<-newparm
        parms.pchange<-parms
        parmsp<-as.matrix(parms.pchange[,9])
        locsub<-which(parmsp==0)
        if(length(locsub)!=0){
          pmin<-min(parmsp[parmsp!=0])
          subvalue<-10^(1.1*log10(pmin))
          parms.pchange[locsub,9]<-subvalue
        }else{
          parms.pchange<-parms
        }

        if(inputform==1){
          #output result1 using mrMLM numeric format
          parmsShow<-parms
          tempparms<-parms[,3:9]
          tempparms[,7]<--log10(tempparms[,7])
          tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
          tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
          parmsShow<-cbind(genRaw[-1,1],parms[,1:2],tempparms,genRaw[-1,4])
          colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (FASTmrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (FASTmrMLM)'","Genotype for code 1")

        }



        ####### ######## ######


        p<-as.vector(parms1[,8])
        ans<-p.adjust(p, method = "bonferroni", n = length(p))

        rm(gen)
        gc()

        rm(XX1)
        XX1 <- t(gen_bim[,c(1,4)])

        ##########p is parameter########
        sigg<-as.vector(which(p<=svpal))
        le1<-length(sigg)

        if(le1!=0){

          if (length(which(ans<0.05))!=0)
          {
            siggbh<-which(ans<0.05)
            nnn1<-cbind(XX1[1,],XX1[2,])
            setloci<-siggbh
            setposi<-c(XX1[2,siggbh])
            num<-dim(nnn1)[1]
            endresult<-numeric()
            for (t in 1:length(siggbh))
            {
              for (i in 1:num){
                temp<-numeric()
                if ((XX1[1,i]==XX1[1,(setloci[t])])&&(abs(nnn1[i,2]-setposi[t])<=waving))
                {
                  temp<-cbind(matrix(nnn1[i,],1,),i)
                  endresult<-rbind(endresult,temp)
                }
              }
            }
            end<-as.vector(endresult[,3])
            sigg2<-sigg[!sigg%in% end]
            sigg1<-sort(c(siggbh,sigg2))
          }else{
            sigg1<-sigg
          }
          if (length(sigg1)>nrow(X1))
          {
            larsres<-lars((genq_BED[match_gen_ID_idex,sigg1]-1), y, type = "lar",trace = FALSE, normalize = TRUE, intercept = TRUE, eps = .Machine$double.eps, use.Gram=FALSE)
            larsc2<-sigg1[which(larsres$entry!=0)]
            if(length(which(larsres$entry>nrow(X1)))!=0)
            {
              ad1<-sigg1[which(larsres$entry>nrow(X1))]
              larsc<-larsc2[!larsc2%in%ad1]
            }else{
              larsc<-larsc2
            }
          }else{
            larsc<-sigg1
          }

          z<-matrix(1,nrow(X1),1)
          z<-matrix()

          if (is.null(psmatrix)==TRUE)
          {
            z<-matrix(1,nrow(X1),1)
          }else{
            z<-cbind(matrix(1,nrow(X1),1),psmatrix)
          }
          le1<-length(larsc)
          xxxnew11<-as.matrix(genq_BED[match_gen_ID_idex,larsc]-1)

          u1<-ebayes_EM(z,xxxnew11,y)
          obj<-u1$u
          result1<-matrix(0,ncol(gen_bed),1)
          for (i in 1:le1)
          {
            result1[(larsc)[i],1]=obj[i]
          }
          Res<- t(as.matrix((rowSums(result1)/ncol(result1))))
          Res1<-as.vector(Res)
          sig1<-which(abs(Res1)>=1e-5)
          le2<-length(which(abs(Res1)>=1e-5))

          if(le2!=0){
            bbo<-matrix(0,le2,1)
            for (i in 1:le2){
              bbo[i,]=Res1[sig1[i]]
            }
            xxxx<-as.matrix(genq_BED[match_gen_ID_idex,sig1]-1)
            yn<-as.matrix(y)
            xxn<-z
            lod<-likelihood(xxn,xxxx,yn,bbo)

            her1<-vector(length=le2)
            for (i in 1:le2){
              p1<-length(as.vector(which((genq_BED[match_gen_ID_idex,sig1[i]]-1)==1)))/length(genq_BED[match_gen_ID_idex,sig1[i]])
              p2<-1-p1
              her1[i]=((p1+p2)-(p1-p2)^2)*(Res1[sig1[i]])^2
            }

            if(var(y)>=sum(her1)+u1$sigma2){
              her<-(her1/as.vector(var(y)))*100

            }else{
              her<-(her1/(sum(her1)+u1$sigma2))*100
            }

            slod<-cbind(sig1,lod,her)

            if(length(which(slod[,2]>=svlod))>=1){

              if(length(which(slod[,2]>=svlod))==1){
                sslod<-t(as.matrix(slod[which(slod[,2]>=svlod),]))
                sig1<-slod[which(slod[,2]>=svlod),1]
              }else if(length(which(slod[,2]>=svlod))>1){
                sslod<-slod[which(slod[,2]>=svlod),]
                sig1<-sslod[,1]
              }
              xxxx<-as.matrix(genq_BED[match_gen_ID_idex,sig1]-1)
              lod<-sslod[,2]
              her<-sslod[,3]

              ii<-as.vector(sig1)
              qqq<-matrix(0,nrow=length(ii),ncol=6)
              qqq[,1]=as.matrix(ii)
              for (j in 1:length(ii)){
                qqq[j,2]=XX1[1,ii[j]]
                qqq[j,3]=XX1[2,ii[j]]
                qqq[j,4]=result1[ii[j],]

                qqq[j,5]=lod[j]
                qqq[j,6]=her[j]
              }

              rm(XX1,X1,genq_BED)
              gc()

              id<-which(qqq[,5]==0)

              if(length(id)!=dim(qqq)[1]){

                if(length(id)!=0){
                  qqq1<-qqq[-id,]
                }else{
                  qqq1<-qqq
                }
                xxmaf<-t(xxxx)
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
                  result<-as.matrix(qqq1[,-1])
                  vees<-matrix("",nrow = nrow(result),1)
                  pees<-matrix("",nrow = nrow(result),1)
                  pees[1,1]<-pee
                  vees[1,1]<-vee

                }else{
                  result<-t(as.matrix(qqq1[,-1]))
                  pees<-as.matrix(pee)
                  vees<-as.matrix(vee)
                }


                if(nrow(qqq1)>1){
                  result<-as.matrix(qqq1[,-1])
                  result<-result
                  temp<-as.matrix(result[,3:5])
                  temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
                  temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
                  wan<-cbind(result[,1:2],temp)
                  snp<-parmsShow[,11]

                }else{
                  result<-t(as.matrix(qqq1[,-1]))
                  result<-result
                  temp<-t(as.matrix(result[,3:5]))
                  temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
                  temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
                  wan<-cbind(t(as.matrix(result[,1:2])),temp)
                  snp<-parmsShow[,11]

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
                log10P <- as.matrix(round(-log10(pchisq(lodscore1*4.605,1,lower.tail = F)),4))
                if(nrow(tempwan)>1){
                  tempwan1 <- cbind(tempwan[,1:5],log10P,tempwan[,6:10])
                }else{
                  tempwan1 <- cbind(t(as.matrix(tempwan[,1:5])),log10P,t(as.matrix(tempwan[,6:10])))
                }
                wan <- tempwan1
                colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen (total)")
                wan<-as.data.frame(wan)
              }
            }
          }
        }
        if(is.null(parmsShow)==FALSE){
          parmsShow<-parmsShow[,-c(4,5,6,8,9,12)]
        }
        output<-list(result1=parmsShow,result2=wan)
        return(output)



      }






    }
  }
  pKWmEB_2.0<-function(gen,phe,genRaw,kk,fin_block = FALSE,read_ll=NULL,read_genRaw=NULL,block_i,match_gen_ID_idex=NULL){

    inputform<-Genformat
    pheRAW<-phe

    if(is.null(kk)){
      if(is.null(gen)==TRUE)
      {
        warning("Please input correct genotypic dataset !")
      }else{
        envgen<-gen[,3:ncol(gen)]
        envgen<-t(envgen)
        m<-ncol(envgen)
        n<-nrow(envgen)
        #kk1<-matrix(0,n,n)
        # for(k in 1:m){
        #   z<-as.matrix(envgen[,k])
        #   kk1<-kk1+z%*%t(z)
        # }
        kk1<-mrMLM::multiplication_speed(envgen,t(envgen))
        cc<-mean(diag(kk1))
        kk1<-kk1/cc
        kk<-as.matrix(kk1)
      }
      rm(envgen,kk1)
      gc()

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
          p<-pchisq(f,1,lower.tail = F)
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

      parmsShow<-NULL
      wan<-NULL
      parms.pchange<-NULL
      parmsm<-NULL

      K.data <- kk
      Y.data <- phe
      rawgen <- gen
      rawphe <- Y.data

      gene.data<-rawgen[,3:ncol(rawgen)]
      nsample <- ncol(gene.data)

      fix <- matrix(1,nsample,1)
      sam <- nsample
      Y.data <- matrix(Y.data,nsample,1)
      n<-dim(Y.data)[1]
      W.orig<-matrix(1,n,1)
      W <- W.orig
      K <- K.data
      YY <- Y.data
      rm(K.data)
      gc()
      p_value <- svpal
      ffpptotal <- numeric()
      gglartotal <- numeric()
      pvaluetotal <- numeric()

      #for(ii in 1:1){
      ii<-1
      remle2<-emma.REMLE(YY[,ii], W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
      remle1.B1<-emma.maineffects.B(Z=NULL,K,remle2$delta)
      C2<-remle1.B1$mC

      rm(K,remle1.B1)
      gc()

      Y_c <- C2%*%YY[,ii]
      W_c <- C2%*%W
      G_c <- C2%*%t(gene.data)

      GGG <- t(G_c)

      rm(C2,G_c)
      gc()

      allrowmean <- rowMeans(GGG)
      nnG <- nrow(GGG)

      for(jj in 1:nnG)
      {
        GGG[jj,which(GGG[jj,]>=allrowmean[jj])] <- 1
        GGG[jj,which(GGG[jj,]<allrowmean[jj])] <- -1
      }

      gentran <- GGG

      rm(GGG)
      gc()

      phetran <- Y_c
      nn <- dim(gentran)[1]
      bb<-numeric()
      cc <- numeric()
      ff <- numeric()

      newphe <- cbind(matrix(c(1:sam),,1),phetran)
      ph <- unique(newphe[,2])
      newph <- newphe[match(ph,newphe[,2],0L),]
      newy <- newph[,2]
      sob <- newph[,1]

      ff<- foreach(i=1:nn)%do%
        {
          temp <- as.matrix(gentran[i,sob])
          temp<-factor(temp)
          loc<-which(as.numeric(levels(temp))==1)

        }

      fff<-unlist(ff)
      sameloc<-which(fff==1)

      if(length(sameloc)!=0){
        gentran1<-gentran[-c(sameloc),]
      }else if(length(sameloc)==0){
        gentran1<-gentran
      }

      nnn<-dim(gentran1)[1]

      rm(gentran)
      gc()

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

      unsameloc=foreach(i=1:nnn, .combine = 'rbind')%dopar%
        {
          requireNamespace("coin")
          requireNamespace("lars")
          temp <- as.matrix(gentran1[i,sob])
          xy <- cbind(temp,newy)
          b <- unique(xy[,1])

          temp <- factor(temp)
          snp <- data.frame(newy,temp)
          kw <- kruskal_test(newy~temp, data = snp,distribution = "asymptotic")
          kw <- pvalue(kw)
          aa <- kw[1]
        }

      stopCluster(cl)
      a<-matrix(0,nrow = nn,ncol=1)
      a[c(sameloc)]<-1
      a[which(a[]==0)]<-unsameloc
      bb<-a

      rm(gentran1)
      gc()

      kk <- matrix(seq(1:nn),nn,1)
      bb <- matrix(bb,nn,1)
      cc <- cbind(ii,kk,bb)
      pvaluetotal <- cc[,2:3]
      ff <- cc[which(cc[,3] < p_value),]
      ffpptotal <- ff
      pvaluetotal <- pvaluetotal


      if(inputform==1){
        #output result1 using mrMLM numeric format
        parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
        tempparms<-parmsShow
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        kong<-matrix("",nrow(tempparms),1)
        parmsShow<-data.frame(genRaw[-1,1],gen[,1:2],kong,tempparms,genRaw[-1,4])
        colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","SNP effect (pKWmEB)","'-log10(P) (pKWmEB)'","Genotype for code 1")

      }


      ###### ###### ###### ###### ###### ###### ###### ######
      ###### ###### ###### ###### ###### ###### ###### ######

      if(fin_block==FALSE){
        wan <- NULL
        output<-list(result1=parmsShow,result2=wan,result3=cc)
        return(output)
      }else{




        cc <- rbind(read_ll,cc)
        cc[,2] <- c(1:nrow(cc))
        genRaw <- rbind(read_genRaw,genRaw[-1,])

        if(length(c(which(genRaw[,1]=="rs#")))!=1){
          genRaw <- genRaw[-c(which(genRaw[,1]=="rs#")[-1]),]
        }

        pvaluetotal <- cc[,2:3]
        ff <- cc[which(cc[,3] < p_value),]
        ffpptotal <- ff
        pvaluetotal <- pvaluetotal


        if(inputform==1){
          #output result1 using mrMLM numeric format
          parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
          tempparms<-parmsShow
          tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
          tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
          kong<-matrix("",nrow(tempparms),1)
          parmsShow<-data.frame(genRaw[-1,1],gen_bim[,c(1,4)],kong,tempparms,genRaw[-1,4])
          colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","SNP effect (pKWmEB)","'-log10(P) (pKWmEB)'","Genotype for code 1")

        }

        ############lars###########################
        gg <- numeric()
        nchoice <- ff[,2]
        genchoice <- t(gen_bed[match_gen_ID_idex,nchoice]-1)
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
          optgen <- t(gen_bed[match_gen_ID_idex,c(optloci)]-1)
          newphebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
          bbeff <- ebayes_EM(fix,t(optgen),newphebayes)
          lod <- likelihood(fix,t(optgen),newphebayes,bbeff$u)
          optlod <- which(lod>svmlod)
          if(length(optlod)>0){
            locich <- optloci[optlod]
            ggbayes <- cbind(ii,locich,matrix(gen_bim[locich,c(1,4)],,2),bbeff$u[optlod],lod[optlod],bbeff$sigma2)
          }
          gglartotal <- ggbayes
          rm(rawgen)
          gc()


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
          optgen50 <- t(gen_bed[match_gen_ID_idex,c(optloci50)]-1)
          phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
          bbeff50 <- ebayes_EM(fix,t(optgen50),phebayes)
          lod50 <- likelihood(fix,t(optgen50),phebayes,bbeff50$u)

          optlod50 <- which(lod50>svmlod)
          if(length(optlod50)>0){
            locich50 <- optloci50[c(optlod50)]
            ggbayes50 <- cbind(ii,locich50,gen_bim[locich50,c(1,4)],bbeff50$u[optlod50],lod50[optlod50],bbeff50$sigma2)
            hhbayes50 <- rbind(hhbayes50,as.matrix(ggbayes50))
          }

          ##################choose 100 number variable from lars#####################
          optgen100 <- t(gen_bed[match_gen_ID_idex,c(optloci100)]-1)
          phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
          bbeff100 <- ebayes_EM(fix,t(optgen100),phebayes)
          lod100 <- likelihood(fix,t(optgen100),phebayes,bbeff100$u)

          optlod100 <- which(lod100>svmlod)
          if(length(optlod100)>0){
            locich100 <- optloci100[optlod100]
            ggbayes100 <- cbind(ii,locich100,gen_bim[locich100,c(1,4)],bbeff100$u[optlod100],lod100[optlod100],bbeff100$sigma2)
            hhbayes100 <- rbind(hhbayes100,as.matrix(ggbayes100))
          }

          ##################choose 150 number variable from lars#####################
          optgen150 <- t(gen_bed[match_gen_ID_idex,c(optloci150)]-1)
          phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
          bbeff150 <- ebayes_EM(fix,t(optgen150),phebayes)
          lod150 <- likelihood(fix,t(optgen150),phebayes,bbeff150$u)

          optlod150 <- which(lod150>svmlod)
          if(length(optlod150)>0){
            locich150 <- optloci150[optlod150]
            ggbayes150 <- cbind(ii,locich150,gen_bim[locich150,c(1,4)],bbeff150$u[optlod150],lod150[optlod150],bbeff150$sigma2)
            hhbayes150 <- rbind(hhbayes150,as.matrix(ggbayes150))
          }

          rm(rawgen)
          gc()

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
            xx1 <- as.matrix(t(gen_bed[match_gen_ID_idex,ggbayes50[,2]]-1))
            lmres1 <- lm(phebayes~xx1)
            aic1 <- AIC(lmres1)
          }
          if(length(optlod100)==1)
          {
            xx2 <- as.matrix(t(gen_bed[match_gen_ID_idex,ggbayes100[,2]]-1))
            lmres2 <- lm(phebayes~xx2)
            aic2 <- AIC(lmres2)
          }
          if(length(optlod150)==1)
          {
            xx3 <- as.matrix(t(gen_bed[match_gen_ID_idex,ggbayes150[,2]]-1))
            lmres3 <- lm(phebayes~xx3)
            aic3 <- AIC(lmres3)
          }

          if(length(optlod50)>1)
          {
            xx1 <- gen_bed[match_gen_ID_idex,unlist(ggbayes50[,2])]-1
            lmres1 <- lm(phebayes~xx1)
            aic1 <- AIC(lmres1)
          }
          if(length(optlod100)>1)
          {
            xx2 <- gen_bed[match_gen_ID_idex,unlist(ggbayes100[,2])]-1
            lmres2 <- lm(phebayes~xx2)
            aic2 <- AIC(lmres2)
          }
          if(length(optlod150)>1)
          {
            xx3 <- gen_bed[match_gen_ID_idex,unlist(ggbayes150[,2])]-1
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
        #}

        gglartotal <- gglartotal
        #
        # if(inputform==1){
        #   #output result1 using mrMLM numeric format
        #   parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
        #   tempparms<-parmsShow
        #   tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        #   tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        #   kong<-matrix("",nrow(tempparms),1)
        #   parmsShow<-data.frame(genRaw[-1,1],gen[,1:2],kong,tempparms,genRaw[-1,4])
        #   colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","SNP effect (pKWmEB)","'-log10(P) (pKWmEB)'","Genotype for code 1")
        #
        # }

        finalres <- gglartotal

        if(length(finalres)!=0){

          if(length(finalres[,2])>1){

            if((flagps==1)||(exists("psmatrix")==FALSE))
            {
              ex<-cbind(fix,(gen_bed[match_gen_ID_idex,finalres[,2]]-1))
            }else if(flagps==0)
            {
              ex<-cbind(cbind(fix,psmatrix),(gen_bed[match_gen_ID_idex,finalres[,2]]-1))
            }

          }else{

            if((flagps==1)||(exists("psmatrix")==FALSE))
            {
              ex<-cbind(fix,as.matrix(gen_bed[match_gen_ID_idex,unlist(finalres[,2])]-1))
            }else if(flagps==0)
            {
              ex<-cbind(cbind(fix,psmatrix),as.matrix(t(gen_bed[match_gen_ID_idex,unlist(finalres[,2])]-1)))
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

          gc()

          xxxx<-as.matrix(gen_bed[match_gen_ID_idex,unlist(finalres[,2])][3:nrow(gen_bed),]-1)

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

          eeff <- unlist(finalres[,5])
          lo <- unlist(finalres[,6])
          eeff[which(abs(eeff)>=1e-4)] <- round(eeff[which(abs(eeff)>=1e-4)],4)
          eeff[which(abs(eeff)<1e-4)] <- as.numeric(sprintf("%.4e",eeff[which(abs(eeff)<1e-4)]))
          lo[which(abs(lo)>=1e-4)] <- round(lo[which(abs(lo)>=1e-4)],4)
          lo[which(abs(lo)<1e-4)] <- as.numeric(sprintf("%.4e",lo[which(abs(lo)<1e-4)]))
          her[which(abs(her)>=1e-4)] <- round(her[which(abs(her)>=1e-4)],4)
          her[which(abs(her)<1e-4)] <- as.numeric(sprintf("%.4e",her[which(abs(her)<1e-4)]))
          needrs <- genRaw[-1,1]
          needrs <- as.matrix(needrs[unlist(finalres[,2])])
          needgenofor <- as.character()
          if(inputform==1)
          {
            needgenofor <- genRaw[-1,4]
            needgenofor <- as.matrix(needgenofor[unlist(finalres[,2])])
          }
          if(inputform==2)
          {
            needgenofor <- outATCG
            needgenofor <- as.matrix(needgenofor[unlist(finalres[,2])])
          }
          if(inputform==3)
          {
            needgenofor <- outATCG
            needgenofor <- as.matrix(needgenofor[unlist(finalres[,2])])
          }

          phevartotal<-var(pheRAW)
          if(finalres[1,7]>=1e-4){finalres[1,7]<-round(finalres[1,7],4)}
          if(finalres[1,7]<1e-4){finalres[1,7]<-as.numeric(sprintf("%.4e",finalres[1,7]))}
          if(phevartotal>=1e-4){phevartotal<-round(phevartotal,4)}
          if(phevartotal<1e-4){phevartotal<-as.numeric(sprintf("%.4e",phevartotal))}
          tempvar <- dim(as.matrix(lo))[1]
          if(tempvar==1)
          {
            wan<-data.frame(needrs,t(as.matrix(gen_bim[unlist(finalres[,2]),c(1,4)])),as.matrix(eeff),as.matrix(lo),her,maf,needgenofor,as.matrix(finalres[,7]),phevartotal)
          }else if(tempvar>1)
          {
            wan<-data.frame(needrs,gen_bim[unlist(finalres[,2]),c(1,4)],eeff,lo,her,maf,needgenofor)
            wan<-wan[order(wan[,2]),]
            wan<-data.frame(wan,rbind(finalres[1,7],as.matrix(rep("",(tempvar-1))),use.names=FALSE),rbind(phevartotal,as.matrix(rep("",(tempvar-1)))))
          }

          tempwan <- wan
          lodscore1 <- as.numeric(tempwan[,5])
          log10P <- as.matrix(round(-log10(pchisq(lodscore1*4.605,1,lower.tail = F)),4))
          tempwan1 <- cbind(tempwan[,1:5],log10P,tempwan[,6:10])
          wan <- tempwan1

          colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen(total)")
          wan<-as.data.frame(wan)
        }#change20190125
        output<-list(result1=parmsShow,result2=wan)
        return(output)


      }
    }
  }
  pLARmEB_2.0<-function(phe,match_gen_ID_idex=NULL,CriLOD=NULL){

    lodvalue<-CriLOD
    gene.data<-1
    genRaw <- as.matrix(rbind(t(c("rs#","chrom","pos","genotype for code 1")),gen_bim[,c(2,1,4,5)],use.names=FALSE))

    gc()

    inputform<-Genformat

    if(is.null(psmatrix)){
      flagps<-1
    }else{
      flagps<-0
    }

    if(is.null(lodvalue)==TRUE||is.null(lars1)==TRUE){
      warning("Please set parameter!")
    }
    if(lodvalue<0)
    {
      warning("Please input critical LOD score: > 0 !")
    }
    if(lars1<0||lars1>=nrow(phe))
    {
      warning("Please input the number of most relevant variables select by LARS: >0 and less than numbers of sample!")
    }
    if(is.null(gene.data)==TRUE)
    {
      warning("Please input correct genotypic data !")

    }
    if(is.null(phe)==TRUE)
    {
      warning("Please input correct phenotypic data !")
    }
    if((is.null(gene.data)==FALSE)&&(is.null(phe)==FALSE)&&(nrow(gen_bed)!=(nrow(phe))))
    {
      warning("Sample size in genotypic dataset doesn't equal to the sample size in phenotypic dataset !")
    }

    if((is.null(gene.data)==FALSE)&&(is.null(phe)==FALSE)&&((nrow(gen_bed)==(nrow(phe))))&&(lodvalue>=0)&&(lars1>0))
    {

      wan<-NULL
      result<-NULL

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
          p<-pchisq(f,1,lower.tail = F)
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

      lars <-  function(x, y, type = c("lasso", "lar", "forward.stagewise","stepwise"), trace = FALSE,
                        normalize=TRUE, intercept=TRUE, Gram,
                        eps = .Machine$double.eps,  max.steps, use.Gram = TRUE)
      {

        call <- match.call()
        type <- match.arg(type)
        TYPE <- switch(type,
                       lasso = "LASSO",
                       lar = "LAR",
                       forward.stagewise = "Forward Stagewise",
                       stepwise = "Forward Stepwise")
        if(trace)
          cat(paste(TYPE, "sequence\n"))

        nm <- dim(x)
        n <- nm[1]
        m <- nm[2]
        im <- inactive <- seq(m)
        one <- rep(1, n)
        vn <- dimnames(x)[[2]]
        ### Center x and y, and scale x, and save the means and sds
        if(intercept){
          meanx <- drop(one %*% x)/n
          x <- scale(x, meanx, FALSE)  # centers x
          mu <- mean(y)
          y <- drop(y - mu)
        }
        else {
          meanx <- rep(0,m)
          mu <- 0
          y <- drop(y)
        }
        if(normalize){
          normx <- sqrt(drop(one %*% (x^2)))
          nosignal<-normx/sqrt(n) < eps
          if(any(nosignal))# ignore variables with too small a variance
          {
            ignores<-im[nosignal]
            inactive<-im[-ignores]
            normx[nosignal]<-eps*sqrt(n)
            if(trace)
              cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")  #
          }
          else ignores <- NULL #singularities; augmented later as well
          names(normx) <- NULL
          x <- scale(x, FALSE, normx)	# scales x
        }
        else {
          normx <- rep(1,m)
          ignores <- NULL
        }
        if(use.Gram & missing(Gram)) {
          if(m > 500 && n < m)
            cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n"
            )
          if(trace)
            cat("Computing X'X .....\n")
          Gram <- t(x) %*% x	#Time saving
        }
        Cvec <- drop(t(y) %*% x)
        ssy <- sum(y^2)	### Some initializations
        residuals <- y
        if(missing(max.steps))
          max.steps <- 8*min(m, n-intercept)
        beta <- matrix(0, max.steps + 1, m)	# beta starts at 0
        lambda=double(max.steps)
        Gamrat <- NULL
        arc.length <- NULL
        R2 <- 1
        RSS <- ssy
        first.in <- integer(m)
        active <- NULL	# maintains active set
        actions <- as.list(seq(max.steps))

        drops <- FALSE
        Sign <- NULL
        R <- NULL	###

        k <- 0
        while((k < max.steps) & (length(active) < min(m - length(ignores),n-intercept)) )
        {
          action <- NULL
          C <- Cvec[inactive]	#

          Cmax <- max(abs(C))
          if(Cmax<eps*100){
            if(trace)cat("Max |corr| = 0; exiting...\n")
            break
          }
          k <- k + 1
          lambda[k]=Cmax

          if(!any(drops)) {
            new <- abs(C) >= Cmax - eps
            C <- C[!new]	# for later
            new <- inactive[new]	# Get index numbers

            for(inew in new) {
              if(use.Gram) {
                R <- updateR(Gram[inew, inew], R, drop(Gram[
                  inew, active]), Gram = TRUE,eps=eps)
              }
              else {
                R <- updateR(x[, inew], R, x[, active], Gram
                             = FALSE,eps=eps)
              }
              if(attr(R, "rank") == length(active)) {

                nR <- seq(length(active))
                R <- R[nR, nR, drop = FALSE]
                attr(R, "rank") <- length(active)
                ignores <- c(ignores, inew)
                action <- c(action,  - inew)
                if(trace)
                  cat("LARS Step", k, ":\t Variable", inew,
                      "\tcollinear; dropped for good\n")	#
              }
              else {
                if(first.in[inew] == 0)
                  first.in[inew] <- k
                active <- c(active, inew)
                Sign <- c(Sign, sign(Cvec[inew]))
                action <- c(action, inew)
                if(trace)
                  cat("LARS Step", k, ":\t Variable", inew,
                      "\tadded\n")
              }
            }
          }
          else action <-  - dropid
          Gi1 <- backsolve(R, backsolvet(R, Sign))

          dropouts<-NULL
          if(type == "forward.stagewise") {
            directions <- Gi1 * Sign
            if(!all(directions > 0)) {
              if(use.Gram) {
                nnls.object <- nnls.lars(active, Sign, R,
                                         directions, Gram[active, active], trace =
                                           trace, use.Gram = TRUE,eps=eps)
              }
              else {
                nnls.object <- nnls.lars(active, Sign, R,
                                         directions, x[, active], trace = trace,
                                         use.Gram = FALSE,eps=eps)
              }
              positive <- nnls.object$positive
              dropouts <-active[-positive]
              action <- c(action, -dropouts)
              active <- nnls.object$active
              Sign <- Sign[positive]
              Gi1 <- nnls.object$beta[positive] * Sign
              R <- nnls.object$R
              C <- Cvec[ - c(active, ignores)]
            }
          }
          A <- 1/sqrt(sum(Gi1 * Sign))
          w <- A * Gi1	# note that w has the right signs
          if(!use.Gram) u <- drop(x[, active, drop = FALSE] %*% w)	###

          if( (length(active) >=  min(n-intercept, m - length(ignores) ) )|type=="stepwise") {
            gamhat <- Cmax/A
          }
          else {
            if(use.Gram) {
              a <- drop(w %*% Gram[active,  - c(active,ignores), drop = FALSE])
            }
            else {
              a <- drop(u %*% x[,  - c(active, ignores), drop=FALSE])
            }
            gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))

            gamhat <- min(gam[gam > eps], Cmax/A)
          }
          if(type == "lasso") {
            dropid <- NULL
            b1 <- beta[k, active]	# beta starts at 0
            z1 <-  - b1/w
            zmin <- min(z1[z1 > eps], gamhat)
            if(zmin < gamhat) {
              gamhat <- zmin
              drops <- z1 == zmin
            }
            else drops <- FALSE
          }
          beta[k + 1,  ] <- beta[k,  ]
          beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
          if(use.Gram) {
            Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% w
          }
          else {
            residuals <- residuals - gamhat * u
            Cvec <- drop(t(residuals) %*% x)
          }
          Gamrat <- c(Gamrat, gamhat/(Cmax/A))
          arc.length <- c(arc.length, gamhat)
          if(type == "lasso" && any(drops)) {
            dropid <- seq(drops)[drops]
            for(id in rev(dropid)) {
              if(trace)
                cat("Lasso Step", k+1, ":\t Variable", active[
                  id], "\tdropped\n")
              R <- downdateR(R, id)
            }
            dropid <- active[drops]
            beta[k+1,dropid]<-0
            active <- active[!drops]
            Sign <- Sign[!drops]
          }
          if(!is.null(vn))
            names(action) <- vn[abs(action)]
          actions[[k]] <- action
          inactive <- im[ - c(active, ignores)]
          if(type=="stepwise")Sign=Sign*0
        }
        beta <- beta[seq(k + 1), ,drop=FALSE ]
        lambda=lambda[seq(k)]
        dimnames(beta) <- list(paste(0:k), vn)
        if(trace)
          cat("Computing residuals, RSS etc .....\n")
        residuals <- y - x %*% t(beta)
        beta <- scale(beta, FALSE, normx)
        RSS <- apply(residuals^2, 2, sum)
        R2 <- 1 - RSS/RSS[1]
        actions=actions[seq(k)]
        netdf=sapply(actions,function(x)sum(sign(x)))
        df=cumsum(netdf)### This takes into account drops
        if(intercept)df=c(Intercept=1,df+1)
        else df=c(Null=0,df)
        rss.big=rev(RSS)[1]
        df.big=n-rev(df)[1]
        if(rss.big<eps|df.big<eps)sigma2=NaN
        else
          sigma2=rss.big/df.big
        Cp <- RSS/sigma2 - n + 2 * df
        attr(Cp,"sigma2")=sigma2
        attr(Cp,"n")=n
        object <- list(call = call, type = TYPE, df=df, lambda=lambda,R2 = R2, RSS = RSS, Cp = Cp,
                       actions = actions[seq(k)], entry = first.in, Gamrat = Gamrat,
                       arc.length = arc.length, Gram = if(use.Gram) Gram else NULL,
                       beta = beta, mu = mu, normx = normx, meanx = meanx)
        class(object) <- "lars"
        object
      }

      Y.data<-as.matrix(phe)
      if(is.null(psmatrix)==FALSE){
        psmatrix<-as.matrix(psmatrix)
      }
      nsam <-nrow(gen_bed)
      chrnum<-nrow(unique(gen_bim[,1]))

      W.orig<-matrix(1,nsam,1)
      if(is.null(psmatrix)==FALSE){
        W1 <-cbind(W.orig,psmatrix)
      }else{
        W1<-W.orig
      }

      kk<-list(NULL)
      cc<-list(NULL)
      kktotal<-matrix(0,nsam,nsam)

      for(i in 1:chrnum){
        xot <-gen_bed[match_gen_ID_idex,which(gen_bim==i)]-1
        #kk[[i]]<-mrMLM::multiplication_speed(xot,t(xot))
        kk[[i]]<-xot%*%t(xot)
        cc[[i]]<-mean(diag(kk[[i]]))
        kktotal<-kktotal+kk[[i]]
        rm(xot)
      }
      gc()
      larsres <- numeric(0)
      for(i in 1:chrnum){

        xx1 <- t(gen_bed[match_gen_ID_idex,which(gen_bim[,1]==i)]-1)
        YY1 <- matrix(Y.data,,1)

        K1 <- (kktotal-kk[[i]])/(sum(unlist(cc))-as.numeric(cc[i]))

        repl<-numeric()
        if(Bootstrap==TRUE){

          res1<-foreach(repl=1:5,.multicombine=TRUE,.combine='cbind')%do%{

            if(repl==1){
              YY<-YY1
              xx<-xx1
              K<-K1
              W<-W1
            }else{
              s<-srswr(nrow(YY1),nrow(YY1))
              ind<-(1:nrow(YY1))[s!=0]
              n<-s[s!=0]
              ind<-rep(ind,times=n)
              YY<-as.matrix(YY1[ind,])
              xx<-xx1[,ind]

              K <- K1[ind,ind]
              W<-as.matrix(W1[ind,])
            }

            remle2<-emma.REMLE(YY, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
            remle1.B1<-emma.maineffects.B(Z=NULL,K,remle2$delta)
            rm(K)
            gc()
            C2<-remle1.B1$mC
            Y_c <- C2%*%YY
            W_c <- C2%*%W
            G_c <- C2%*%t(xx)
            GGG <- t(G_c)
            rm(G_c)
            gc()
            ylars <- as.matrix(Y_c)
            xlars <- cbind(W_c,t(GGG))
            rm(GGG)
            gc()
            LAR <- lars(xlars,ylars,type="lar",use.Gram=F,max.steps=lars1)
            rm(xlars)
            gc()
            LAR$beta[nrow(LAR$beta),]
          }

        }else if(Bootstrap==FALSE){
          res1 <- numeric()
          remle2<-emma.REMLE(YY1, W1, K1, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
          remle1.B1<-emma.maineffects.B(Z=NULL,K1,remle2$delta)
          rm(K1)
          gc()
          C2<-remle1.B1$mC
          Y_c <- C2%*%YY1
          W_c <- C2%*%W1
          G_c <- C2%*%t(xx1)
          rm(xx1)
          gc()
          GGG <- t(G_c)
          rm(G_c)
          gc()
          ylars <- as.matrix(Y_c)
          xlars <- cbind(W_c,t(GGG))
          rm(GGG)
          gc()
          LAR <- lars(xlars,ylars,type="lar",use.Gram=F,max.steps=lars1)
          rm(xlars)
          gc()
          res1<-cbind(res1,LAR$beta[nrow(LAR$beta),])
        }

        if(is.null(psmatrix)==FALSE){
          rr <- as.matrix(res1[-c(1:(ncol(psmatrix)+1)),])
        }else{
          rr <- as.matrix(res1[-1,])
        }
        larsres <- rbind(larsres,rr)
      }


      rm(kk,kktotal)
      gc()

      if(Bootstrap==TRUE){
        count <- matrix(rep(0,nrow(larsres)),nrow(larsres),1)

        ttt <- numeric()
        for(ii in 1:nrow(larsres))
        {
          tt <- 0
          for(jj in 1:ncol(larsres))
          {
            if ((abs(larsres[ii,jj]))>0){tt <- tt+1}
          }
          count[ii] <-tt
        }
        larsres <-cbind(larsres,count)

        for(ii in 1:nrow(larsres))
        {
          if(larsres[ii,ncol(larsres)]>=3){ttt <- cbind(ttt,ii)}
        }

        countnum <- ttt

      }else{

        countnum <- numeric()
        for(ii in 1:nrow(larsres))
        {
          if ((abs(larsres[ii]))>0){countnum <- cbind(countnum,ii)}
        }
      }
      if(ncol(countnum)>nrow(phe)){

        if(length(countnum)==1){
          xx2 <- matrix((gen_bed[match_gen_ID_idex,c(countnum)]-1),1,)

        }else{
          xx2 <- as.matrix(t(gen_bed[match_gen_ID_idex,c(countnum)]-1))
        }
        YY2 <- matrix(Y.data,,1)

        ylars <- as.matrix(YY2)
        xlars <- cbind(W1,t(xx2))
        LAR <- lars(xlars,ylars,type="lar",use.Gram=F)

        res1<-as.matrix(LAR$beta[nrow(LAR$beta),])

        rm(xlars,xx2)
        gc()

        if(is.null(psmatrix)==FALSE){
          rr <- as.matrix(res1[-c(1:(ncol(psmatrix)+1)),])
        }else{
          rr <- as.matrix(res1[-1,])
        }

        ct <- numeric()
        for(ii in 1:nrow(rr))
        {
          if ((abs(rr[ii]))>0){ct <- cbind(ct,ii)}
        }

        inct<-c(ct)
        countnum<-countnum[,inct]

      }

      if(length(countnum)==1){
        xeb <- cbind(gen_bim[c(countnum),c(1,4)],matrix((gen_bed[match_gen_ID_idex,c(countnum)]-1),1,))
        ebrow <-matrix(xeb[,1:2],,2)
        xeb1<-matrix(xeb[,3:ncol(xeb)],1,)
        xxeb <- as.matrix(t(xeb1))
        nmak <- ncol(xxeb)

      }else{
        xeb <- cbind(gen_bim[c(countnum),c(1,4)],as.matrix(t(gen_bed[match_gen_ID_idex,c(countnum)]-1)))
        ebrow <-as.matrix(xeb[,1:2])
        xeb1<-as.matrix(xeb)
        xxeb <- as.matrix(t(xeb1[,-c(1:2)]))
        nmak <- ncol(xxeb)
      }
      bayeslodres <- numeric()

      genname<-gen_bim[,c(1,4)]

      rm(xeb,gene.data)
      gc()

      yeb <- as.matrix(phe)

      if(is.null(psmatrix)==FALSE){
        u1<-ebayes_EM(cbind(matrix(1,nrow(xxeb),1),psmatrix),xxeb,yeb)
        xb<-u1$u
      }else{
        u1<-ebayes_EM(matrix(1,nrow(xxeb),1),xxeb,yeb)
        xb<-u1$u
      }
      xb<-as.matrix(xb)
      if(is.null(psmatrix)==FALSE){
        temp<-cbind(matrix(1,nrow(xxeb),1),psmatrix)
      }else{
        temp<-matrix(1,nrow(xxeb),1)
      }

      lodres<-likelihood(temp,xxeb,yeb,xb)
      lodres<-as.matrix(lodres)
      #### compute heredity#######
      ch_er <- as.numeric()
      ch_x <- cbind(matrix(1,nrow(xxeb),1),xxeb)

      ch_bb <- rbind(mean(yeb),as.matrix(xb))

      rm(xxeb)
      gc()


      for(i in 1:(ncol(ch_x)-1))
      {
        ch_xi <- ch_x[,(1+i)]
        as1 <- length(which(ch_xi==1))/nrow(ch_x)
        as2 <- 1-as1
        ch_er <- rbind(ch_er,(1-(as1-as2)*(as1-as2))*ch_bb[i+1]*ch_bb[i+1])
      }
      ch_v0 <- (1/(nrow(ch_x)-1))*(t(yeb-ch_x%*%ch_bb)%*%(yeb-ch_x%*%ch_bb))

      rm(ch_x)
      gc()


      if(var(yeb)>=sum(ch_er)+ch_v0){
        hered <- (ch_er/as.vector(var(yeb)))*100
      }else{
        hered <- (ch_er/as.numeric(sum(ch_er)+ch_v0))*100
      }

      bayeslodres<-cbind(ebrow,xb,lodres,hered)


      lodid<-which(bayeslodres[,4]>lodvalue)
      if(length(lodid)!=0){

        if(length(lodid)==1){
          lastres<-matrix(bayeslodres[lodid,],1,)
          xeb2<-matrix(xeb1[lodid,],1,)
        }else{
          lastres<-bayeslodres[lodid,]
          xeb2<-as.matrix(xeb1[lodid,])
        }

        rm(xeb1)
        gc()

        xxmaf<- xeb2
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
        pee<-round(var(yeb),4)

        vees<-matrix("",nrow = nrow(lastres),1)
        pees<-matrix("",nrow = nrow(lastres),1)
        pees[1,1]<-pee
        vees[1,1]<-vee
        result<-lastres
        result<-result

        if(nrow(result)>1){
          temp<-as.matrix(result[,3:5])
          temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
          temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
          wan<-cbind(result[,1:2],temp)
        }else{
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
        log10P <- as.matrix(round(-log10(pchisq(lodscore1*4.605,1,lower.tail = F)),4))
        if(nrow(tempwan)>1){
          tempwan1 <- cbind(tempwan[,1:5],log10P,tempwan[,6:10])
        }else{
          tempwan1 <- cbind(t(as.matrix(tempwan[,1:5])),log10P,t(as.matrix(tempwan[,6:10])))
        }
        wan <- tempwan1
        colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen (total)")
        wan<-as.data.frame(wan)
      }
      output<-list(result=wan)
    }
    return(output)
  }

  gen_bed <- BEDMatrix(paste(fileGen,".bed",sep=""))
  print("Running mrMLM programs with low RAM consumption, please be patient...")
  gen_bim <- fread(paste(fileGen,".bim",sep=""))
  gen_fam <- fread(paste(fileGen,".fam",sep=""))
  genRaw_dup_TF <- duplicated(gen_bim[,c(1,4)])
  if(sum(genRaw_dup_TF)!=0){
    gen_bim[genRaw_dup_TF,4] <- gen_bim[genRaw_dup_TF,4]+1
  }
  phy <- as.matrix(fread(filePhe,header=FALSE))
  PheName <- as.matrix(phy[1,-1])
  phy_match_list <- phy_match(gen_fam,phy)
  match_gen_ID_idex <- phy_match_list[[2]]
  samename_genphy <- phy_match_list[[1]][-1,1]
  phy <- as.matrix(apply(phy_match_list[[1]][-1,-1],2,as.numeric))


  # psmatrix
  psmatrix <- NULL
  if(!is.null(filePS)){
    filePS<-fread(filePS,header = FALSE,stringsAsFactors=T)
    filePS<-as.matrix(filePS)

    nnpprow<-dim(filePS)[1]
    nnppcol<-dim(filePS)[2]
    filePS[1,2:nnppcol]<-"  "
    psmatrixPre<-filePS[3:nnpprow,]
    namePop<-as.matrix(psmatrixPre[,1])
    sameGenPop<-intersect(samename_genphy,namePop)
    locPop<-match(sameGenPop,namePop)
    selectpsmatrixq<-psmatrixPre[locPop,-1]
    if(PopStrType=="Q"){
      selectpsmatrix<-matrix(as.numeric(selectpsmatrixq),nrow = length(locPop))
      coldelet<-which.min(apply(selectpsmatrix,2,sum))
      psmatrix<-as.matrix(selectpsmatrix[,-coldelet])
    }else if(PopStrType=="PCA"){
      psmatrix<-matrix(as.numeric(selectpsmatrixq),nrow = length(locPop))
    }else if(PopStrType=="EvolPopStr"){
      otrait_ind<-sort(unique(selectpsmatrixq))
      pop_col<-length(otrait_ind)-1
      pop_each<-numeric()
      for(j in 1:length(selectpsmatrixq)){
        if(selectpsmatrixq[j]==otrait_ind[1]){
          pop_0<-matrix(-1,1,pop_col)
        }else{
          pop_0<-matrix(0,1,pop_col)
          popnum_loc<-which(otrait_ind[]==selectpsmatrixq[j])
          pop_0[1,popnum_loc-1]<-1
        }
        pop_each<-rbind(pop_each,pop_0)
      }

      psmatrix=pop_each
    }
  }

  # fileCovphy
  covmatrixRaw<-NULL
  if(!is.null(fileCov)){
    covmatrixRaw<-fread(fileCov,header = FALSE,stringsAsFactors=T)
    covmatrixRaw<-as.matrix(covmatrixRaw)
  }
  if(is.null(covmatrixRaw)){
    phy<-phy
  }else{
    nncovrow<-nrow(covmatrixRaw)
    covmatrixPre<-covmatrixRaw[3:nncovrow,]
    namecov<-as.matrix(covmatrixPre[,1])
    sameGencov<-intersect(samename_genphy,namecov)
    loccov<-match(sameGencov,namecov)
    selectcovmatrixq<-covmatrixPre[loccov,-1]

    covname<-covmatrixRaw[2,-1]
    label<-substr(covname,1,3)
    if(("Cat"%in%label)&&("Con"%in%label)){
      cat_loc<-as.numeric(which(label=="Cat"))
      con_loc<-as.numeric(which(label=="Con"))
      selectcovmatrixqq<-selectcovmatrixq
      selectcovmatrixq<-selectcovmatrixq[,cat_loc]
      covnum<-t(selectcovmatrixq)
      yygg1<-numeric()
      for(i in 1:nrow(covnum)){
        otrait_ind<-sort(unique(covnum[i,]))
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
      yygg1<-cbind(yygg1,as.matrix(selectcovmatrixqq[,con_loc]))
    }else if(all(label=="Cat")){
      covnum<-t(selectcovmatrixq)
      yygg1<-numeric()
      for(i in 1:nrow(covnum)){
        otrait_ind<-sort(unique(covnum[i,]))
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
    }else if(all(label=="Con")){
      yygg1<-selectcovmatrixq
    }

    W.orig<-matrix(1,nrow(phy),1)
    xenvir<-cbind(W.orig,yygg1)
    xenvir<-apply(xenvir,2,as.numeric)
    beta<-solve(t(xenvir)%*%xenvir)%*%t(xenvir)%*%phy
    phy<-phy-xenvir%*%beta+W.orig
  }

  # svmlod
  if(is.null(svmlod)){svmlod <- 3}
  # dir
  if(is.null(dir)){dir <- getwd()}

  screen_PC<-function(reMR,phe_num){
    if(nrow(reMR)>=200){
      reMR4<-as.matrix(reMR[,4])
      datashuz1<-gen_bim[,2]
      calculate_gene<-gen_bed[match_gen_ID_idex,which(datashuz1%in%reMR4)]-1
      gene_shuzhi<-apply(calculate_gene,2,as.numeric)
      larsres<-lars(gene_shuzhi,phe_num,type = "lar",trace = FALSE,use.Gram=FALSE,max.steps=200)
      X<-gene_shuzhi[,which(larsres$beta[nrow(larsres$beta),]!=0)]
      MR200<-reMR[which(larsres$beta[nrow(larsres$beta),]!=0),]
      z<-cbind(matrix(1,nrow(gene_shuzhi),1),psmatrix)
      u1<-try({sblgwas(z,phe_num,X,t = -4,max.iter = 200,min.err = 1e-8)},silent=TRUE)
      if('try-error' %in% class(u1)){
        u1<-try({sblgwas(z,phe_num,X,t = -2,max.iter = 200,min.err = 1e-8)},silent=TRUE)
      }
      reMRshai<-MR200[which(u1$blup$p_wald<=0.01),]
      ind1<-which(larsres$beta[nrow(larsres$beta),]!=0)
      indz<-ind1[which(u1$blup$p_wald<=0.01)]
    }else if(nrow(reMR)<200){
      reMR4<-as.matrix(reMR[,4])
      datashuz1<-as.matrix(gen_bim[,2])
      calculate_gene<-gen_bed[match_gen_ID_idex,which(datashuz1%in%reMR4)]-1
      gene_shuzhi<-apply(calculate_gene,2,as.numeric)
      X<-gene_shuzhi
      z<-cbind(matrix(1,nrow(gene_shuzhi),1),psmatrix)
      u1<-try({sblgwas(z,phe_num,X,t = -4,max.iter = 200,min.err = 1e-8)},silent=TRUE)
      if('try-error' %in% class(u1)){
        u1<-try({sblgwas(z,phe_num,X,t = -2,max.iter = 200,min.err = 1e-8)},silent=TRUE)
      }
      reMRshai<-reMR[which(u1$blup$p_wald<=0.01),]

      indz<-which(u1$blup$p_wald<=0.01)
    }
    reMR<-cbind(reMRshai[,1:12],reMR[1:nrow(reMRshai),13:14])
    result<-list(reMR,indz)
    return(result)
  }

  trait_i <- 1
  for(trait_i in 1:trait){
    i <- trait_i
    reMR<-NULL;reFMR<-NULL;reFME<-NULL;rePLA<-NULL;rePKW<-NULL;reISIS<-NULL
    re1MR<-NULL;re1FMR<-NULL;re1FME<-NULL;re1PLA<-NULL;re1PKW<-NULL;re1ISIS<-NULL
    remanMR<-NULL;reqqMR<-NULL;remanFMR<-NULL;reqqFMR<-NULL;remanFME<-NULL;reqqFME<-NULL;
    replPLA<-NULL;remanPKW<-NULL;reqqPKW<-NULL; replISIS<-NULL;metaresult<-NULL;result_output<-NULL

    TRY1<-try({

      if("mrMLM"%in%method){
        outMR <- mrMLMFun.PC(gen_bed,gen_bim,gen_fam,phy=phy,phy_match_list,block_m=BLOCK_M,trait=trait_i)
        if(is.null(outMR$result2)==FALSE){
          me<-matrix("mrMLM",nrow(outMR$result2),1)
          tr<-matrix(trait_i,nrow(outMR$result2),1)
          trna<-matrix(PheName[trait_i,],nrow(outMR$result2),1)
          colnames(me)<-"Method"
          colnames(tr)<-"Trait ID"
          colnames(trna)<-"Trait name"
          reMR<-cbind(tr,trna,me,as.matrix(outMR$result2))
          if(nrow(reMR)>50){
            reMR<-screen_PC(reMR,phy[,trait_i])[[1]]
          }
        }
        me1<-matrix("mrMLM",nrow(outMR$result1),1)
        tr1<-matrix(trait_i,nrow(outMR$result1),1)
        tr1na<-matrix(PheName[trait_i,],nrow(outMR$result1),1)
        colnames(me1)<-"Method"
        colnames(tr1)<-"Trait ID"
        colnames(tr1na)<-"Trait name"
        re1MR<-cbind(tr1,tr1na,me1,as.matrix(outMR$result1))
      }
    },silent=FALSE)


    if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){
      TRY2<-try({

        if("FASTmrMLM"%in%method){
          outFMR <- FASTmrMLM.PC(gen_bed,gen_bim,gen_fam,phy=phy,phy_match_list,block_m=BLOCK_M,trait=trait_i)
          if(is.null(outFMR$result2)==FALSE){
            me<-matrix("FASTmrMLM",nrow(outFMR$result2),1)
            tr<-matrix(trait_i,nrow(outFMR$result2),1)
            trna<-matrix(PheName[trait_i,],nrow(outFMR$result2),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            reFMR<-cbind(tr,trna,me,as.matrix(outFMR$result2))
            if(nrow(reFMR)>50){
              reFMR<-screen_PC(reFMR,phy[,trait_i])[[1]]
            }
          }

          me1<-matrix("FASTmrMLM",nrow(outFMR$result1),1)
          tr1<-matrix(trait_i,nrow(outFMR$result1),1)
          tr1na<-matrix(PheName[trait_i,],nrow(outFMR$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1FMR<-cbind(tr1,tr1na,me1,as.matrix(outFMR$result1))
        }
      },silent=FALSE)
    }


    if ('try-error' %in% class(TRY2)|| !('try-error' %in% class(TRY2))){

      TRY3<-try({

        if("FASTmrEMMA"%in%method){
          outFME <- FASTmrEMMA.PC(gen_bed,gen_bim,gen_fam,phy=phy,phy_match_list,block_m=BLOCK_M,trait=trait_i,Likelihood="REML")

          if(is.null(outFME$result2)==FALSE){
            me<-matrix("FASTmrEMMA",nrow(outFME$result2),1)
            tr<-matrix(trait_i,nrow(outFME$result2),1)
            trna<-matrix(PheName[trait_i,],nrow(outFME$result2),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            reFME<-cbind(tr,trna,me,as.matrix(outFME$result2))
            if(nrow(reFME)>50){
              reFME<-screen_PC(reFME,phy[,trait_i])[[1]]
            }
          }
          me1<-matrix("FASTmrEMMA",nrow(outFME$result1),1)
          tr1<-matrix(trait_i,nrow(outFME$result1),1)
          tr1na<-matrix(PheName[trait_i,],nrow(outFME$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1FME<-cbind(tr1,tr1na,me1,as.matrix(outFME$result1))
        }
      },silent=FALSE)

    }


    if ('try-error' %in% class(TRY3)|| !('try-error' %in% class(TRY3))){

      TRY4<-try({

        if("pLARmEB"%in%method){
          outPLA <- pLARmEB.PC(gen_bed,gen_bim,gen_fam,phy=phy,phy_match_list,trait=trait_i)
          if(is.null(outPLA$result)==FALSE){
            me<-matrix("pLARmEB",nrow(outPLA$result),1)
            tr<-matrix(trait_i,nrow(outPLA$result),1)
            trna<-matrix(PheName[trait_i,],nrow(outPLA$result),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            rePLA<-cbind(tr,trna,me,as.matrix(outPLA$result))
            replPLA<-outPLA$plot
            if(nrow(rePLA)>50){
              rePLAQ<-screen_PC(rePLA,phy[,trait_i])
              rePLA<-rePLAQ[[1]]
            }
          }
        }
      },silent=FALSE)

    }


    if ('try-error' %in% class(TRY4)|| !('try-error' %in% class(TRY4))){

      TRY5<-try({

        if("pKWmEB"%in%method){
          outPKW <- pKWmEB.PC(gen_bed,gen_bim,gen_fam,phy=phy,phy_match_list,block_m=BLOCK_M,trait=trait_i)

          if(is.null(outPKW$result2)==FALSE){
            me<-matrix("pKWmEB",nrow(outPKW$result2),1)
            tr<-matrix(trait_i,nrow(outPKW$result2),1)
            trna<-matrix(PheName[trait_i,],nrow(outPKW$result2),1)
            colnames(me)<-"Method"
            colnames(tr)<-"Trait ID"
            colnames(trna)<-"Trait name"
            rePKW<-cbind(tr,trna,me,as.matrix(outPKW$result2))
            if(nrow(rePKW)>50){
              rePKW<-screen_PC(rePKW,phy[,trait_i])[[1]]
            }
          }
          me1<-matrix("pKWmEB",nrow(outPKW$result1),1)
          tr1<-matrix(trait_i,nrow(outPKW$result1),1)
          tr1na<-matrix(PheName[trait_i,],nrow(outPKW$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1PKW<-cbind(tr1,tr1na,me1,as.matrix(outPKW$result1))
        }
      },silent=FALSE)
    }



    if ('try-error' %in% class(TRY5)|| !('try-error' %in% class(TRY5))){
      TRY7<-try({
        output1qq<-list(re1MR,re1FMR,re1FME,re1PKW)
        output1q<-do.call(rbind,output1qq)

        if(isFALSE(all(lengths(output1qq)==0))){
          eff<-numeric()
          logp<-numeric()
          for(bb in c(which(lengths(output1qq)!=0))){
            eff_every<-as.matrix(output1qq[[bb]][,7])
            colnames(eff_every)<-colnames(output1qq[[bb]])[7]
            eff<-cbind(eff,eff_every)

            logp_every<-as.matrix(output1qq[[bb]][,8])
            colnames(logp_every)<-colnames(output1qq[[bb]])[8]
            logp<-cbind(logp,logp_every)
          }
          gencode1<-as.matrix(output1qq[[which(lengths(output1qq)!=0)[1]]][,9])
          colnames(gencode1)<-colnames(output1q)[[9]]

          output1<-cbind(output1qq[[which(lengths(output1qq)!=0)[1]]][,c(1,2,4,5,6)],eff,logp,gencode1)
          if("SNP effect (pKWmEB)"%in%colnames(output1)){
            output1<-output1[,-c(which(colnames(output1)%in%"SNP effect (pKWmEB)"))]
          }
        }else{
          output1<-output1q
        }

        write.table(output1,paste(dir,"/",trait_i,"_intermediate result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)

      },silent=FALSE)
    }

    if ('try-error' %in% class(TRY7)|| !('try-error' %in% class(TRY7))){
      TRY8<-try({

        output<-list(reMR,reFMR,reFME,rePLA,rePKW)
        output<-do.call(rbind,output)
        write.table(output,paste(dir,"/",trait_i,"_Final result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)

      },silent=FALSE)
    }



    if ('try-error' %in% class(TRY8)|| !('try-error' %in% class(TRY8))){
      TRY9<-try({

        if(DrawPlot==TRUE){


          if(isFALSE(all(lengths(output1qq)==0))){

            manwidth<-28000;manhei<-7000;manwordre<-60;manfigurere<-600
            qqwidth<-10000;qqhei<-10000;qqwordre<-60;qqfigurere<-600

            if(Plotformat1=="*.png"){
              png(paste(dir,"/",i,"_Manhattan plot.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=CriLOD)
              dev.off()

              png(paste(dir,"/",i,"_qq plot.png",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
              QQ_mrMLM(data_in=as.matrix(output1q))
              dev.off()

            }else if(Plotformat1=="*.tiff"){
              tiff(paste(dir,"/",i,"_Manhattan plot.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=CriLOD)
              dev.off()

              tiff(paste(dir,"/",i,"_qq plot.tiff",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
              QQ_mrMLM(data_in=as.matrix(output1q))
              dev.off()

            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_Manhattan plot.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=CriLOD)
              dev.off()

              jpeg(paste(dir,"/",i,"_qq plot.jpeg",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
              QQ_mrMLM(data_in=as.matrix(output1q))
              dev.off()

            }else if(Plotformat1=="*.pdf"){
              pdf(paste(dir,"/",i,"_Manhattan plot.pdf",sep=""),width=16,height=4,pointsize = 20)
              manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),CoorLwd=2,lodline=CriLOD)
              dev.off()

              pdf(paste(dir,"/",i,"_qq plot.pdf",sep=""),pointsize = 25)
              QQ_mrMLM(data_in=as.matrix(output1q),CoorLwd=2)
              dev.off()
            }

          }else{
            warning("Draw plot need intermediate result of mrMLM, FASTmrMLM, FASTmrEMMA or pKWmEB!")
          }
        }

      },silent=FALSE)
    }



  }


}else{
  screen<-function(reMR,rawgen,gen_num,phe_num,ps_num){
    if(nrow(reMR)>=200){
      reMR4<-as.matrix(reMR[,4])
      datashuz1<-rawgen[-1,1]
      calculate_gene<-t(gen_num[which(datashuz1%in%reMR4),-c(1,2)])
      gene_shuzhi<-apply(calculate_gene,2,as.numeric)
      larsres<-lars(gene_shuzhi,phe_num,type = "lar",trace = FALSE,use.Gram=FALSE,max.steps=200)
      X<-gene_shuzhi[,which(larsres$beta[nrow(larsres$beta),]!=0)]
      MR200<-reMR[which(larsres$beta[nrow(larsres$beta),]!=0),]
      z<-cbind(matrix(1,nrow(gene_shuzhi),1),ps_num)
      u1<-try({sblgwas(z,phe_num,X,t = -4,max.iter = 200,min.err = 1e-8)},silent=TRUE)
      if('try-error' %in% class(u1)){
        u1<-try({sblgwas(z,phe_num,X,t = -2,max.iter = 200,min.err = 1e-8)},silent=TRUE)
      }
      reMRshai<-MR200[which(u1$blup$p_wald<=0.01),]
      ind1<-which(larsres$beta[nrow(larsres$beta),]!=0)
      indz<-ind1[which(u1$blup$p_wald<=0.01)]
    }else if(nrow(reMR)<200){
      reMR4<-as.matrix(reMR[,4])
      datashuz1<-rawgen[-1,1]
      calculate_gene<-t(gen_num[which(datashuz1%in%reMR4),-c(1,2)])
      gene_shuzhi<-apply(calculate_gene,2,as.numeric)
      X<-gene_shuzhi
      z<-cbind(matrix(1,nrow(gene_shuzhi),1),ps_num)
      u1<-try({sblgwas(z,phe_num,X,t = -4,max.iter = 200,min.err = 1e-8)},silent=TRUE)
      if('try-error' %in% class(u1)){
        u1<-try({sblgwas(z,phe_num,X,t = -2,max.iter = 200,min.err = 1e-8)},silent=TRUE)
      }
      reMRshai<-reMR[which(u1$blup$p_wald<=0.01),]

      indz<-which(u1$blup$p_wald<=0.01)
    }
    reMR<-cbind(reMRshai[,1:12],reMR[1:nrow(reMRshai),13:14])
    result<-list(reMR,indz)
    return(result)
  }

  svrad<-SearchRadius;svmlod<-CriLOD;lars1<-SelectVariable

  if(Genformat=="Num"){Genformat<-1}else if(Genformat=="Cha"){Genformat<-2}else if(Genformat=="Hmp"){Genformat<-3}

  Plotformat1<-paste("*.",Plotformat,sep="");Plotformat2<-paste("*.",Plotformat,sep="")

  readraw<-ReadData(fileGen,filePhe,fileKin,filePS,fileCov,Genformat)

  PheName<-readraw$phename
  CLO<-readraw$CLO

  print("Running in progress, please be patient...")


  for (i in trait){

    InputData<-inputData(readraw,Genformat,method,i,PopStrType)

    reMR<-NULL;reFMR<-NULL;reFME<-NULL;rePLA<-NULL;rePKW<-NULL;reISIS<-NULL
    re1MR<-NULL;re1FMR<-NULL;re1FME<-NULL;re1PLA<-NULL;re1PKW<-NULL;re1ISIS<-NULL
    remanMR<-NULL;reqqMR<-NULL;remanFMR<-NULL;reqqFMR<-NULL;remanFME<-NULL;reqqFME<-NULL;
    replPLA<-NULL;remanPKW<-NULL;reqqPKW<-NULL; replISIS<-NULL;metaresult<-NULL;result_output<-NULL


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
          if(nrow(reMR)>50){
            reMR<-screen(reMR,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)[[1]]
          }
        }
        me1<-matrix("mrMLM",nrow(outMR$result1),1)
        tr1<-matrix(i,nrow(outMR$result1),1)
        tr1na<-matrix(PheName[i,],nrow(outMR$result1),1)
        colnames(me1)<-"Method"
        colnames(tr1)<-"Trait ID"
        colnames(tr1na)<-"Trait name"
        re1MR<-cbind(tr1,tr1na,me1,as.matrix(outMR$result1))
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
            if(nrow(reFMR)>50){
              reFMR<-screen(reFMR,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)[[1]]
            }
          }

          me1<-matrix("FASTmrMLM",nrow(outFMR$result1),1)
          tr1<-matrix(i,nrow(outFMR$result1),1)
          tr1na<-matrix(PheName[i,],nrow(outFMR$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1FMR<-cbind(tr1,tr1na,me1,as.matrix(outFMR$result1))
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
            if(nrow(reFME)>50){
              reFME<-screen(reFME,InputData$doFME$genRaw,InputData$doFME$gen,InputData$doFME$phe,InputData$doFME$psmatrix)[[1]]
            }
          }
          me1<-matrix("FASTmrEMMA",nrow(outFME$result1),1)
          tr1<-matrix(i,nrow(outFME$result1),1)
          tr1na<-matrix(PheName[i,],nrow(outFME$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1FME<-cbind(tr1,tr1na,me1,as.matrix(outFME$result1))
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
            replPLA<-outPLA$plot
            if(nrow(rePLA)>50){
              rePLAQ<-screen(rePLA,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)
              rePLA<-rePLAQ[[1]]
            }
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
            if(nrow(rePKW)>50){
              rePKW<-screen(rePKW,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)[[1]]
            }
          }
          me1<-matrix("pKWmEB",nrow(outPKW$result1),1)
          tr1<-matrix(i,nrow(outPKW$result1),1)
          tr1na<-matrix(PheName[i,],nrow(outPKW$result1),1)
          colnames(me1)<-"Method"
          colnames(tr1)<-"Trait ID"
          colnames(tr1na)<-"Trait name"
          re1PKW<-cbind(tr1,tr1na,me1,as.matrix(outPKW$result1))
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
            replISIS<-outISIS$plot
            if(nrow(reISIS)>50){
              reISISQ<-screen(reISIS,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)
              reISIS<-reISISQ[[1]]
            }
          }
        }
      },silent=FALSE)
    }

    if ('try-error' %in% class(TRY6)|| !('try-error' %in% class(TRY6))){
      TRY7<-try({
        output1qq<-list(re1MR,re1FMR,re1FME,re1PKW)
        output1q<-do.call(rbind,output1qq)

        if(isFALSE(all(lengths(output1qq)==0))){
          eff<-numeric()
          logp<-numeric()
          for(bb in c(which(lengths(output1qq)!=0))){
            eff_every<-as.matrix(output1qq[[bb]][,7])
            colnames(eff_every)<-colnames(output1qq[[bb]])[7]
            eff<-cbind(eff,eff_every)

            logp_every<-as.matrix(output1qq[[bb]][,8])
            colnames(logp_every)<-colnames(output1qq[[bb]])[8]
            logp<-cbind(logp,logp_every)
          }
          gencode1<-as.matrix(output1qq[[which(lengths(output1qq)!=0)[1]]][,9])
          colnames(gencode1)<-colnames(output1q)[[9]]

          output1<-cbind(output1qq[[which(lengths(output1qq)!=0)[1]]][,c(1,2,4,5,6)],eff,logp,gencode1)
          if("SNP effect (pKWmEB)"%in%colnames(output1)){
            output1<-output1[,-c(which(colnames(output1)%in%"SNP effect (pKWmEB)"))]
          }
        }else{
          output1<-output1q
        }

        write.table(output1,paste(dir,"/",i,"_intermediate result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)

      },silent=FALSE)
    }

    if ('try-error' %in% class(TRY7)|| !('try-error' %in% class(TRY7))){
      TRY8<-try({

        output<-list(reMR,reFMR,reFME,rePLA,rePKW,reISIS)
        output<-do.call(rbind,output)
        write.table(output,paste(dir,"/",i,"_Final result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)

      },silent=FALSE)
    }



    if ('try-error' %in% class(TRY8)|| !('try-error' %in% class(TRY8))){
      TRY9<-try({

        if(DrawPlot==TRUE){


          if(isFALSE(all(lengths(output1qq)==0))){

            manwidth<-28000;manhei<-7000;manwordre<-60;manfigurere<-600
            qqwidth<-10000;qqhei<-10000;qqwordre<-60;qqfigurere<-600

            if(Plotformat1=="*.png"){
              png(paste(dir,"/",i,"_Manhattan plot.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=CriLOD)
              dev.off()

              png(paste(dir,"/",i,"_qq plot.png",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
              QQ_mrMLM(data_in=as.matrix(output1q))
              dev.off()

            }else if(Plotformat1=="*.tiff"){
              tiff(paste(dir,"/",i,"_Manhattan plot.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=CriLOD)
              dev.off()

              tiff(paste(dir,"/",i,"_qq plot.tiff",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
              QQ_mrMLM(data_in=as.matrix(output1q))
              dev.off()

            }else if(Plotformat1=="*.jpeg"){
              jpeg(paste(dir,"/",i,"_Manhattan plot.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=CriLOD)
              dev.off()

              jpeg(paste(dir,"/",i,"_qq plot.jpeg",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
              QQ_mrMLM(data_in=as.matrix(output1q))
              dev.off()

            }else if(Plotformat1=="*.pdf"){
              pdf(paste(dir,"/",i,"_Manhattan plot.pdf",sep=""),width=16,height=4,pointsize = 20)
              manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),CoorLwd=2,lodline=CriLOD)
              dev.off()

              pdf(paste(dir,"/",i,"_qq plot.pdf",sep=""),pointsize = 25)
              QQ_mrMLM(data_in=as.matrix(output1q),CoorLwd=2)
              dev.off()
            }

          }else{
            warning("Draw plot need intermediate result of mrMLM, FASTmrMLM, FASTmrEMMA or pKWmEB!")
          }
        }

      },silent=FALSE)
    }

  }
}


}


ReadData<-function(fileGen=NULL,filePhe=NULL,fileKin=NULL,filePS=NULL,fileCov=NULL,Genformat=NULL){
  kkRaw<-NULL
  psmatrixRaw<-NULL
  covmatrixRaw<-NULL
  inputform<-Genformat
  CLO<-NULL
  if(!is.null(fileGen)){
    if(is.character(fileGen)==TRUE){
      genRaw<-fread(fileGen,header = FALSE,stringsAsFactors=T)
      genRaw_dup_TF <- duplicated(genRaw[,c(2,3)])
      if(sum(genRaw_dup_TF)!=0){
        if(sum(duplicated(genRaw))!=0){
          genRaw <- genRaw[!genRaw_dup_TF,]
        }else{
          genRaw_3 <- as.character(genRaw[,3])
          genRaw_3[genRaw_dup_TF] <- as.character(as.numeric(genRaw_3[genRaw_dup_TF])+1)
          genRaw[,3] <- as.factor(genRaw_3)
        }
      }
    }else{
      genRaw<-fileGen
      genRaw_dup_TF <- duplicated(genRaw[,c(2,3)])
      if(sum(genRaw_dup_TF)!=0){
        if(sum(duplicated(genRaw))!=0){
          genRaw <- genRaw[!genRaw_dup_TF,]
        }else{
          genRaw_3 <- as.character(genRaw[,3])
          genRaw_3[genRaw_dup_TF] <- as.character(as.numeric(genRaw_3[genRaw_dup_TF])+seq(1,sum(genRaw_dup_TF),1))
          genRaw[,3] <- as.factor(genRaw_3)
        }
      }
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
    hapName<- c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode")
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
  if(!is.null(fileCov)){
    covmatrixRaw<-fread(fileCov,header = FALSE,stringsAsFactors=T)
    covmatrixRaw<-as.matrix(covmatrixRaw)
  }
  phename<-as.matrix(pheRaw1q[1,2:ncol(pheRaw1q)])
  output<-list(genRaw=genRaw,pheRaw1q=pheRaw1q,kkRaw=kkRaw,psmatrixRaw=psmatrixRaw,covmatrixRaw=covmatrixRaw,phename=phename,CLO=CLO)
  return(output)
}

DoData<-function(genRaw=NULL,Genformat=NULL,pheRaw1q=NULL,kkRaw=NULL,psmatrixRaw=NULL,
                 covmatrixRaw=NULL,trait=NULL,type=NULL,PopStrType=NULL){
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
    gen<-as.matrix(needGen[-1,])
    gen<-matrix(as.numeric(gen),nrow=nrow(gen))
    rm(newGen,needGen)
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
    gen<-as.matrix(needGen[-1,])
    gen<-matrix(as.numeric(gen),nrow=nrow(gen))
    rm(newGen,needGen)
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
    hapName<-matrix(c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode"),1,)
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
    gen<-as.matrix(needGen[-1,])
    gen<-matrix(as.numeric(gen),nrow=nrow(gen))
    rm(newGen,needGen)
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
    selectpsmatrixq<-psmatrixPre[locPop,-1]
    if(PopStrType=="Q"){
      selectpsmatrix<-matrix(as.numeric(selectpsmatrixq),nrow = length(locPop))
      coldelet<-which.min(apply(selectpsmatrix,2,sum))
      psmatrix<-as.matrix(selectpsmatrix[,-coldelet])
    }else if(PopStrType=="PCA"){
      psmatrix<-matrix(as.numeric(selectpsmatrixq),nrow = length(locPop))
    }else if(PopStrType=="EvolPopStr"){
      otrait_ind<-sort(unique(selectpsmatrixq))
      pop_col<-length(otrait_ind)-1
      pop_each<-numeric()
      for(j in 1:length(selectpsmatrixq)){
        if(selectpsmatrixq[j]==otrait_ind[1]){
          pop_0<-matrix(-1,1,pop_col)
        }else{
          pop_0<-matrix(0,1,pop_col)
          popnum_loc<-which(otrait_ind[]==selectpsmatrixq[j])
          pop_0[1,popnum_loc-1]<-1
        }
        pop_each<-rbind(pop_each,pop_0)
      }

      psmatrix=pop_each
    }

  }
  if(is.null(covmatrixRaw)){
    phe<-phe
  }else{
    nncovrow<-nrow(covmatrixRaw)
    covmatrixPre<-covmatrixRaw[3:nncovrow,]
    namecov<-as.matrix(covmatrixPre[,1])
    sameGencov<-intersect(sameName,namecov)
    loccov<-match(sameGencov,namecov)
    selectcovmatrixq<-covmatrixPre[loccov,-1]

    covname<-covmatrixRaw[2,-1]
    label<-substr(covname,1,3)
    if(("Cat"%in%label)&&("Con"%in%label)){
      cat_loc<-as.numeric(which(label=="Cat"))
      con_loc<-as.numeric(which(label=="Con"))
      selectcovmatrixqq<-selectcovmatrixq
      selectcovmatrixq<-selectcovmatrixq[,cat_loc]
      covnum<-t(selectcovmatrixq)
      yygg1<-numeric()
      for(i in 1:nrow(covnum)){
        otrait_ind<-sort(unique(covnum[i,]))
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
      yygg1<-cbind(yygg1,as.matrix(selectcovmatrixqq[,con_loc]))
    }else if(all(label=="Cat")){
      covnum<-t(selectcovmatrixq)
      yygg1<-numeric()
      for(i in 1:nrow(covnum)){
        otrait_ind<-sort(unique(covnum[i,]))
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
    }else if(all(label=="Con")){
      yygg1<-selectcovmatrixq
    }

    W.orig<-matrix(1,nrow(phe),1)
    xenvir<-cbind(W.orig,yygg1)
    xenvir<-apply(xenvir,2,as.numeric)
    beta<-solve(t(xenvir)%*%xenvir)%*%t(xenvir)%*%phe
    phe<-phe-xenvir%*%beta+W.orig
  }
  genRaw<-genRaw[,1:12]
  doresult<-list(gen=gen,phe=phe,outATCG=outATCG,genRaw=genRaw,kk=kk,psmatrix=psmatrix)
  return(doresult)
}

inputData<-function(readraw,Genformat=NULL,method=NULL,trait=NULL,PopStrType=NULL){

  doMR<-NULL;doFME<-NULL

  if("mrMLM"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType)
  }

  if("FASTmrMLM"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType)
  }

  if("FASTmrEMMA"%in%method){
    doFME<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=1,PopStrType)
  }

  if("pLARmEB"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType)
  }
  if("pKWmEB"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType)
  }

  if("ISIS EM-BLASSO"%in%method){
    doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType)
  }

  output<-list(doMR=doMR,doFME=doFME)
  return(output)

}


















