\name{mrMLM}
\alias{mrMLM}
\alias{mrMLM}
\docType{package}
\title{Multi-Locus Random-SNP-Effect Mixed Linear Model Tools for Genome-Wide Association Study 
}
\description{
Conduct multi-locus genome-wide association study under the framework of multi-locus random-SNP-effect mixed linear model (mrMLM). First, each marker on the genome is scanned. Bonferroni correction is replaced by a less stringent selection criterion for significant test. Then, all the markers that are potentially associated with the trait are included in a multi-locus genetic model, their effects are estimated by empirical Bayes and all the nonzero effects were further identified by likelihood ratio test for ture QTL.
}
\details{
\tabular{ll}{
Package: \tab mrMLM\cr
Type: \tab Package\cr
Version: \tab 4.0\cr
Date: \tab 2019-10-21\cr
Depends: \tab ggplt2,lars\cr
Imports: \tab methods,qqman,foreach,ncvreg,coin,sampling,data.table,doParallel
\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\usage{
mrMLM(fileGen,filePhe,fileKin,filePS,PopStrType,fileCov,Genformat,
method,Likelihood,trait,SearchRadius,CriLOD,SelectVariable,Bootstrap,
DrawPlot,Plotformat,Resolution,dir)
}
\arguments{
	\item{fileGen}{File path and name in your computer of Genotype.}
	\item{filePhe}{File path and name in your computer of Phenotype.}
	\item{fileKin}{File path and name in your computer of Kinship.}
	\item{filePS}{File path and name in your computer of Population Structure.}
      \item{PopStrType}{The type of population structure.}                    
      \item{fileCov}{File path and name in your computer of covariate.} 
	\item{Genformat}{Format for genotypic codes, Num (number), Cha (character) and Hmp (Hapmap).}
	\item{method}{Six multi-locus GWAS methods. Users may select one to six methods, including mrMLM, FASTmrMLM, FASTmrEMMA, pLARmEB, pKWmEB and ISIS EM-BLASSO.}
	\item{Likelihood}{This parameter is only for FASTmrEMMA, including restricted maximum likelihood (REML) and maximum likelihood (ML).}
      \item{trait}{Traits analyzed from number 1 to number 2.}
	\item{SearchRadius}{This parameter is only for mrMLM and FASTmrMLM, indicating Search Radius in search of potentially associated QTN.}
	\item{CriLOD}{Critical LOD score for significant QTN.} 
      \item{SelectVariable}{This parameter is only for pLARmEB. SelectVariable=50 indicates that 50 potentially associated variables are selected from each chromosome. Users may change this number in real data analysis in order to obtain the best results as final results.}
      \item{Bootstrap}{This parameter is only for pLARmEB, including FASLE and TRUE, Bootstrap=FALSE indicates the analysis of only real dataset, Bootstrap=TRUE indicates the analysis of both real dataset and four resampling datasets.}
       \item{DrawPlot}{This parameter is for all the six methods, including FALSE and TRUE, DrawPlot=FALSE indicates no figure output, DrawPlot=TRUE indicates the output of the Manhattan, QQ and LOD score against genome position figures.}
     \item{Plotformat}{This parameter is for all the figure files, including *.jpeg, *.png, *.tiff and *.pdf.}
     \item{Resolution}{This parameter is for all the figure files, including Low and High.}
     \item{dir}{This parameter is for the save path.} 
}
\author{
Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
1. Zhang YM, Mao Y, Xie C, Smith H, Luo L, Xu S*. Genetics 2005,169:2267-2275.
2. Wang SB, Feng JY, Ren WL, Huang B, Zhou L, Wen YJ, et al. Sci Rep 2016,6:19444.
3. Tamba CL, Ni YL, Zhang YM*. PLoS Comput Biol 2017,13(1):e1005357.
4. Zhang J, Feng JY, Ni YL, Wen YJ, Niu Y, Tamba CL, et al. Heredity 2018,118(6):517-524.
5. Ren WL, Wen YJ, Dunwell JM, Zhang YM*. Heredity 2018,120(3): 208-218.
6. Wen YJ, Zhang H, Ni YL, Huang B, Zhang J, Feng JY, et al. Brief Bioinform 2018,19(4):700-712.
7. Tamba CL, Zhang YM. bioRxiv,preprint first posted online Jun. 7, 2018, doi:https://doi.org/10.1101/341784.
8. Zhang YW, Tamba CL, Wen YJ, Li P, Ren WL, Ni YL, et al. Genomics, Proteomics & Bioinformatics, resubmission.
}
\examples{
Ge1=data(Genotype)
Ph1=data(Phenotype)
mrMLM(fileGen=Genotype,filePhe=Phenotype,Genformat="Num",method=c("FASTmrMLM"),
trait=1,CriLOD=3,dir=tempdir())
}