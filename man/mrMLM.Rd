\name{mrMLM}
\alias{mrMLM}
\alias{mrMLM}
\docType{package}
\title{Multi-Locus Random-SNP-Effect Mixed Linear Model Tools for Genome-Wide Association Study 
}
\description{
Conduct multi-locus genome-wide association study under the framework of random-SNP-effect mixed linear model (mrMLM). First, each marker on the genome is scanned. Bonferroni correction is replaced by a less stringent selection criterion for significant test. Then, all the markers that are potentially associated with the trait are included in a multi-locus model, their effects are estimated by empirical Bayes and true QTNs are identified by likelihood ratio test.

}
\details{
\tabular{ll}{
Package: \tab mrMLM\cr
Type: \tab Package\cr
Version: \tab 3.1\cr
Date: \tab 2018-8-4\cr
Depends: \tab MASS,data.table,doParallel,foreach\cr
Imports: \tab methods,openxlsx,stringr,qqman,ggplot2,lars,ncvreg,coin\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\usage{
mrMLM(fileGen,filePhe,fileKin,filePS,Genformat,method,Likelihood,trait,
SearchRadius,CriLOD,SelectVariable,Bootstrap,DrawPlot,Plotformat,Resolution,dir)
}
\arguments{
	\item{fileGen}{File path and name in your computer of Genotype.}
	\item{filePhe}{File path and name in your computer of Phenotype.}
	\item{fileKin}{File path and name in your computer of Kinship.}
	\item{filePS}{File path and name in your computer of Population Structure.}
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
Zhang Ya-Wen, Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
Methodological implementation of mixed linear models in multi-locus genome-wide association studies.
Yang-Jun Wen,Hanwen Zhang,Yuan-Li Ni,Bo Huang,Jin Zhang,Jian-Ying Feng,Shi-Bo Wang,Jim M.Dunwell,Yuan-Ming Zhang,Rongling Wu.
}
\keyword{
 R
mrMLM
package
} 

\examples{
G1=data(Gen)
P1=data(Phe)
mrMLM(fileGen=Gen,filePhe=Phe,fileKin=NULL,filePS=NULL,Genformat="Num",
method=c("ISIS EM-BLASSO"),Likelihood="REML",trait=1:1,SearchRadius=20,
CriLOD=3,SelectVariable=50,Bootstrap=FALSE,
DrawPlot=FALSE,Plotformat="jpeg",Resolution="Low",dir=tempdir())
}