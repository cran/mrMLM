\name{mrMLM-package}
\alias{mrMLM-package}
\alias{mrMLM}
\docType{package}
\title{
Multilocus Random Mixed Linear Model for Genome-Wide Association Study and Linkage Analysis
}
\description{
Conduct multi-locus random-SNP-effect mixed linear model analysis for genome-wide association study and linkage analysis. First, each position (or marker) on the genome is scanned using random-SNP-effect mixed linear model method. Then, all the markers (or putative QTL) selected from the first step are included in the multi-locus model and their effects are estimated by an EM empirical Bayes method. Bonferroni correction is replaced by a less stringent selection criterion for significance test.
}
\details{

\tabular{ll}{
Package: \tab mrMLM\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2016-06-28\cr
Depends: \tab gWidgets,gWidgetsRGtk2,RGtk2Extras\cr
Imports: \tab RGtk2,methods,openxlsx,stringr,ggplot2,qqman\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
Users can use library(mrMLM) to start the GUI and use 'mrMLM()' to restart the programe.Please Refer to the documentation section of the project homepage for the usage of mrMLM.

Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>\cr
Author:Wenlong Ren,Yuanli Ni,Shibo Wang,Bo Huang and Yuanming Zhang

}
\references{
1. Wang S, Feng J, Ren W and etc. 2016.Improving power and accuracy of genome-wide association studies via a multi-locus linear model methodology. Scientific Reports 6: 19444.
2. Wang Shi-Bo, Wen Yang-Jun, Ren Wen-Long, Ni Yuan-Li, Zhang Jin, Feng Jian-Ying, Zhang Yuan-Ming.Mapping small-effect and linked quantitative trait loci for complex traits in backcross or DH populations via a multi-locus GWAS methodology. Scientific Reports 2016, 6: 29951.
}

\keyword{
R
mrMLM
package 
}
\examples{
library(mrMLM)
}
