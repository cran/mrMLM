\name{mrMLM-package}
\alias{mrMLM-package}
\alias{mrMLM}
\docType{package}
\title{
Multi-Locus Random-SNP-Effect Mixed Linear Model for Genome-Wide Association Studies and Linkage Analyses
}
\description{
Conduct multi-locus GWAS and multi-locus QTL mapping under the framework of random-SNP-effect mixed linear model (mrMLM). First, each position (or marker) on the genome is scanned by mrMLM algorithm. Bonferroni correction is replaced by a less stringent selection criterion for significant test. Then, all the markers (or QTL) that are potentially associated with the trait are included in a multi-locus model, their effects are estimated by empirical Bayes and true QTN or QTL are identified by likelihood ratio test.
}
\details{
\tabular{ll}{
Package: \tab mrMLM\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2016-12-19\cr
Depends: \tab gWidgets,gWidgetsRGtk2,RGtk2Extras,MASS\cr
Imports: \tab RGtk2,methods,openxlsx,stringr,ggplot2,qqman\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
Users can use library(mrMLM) to start the GUI and use 'mrMLM()' to restart the programe.Please Refer to the documentation section of the project homepage for the usage of mrMLM.

Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>\cr
Author: Wenlong Ren,Yuanli Ni,Shibo Wang,Yangjun Wen,Bo Huang and Yuanming Zhang
}
\references{
1.	Wang Shi-Bo, Feng Jian-Ying, Ren Wen-Long, Huang Bo, Zhou Ling, Wen Yang-Jun, Zhang Jin, Jim M. Dunwell, Xu Shizhong*, Zhang Yuan-Ming*. Improving power and accuracy of genome-wide association studies via a multi-locus mixed linear model methodology. Scientific Reports 2016, 6: 19444.\cr
2.	Wang Shi-Bo, Wen Yang-Jun, Ren Wen-Long, Ni Yuan-Li, Zhang Jin, Feng Jian-Ying, Zhang Yuan-Ming*. Mapping small-effect and linked quantitative trait loci for complex traits in backcross or DH populations via a multi-locus GWAS methodology. Scientific Reports 2016, 6: 29951.\cr
3.	Wen Yang-Jun, Zhang Hanwen, Ni Yuan-Li, Huang Bo, Zhang Jin, Feng Jian-Ying, Wang Shi-Bo, Jim M. Dunwell, Zhang Yuan-Ming*, Wu Rongling*. Methodological implementation of mixed linear models in multi-locus genome-wide association studies. Briefings in Bioinformatics 2016, DOI: 10.1093/bib/bbw145\cr
4.  Zhang Yuan-Ming et al. Mapping quantitative trait loci using naturally occurring genetic variance among commercial inbred line of maize (Zea mays L.). Genetics 2005, 169:2267-2275.\cr
}

\keyword{ 
R
mrMLM
package 
}
\examples{
library(mrMLM)
}
