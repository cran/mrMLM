\name{mrMLM-package}
\alias{mrMLM-package}
\alias{mrMLM}
\docType{package}
\title{
Multi-Locus Random-SNP-Effect Mixed Linear Model for Genome-Wide Association Studies and Linkage Analyses
}
\description{
 Conduct multi-locus GWAS and multi-locus QTL mapping under the framework of random-SNP-effect mixed linear model (mrMLM). First, each position (or marker) on the genome is scanned by mrMLM algorithm.
 Bonferroni correction is replaced by a less stringent selection criterion for significant test. Then, all the markers (or QTL) that are potentially associated with the trait are included in a multi-locus model, their effects are estimated by empirical Bayes and true QTN or QTL are identified by likelihood ratio test.

}
\details{

\tabular{ll}{
Package: \tab mrMLM\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2017-02-08\cr
Depends: \tab gWidgets,gWidgetsRGtk2,RGtk2Extras,MASS\cr
Imports: \tab RGtk2,methods,openxlsx,stringr,qqman,ggplot2,lars,ncvreg\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
Users can use library(mrMLM) to start the GUI and use 'mrMLM()' to restart the programe.Please Refer to the documentation section of the project homepage for the usage of FASTmrEMMA.
}
\author{
Ren Wen-Long,Ni Yuan-Li,Wen Yang-Jun,Zhang Jin,Cox Lwaka Tamba,Wang Shi-Bo,Huang Bo,Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
Improving power and accuracy of genome-wide association studies via a multi-locus mixed linear model methodology.
Wang Shi-Bo,Feng Jian-Ying,Ren Wen-Long,Huang Bo,Zhou Ling,Wen Yang-Jun,Zhang Jin,Jim M. Dunwell,Xu Shizhong*,Zhang Yuan-Ming*
}

\keyword{
 R
mrMLM
package
}

\examples{
library(mrMLM)
}
