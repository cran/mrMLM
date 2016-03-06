\name{mrMLM-package}
\alias{mrMLM-package}
\alias{mrMLM}
\docType{package}
\title{
Multilocus Random Mixed Linear Model for Genome-Wide Association Study
}
\description{
Conduct multi-locus random-SNP-effect mixed linear model analysis for genome-wide association study. First, each marker is scanned using random-SNP-effect mixed linear model method. Then, all the markers selected from the first step are included in the multi-locus model and their effects are estimated by an EM empirical Bayes method. No multiple test correction is needed.
}
\details{
\tabular{ll}{
Package: \tab mrMLM\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-03-05\cr
Depends: \tab gWidgets, gWidgetsRGtk2, RGtk2Extras, qqman\cr
Imports: \tab RGtk2,methods\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
Users can use library(mrMLM) to start the GUI and use 'mrMLM()' to restart the programe.Please Refer to the documentation section of the project homepage for the usage of mrMLM.
}
\author{
Wenlong Ren,Shibo Wang and Yuanming Zhang\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
Wang S, Feng J, Ren W and etc. 2016.Improving power and accuracy of genome-wide association studies via a multi-locus linear model methodology. Scientific Reports 6: 19444.
}
\keyword{
R
mrMLM
package 
}
\examples{
library(mrMLM)
}
