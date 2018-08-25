\name{inputData}
\alias{inputData}
\title{Input data which have been transformed}
\description{
Input all the dataset which have been transformed}
\usage{
inputData(readraw,Genformat,method,trait)
}
\arguments{
	\item{readraw}{genotype matrix.}
	\item{Genformat}{genotype format.}
	\item{method}{which method to analysis.}
	\item{trait}{which trait to analysis.}
      }
\author{
Zhang Ya-Wen, Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(Gen)
P1=data(Phe)
Readraw=ReadData(fileGen=Gen,filePhe=Phe,fileKin=NULL,filePS =NULL,
Genformat=1)
result=inputData(readraw=Readraw,Genformat=3,method="mrMLM",trait=1)
}



