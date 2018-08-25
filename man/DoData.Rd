\name{DoData}
\alias{DoData}
\title{process raw data}
\description{
process raw data for later use}
\usage{
DoData(genRaw,Genformat,pheRaw1q,kkRaw,psmatrixRaw,trait,type)
}
\arguments{
	\item{genRaw}{raw genotype matrix.}
	\item{Genformat}{genotype format.}
	\item{pheRaw1q}{raw phenotype matrix.}
	\item{kkRaw}{raw kinship matrix.}
      \item{psmatrixRaw}{raw population structure matrix.}
      \item{trait}{which trait to analysis.} 
      \item{type}{which type to transform.}
      }
\author{
Zhang Ya-Wen, Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(Gen)
P1=data(Phe)
readraw=ReadData(fileGen=Gen,filePhe=Phe,fileKin=NULL,filePS =NULL,
Genformat=1)
result=DoData(readraw$genRaw,Genformat=3,readraw$pheRaw1q,readraw$kkRaw,
readraw$psmatrixRaw,trait=1,type=2)
}



