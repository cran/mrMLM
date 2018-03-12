\name{ISIS}
\alias{ISIS}
\title{To perform GWAS with ISIS EM-BLASSO method}
\description{
Iterative Sure Independence Screening EM-Bayesian LASSO}
\usage{
ISIS(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,CLO)
}
\arguments{
	\item{gen}{genotype matrix.}
	\item{phe}{phenotype matrix.}
	\item{outATCG}{genotype for code 1.}
	\item{genRaw}{raw genotype.}
	\item{kk}{kinship matrix.}
	\item{psmatrix}{population structure matrix.}
	\item{svpal}{Critical P-value for selecting variable.}
	\item{svmlod}{Critical LOD score for significant QTN.}
	\item{Genformat}{Format for genotypic codes.}
      \item{CLO}{number of CPU.}
	}
\author{
Zhang Ya-Wen, Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}


