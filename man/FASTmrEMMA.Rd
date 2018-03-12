\name{FASTmrEMMA}
\alias{FASTmrEMMA}
\title{To perform GWAS with FASTmrEMMA method}
\description{
FAST multi-locus random-SNP-effect EMMA}
\usage{
FASTmrEMMA(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,Likelihood,CLO)
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
      \item{Likelihood}{restricted maximum likelihood (REML) and maximum likelihood (ML).}
      \item{CLO}{number of CPU.}
}
\author{
Zhang Ya-Wen, Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}


