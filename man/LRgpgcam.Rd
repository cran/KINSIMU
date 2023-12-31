% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRgpgcam.R
\name{LRgpgcam}
\alias{LRgpgcam}
\title{LR in grandparent-child identification with reference}
\usage{
LRgpgcam(A, C, GP, af, rare = NULL, allelename = FALSE, M = NULL)
}
\arguments{
\item{A}{Genotype data of the alleged uncle/aunt, should be data.frame with 2 columns and ss rows, where ss stand for sample size;}

\item{C}{Genotype data of the child, with the same form with \code{A}}

\item{GP}{Genotype data of the alleged grandparent of C, with the same form with \code{A}}

\item{af}{name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names}

\item{rare}{frequency of rare allele on the locus}

\item{allelename}{if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the af matrix}

\item{M}{Genotype data of the other parent of C, which can be \code{NULL} or as that of \code{A}}
}
\value{
a data frame with one column and ss rows, containing log10 value of the CLR of each case
}
\description{
LR when a child (C) is alleged to be a grand-child of another individual (GP), and an offspring(A) of the alleged grandparent participated, with or without the assistant of B's other parent (M). Hp is that, C is an offspring of A's full-sibling; and Hd that, C is unrelated to GP and A. Inbreeding is not considered.
}
\details{
There might be no allele sharing between GP and A, or between M and B, if so, the related part in the LR calculation would be treated as 0, which can be further optimized in future version.
}
\examples{
# Construct pedi data.frames for two types of pedigrees
pedi1 <- data.frame(Person=c("GF","GM","F","A","M","C"),
Father=c("RI","RI","GF","GF","RI","F"),
Mother=c("RI","RI","GM","GM","RI","M"))
pedi2 <- data.frame(Person=c("GF","GM","F","A","M","C"),
Father=c("RI","RI","RI","GF","RI","F"),
Mother=c("RI","RI","RI","GM","RI","M"))
LR_1=LR_2=data.frame(Log10CLR=rep(0,10000))
for (i in 1:42) { 
# Simulate 10000 group of pedigrees where the Hp is true
  Genotype <- pedisimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,pedi = pedi1)
  LR_1 <-LR_1+ LRgpgcam(A=Genotype[,7:8],C=Genotype[,11:12],GP=Genotype[,1:2],M=Genotype[,9:10],
                        af=FortytwoSTR$afmatrix[[i]],rare=FortytwoSTR$rare[i][1,1])
#Simulate 10000 group of false pedigrees, i.e., P and C is unrelated to GP and A
Genotype <- pedisimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,pedi = pedi2)
  LR_2 <-LR_2+ LRgpgcam(A=Genotype[,7:8],C=Genotype[,11:12],GP=Genotype[,1:2],M=Genotype[,9:10],
                        af=FortytwoSTR$afmatrix[[i]],rare=FortytwoSTR$rare[i][1,1])
}
xmin<-floor(min(min(LR_1$Log10CLR),min(LR_2$Log10CLR)))
xmax<-ceiling(max(max(LR_1$Log10CLR),max(LR_2$Log10CLR)))
par(mfrow = c(1, 2))
hist(LR_2$Log10CLR,xlab = expression(log[10]~CLR),main = "False pedigree",
     xlim = c(xmin,xmax), col = "red")
hist(LR_1$Log10CLR,xlab = expression(log[10]~CLR),main = "True cases",
     xlim = c(xmin,xmax), col = "blue")  


}
