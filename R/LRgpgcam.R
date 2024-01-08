#' @title LR in grandparent-child identification with reference
#' @description
#' LR when a child (C) is alleged to be a grand-child of another individual (GP), and an offspring(A) of the alleged grandparent participated, with or without the assistant of B's other parent (M). Hp is that, C is an offspring of A's full-sibling; and Hd that, C is unrelated to GP and A. Inbreeding is not considered.
#'
#' @param A Genotype data of the alleged uncle/aunt, should be data.frame with 2 columns and ss rows, where ss stand for sample size;
#' @param C Genotype data of the child, with the same form with \code{A}
#' @param GP Genotype data of the alleged grandparent of C, with the same form with \code{A}
#' @param M Genotype data of the other parent of C, which can be \code{NULL} or as that of \code{A}
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names
#' @param rare frequency of rare allele on the locus
#' @param allelename if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the af matrix
#'
#' @return a data frame with one column and ss rows, containing log10 value of the CLR of each case
#' @details There might be no allele sharing between GP and A, or between M and B, if so, the related part in the LR calculation would be treated as 0, which can be further optimized in future version.
#'
#'
#' @export
#' @examples
#' # Construct pedi data.frames for two types of pedigrees
#' pedi1 <- data.frame(Person=c("GF","GM","F","A","M","C"),
#' Father=c("RI","RI","GF","GF","RI","F"),
#' Mother=c("RI","RI","GM","GM","RI","M"))
#' pedi2 <- data.frame(Person=c("GF","GM","F","A","M","C"),
#' Father=c("RI","RI","RI","GF","RI","F"),
#' Mother=c("RI","RI","RI","GM","RI","M"))
#' LR_1=LR_2=data.frame(Log10CLR=rep(0,10000))
#' for (i in 1:42) { 
#' # Simulate 10000 group of pedigrees where the Hp is true
#'   Genotype <- pedisimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,pedi = pedi1)
#'   LR_1 <-LR_1+ LRgpgcam(A=Genotype[,7:8],C=Genotype[,11:12],GP=Genotype[,1:2],M=Genotype[,9:10],
#'                         af=FortytwoSTR$afmatrix[[i]],rare=FortytwoSTR$rare[i][1,1])
#'#Simulate 10000 group of false pedigrees, i.e., P and C is unrelated to GP and A
#' Genotype <- pedisimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,pedi = pedi2)
#'   LR_2 <-LR_2+ LRgpgcam(A=Genotype[,7:8],C=Genotype[,11:12],GP=Genotype[,1:2],M=Genotype[,9:10],
#'                         af=FortytwoSTR$afmatrix[[i]],rare=FortytwoSTR$rare[i][1,1])
#' }
# histograms of CLR distributions in the two groups
#' xmin<-floor(min(min(LR_1$Log10CLR),min(LR_2$Log10CLR)))
#' xmax<-ceiling(max(max(LR_1$Log10CLR),max(LR_2$Log10CLR)))
#' par(mfrow = c(1, 2))
#' hist(LR_2$Log10CLR,xlab = expression(log[10]~CLR),main = "False pedigree",
#'      xlim = c(xmin,xmax), col = "red")
#' hist(LR_1$Log10CLR,xlab = expression(log[10]~CLR),main = "True cases",
#'      xlim = c(xmin,xmax), col = "blue")  

#' 
#' 

LRgpgcam<-function(A,C,GP,af,rare=NULL,allelename=FALSE,M=NULL){
  if (ncol(A)!=2 || ncol(C)!=2 || ncol(GP)!=2 || nrow(A)!=nrow(C) || nrow(A)!=nrow(GP) || nrow(C)!=nrow(GP)) {
    stop("false in individual data")
  }
  if (isTRUE(allelename)) {
    pa<-af[as.character(A[,1]),]
    pb<-af[as.character(A[,2]),]
    pc<-af[as.character(C[,1]),]
    pd<-af[as.character(C[,2]),]
  } else {
    pa<-af[A[,1],]
    pb<-af[A[,2],]
    pc<-af[C[,1],]
    pd<-af[C[,2],]
  }
  pa<-as.numeric(pa)
  pb<-as.numeric(pb)
  pc<-as.numeric(pc)
  pd<-as.numeric(pd)
  if (any(is.na(pa)) || any(is.na(pb)) || any(is.na(pc)) || any(is.na(pd))) {
    if (is.null(rare)) {
      stop("please input frequency data of rare alleles")
    }
    pa[is.na(pa)]<-rare
    pb[is.na(pb)]<-rare
    pc[is.na(pc)]<-rare
    pd[is.na(pd)]<-rare
  }
  if (isFALSE(is.null(M))) {
    if (ncol(M)!=2 || nrow(M)!=nrow(GP)) {
      stop("false in indiviudal data")
    }
    dmc<-as.double(M[,1]==C[,1])+as.double(M[,2]==C[,1])
    dmd<-as.double(M[,1]==C[,2])+as.double(M[,2]==C[,2])
  } else {
    dmc<-pc
    dmd<-pd
  }
  dgpa<-as.double(GP[,1]==A[,1])+as.double(GP[,2]==A[,1])
  dgpb<-as.double(GP[,1]==A[,2])+as.double(GP[,2]==A[,2])
  dgpc<-as.double(GP[,1]==C[,1])+as.double(GP[,2]==C[,1])
  dgpd<-as.double(GP[,1]==C[,2])+as.double(GP[,2]==C[,2])
  ac<-as.double(A[,1]==C[,1])
  ad<-as.double(A[,1]==C[,2])
  bc<-as.double(A[,2]==C[,1])
  bd<-as.double(A[,2]==C[,2])

  LR1<-(dmc*dgpd+dmd*dgpc)/(4*(pc*dmd+pd*dmc))
  LR2<-(pa*dgpb*ac*dmd+dgpa*pb*bc*dmd+pa*dgpb*ad*dmc+dgpa*pb*bd*dmc)/
    (4*(pa*dgpb+pb*dgpa)*(pc*dmd+pd*dmc))
  LR1[is.na(LR1)]=0
  LR2[is.na(LR2)]=0
  LR1[is.infinite(LR1)]=0
  LR2[is.infinite(LR2)]=0
  results=data.frame(Log10CLR=log10(1/4+LR1+LR2))

  return(results)
}
