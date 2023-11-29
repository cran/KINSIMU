#' @title LR in half-sibling identification with the identical parent participated
#' @description
#' LR when a pair of siblings and one of their identical parent participated, Hp and Hd of which being "the other parents of the two siblings being specific related" and "the other parents of them being unrelated". Inbreeding factors are not taken into consideration.
#'
#' @param A Genotype data of the first sibling, should be data.frame with 2 columns and ss rows, where ss stand for sample size;
#' @param B Genotype data of the second sibling, with the same form with \code{A}
#' @param P Genotype data of the identical parent of \code{A} and \code{B}, with the same form with them
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names
#' @param rare frequency of rare allele on the locus
#' @param allelename if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the af matrix
#' @param phi kinship coefficient of the other parents of the two siblings under Hp, with default of 1/2, i.e, being identical or MZ
#'
#' @return a data frame with one column and ss rows, containing log10 value of the CLR of each case
#' @details Mutation might be found between P with A or B, if so, LR would be output as 1-phi, which can be further optimized in the future version.
#' @export
#' @examples
#' # Simulate 10000 groups of A/B/P where A is full sibiling of B
#' pedi <- data.frame(Person=c("F","M","A","B"),
#' Father=c("RI","RI","F","F"),
#' Mother=c("RI","RI","M","M"))
#' Genotype=pedisimu(af = FortytwoSTR$afmatrix[[1]],ss = 10000,pedi = pedi)
#' #Calculation
#' LR_1=LRhsip(A=Genotype[,5:6],B=Genotype[,7:8],P=Genotype[,3:4],
#' af = FortytwoSTR$afmatrix[[1]],rare=FortytwoSTR$rare[1])
#'
#' # Simulate 10000 groups of A/B/P where A is half sibling of B, i.e., the true phi=0
#' pedi <- data.frame(Person=c("M","A","B"),
#' Father=c("RI","RI","RI"),
#' Mother=c("RI","M","M"))
#' Genotype=pedisimu(af = FortytwoSTR$afmatrix[[1]],ss = 10000,pedi = pedi)
#' #Calculation
#' LR_2=LRhsip(A=Genotype[,3:4],B=Genotype[,5:6],P=Genotype[,1:2],
#' af = FortytwoSTR$afmatrix[[1]],rare=FortytwoSTR$rare[1])
#'

LRhsip<-function(A,B,P,af,rare=NULL,allelename=FALSE,phi=0.5){
  if (ncol(A)!=2 || ncol(B)!=2 || ncol(P)!=2 || nrow(A)!=nrow(B) || nrow(A)!=nrow(P) || nrow(B)!=nrow(P)) {
    stop("false in individual data")
  }
  if (phi>1/2) {
    stop("false in phi value")
  }
  dma<-as.double(P[,1]==A[,1])+as.double(P[,2]==A[,1])
  dmb<-as.double(P[,1]==A[,2])+as.double(P[,2]==A[,2])
  dmc<-as.double(P[,1]==B[,1])+as.double(P[,2]==B[,1])
  dmd<-as.double(P[,1]==B[,2])+as.double(P[,2]==B[,2])
  ac<-as.double(A[,1]==B[,1])
  ad<-as.double(A[,2]==B[,1])
  bc<-as.double(A[,1]==B[,2])
  bd<-as.double(A[,2]==B[,2])
  if (isTRUE(allelename)) {
    pa<-af[as.character(A[,1]),]
    pb<-af[as.character(A[,2]),]
    pc<-af[as.character(B[,1]),]
    pd<-af[as.character(B[,2]),]
    if (any(is.na(pa)) || any(is.na(pb)) || any(is.na(pc)) || any(is.na(pd))) {
      if (is.null(rare)) {
        stop("please input frequency data of rare alleles")
      }
      pa[is.na(pa)]<-rare
      pb[is.na(pb)]<-rare
      pc[is.na(pc)]<-rare
      pd[is.na(pd)]<-rare
    }
    pa<-as.numeric(pa)
    pb<-as.numeric(pb)
    pc<-as.numeric(pc)
    pd<-as.numeric(pd)
  } else {
    pa<-af[A[,1],]
    pb<-af[A[,2],]
    pc<-af[B[,1],]
    pd<-af[B[,2],]
  }

  LR1<-(pa*dmb*(ac*dmd+ad*dmc)+pb*dma*(bc*dmd+bd*dmc))/((dmc*pd+dmd*pc)*(dma*pb+dmb*pa))
  LR1[is.na(LR1)]=0
  results=data.frame(Log10CLR=log10(1-phi+phi*LR1))
  return(results)

}
