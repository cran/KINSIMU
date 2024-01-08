#' @title Incest index
#' @description
#' Calculation of the ratio for a parent-child pair between the probability that the child's other parent is a relative of the present parent to the probability that the child's parents are unrelated
#'
#' @param Parent Genotypes of individual A of each case, which should be data.frame with 2 columns and ss rows, where ss stand for sample size
#' @param Child Genotypes of individual B of each case, which should be data.frame with 2 columns and ss rows, where ss stand for sample size
#' @param af name of allele frequency matrix, a data.frame of 1 column with the allele name being row names
#' @param rare frequency of rare allele
#' @param allelename if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the af matrix
#' @param phi kinship coefficent between the parents under Hp (that under Hd euquals to 0),with a defualt of 0.25, i.e., father-daughter incest or full-sibling incest
#' @return a data.frame containing 2 columns: Ngs (the genotype similarity score, 1 if both of the two alleles in child's genotype can be inherited from the mother) and Log10LR for each simulation
#' @details The premise of using this function should be the confirmation of the parent-child relationship between the two individuals, if there is no sharing alleles between them, 1-2phi would be regarded as the output
#'
#' @examples
#' # Simulate 10,000 mother-child pairs with father-daughter incest with pedisimu() function
#' # based on the 42 STRs in "FortytwoSTR" data.
#' pedi<-data.frame(Person=c("F","M","C"),Father=c("RI","F","F"),Mother=c("RI","RI","M"))
#' II_1=II_2=data.frame(Ngs=rep(0,10000),IIphi=rep(0,10000))
#' for (i in 1:42) {
#' Genotype1<-pedisimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,pedi = pedi)
#' # Calculate II for each case.
#' II_1<-II_1+IICAL(Parent = Genotype1[,3:4],Child = Genotype1[,5:6],af=FortytwoSTR$afmatrix[[i]],
#' rare=FortytwoSTR$rare[i][1,1],phi=0.25)
#' #Simulate 10,000 non-inbred mother-child pairs
#' Genotype2<-pairsimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,delta = c(0,1,0),allelename = FALSE)
#' II_2<-II_2+IICAL(Parent = Genotype2[,1:2],Child = Genotype2[,3:4],af=FortytwoSTR$afmatrix[[i]],
#' rare=FortytwoSTR$rare[i][1,1],phi=0.25)
#' }
#' # histograms of CII distributions in the two groups
#' xmin<-floor(min(min(II_1$IIphi),min(II_2$IIphi)))
#' xmax<-ceiling(max(max(II_1$IIphi),max(II_2$IIphi)))
#' par(mfrow = c(1,2))
#' hist(II_2$IIphi,xlab = expression(log[10]~CII),main = "Non-inbreed cases",
#' xlim = c(xmin,xmax), col = "red")
#' hist(II_1$IIphi,xlab = expression(log[10]~CII),main = "Inbreed cases",
#' xlim = c(xmin,xmax), col = "blue")
#' @export
#'

IICAL<-function(Parent,Child,af,rare=NULL,allelename=FALSE,phi=0.25){
  if (ncol(Parent)!=2 || ncol(Child)!=2 || nrow(Parent)!=nrow(Child)) {
    stop(paste("false in individual data"))
  }
  if (phi>0.25) {
    stop(paste("Wrong phi input"))
  }
  colnames(Parent)=c("P","M")
  colnames(Child)=c("P","M")
  if (isTRUE(allelename)) {
    pc<-af[as.character(Child$P),]
    pd<-af[as.character(Child$M),]
  } else {
    pc<-af[Child$P,]
    pd<-af[Child$M,]
  }
  if (any(is.na(pc)) || any(is.na(pd))) {
    if (is.null(rare)) {
      stop("please input frequency data of rare alleles")
    }
    pc[is.na(pc)]<-rare
    pd[is.na(pd)]<-rare
  }
  pc<-as.numeric(pc)
  pd<-as.numeric(pd)
  dc<-as.double(Parent$P==Child$P)+as.double(Parent$M==Child$P)
  dd<-as.double(Parent$P==Child$M)+as.double(Parent$M==Child$M)
  Ngs<-1-as.double(dc*dd==0)
  IIphi<-log10(2*phi*dc*dd/(dc*pd+dd*pc)+1-2*phi)
  IIphi[is.na(IIphi)]=log10(1-2*phi)
  IIphi[is.infinite(IIphi)]=log10(1-2*phi)
  results<-data.frame(Ngs=Ngs,IIphi=IIphi)
  return(results)
}
