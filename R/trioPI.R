#' @title LR in standard trio cases
#' @description Calculating LR in cases where 3 participants being available, a child, his/her biological mother (or father) whose parentage is confirmed, and a male (or female) who is unrelated to the confirmed parent and alleged to be specific relative of the child, usually the father. Null hypothesis, i.e., that the alleged participant is unrelated to the child, is taken as Hd.
#'
#' @param AR Genotype of the alleged relative, which should be data.frame with 2 columns and ss rows, where ss stand for sample size;
#' @param C Genotype of the child, see in \code{AR} for the data form
#' @param TP Genotype of the confirmed parent of the child, see in \code{AR} for the data form
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names
#' @param rare frequency of rare allele on the locus
#' @param muAtoC mutation rate from \code{AR} to \code{C} if AR is alleged to be the child's parent, with a default of 0.002, please note that mistakes would be introduced if the mutation rate is larger than 0.2
#' @param muTtoC mutation rate from \code{TP} to \code{C}, with a default of 0.002/3.5, please note that mistakes would be introduced if the mutation rate is larger than 0.2
#' @param allelename if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the af matrix
#' @param kappa1 kappa_1 of the alleged relationship between AR and C, with a default of 1, meaning the AR is alleged to be the other parent of C
#'
#'
#' @details If any one of TP or AR cannot provide any of the C's alleles through integer steps, \code{muAtoC} would be output when calculating PI, and under other situations, the case with minimum mutation steps, (AR to C)+(TP to C) under Hp and (TP to C) under Hd, would be considered. Much more is required for further discussion in PI calculation involving mutations.
#'
#' @return a data.frame containing \code{ss} rows and 1 column, containing the log10 of LR for each group
#' @export
#' @examples
#'  # Three types of pedigrees are simulated: pedi1: father-mother-child; 
#'  # pedi2: random male-mother-child; and pedi3: uncle-mother-child
#'  pedi1 <- data.frame(Person=c("F","M","C"),Father=c("RI","RI","F"),Mother=c("RI","RI","M"))
#'  pedi2 <- data.frame(Person=c("F","M","C"),Father=c("RI","RI","RI"),Mother=c("RI","RI","M"))
#'  pedi3 <- data.frame(Person=c("GF","GM","AR","F","M","C"),
#'  Father=c("RI","RI","GF","GF","RI","F"),
#'  Mother=c("RI","RI","GM","GM","RI","M"))
#'  # Two types of LRs are calculated: PI_1 and PI_2: Paternity index in trio cases; 
#'  # AI_1 and AI_2: Avuncular index in trio cases.
#'  PI_1=PI_2=AI_1=AI_2=data.frame(Log10CLR=rep(0,10000))
#'  # Simulation are carried out based on the frequency data of the 42 STRs in FortytwoSTR dataset
#'  # Setting sample size as 10,000
#'  for (i in 1:42) {
#'    Genotype1<-pedisimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,pedi = pedi1)
#'    PI_1<-PI_1+trioPI(AR=Genotype1[,1:2],TP=Genotype1[,3:4],
#'                 C=Genotype1[,5:6],af=FortytwoSTR$afmatrix[[i]],
#'                 rare=FortytwoSTR$rare[i])
#'    Genotype2<-pedisimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,pedi = pedi2)
#'    PI_2<-PI_2+trioPI(AR=Genotype2[,1:2],TP=Genotype2[,3:4],
#'                 C=Genotype2[,5:6],af=FortytwoSTR$afmatrix[[i]],
#'                 rare=FortytwoSTR$rare[i])
#'    AI_2<-AI_2+trioPI(AR=Genotype2[,1:2],TP=Genotype2[,3:4],
#'                 C=Genotype2[,5:6],af=FortytwoSTR$afmatrix[[i]],
#'                 rare=FortytwoSTR$rare[i],kappa1=0.5)             
#'    Genotype3<-pedisimu(af = FortytwoSTR$afmatrix[[i]],ss = 10000,pedi = pedi3)
#'    AI_1<-AI_1+trioPI(AR=Genotype3[,5:6],TP=Genotype3[,9:10],
#'                 C=Genotype3[,11:12],af=FortytwoSTR$afmatrix[[i]],
#'                 rare=FortytwoSTR$rare[i],kappa1=0.5)
#'  }
#'  #histogram of the final results
#'  xmin1<-floor(min(min(PI_1$Log10CLR),min(PI_2$Log10CLR)))
#'  xmax1<-ceiling(max(max(PI_1$Log10CLR),max(PI_2$Log10CLR)))
#'  xmin2<-floor(min(min(AI_1$Log10CLR),min(AI_2$Log10CLR)))
#'  xmax2<-ceiling(max(max(AI_1$Log10CLR),max(AI_2$Log10CLR)))
#'  par(mfrow = c(2, 2))
#'  hist(PI_1$Log10CLR,xlab = expression(log[10]~CPI),main = "True parentage cases",
#'       xlim = c(xmin1,xmax1), col = "blue")
#'  hist(AI_1$Log10CLR,xlab = expression(log[10]~CAI),main = "True avuncular cases",
#'       xlim = c(xmin2,xmax2), col = "blue")  
#'  hist(PI_2$Log10CLR,xlab = expression(log[10]~CPI),main = "False pedigree in parentage cases",
#'       xlim = c(xmin1,xmax1), col = "red")
#'  hist(AI_2$Log10CLR,xlab = expression(log[10]~CAI),main = "False pedigree in avuncular cases",
#'       xlim = c(xmin2,xmax2), col = "red")

#'


#'

trioPI<-function(AR,C,TP,af,rare=NULL,allelename=FALSE,muAtoC=0.002,muTtoC=0.002/3.5,kappa1=1){
  if (ncol(AR)!=2 || ncol(C)!=2 || ncol(TP)!=2 || ncol(AR)!=ncol(C) || ncol(AR)!=ncol(TP) || ncol(C)!=ncol(TP)) {
    stop(paste("false in individual data"))
  }
  if (kappa1>1) {
    stop(paste("kappa1 > 1"))
  }
  if (isTRUE(allelename)) {
    pc<-af[as.character(C[,1]),]
    pd<-af[as.character(C[,2]),]
  } else{
    pc<-af[C[,1],]
    pd<-af[C[,2],]
    an<-as.data.frame(as.numeric(row.names(af)))
    AR[,1]<-an[AR[,1],]
    AR[,2]<-an[AR[,2],]
    C[,1]<-an[C[,1],]
    C[,2]<-an[C[,2],]
    TP[,1]<-an[TP[,1],]
    TP[,2]<-an[TP[,2],]
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
  dmc<-as.double(TP[,1]==C[,1])+as.double(TP[,2]==C[,1])
  dmd<-as.double(TP[,1]==C[,2])+as.double(TP[,2]==C[,2])
  dfc<-as.double(AR[,1]==C[,1])+as.double(AR[,2]==C[,1])
  dfd<-as.double(AR[,1]==C[,2])+as.double(AR[,2]==C[,2])
  if (kappa1<1) {
    if (isTRUE(all(dmc+dmd>0))) {
      results<-log10(1-kappa1+kappa1*(dfc*dmd+dfd*dmc)/(2*(dmc*pd+dmd*pc)))
    } else {
      ac<-round(abs(TP[,1]-C[,1])+abs(TP[,1]-C[,1])%%1*100000,1)
      ad<-round(abs(TP[,1]-C[,2])+abs(TP[,1]-C[,2])%%1*100000,1)
      bc<-round(abs(TP[,2]-C[,1])+abs(TP[,2]-C[,1])%%1*100000,1)
      bd<-round(abs(TP[,2]-C[,2])+abs(TP[,2]-C[,2])%%1*100000,1)
      msteps<-pmin(ac,ad,bc,bd)
      dmc<-as.double(ac==msteps)+as.double(bc==msteps)
      dmd<-as.double(ad==msteps)+as.double(bd==msteps)
      mstepc<-pmin(ac,bc)
      mstepd<-pmin(ad,bd)
      dmc2<-1+as.double(ac==bc)
      dmd2<-1+as.double(ad==bd)
      results<-data.frame(Log10CLR=log10(1-kappa1+
                       kappa1*as.double(msteps==0)*(dfc*dmd+dmc*dfd)/(2*(dmc*pd+dmd*pc))+
                       kappa1*as.double(msteps>0)*(dfc*dmd2*10^(msteps-mstepd)+dfd*dmc2*10^(msteps-mstepc))/(2*(dmc*pd+dmd*pc))))
    }
  } else {
    if (isTRUE(all(dmc*dfd+dmd*dfc>0))) {
      results<-data.frame(Log10CLR=log10((dfc*dmd+dfd*dmc)/(2*(dmc*pd+dmd*pc))))
    } else {
        ac<-round(abs(AR[,1]-C[,1])+abs(AR[,1]-C[,1])%%1*100000,1)
        ad<-round(abs(AR[,1]-C[,2])+abs(AR[,1]-C[,2])%%1*100000,1)
        bc<-round(abs(AR[,2]-C[,1])+abs(AR[,2]-C[,1])%%1*100000,1)
        bd<-round(abs(AR[,2]-C[,2])+abs(AR[,2]-C[,2])%%1*100000,1)
        ec<-round(abs(TP[,1]-C[,1])+abs(TP[,1]-C[,1])%%1*100000,1)
        ed<-round(abs(TP[,1]-C[,2])+abs(TP[,1]-C[,2])%%1*100000,1)
        fc<-round(abs(TP[,2]-C[,1])+abs(TP[,2]-C[,1])%%1*100000,1)
        fd<-round(abs(TP[,2]-C[,2])+abs(TP[,2]-C[,2])%%1*100000,1)
        f1m1<-ac+ed
        f1m2<-ac+fd
        f2m1<-bc+ed
        f2m2<-bc+fd
        m1f1<-ec+ad
        m1f2<-ec+bd
        m2f1<-fc+ad
        m2f2<-fc+bd
        ac2=ac
        ad2=ad
        bc2=bc
        bd2=bd
        ec2=ec
        ed2=ed
        fc2=fc
        fd2=fd
        
        ac2[ac==0]<-1-log10(2/muAtoC)
        ad2[ad==0]<-1-log10(2/muAtoC)
        bc2[bc==0]<-1-log10(2/muAtoC)
        bd2[bd==0]<-1-log10(2/muAtoC)
        ec2[ec==0]<-1-log10(2/muTtoC)
        ed2[ed==0]<-1-log10(2/muTtoC)
        fc2[fc==0]<-1-log10(2/muTtoC)
        fd2[fd==0]<-1-log10(2/muTtoC) # transforming 0 to make no mutation case suitable for mutation formula, less than 0 if mutation rate less than 0.2
        # Thus, if mutation rate is larger than 0.2, the formula would be misleading

        minall<-pmin(f1m1,f1m2,f2m1,f2m2,m1f1,m1f2,m2f1,m2f2)
        minf<-pmin(ac,ad,bc,bd)
        minm<-pmin(ec,ed,fc,fd)

        x<-as.double(f1m1==minall)*muAtoC*muTtoC*10^(2-ac2-ed2)/4+
          as.double(f1m2==minall)*muAtoC*muTtoC*10^(2-ac2-fd2)/4+
          as.double(f2m1==minall)*muAtoC*muTtoC*10^(2-bc2-ed2)/4+
          as.double(f2m2==minall)*muAtoC*muTtoC*10^(2-bc2-fd2)/4+
          as.double(m1f1==minall)*muAtoC*muTtoC*10^(2-ec2-ad2)/4+
          as.double(m1f2==minall)*muAtoC*muTtoC*10^(2-ec2-bd2)/4+
          as.double(m2f1==minall)*muAtoC*muTtoC*10^(2-fc2-ad2)/4+
          as.double(m2f2==minall)*muAtoC*muTtoC*10^(2-fc2-bd2)/4
        y<-as.double(ec==minm)*muTtoC*10^(1-ec2)+
          as.double(ed==minm)*muTtoC*10^(1-ed2)+
          as.double(fc==minm)*muTtoC*10^(1-fc2)+
          as.double(fd==minm)*muTtoC*10^(1-fd2)
        results<-data.frame(Log10CLR=log10(as.double(minall+minf+minm>10000)*muAtoC+
                         as.double(minf+minm+minall<10000)*x/y))
        }
    }
  return(results)
}


