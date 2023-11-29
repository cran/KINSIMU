#' @title CLR for a single case
#' @description
#' CLR calculation for a single case, where two individual participated and labeled as A and B
#' @param AB Genotypes of two individuals in a case, which should be data.frame with 4 columns (2 for each individual) and \code{nl} rows, where \code{nl} stand for number of loci. The row names should be the name of each locus;
#' @param afmatrix name of allele frequency list, which can be loaded with "EvaluatePanel" function
#' @param rare a data frame containing the frequency of rare allele on the locus, with 1 row and multiple columns, each column for a marker;
#' @param allelename if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the afmatrix
#' @param stepPI If TRUE, empirical decreasing model of STR mutation would be taken when paternity index is needed to be calculated, otherwise, mutation rate would be taken as PI if IBS=0 between an alleged PC pair.
#' @param adelta3 distributions of the IBD coefficient of the outbred plaintiff's hypotheses in LR calculation, which should be a data.frame with 3 columns and x rows, where x stood for the number of such LR being calculated. The names of columns should be "k0", "k1" and "k2", and those of rows the name of LRs
#' @param adelta9 distributions of the Jacquard coefficient of the inbred plaintiff's hypotheses in LR calculation, which should be a data.frame with 9 columns and x rows, where x stood for the number of such LR being calculated. The names of columns should be "D1-D9", and those of rows the name of LRs
#' @param mu mutation rate when paternity index is needed to be calculated, defualts to 0.002.
#'
#' @return a list of two data.frames: (i) results_on_each_marker: multiple types of parameters calculated on each marker, including IBS and multiple log10 of LRs; (ii) total_results_of_the _case: the CIBS and log10CLR results for the whole case.
#' @export
#'
#' @examples
#' # example code
#' AB<-data.frame(a=rep(0,42),b=rep(0,42),c=rep(0,42),d=rep(0,42))
#' for (i in 1:42) {
#'   temp<-pairsimu(af = FortytwoSTR$afmatrix[[i]],ss = 1,delta = c(0,1,0),allelename = FALSE)
#'   AB[i,]=temp
#'   rownames(AB)[i]=names(FortytwoSTR$afmatrix)[i]
#' }
#' adelta3<-data.frame(k0=c(0,0.25,0.5),k1=c(1,0.5,0.5),k2=c(0,0.25,0),row.names = c("PC","FS","HS"))
#' adelta9<-data.frame(D1=0,D2=0,D3=0,D4=0,D5=0.25,D6=0,D7=0.25,D8=0.5,D9=0,row.names = "FIMCpair")
#' results<-logLR(AB=AB,afmatrix=FortytwoSTR$afmatrix,rare=FortytwoSTR$rare,stepPI=TRUE,
#' adelta3=adelta3,adelta9=adelta9)

logLR<-function(AB,afmatrix=NULL,rare=NULL,allelename=FALSE,stepPI=FALSE,
                adelta3=NULL,adelta9=NULL,mu=0.002){
  if (ncol(AB)!=4) {
    stop(paste("false in individual data"))
  }
  if ((is.null(afmatrix) || is.null(rare))) {
    if (isTRUE(stepPI) || nrow(as.data.frame(adelta3))+nrow(as.data.frame(adelta9))>0) {
      stop("Please input the frequency data of regular and rare alleles")
    }
  }
  if (any(duplicated(rownames(AB)))) {
    stop("Duplicate locus names in AB data frame")
  }
  if (any(duplicated(names(afmatrix)))) {
    stop("Duplicate locus names in afmatrix list")
  }
  if (any(duplicated(colnames(rare)))) {
    stop("Duplicate locus names in rare data frame")
  }
  if (!isTRUE(all(rownames(AB) %in% names(afmatrix))) || !isTRUE(all(rownames(AB) %in% colnames(rare)))) {
    stop("There are markers in AB data frame not included in afmatrix or rare data frame")
  }
  nl<-nrow(AB)
  para<-as.data.frame(matrix(0,nrow = nl,ncol = 4))
  colnames(para)<-c("ac","ad","bc","bd")
  para$ac<-as.double(AB[,1]==AB[,3])
  para$ad<-as.double(AB[,1]==AB[,4])
  para$bc<-as.double(AB[,2]==AB[,3])
  para$bd<-as.double(AB[,2]==AB[,4])
  markerresults<-data.frame(IBS=as.double(para$ac+para$ad+para$bc+para$bd>0)+
                              as.double(para$ac*para$bd+para$ad*para$bc>0))
  if (is.null(adelta3) & is.null(adelta9)) {
    warning("No alleged hypothesis was input, and only IBS results are output")
  } else {
    #pc and pd
    if (isTRUE(allelename)) {
      for (i in 1:nl) {
        para$pc[i]<-afmatrix[[rownames(AB)[i]]][as.character(AB[i,3]),]
        para$pd[i]<-afmatrix[[rownames(AB)[i]]][as.character(AB[i,4]),]
        if(is.na(para$pc[i])){
          para$pc[i]<-rare[[rownames(AB)[i]]]
        }
        if(is.na(para$pd[i])){
          para$pd[i]<-rare[[rownames(AB)[i]]]
        }
      }
    } else{
      for (i in 1:nl) {
        para$pc[i]<-afmatrix[[rownames(AB)[i]]][AB[i,3],]
        para$pd[i]<-afmatrix[[rownames(AB)[i]]][AB[i,4],]
        if(is.na(para$pc[i])){
          para$pc[i]<-rare[[rownames(AB)[i]]]
        }
        if(is.na(para$pd[i])){
          para$pd[i]<-rare[[rownames(AB)[i]]]
        }
      }
    }
    para$pc=as.numeric(para$pc)
    para$pd=as.numeric(para$pd)
    para$LRid<-(para$ac*para$bd+para$ad*para$bc)/(2*para$pc*para$pd)
    para$PInomu<-((para$ac+para$bc)/para$pc+(para$ad+para$bd)/para$pd)/4
    if (isTRUE(stepPI)) {#stepwise mutation model
      if (isTRUE(allelename)){
        #allele name
        para$acm=abs(AB[,1]-AB[,3])+abs(AB[,1]-AB[,3])%%1*100000
        para$adm=abs(AB[,1]-AB[,4])+abs(AB[,1]-AB[,4])%%1*100000
        para$bcm=abs(AB[,2]-AB[,3])+abs(AB[,2]-AB[,3])%%1*100000
        para$bdm=abs(AB[,2]-AB[,4])+abs(AB[,2]-AB[,4])%%1*100000
      } else {
        # allele position
        for (j in 1:nl) {
          an<-as.data.frame(as.numeric(row.names(afmatrix[[rownames(AB)[j]]])))
          para$acm[j]=abs(an[AB[j,1],]-an[AB[j,3],])+abs(an[AB[j,1],]-an[AB[j,3],])%%1*100000
          para$adm[j]=abs(an[AB[j,1],]-an[AB[j,4],])+abs(an[AB[j,1],]-an[AB[j,4],])%%1*100000
          para$bcm[j]=abs(an[AB[j,2],]-an[AB[j,3],])+abs(an[AB[j,2],]-an[AB[j,3],])%%1*100000
          para$bdm[j]=abs(an[AB[j,2],]-an[AB[j,4],])+abs(an[AB[j,2],]-an[AB[j,4],])%%1*100000
        }
      }
      para$steps<-pmin(para$acm,para$adm,para$bcm,para$bdm)
      para$dc<-as.double(para$acm==para$steps)+as.double(para$bcm==para$steps)
      para$dd<-as.double(para$adm==para$steps)+as.double(para$bdm==para$steps)
      para$PImu<-as.double(para$steps==0)*(para$dc/para$pc+para$dd/para$pd)/4+ # ignore mutation if IBS>0
        as.double(para$steps>10000)*mu+ # take mutation rate as PI if there is no integer step of mutation
        as.double(para$steps>0 & para$steps<10000)*mu*10^(1-para$steps)*
        (para$dc/para$pc+para$dd/para$pd)/8 # step wise model of mutation
    } else {
      para$PImu=para$PInomu
      para$PImu[para$PInomu==0]<-mu
    }
    if (is.null(adelta9)) {
      if (any(apply(adelta3,1,sum)!=1) || any(adelta3<0)) {
        stop("Error in alleged IBD coefficient")
      }
      for (j in 1:nrow(adelta3)) {
        markerresults[,j+1]<- log10(adelta3$k2[j]*para$LRid+
                                      as.numeric(adelta3$k1[j]==1)*adelta3$k1[j]*para$PImu+
                                      as.numeric(adelta3$k1[j]<1)*adelta3$k1[j]*para$PInomu+
                                      adelta3$k0[j])
      }
      colnames(markerresults)[2:(nrow(adelta3)+1)]=paste("LR_",rownames(adelta3),sep = "")
    } else if(is.null(adelta3)) {
      if (any(apply(adelta9,1,sum)!=1) || any(adelta9<0)) {
        stop("Error in alleged Jacquard coefficient")
      }
      if (sum(adelta9$D1,adelta9$D2,adelta9$D3,adelta9$D4,adelta9$D5,adelta9$D6)==0){
        # take as adelta3
        for (j in 1:nrow(adelta9)) {
          markerresults[,j+1]<- log10(adelta9$D7[j]*para$LRid+
                                        as.numeric(adelta9$D8[j]==1)*adelta9$D8[j]*para$PImu+
                                        as.numeric(adelta9$D8[j]<1)*adelta9$D8[j]*para$PInomu+
                                        adelta9$D9[j])
        }
      } else {
        #calculate of F1-F6
        para$ab<-as.double(AB[,1]==AB[,2])
        para$cd<-as.double(AB[,4]==AB[,3])
        if (isTRUE(allelename)) {
          for (i in 1:nl) {
            para$pa[i]<-afmatrix[[rownames(AB)[i]]][as.character(AB[i,1]),]
            if(is.na(para$pa[i])){
              para$pa[i]<-rare[[rownames(AB)[i]]]
            }
          }
        } else{
          for (i in 1:nl) {
            para$pa[i]<-afmatrix[[rownames(AB)[i]]][AB[i,1],]
            if(is.na(para$pa[i])){
              para$pa[i]<-rare[[rownames(AB)[i]]]
            }
          }
        }
        para$pa=as.numeric(para$pa)
        para$FD1<-para$ab*para$ac*para$ad/(para$pa^3)
        para$FD2<-para$ab*para$cd/(para$pa*para$pc)
        para$FD3<-para$ab*(para$ac+para$ad)/(2*para$pa^2)
        para$FD4<-para$ab/para$pa
        para$FD5<-(para$ac+para$bc)*para$cd/(2*para$pc^2)
        para$FD6<-para$cd/para$pc
        for (j in 1:nrow(adelta9)) {
          markerresults[,j+1]<-log10(adelta9$D1[j]*para$FD1+
                                 adelta9$D2[j]*para$FD2+
                                 adelta9$D3[j]*para$FD3+
                                 adelta9$D4[j]*para$FD4+
                                 adelta9$D5[j]*para$FD5+
                                 adelta9$D6[j]*para$FD6+
                                 adelta9$D7[j]*para$LRid+
                                 as.numeric(adelta9$D9[j]==0)*adelta9$D8[j]*para$PImu+
                                 as.numeric(adelta9$D9[j]>0)*adelta9$D8[j]*para$PInomu+
                                 adelta9$D9[j])
          }
      }
      colnames(markerresults)[2:(nrow(adelta9)+1)]=paste("LR_",rownames(adelta9),sep = "")
    } else {
      if (any(apply(adelta9,1,sum)!=1) || any(adelta9<0)) {
        stop("Error in alleged Jacquard coefficient")
      }
      if (any(apply(adelta3,1,sum)!=1) || any(adelta3<0)) {
        stop("Error in alleged IBD coefficient")
      }
      for (j in 1:nrow(adelta3)) {
        markerresults[,j+1]<- log10(adelta3$k2[j]*para$LRid+
                                      as.numeric(adelta3$k1[j]==1)*adelta3$k1[j]*para$PImu+
                                      as.numeric(adelta3$k1[j]<1)*adelta3$k1[j]*para$PInomu+
                                      adelta3$k0[j])
      }
      colnames(markerresults)[2:(nrow(adelta3)+1)]=paste("LR_",rownames(adelta3),sep = "")
      if (sum(adelta9$D1,adelta9$D2,adelta9$D3,adelta9$D4,adelta9$D5,adelta9$D6)==0){
        # take as adelta3
        for (j in 1:nrow(adelta9)) {
          markerresults[,j+nrow(adelta3)+1]<- log10(adelta9$D7[j]*para$LRid+
                                        as.numeric(adelta9$D8[j]==1)*adelta9$D8[j]*para$PImu+
                                        as.numeric(adelta9$D8[j]<1)*adelta9$D8[j]*para$PInomu+
                                        adelta9$D9[j])
        }
      } else {
        #calculate of F1-F6
        para$ab<-as.double(AB[,1]==AB[,2])
        para$cd<-as.double(AB[,4]==AB[,3])
        if (isTRUE(allelename)) {
          for (i in 1:nl) {
            para$pa[i]<-afmatrix[[rownames(AB)[i]]][as.character(AB[i,1]),]
            if(is.na(para$pa[i])){
              para$pa[i]<-rare[[rownames(AB)[i]]]
            }
          }
        } else{
          for (i in 1:nl) {
            para$pa[i]<-afmatrix[[rownames(AB)[i]]][AB[i,1],]
            if(is.na(para$pa[i])){
              para$pa[i]<-rare[[rownames(AB)[i]]]
            }
          }
        }
        para$pa=as.numeric(para$pa)
        para$FD1<-para$ab*para$ac*para$ad/(para$pa^3)
        para$FD2<-para$ab*para$cd/(para$pa*para$pc)
        para$FD3<-para$ab*(para$ac+para$ad)/(2*para$pa^2)
        para$FD4<-para$ab/para$pa
        para$FD5<-(para$ac+para$bc)*para$cd/(2*para$pc^2)
        para$FD6<-para$cd/para$pc
        for (j in 1:nrow(adelta9)) {
          markerresults[,j+nrow(adelta3)+1]<-log10(adelta9$D1[j]*para$FD1+
                                 adelta9$D2[j]*para$FD2+
                                 adelta9$D3[j]*para$FD3+
                                 adelta9$D4[j]*para$FD4+
                                 adelta9$D5[j]*para$FD5+
                                 adelta9$D6[j]*para$FD6+
                                 adelta9$D7[j]*para$LRid+
                                 as.numeric(adelta9$D9[j]==0)*adelta9$D8[j]*para$PImu+
                                 as.numeric(adelta9$D9[j]>0)*adelta9$D8[j]*para$PInomu+
                                 adelta9$D9[j])
        }
      }
      colnames(markerresults)[(nrow(adelta3)+2):(nrow(adelta3)+nrow(adelta9)+1)]=
        paste("LR_",rownames(adelta9),sep = "")
    }
  }
  caseresults<-as.data.frame(colSums(markerresults))
  colnames(caseresults)<-"Results"
  results<-list(markerresults,caseresults)
  names(results)<-c("results_on_each_marker","total_results_of_the_case")
  return(results)
}
