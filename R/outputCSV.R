#' @title Output frequency data into ISFG format
#' @param data frequency data, a list containing 4 data frames.
#' @param strpath name of the pathway of the output csv file
#'
#' @return a .csv file in ISFG formate
#' @export
#'
#' @examples
#' path<-tempdir()
#' outputCSV(FortytwoSTR,file.path(path,'data.csv'))
#' @importFrom utils write.csv

outputCSV<-function(data,strpath){
  if (!is.list(data) || length(data)!=4 ||
      isFALSE(all(names(data) %in% c("afmatrix","rare","indicators","panelparas"))) ||
      nrow(data$indicators)!=13) {
    stop("please input correct data")
  }
  nl<-length(data$afmatrix)
  allelenames<-as.data.frame(matrix(0,nrow = 0,ncol = 1))
  for (i in 1:nl) {
    temp<-as.data.frame(row.names(data$afmatrix[[i]]))
    allelenames<-rbind(allelenames,temp)
  }
  allelenames<-as.numeric(unique(allelenames)[,1])
  allelenames<-sort(allelenames)

  results<-as.data.frame(matrix(0,nrow = length(allelenames)+1,ncol = nl+1))
  colnames(results)=c("Allele", colnames(data$rare))
  results[1:length(allelenames),1]=allelenames
  results[(length(allelenames)+1),1]="N"
  results[(length(allelenames)+1),2:(nl+1)]=data$indicators[13,]
  for (i in 1:nl) {
    freq<-data$afmatrix[[i]]
    results[(1:length(allelenames)),(i+1)]=freq[as.character(results[1:length(allelenames),1]),]
  }
  results[is.na(results)]=0
  write.csv(results,file = strpath,row.names = FALSE)
}
