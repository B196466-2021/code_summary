#!/usr/bin/env Rscript

library(dplyr)
library(tibble)
library(tidyr)
library(oligo)
library(ggplot2)
library(devtools)
library(magrittr)
library(purrr)
library(FactoMineR)
library(factoextra)
library(readxl)
library(conflicted)
library(R.utils)


batch="GSE25065_RAW"
platform="96"
dir="C:/Users/11042/Desktop/BRCA_neoadjuvant/GSE25065/GSE25065_RAW"

geo_process<-function(batch,dir,gpl){
  t<-paste0(dir,batch)
  cel<-list.celfiles(t,full.names=TRUE)
  raw<-read.celfiles(cel)
  gc()
  eset<-exprs(rma(raw))
  gc()
  f=paste0(substr(batch,0,10),"_",gpl,"_RMA_log2intensity_",dim(eset)[2],"_",dim(eset)[1],".txt")
  data.table::fwrite(as.data.frame(eset),file=paste0(substr(t,0,nchar(t)-3),f),sep = '\t',row.names = T,col.names = T)
}
geo_process(batch,dir,platform)
