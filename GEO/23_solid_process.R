library(affy)
#library(dplyr)
library(tibble)
library(tidyr)
#library(oligo)
library(ggplot2)
library(sva)
library(devtools)
library(ggloop)
library(magrittr)
library(purrr)
library(FactoMineR)
library(factoextra)
library(readxl)
library(conflicted)
##Step1:下载GEO,手动下载raw cel数据
##Step2：用oligo处理cel(affy不好用)，RMA标准化，要整批次处理
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

batch="GSE196723"
dir="D:/23 solid tumors from GEO/4_Breast Cancer/affy/"

##bioconductor无法下载GSE131418注释文件的解决办法
BiocManager::install("makecdfenv")
library("makecdfenv")
make.cdf.package("GPL23985_ADXECv1a520743.cdf.gz",'adxecv1a520743cdf',species = "Homo sapiens",compress = T)
install.packages("adxecv1a520743cdf",repos = NULL,type = "source")
library(adxecv1a520743cdf)
##Step3：注释基因symbol,需要手动下载注释文件
##Step4：python 完成merge
exp1<-data.table::fread("merge_hta_493_32670.txt")
exp1<-data.frame(exp1)
exp1.1<-aggregate(exp1[,-1],by=list(exp1[,1]),mean)
f<-data.frame(read_xlsx("batch.xlsx"))
x1<-prcomp(t(exp1.1[,-1]),scale. = TRUE)
df1<-as.data.frame(x1[["x"]][,1:2])
df1$group<-f$group
ggplot(df1,aes(x=PC1,y=PC2,color=group))+geom_point()
row.names(exp1.1)<-exp1.1$Group.1
combat_exp<-ComBat(dat=exp1.1[,-1], batch=f$group,par.prior=TRUE, prior.plots=FALSE)
x2<-prcomp(t(combat_exp),scale. = TRUE)
df2<-as.data.frame(x2[["x"]][,1:2])
df2$group<-f$group
ggplot(df2,aes(x=PC1,y=PC2,color=group))+geom_point()
data.table::fwrite(as.data.frame(combat_exp),file = 'merge_hta_combat_30905_491.txt',sep = '\t',row.names = T)