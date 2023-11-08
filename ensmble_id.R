#! /usr/bin/Rscript

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
gene<-as.character(gene_list$Symbol)
ensembls <- mapIds(org.Hs.eg.db, keys = gene, keytype = "SYMBOL", column="ENSEMBL")
ensembls
ensembls<-as.data.frame(ensembls)