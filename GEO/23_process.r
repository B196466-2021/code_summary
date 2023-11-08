#! /usr/bin/Rscript

library(affy)
library(dplyr)
library(tibble)
library(tidyr)
library(oligo)
library(ggplot2)
library(sva)
library(magrittr)
library(purrr)
library(FactoMineR)
library(factoextra)
library(readxl)
library(conflicted)
library(R.utils)

batch="GSE20181_RAW"
platform="96"

# Set the folder path to extract
folder_path <- my_path

# Get all. gz files in the folder
gz_files <- list.files(path = folder_path, pattern = "\\.gz$")

# Loop through each. gz file in the folder and extract them
for (gz_file in gz_files) {
  # Construct the extracted file name
  extracted_file <- sub("\\.gz$", "", gz_file)
  
  # unzip
  gunzip(file = file.path(folder_path, gz_file), destname = file.path(folder_path, extracted_file))
}

#files <- dir(pattern="gz$")
#sapply(files, gunzip)

##Step 1: Download GEO and manually download raw cel data
##Step 2: Use oligo to process cel (affy is not easy to use), standardize RMA, and process the entire batch
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

GPL<-"GPL16570"
expr_df <- read.table("C:/Users/11042/Desktop/GSE136157_series_matrix_no_anno.txt",header=TRUE) 


class(expr_df) # dataframe
dim(expr_df)
expr_df[1:3,]

sample_names <- colnames(expr_df)[-1]
probe_ids <- expr_df$ID_REF

library(GEOquery)
gpl <- getGEO(GPL, destdir=".") 
colnames(Table(gpl))
head(Table(gpl)[,c(1,10)])
ID2symbol=Table(gpl)[,c(1,10)]

id_table <- read.table("C:/Users/11042/Desktop/GPL16570-1802_no_anno.txt",header=TRUE,sep = "\t",comment.char = "#",fill=TRUE)
id_table[1:11,1:2]
colnames(id_table)
 
# install.packages("data.table")
require(data.table)
 
probe2symbol <- id_table[,c("ID","gene_assignment")]
 
#Obtaining a gene symbol, one probe may correspond to multiple genes.
symbol <-tstrsplit(id_table$gene_assignment, "//", fixed=TRUE)[[2]]
#Remove spaces
symbol<- trimws(symbol, which = c("both", "left", "right"),whitespace = "[ \t\r\n]")
probe2symbol["symbol"] <- symbol
#Remove gene_ Assignment column
probe2symbol <-  probe2symbol[,c("ID","symbol")]
 
head(probe2symbol)
#Change the ID column name to probe_ ID
colnames(probe2symbol) <- c("probe_id","symbol")

library(mogene10sttranscriptcluster.db)
## Bimap interface:
x <- mogene10sttranscriptclusterSYMBOL
# Get the probe identifiers that are mapped to a gene symbol 
mapped_probes <- mappedkeys(x)
# Convert to dataframe
probe2symbol2 <- as.data.frame(x[mapped_probes])

merged_expr_df <- merge(x = expr_df, y = probe2symbol, by.x = "ID_REF",
                        by.y = "probe_id", all.x= TRUE)
#Remove probe_ Line with ID without corresponding gene symbol
filt_expr_df <- merged_expr_df[complete.cases(merged_expr_df),]
#Remove ID_ REF column (probe id)
filt_expr_df <- subset(filt_expr_df, select = -ID_REF)
#Take the average or maximum value of all probes for each gene as the expression level of the gene
m_df <- aggregate(.~symbol,data=filt_expr_df,mean)
m_df <- aggregate(.~symbol,data=filt_expr_df,max)
 
rownames(m_df) <- m_df$symbol
 
m_df <- subset(m_df, select = -symbol) #去掉symbol列
exprSet <- as.matrix(m_df)
head(exprSet)
write.table(exprSet,file="GSE136157.txt",sep="\t")