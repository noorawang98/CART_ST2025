library(CIBERSORT)
library(dplyr)
library(tidyverse)
# # sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
# # mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
# setwd('~/LJX/st/rds')
# #=====preprocessed scaled TCGA.fpkm==========
# tcga = readr::read_csv('~/TCGA/tcga_sel.csv')
# minvalue = min(as.matrix(tcga[,-1]),na.rm=T)
# if(minvalue<0){
#   offset = abs(minvalue)
#   tcga[,-1] =tcga[,-1]+offset
# }
# 
# input = tcga%>%mutate(across(-gene,~log2(.+1)))
# 
# input = input%>%group_by(gene)%>%summarise(across(everything(),mean,na.rm=T)) %>% ungroup
# 
# write.table(input,'/home/wr/TCGA/tcga_sel.txt',sep='\t',row.names = F,quote=F)
#============= run CIBERSORT==============
mixture_file = '/home/wr/TCGA/tcga_sel.txt'
sig_matrix = '/home/wr/LJX/st/CIBERSORTx_Job22_ST_NR_mm2hsa.txt'
results <- cibersort(sig_matrix, mixture_file)
write.csv(results,'CIBSERORT_st_TCGA.csv')