setwd('~/TCGA')
library(dplyr)
tcga = readr::read_tsv('./tcga_RSEM_gene_fpkm.txt')
colnames(tcga)=c('id',colnames(tcga)[-1])
samplenames = colnames(tcga)[-1]
id2gene= readr::read_tsv('./ID_CONVERT_probeMap_gencode.v23.annotation.gene.probemap')
tcga = left_join(id2gene,tcga,'id')
tcga = na.omit(tcga)
tcga = tcga[,c('gene',samplenames)]
readr::write_csv(tcga,'tcga_RSEM_genename_fpkm.csv')

tcga_merge= aggregate(.~gene,tcga,max)
features = readr::read_tsv('../LJX/st/ST_NR_scRNA_main.txt')
readr::write_csv(tcga,'tcga_RSEM_genename_fpkm.csv')

