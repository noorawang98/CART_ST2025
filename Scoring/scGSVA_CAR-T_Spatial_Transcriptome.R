#-------ECM associated spatial distribution-------------
setwd('/home/wr/LJX/st')
library(Seurat)

L13_integrated = readRDS('./L13_integrated.RDS')

library(GSVA)
library(tidyverse)
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
GOID <- c("GO:0030198")

# GOgeneID <- get(GOID, org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()

GOgeneID <- get(GOID, org.Mm.egGO2ALLEGS) %>% mget(org.Mm.egSYMBOL)%>%unlist()

gsea.gset = read.gmt('~/LJX/mm/GMT_ECM.Mm/ECMgenesets.v2024.1.Mm.gmt')

gset = split(gsea.gset$gene,gsea.gset$term)
# list('ECM organization(GO:0030198)'=names(GOgeneID))

gset$`ECM organization(GO:0030198)`=names(GOgeneID)

mat= as(object = L13_integrated[["RNA"]], Class = "Assay")
mat = as.matrix(L13_integrated@assays$RNA@counts)
gsva.mat = gsva(mat,gset)
write.csv(gsva.mat,'ECM_GSVA_L13.csv')
