library(Seurat)
library(dplyr)

setwd("~/LJX/sc_in")

immune.combined =readRDS("immune.combined.RDS")

sampleinfo = immune.combined@meta.data

sampleinfo = sampleinfo[!is.na(sampleinfo$sample),]

immune.combined=subset(immune.combined,cells = rownames(sampleinfo))

immune.combined$sample = factor(immune.combined$sample,levels=c("N9","R5"),
                                labels = c("Non-Responder","Responder"))

DefaultAssay(immune.combined) <- "RNA"

s.gene = readr::read_csv("~/LJX/mm/mm.s.cellcyclegenes.csv")
s.genes= s.gene$x
g2m.gene =readr::read_csv("~/LJX/mm/mm.g2m.cellcyclegenes.csv")
g2m.genes = g2m.gene$x


immune.combined <- CellCycleScoring(immune.combined, 
                                    s.features = s.genes, 
                                    g2m.features = g2m.genes, 
                                    set.ident = TRUE)

write.csv(immune.combined@meta.data,"immune.combined.metadata.csv")
immune.combined = RunTSNE(immune.combined)
saveRDS(immune.combined,"immune.combined.RDS")
