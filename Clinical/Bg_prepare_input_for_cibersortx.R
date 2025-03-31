library(dplyr)
library(tidyverse)
library(tibble)
# setwd('/home/Data/卫健委项目/公共数据/GSE192742/')
setwd('/home/Data/卫健委项目/公共数据/GSE192742/scRNA/')
# # library(rhdf5)
# # h5filelist = grep('.h5',dir('./'),value=T)
# txtfilelist = grep('.txt',dir('./'),value=T)
# 
# metadata = read.csv('./annot_mouseStStAll.csv')
# cell.immune = metadata
# cell.immune$barcode=cell.immune$cell
# 
# cell.immune = separate(cell.immune,col = 'barcode',into = c('barcode_seq','barcode_sample'),sep = '-')
# cell.immune$barcodeid = paste0(cell.immune$barcode_seq,'-1')
# cell.immune = cell.immune[grep('scRnaSeq',cell.immune$typeSample),]
# 
# # df.h5filelist=lapply(h5filelist,function(i)dior::read_h5(i))
# # lapply(c(1,4,5),function(i){
# #   df = readr::read_tsv(txtfilelist[i])
# #   df = aggregate(.~...1,data=df,mean)
# #   genes = df$...1
# #   cells = intersect(cell.immune$barcodeid,colnames(df))
# #   df=df[,cells]%>%as.data.frame()
# #   rownames(df)=genes
# #   df= rownames_to_column(df,var='Gene')
# #   readr::write_tsv(df,txtfilelist[i])
# # })
# 
# library(Seurat)
# df.mat=readr::write_tsv("../MICA.exp.txt")
# cells = intersect(cells,cell.immune$barcodeid)
# cell.meta = cell.immune[cell.immune$barcodeid%in%cells,]
# cell.meta = data.frame('barcode'=cell.meta$barcodeid,'Celltype'=cell.meta$annot)
# readr::write_tsv(cell.meta,"../MICA.meta.txt")
# df.mat = data.frame(df.mat,row.names = 1)
# colnames(df.mat)=gsub('\\.','-',colnames(df.mat))
# cell.meta.noreplicate =cell.meta[!duplicated(cell.meta$barcode),]
# rownames(cell.meta.noreplicate)=cell.meta.noreplicate$barcode
# cell.meta.noreplicate = cell.meta.noreplicate[colnames(df.mat),]
# all.equal(rownames(cell.meta.noreplicate),colnames(df.mat))
# obj = CreateSeuratObject(df.mat,meta.data = cell.meta.noreplicate)
# obj = obj %>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures()
# 
# cell.meta = cell.immune[cell.immune$barcodeid%in%cell.meta$barcode,]
# cell.meta  =cell.meta[!duplicated(cell.meta$barcodeid),]
# rownames(cell.meta) = cell.meta$barcodeid
# cell.meta = cell.meta[rownames(cell.meta.noreplicate),]
# obj = AddMetaData(obj,cell.meta)
# obj = RunPCA(obj)%>%FindNeighbors(dims=1:30)%>%RunUMAP(dims=1:30)
# #check the batch effect
# DimPlot(obj,group.by = "sample",label = T,repel = T)+NoLegend()
# DimPlot(obj,group.by = "Celltype",label = T,repel = T)+NoLegend()
# 
# df.mat = as.matrix(df.mat)
# features = VariableFeatures(obj,n=2000)
# df.mat.hvg = df.mat[features,]%>%as.matrix()
# # write.table(df.mat.hvg,'../MICA.scRNA.txt',sep = '\t',row.names=T，col.names=T)

# saveRDS(obj,'../MICA.scRNAobj.rds')

obj = readRDS('../MICA.scRNAobj.rds')
df.mat = GetAssayData(obj,assay = 'RNA',layer = 'counts')
df.mat = as.data.frame(df.mat)
df.mat$Gene=rownames(df.mat)
df.mat = df.mat[,c(ncol(df.mat),1:(ncol(df.mat)-1))]
features =VariableFeatures(obj,n=3000)
df.mat.hvg = df.mat[c('Ccr1','Spp1','Fap',features),]

#mouse2human for tcga 
library(biomaRt)
ensemble = useEnsembl(biomart = 'genes',mirror = 'asia')
# mouse = useDataset(dataset = 'mmusculus_gene_ensembl',mart=ensemble)
# human= useDataset(dataset = 'hsapiens_gene_ensembl',mart=ensemble)
features =rownames(df.mat.hvg)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mus2hsa = getLDS(attributes = c('mgi_symbol'),
                 filters='mgi_symbol',
                 values=features,mart=mouse,
                 attributesL = c('hgnc_symbol'),
                 martL = human,uniqueRows = T)


library()

cell.meta = data.frame(barcode=obj$barcodeid,celltype=obj$Celltype)
cell.meta =data.frame(obj@active.ident)
cell.meta =data.frame(barcode = rownames(cell.meta),celltype=paste0('NR_',cell.meta$obj.active.ident))
cell.meta.immune = cell.meta[!cell.meta$celltype%in%c('Hepatocytes','Fibroblasts','Endothelial cells','Cholangiocytes'),]

df.mat.hvg= df.mat.hvg[,c('Gene',cell.meta$barcode)]
df.mat.hvg = rbind(c('Gene',as.character(obj$Annotaion_MainCelltype)),df.mat.hvg)
df.mat.hvg = rbind(c('Gene',as.character(cell.meta$celltype)),df.mat.hvg)
readr::write_tsv(df.mat.hvg,'../MICA_scRNA_main.txt',col_names = F)
readr::write_tsv(df.mat.hvg,'../ST_NR_scRNA_main.txt',col_names = F)
# readr::write_tsv(df.mat.hvg,'../MICA.exp.txt')
# readr::write_tsv(cell.meta.immune,'../MICA.meta.txt')

markers= c(
  'Alb','Cyp2e1','Apoa1',#Hep
  'Cd3d','Cd3e','Cd4','Cd8a','Cd8b1',#T
  'Il7r','Sell',#Naive
  'Trdc',"Il17a" ,#
           'Foxp3','Il2ra',#Treg
           'Nkg7','Klrk1','Klrb1c',#NK
           'Ncr1','Tbx21','Ifng','Il12rb1',#ILC1
           'Cd19','Cd79b','Cd79a','Ms4a1',#B
           'Adgre1','Cd14','Itgam',#Monocyte
           'Marco','Cd68','Vsig4',#KC
           'Csf1r','C1qa','Il1b','Cd86', #M1
           'Cd163','Arg1','Mrc1',#M2
           'Csf3r','S100a8','S100a9',#Neu
           'Il3ra','Hdc','Fcer1a',#Bas
           'Itgax','Xcr1','Batf3','Irf8','Cd24a',#cDC1
           'Zbtb46','Sirpa','Irf4',#cDC2
           'Siglech' #pDC
           )

DotPlot(obj,markers,cluster.idents = T)+coord_flip()+
  scale_color_gradientn(colours = rev(rcolors$RdBu))

rename = c('0'='KC2-like',
           '1'='M2-like',
           '2'='pDCs',
           '3'='CD4T Cells',
           '4'='KC2-like',
           '5'='CD8T Cells',
           '6'='Mono-derived',
           '7'='ILC1',
           '8'='B Cells',
           '9'='Neutrophils',
           '10'='pDCs',
           '11'='CD4T Cells',
           '12'='Mono-derived',
           '13'='M2-like',
           '14'='cDC1',
           '15'='NK Cells',
           '16'='cDC2',
           '17'='CD8T Cells',
           '18'='Mono-derived',
           '19'='pDCs',
           '20'='Mono-derived',
           '21'='CD4T Cells',
           '22'='Mono-derived',
           '23'='Monocytes',
           '24'='CD8T Cells',
           '25'='CD8T Cells',
           '26'='Neutrophils',
           '27'='cDC1',
           '28'='NK Cells',
           '29'='Mono-like Macrophage',
           '30'='cDC2',
           '31'='B Cells',
           '32'='unknown',
           '33'='CD4T Cells',
           '34'='Mig_DC',
           '35'='M1-like',
           '36'='Neutrophils',
           '37'='Hepatocytes',
           '38'='M2-like',
           '39'='Basophils',
           '40'='KC1-like'
           )
obj@active.ident = obj$seurat_clusters
obj@active.ident = factor(obj$seurat_clusters,levels=0:40)
obj = RenameIdents(obj,rename)
obj = AddMetaData(obj,obj@active.ident,col.name = 'Annotaion_MainCelltype')
DimPlot(obj,label = T)
write.csv(obj@meta.data,'../cell.meta.csv')


