#reanalysis of scRNA
setwd('~/LJX/code/celescope/shell_rna/')
setwd('~/LJX/sc_in/adata_h5ad/filter_adata_20240420/')
devtools::install_github("cellgeni/schard")
devtools::install_github('immunogenomics/harmony', force = TRUE)

library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)


ad_path = 'Control/outs/rna.h5ad'
# adata= schard::h5ad2seurat(ad_path)

library(Seurat)
library(rhdf5)
source('~/LJX/code/schard-main/R/functions.R')
source('~/LJX/code/schard-main/R/h5ad_util.R')
# query = h5ad2seurat('./query.h5ad')
adata = h5ad2seurat('./adata_concat1.h5ad')
control = h5ad2seurat('./../Control.h5ad')

adata_concat = h5ad2seurat('./adata_concat1.h5ad')

VlnPlot(adata_concat,features = c('n_genes_by_counts','pct_counts_mt','doublet_score'),pt.size=0,group.by = 'batch')
adata=subset(adata_concat,subset=(n_genes_by_counts<5000&pct_counts_mt<10&doublet_score<0.3))

adata=subset(adata,subset=(n_genes_by_counts<6000&pct_counts_mt<10&doublet_score<0.3))



VlnPlot(adata,features = c('n_genes_by_counts','pct_counts_mt','doublet_score'),pt.size=0,group.by = 'batch')



adata = h5ad2seurat('./adata_concat_scvi_annotated_2500.h5ad')
adata = h5ad2seurat('./adata_concat_scvi_annotated_2500_noexo.h5ad')

adata_concat = subset(adata_concat, cells = names(adata@active.ident))

adata_concat@reductions = adata@reductions

adata_concat = adata_concat %>% FindNeighbors(reduction='XscVI_',dims = 1:30)%>%RunUMAP(reduction='XscVI_',dims = 1:30)

adata_concat = AddMetaData(adata_concat,adata$majority_voting,'majority_voting')

adata_concat = AddMetaData(adata_concat,adata$broad_celltype,'broad_celltype')

DimPlot(adata_concat,reduction = 'umap',group.by = 'majority_voting')



batchcolor=c('lightgrey','#D3706C','#60C7C4')

VlnPlot(adata,features = c('n_genes_by_counts','pct_counts_mt','doublet_score'),
        group.by = 'batch',pt.size = 0,cols = batchcolor)

# ggsave('Vlnplot_sc_CD45_immune.pdf',p,width =6,height = 3 )

#----add SCANVI obsm to adata---------------
# query = query %>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures()%>%RunPCA()

# adata =subset(adata,cells= rownames(query@meta.data))
# 
# adata = AddMetaData(adata,query$predictions_scanvi,'predictions_scanvi')
# 
# VlnPlot(adata,features = c('nFeature_RNA','nCount_RNA'),
#           group.by = 'batch',pt.size = 0,cols = batchcolor)


#--------RunHarmony to remove batch effect--------------
adata = adata %>%NormalizeData()%>%ScaleData()%>%
  FindVariableFeatures(nfeatures = 2000)%>%RunPCA()

#use harmony remove the batch effect
adata <- RunHarmony(adata, "batch")

adata <- adata %>% 
  RunHarmony("batch", 
             plot_convergence = TRUE, 
             nclust = 50, max_iter = 10, early_stop = F)


p1 <- DimPlot(object = adata, reduction = "harmony", pt.size = .1, group.by = "batch")
p2 <- VlnPlot(object = adata, features = "harmony_1", group.by = "batch",  pt.size = .1)
plot_grid(p1,p2)

DimHeatmap(object =adata, reduction = "harmony", cells = 500, dims = 1:3)

adata = RunUMAP(adata,reduction = 'harmony',dims=1:20)

p1=DimPlot(adata,reduction = 'umap',
           # group.by = 'leiden',
           label = T,repel = T)+theme_map()
p2=DimPlot(adata,reduction = 'umap',group.by = 'batch',cols = batchcolor)+theme_map()
# p3=DimPlot(adata,reduction = 'umap',group.by = 'majority_voting',label = T,repel = T)+theme_map()
# p4=DimPlot(adata,reduction = 'umap',group.by = 'broad_celltype',label = T,repel = T)+theme_map()
plot_grid(p1,p3,p2,p4,ncol = 2)

plot_grid(p1,p2,ncol = 2)

adata <- adata %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 1)
# %>%RunUMAP(dims=1:30,reduction='harmony')
# 
# adata <-adata %>%
#   # RunUMAP(reduction.model = "harmony")%>%
#   RunTSNE(reduction = 'harmony')

# 
# p1 <- DimPlot(adata, reduction = "tsne", group.by = "batch", pt.size = .1)
# p2 <- DimPlot(adata, reduction = "tsne", label = TRUE, pt.size = .1)
# plot_grid(p1, p2)

# DimPlot(adata, reduction = "umap", group.by = "batch", pt.size = .1)

?RunUMAP
# adata=RunUMAP(adata,reduction.model = "harmony")

adata@reductions$XscANVI_ = query@reductions$XscANVI_

adata = adata %>%FindNeighbors(dims=1:30)%>%
  RunUMAP(dims=1:30)%>%FindClusters()

FindNeighbors(dims=1:10,reduction ='XscANVI_' )%>%
  RunUMAP(dims=1:10)%>%FindClusters()


DotPlot(adata,assay = 'RNA',
        # group.by = 'leiden',
        # group.by = 'predictions_scanvi',
        features = c('mCherry-CAR','EGFP', 'Acta2',
                               'Cd19','Cd79a','Cd79b','Ighm', 'Ms4a1',
                     'Hba-a1','Hbb-bt',
                     'Dcn','Col1a1','Vwf',
                               'S100a8','S100a9','Csf1r','Csf3r',
                               'Adgre1','Il1b','Mrc1','C1qc',
                     'Itgax','Xcr1','Cd63','Cst3','Irf8',
                               'Cd3e','Cd8a','Cd8b1','Cd4','Il7r','Mki67','Klrk1','Klrb1c','Klra9'),
        cluster.idents = T)+scale_color_gradient2(low = 'blue',mid = 'lightgrey',high = 'red')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+xlab('')+ylab('')|DimPlot(adata,group.by = 'seurat_clusters',reduction = 'umap',label = T)



DimPlot(adata,group.by ='seurat_clusters',label = T,reduction = 'umap')|DimPlot(adata,group.by ='batch',cols = batchcolor,label = T,reduction = 'umap')

table(adata$batch,adata$seurat_clusters)

p1=DimPlot(adata,group.by = 'batch',reduction = 'umap',cols = batchcolor)
p2 = DimPlot(adata,group.by = 'seurat_clusters',reduction='umap',label = T,repel = T)
plot_grid(p1, p2)


DimPlot(query,group.by = 'leiden',label = T)

# cells =names(adata@active.ident)[!adata@active.ident%in%c('3','24','10')]

# cells =names(adata@active.ident)[!adata@active.ident%in%c('3','22')]


adata = h5ad2seurat('./adata_concat1.h5ad')

adata = subset(adata,cells =cells)

adata = adata %>%NormalizeData()%>%ScaleData()%>%
  FindVariableFeatures(nfeatures = 2000)%>%RunPCA()

adata = adata %>%ScaleData()%>%FindVariableFeatures(nfeatures = 2000)%>%RunPCA()


adata = adata %>%FindNeighbors(dims=1:30)%>%FindClusters(resolution = 1)

# adata = adata%>%FindNeighbors(dims=1:30)%>%
#   RunUMAP(dims=1:30)%>%FindClusters()



DotPlot(adata,features = c('EGFP', 'Acta2','Epcam','Krt18',
                           'Col1a1','Dcn','Vwf','Eng',
                           'Cd19','Cd79a','Ms4a1', 
                           'Cd247','Cd3e','Cd4','Cd8b1',
                           'Ccr7','Il7r','Trdc','Trgc1','Mki67',
                           'Klrk1', 'Klra9',
                           'S100a8','S100a9','Csf3r',
                           'Adgre1','Cd14',
                           'C1qc','C1qb','Csf1r',
                           'Itgax','Xcr1','Clec4n','H2-Aa','Siglech'),cluster.idents = T)+
  scale_color_gradientn(colours =hcl.colors(palette = 'Spectral',n=22,rev=T) )+
  # scale_color_gradient2(low = 'blue',mid = 'lightgrey',high = 'red')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+xlab('')+ylab('')+coord_flip()


P=DotPlot(adata,features =c('Epcam','Krt18',
                          'EGFP', 'Acta2',
                          'Col1a1','Dcn',
                          
                          'Cd19','Cd79a',
                          'Trdc','Trgc1',
                          
                          'Cd3e','Cd2','Cd4','Cd8b1',
                          'Klrk1', 'Klra9',
                          'Lyz2','Cd14'
                          ),cluster.idents = F)+
  scale_color_gradientn(colours =hcl.colors(palette ='RdBu',n=22,rev=T) )+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+xlab('')+ylab('')+coord_flip()

df= P$data

adata.split = SplitObject(adata)

# fibro = adata.split$



p=ggplot(df[df$pct.exp>20&df$avg.exp.scaled>1.1,],
       aes(x=id,y=features.plot,fill=avg.exp.scaled,size=pct.exp))+
  geom_point(shape=21)+scale_fill_gradientn(colors = rcolors::rcolors$OrRd)+
  scale_size_continuous(range = c(4,6))+
  ggpubr::theme_pubr(legend='right')+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+ylab('')
p

ggsave('DotPlot_maintype.pdf',p,width=5,height = 5)

adata@active.ident = factor(adata@active.ident,levels=c('Epithelium','EGFP+ Cells',
                                                        'Fibro/Endo Cells',
                                                        'B Cells','gdT Cells','T Cells','NK/NKT Cells','Myeloid'))

colors=get_color(col = rev(rcolors$amwg),n = 8,show = T)

p=DimPlot(adata,# group.by = 'seurat_clusters',
        reduction = 'umap',label = T,repel=T,cols = colors)+theme_map()

ggsave('UMAP_all.pdf',p,width=6,height=5)

p2= DimPlot(adata,
        # label=T,repel=T,
        reduction='umap',split.by = 'batch',cols = colors,ncol=1)&theme_map()

ggsave('UMAP_Split_by_batch.pdf',p2,height=7,width=4)

myumapplot(adata,reduction = 'umap',colors = colors,title = 'CART treated Mouse')

df = table(query$Maintype,query$batch)
df.notumor = df[grep('Tumor',rownames(df),invert = T),]
df.immune = df.notumor[grep('Stromal',rownames(df.notumor),invert = T),]
df=apply(df.immune,2,function(x)x/sum(x))%>%reshape2::melt()
df=apply(df,2,function(x)x/sum(x))%>%reshape2::melt()


# df$Var1 =factor(df$Var1,levels=c(''))

p=ggplot(df,aes(x=reorder(Var1,value),y=value,fill=Var2))+
  geom_bar(stat='identity',position = 'dodge2')+
  scale_fill_manual(values = batchcolor,name='')+
  ggpubr::theme_pubr(legend = c(0.5,0.4))+coord_flip()+ylab('All Cell Ratio')+
  xlab('')+ scale_y_continuous(expand=c(0,0))

p
ggsave('Barplot_ALL_cellratio.pdf',p,width=3.5,height = 3.5,device = cairo_pdf)



df.immune=apply(df.notumor,2,function(x)x/sum(x))%>%reshape2::melt()
p=ggplot(df.immune,aes(x=reorder(Var2,value),y=value,fill=Var1))+
  geom_bar(stat='identity',position = 'fill')+
  scale_fill_manual(values = cols,name='')+
  ggpubr::theme_pubr(legend = 'right',border = T)+
  # ggpubr::theme_pubr(legend = c(0.5,0.4))+
  ylab('Tumor microenvironment Cell Ratio')+
  xlab('')+ scale_y_continuous(expand=c(0,0))

p
ggsave('Barplot_TME_cellratio.pdf',p,width=3.5,height = 3.5,device = cairo_pdf)



#--------subset adata----------
adata.split = SplitObject(adata)

maintype = names(adata.split)
selcells = c(3,4,5,6)

subdata = lapply(selcells,function(i){
  obj =adata.split[[i]]
  obj = ScaleData(obj)%>%FindVariableFeatures(nfeatures=1000)%>%RunHarmony("batch", 
                                                                           plot_convergence = TRUE, 
                                                                           nclust = 50, max_iter = 10, early_stop = F)%>%
  RunUMAP(reduction = 'harmony',dims=1:20)%>%FindNeighbors(reduction='harmony',dims=1:20)%>%FindClusters(cluster.name = paste0('subtype'))
  return(obj)
  
})


names(subdata) = maintype[selcells]

plt=lapply(1:length(subdata), function(i){
  obj = subdata[[i]]
  DimPlot(obj,label = T,repel = T,reduction='umap',pt.size = 0.5,
          cols = get_color(col = rcolors$amwg,n = length(levels(obj@active.ident))))+
    ggtitle(names(subdata)[i])+
    theme(axis.text = element_blank(),axis.ticks = element_blank())
})


p=CombinePlots(plt,ncol=2)


p


ggsave('UMAP_subtype.pdf',p,width=8,height = 6)

# devtools::install_github('RGLab/MAST') install packages MAST  

marker_dic = lapply(subdata,function(x){
  marker = FindAllMarkers(x,only.pos = T,min.pct = 0.3)
  return(marker)
})



fibroendo = subdata$`Fibro/Endo Cells`

DotPlot(fibroendo,features = c('Pecam1','Thy1','Eng','Vwf','Lyve1','Col1a1','Dcn','Cxcl1','Cxcl10'),cluster.idents = T)


fibrocells = names(fibroendo@active.ident)[!fibroendo@active.ident%in%c('2','3','6')]

fibrocells = names(fibroendo@active.ident)[fibroendo@active.ident%in%c('0','5')]

fibroendo =subset(fibroendo,cells=fibrocells)

fibroendo = fibroendo %>%ScaleData()%>%
  FindVariableFeatures(nfeatures = 1000)%>%
  RunPCA()%>%
  RunHarmony(plot_convergence = TRUE, 'batch',
             nclust = 50, max_iter = 10, early_stop = F)%>%
  RunUMAP(reduction='harmony',dims=1:30)%>%FindNeighbors(dims=1:30)

fibroendo = fibroendo%>%FindClusters()

DimPlot(fibroendo,reduction = 'umap',label = T)|DimPlot(fibroendo,reduction = 'umap',group.by = 'batch')

markers.fibro = FindAllMarkers(fibroendo,only.pos = T,test.use = 'MAST',min.pct = 0.3)


DotPlot(fibroendo,features = c('Col1a1','Col1a2','Dcn','Eng','Vwf','Lyve1',
                               grep('Cxcl',rownames(fibroendo),value = T)),
        cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  coord_flip()|DimPlot(fibroendo,label=T,reduction = 'umap')

DimPlot(fibroendo,label=T,reduction = 'umap')|DimPlot(fibroendo,group.by = 'batch',label=T,reduction = 'umap')



myeloid =subdata$Myeloid


DotPlot(myeloid,features =c('Cd14','Adgre1','C1qa','Csf1r',
                            'Cd86','Il1b','Cd3e','Cd3d',
                            'Arg1','Mrc1','Cd163',
                            'Vcan','Vegfa','Spp1',
                            'Csf3r','S100a8','S100a9',
                            'Itgax','H2-Aa','Xcr1','Siglech'),cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))|plt[[3]]


myeloidcells =names(myeloid@active.ident)[!myeloid@active.ident%in%c('2','5','6')]


rename = c('0'='Mo-Mac',
           '1'='M2-like',
           '2'='2',
           '3'='cDC',
           '4'='M2-like',
           '5'='5',
           '6'='6',
           '7'='pDC',
           '8'='Neu')
myeloid = RenameIdents(myeloid,rename)


bc = subdata$`B Cells`

DotPlot(bc,features =c(
  'EGFP','Cd19','Ms4a1','Cd79a','Cd79b','Ighm','Sdc1','Jchain'),cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))|plt[[2]]

bcells = names(bcells@active.ident)

tc = subdata$`T Cells`

surfacemarker=tc@assays$RNA@data[c('Cd4','Cd8a','Cd8b1'),]


surfacemarker = surfacemarker[surfacemarker$Cd4*surfacemarker$Cd8a*surfacemarker$Cd8b1==0,]

tc = subset(tc,cells= rownames(surfacemarker))

DimPlot(tc,reduction = 'umap')

FeaturePlot(tc,features = c('Cd4','Cd8a'),reduction = 'umap')

markers = list(
  'Mono/Gran'=c('Itgam','Itgae','Cd14','Lyz2','Ctss',
                'Adgre1','C1qb','Csf1r','Fcer1g',
                'S100a8','S100a9',
                'Csf3r','Ly6c2',
                'Spp1','Vegfa',
                'Cd163','Mrc1','Cd86'),
  'DC'=c('Cst3','H2-Aa','Itgax','Xcr1','Cd209a','Clec10a','Siglech'),
  'inflo'=c('EGFP','mCherry-CAR'),
  'Stromal'=c('Acta2','Col1a1','Dcn','Eng','Epcam','Krt18','Kdr','Vwf'),
  'B'=c('Cd19','Cd79a','Cd79b','Ms4a1','Sdc1'),
  'T_surface'=c('Cd3d','Cd3e','Cd3g','Cd4','Cd8b1','Trbc1','Trgc2','Trdc'),
  'Stem/Active'=c('Il2ra','Cd69','Mki67','Ccr7','Lef1'),
  'Progenitor'=c('Kit','Slamf15','Cd34','Ly6d','Ly6a','Mpo','Gata1'),
  'Chemokine'=c('Cxcl1','Cxcl10','Cxcr3','Cxcr4','Ccr1','Ccr2'),
  
  'Cytokine'=c('Il1b','Il4','Il6','Il10','Tnf','Ifng','Gzma','Gzmb','Gzmk'),
  'T_function'=c('Cd2','Cd5','Cd7','C1qbp',
                 'Il7r','Cd44','Sell','Tcf7',
                 'Foxp3','Ctla4'),
  'Exhaustion'=c('Pdcd1','Havcr2','Tight','Fas','Fasl','Tox'),
  'NK'=c('Klrb1c','Klrk1','Ncam1')
  
)

DotPlot(query,group.by='predicted.id',
        markers,cluster.idents = T)+
  theme(axis.text.x=element_text(angle=45,hjust = 1,vjust = 1))+
  scale_color_gradient2(low = 'steelblue',mid = 'lightgrey',high = 'red')


query = query %>%ScaleData()%>%RunUMAP(reduction='Xumap_',dims=1:20)

# query@active.ident =factor(query$predicted.id,levels=unique(query$predicted.id))

DimPlot(query)


T_surface =query@assays$RNA@data[c('Cd4','Cd8b1'),]%>%t()%>%as.data.frame()
surfacemarker2=query@assays$RNA@data[c('Cd3e','Cd3g','Cd3d'),]%>%t()%>%as.data.frame()

query = AddMetaData(query,metadata = surfacemarker2)

query = AddMetaData(query,metadata = T_surface)

metadata = query@meta.data

maincelltype = list(
  'B'=c('B','Plasma'),
  'CD4T'=c('CD4_Th','CD4_Tpex','CD4_Tn/cm','Treg'),
                    'CD8T'=c('CD8_Tpex','CD8_Tem','CD8_Trm','CD8_Tn/cm'),
                    'γδT'=c('γδ_Tex','γδ_Teff'),
                    'NKT'=c('NKT'),
                    'Prolif T'=c('Tcycling'),
                    'Myeloid'=c('MDSC','M1','SPP1+ Mø','M2','Neutrophil','Mo-DC','cDC1','cDC2','pDC'),
  'Stromal'=c('Vein EC/Fibro','Fibro','EpC'),
  'Tumor'=c('Tumor'))


maincelltype = reshape2::melt(maincelltype)

colnames(maincelltype)=c('predicted.id','maintype')

metadata$id = rownames(metadata)

metadata = left_join(metadata,maincelltype,"predicted.id")

metadata = data.frame(metadata,row.names = 'id')

query = AddMetaData(query,metadata)

DimPlot(query,group.by = 'maintype',label = T)|DimPlot(query,label = T,repel = T)+NoLegend()

metadata = query@meta.data

nocells = c(rownames(metadata)[metadata$Cd4>0&metadata$maintype=='CD8T'],
            rownames(metadata)[(metadata$Cd4+metadata$Cd8b1)>0&metadata$maintype=='Stromal'],
            # rownames(metadata)[(metadata$Cd3e+metadata$Cd3d+metadata$Cd3g)>0&metadata$maintype=='Myeloid']
            )%>%unique()
length(nocells)

cells= setdiff(rownames(metadata),nocells)

query = subset(query,cells=cells)


rename=levels(query@active.ident)

names(rename)=rename

rename['Tex']='CD4_Tpex'

rename['M0']='Mo-DC'
rename['Neutrophi']='Neutrophil'


query=RenameIdents(query,rename)


df = DotPlot(adata,
             # group.by = 'subtype',
             # group.by='predicted.id',
             markers,cluster.idents = T)+
  theme(axis.text.x=element_text(angle=45,hjust = 1,vjust = 1))+
  scale_color_gradient2(low = 'steelblue',mid = 'lightgrey',high = 'red')

df = df$data

p1=ggplot(df[df$pct.exp>20&df$avg.exp.scaled>0.75,],aes(x=features.plot,y=id,fill=avg.exp.scaled,size=pct.exp))+
  geom_point(shape=21)+facet_grid(facets = ~feature.groups, scales = "free_x", 
                                  space = "free_x", switch = "y")+
  ggpubr::theme_pubr()+scale_size_continuous(range = c(2,6))+
  scale_fill_viridis_c(option = 'A',direction = -1,alpha = 0.8)+
  # scale_fill_gradient2(low = 'steelblue',mid = 'gold',high = 'red')+
  theme(axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5))+xlab('')+ylab('')

ggsave('2_DotPlot_marker.pdf',p1,width=15,height = 8,device = cairo_pdf)

p2=DimPlot(query,reduction='Xumap_',group.by = 'subtype',label = T, repel = T)+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  xlab('UMAP_1')+ylab('UMAP_2')+ggtitle('All Cell Clusters (33242 cells)')
ggsave('1_Cellcluster_UMAP.pdf',p2,width=8,height = 6,device = cairo_pdf)




p3=DimPlot(query,reduction='Xumap_',group.by = 'Maintype',
           cols = get_color(col = rcolors$amwg,n = length(unique(query$Maintype))),
           label = T, repel = T)+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  xlab('UMAP_1')+ylab('UMAP_2')+ggtitle('Main Cell Clusters (33242 cells)')

cols=get_color(col = rcolors$amwg,n = length(unique(query$Maintype)))
names(cols)=c('B','CD4T','CD8T','Myeloid','NKT','Prolif T','Stromal','Tumor','γδT')
p4=DimPlot(query,reduction='Xumap_',group.by = 'Maintype',split.by = 'batch',ncol = 1,
        cols = cols,
        # label = T, repel = T
        )&theme_map()&ggtitle('')

p4
ggsave('3_Cellcluster_UMAP_maintype.pdf',p3,width=7,height = 6,device = cairo_pdf)
ggsave('4_Cellcluster_UMAP_maintype_split.pdf',p4,width=3,height = 6,device = cairo_pdf)

query$Maintype = factor(query$Maintype,levels=c('B','Tumor','Stromal',
                                                'CD8T','CD4T','NKT','γδT','Prolif T',
                                                'Myeloid'))

df=DotPlot(query,group.by = 'Maintype',features=c(
  
  'Cd19','Cd79a','Ms4a1',
  'EGFP','Acta2','Epcam','Col1a1','Dcn',
  'Cd3e','Cd3d','Cd8b1','Cd8a',
  'Cd4','Il7r',
  'Klrb1c','Klrk1','Gzma',
  
  'Trdc','Trgc1',
  'Mki67',
  #'Cd69','Il2ra',
  'Cd14','Lyz2','H2-Aa'
  
   ),cluster.idents = F)+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

df = df$data

p4=ggplot(df[df$pct.exp>27&df$avg.exp.scaled>0.27,],aes(x=features.plot,y=id,fill=avg.exp.scaled,size=pct.exp))+
  geom_point(shape=21)+
  # facet_grid(facets = ~feature.groups, scales = "free_x", 
  #                                 space = "free_x", switch = "y")+
  ggpubr::theme_pubr()+scale_size_continuous(range = c(2,6))+
  scale_fill_viridis_c(option = 'A',direction = -1,alpha = 0.8)+
  # scale_fill_gradient2(low = 'steelblue',mid = 'gold',high = 'red')+
  theme(axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1))+xlab('')+ylab('')+coord_flip()+theme(legend.position = 'right')

ggsave('./2_DotPlot_marker_maintype.pdf',p4,width = 5,height = 5,device=cairo_pdf)

integrated = CreateSeuratObject(counts = query@assays$RNA@counts,
                                meta.data =  query@meta.data)

integrated= integrated %>%NormalizeData()%>%
  ScaleData()%>%FindVariableFeatures(nfeatures = 1000)%>%
  RunPCA()%>%
  RunHarmony(plot_convergence = TRUE, 'batch',
             nclust = 50, max_iter = 10, early_stop = F)%>%
  FindNeighbors(reduction='harmony',dims=1:30)%>%
  RunUMAP(reduction='harmony',dims=1:30)

??write10xCounts
devtools::install_github("MarioniLab/DropletUtils")
library(DropletUtils)
# DimPlot(query,reduction='Xumap_')

celltype =list('T'=c('CD4T','CD8T','NKT','γδT','Prolif T'),
               'B'=c('B'),
               'Myeloid'=c('Myeloid'),
               'Stromal'=c('Stromal'),
               'Tumor'=c('Tumor'))
celltype = reshape2::melt(celltype)
colnames(celltype)=c('Maintype','Split')

metadata = query@meta.data

metadata =left_join(metadata,celltype,"Maintype")

query = AddMetaData(query,metadata = metadata$Split,col.name = 'Split')

query$Split[is.na(query$Split)] ='Tumor' 

querylist = SplitObject(query,split.by = 'Split')

querylist = lapply(querylist,function(obj){
  obj = obj %>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures(nfeatures=1000)%>%
    RunPCA()%>%FindNeighbors(dims=1:30)%>%RunUMAP(dims=1:30)
  return(obj)
})


#----CART cells--------
# query = adata

cart = subset(query,subset=(`mCherry-CAR`>0))

cart = cart %>%ScaleData()%>%FindVariableFeatures(nfeatures = 1500)%>%RunPCA()%>%FindNeighbors(dims=1:30)%>%
  RunUMAP(dims=1:30)
cart =FindClusters(cart)

saveRDS(cart,'mcherry_cart_from_query.rds')
DimPlot(cart,group.by = 'batch',reduction='umap')

pdf('6_UMAP_CART.pdf',width = 10,height =3 )
DimPlot(cart,reduction='umap',cols = get_color(rcolors$Paired,n = 5),
        label = T,repel = T)+ggtitle('CAR-T cell type')|
  DimPlot(cart,group.by = 'batch',reduction='umap',cols = batchcolor[2:3])|
  FeaturePlot(cart,pt.size = 0.1,features = c('mCherry','Cd3e','Cd8b1','Cd4',
                                'Tcf7','Ccr7','Myb','Sell',
                                'Mki67','Nkg7','Gzmb',
                                'Pdcd1','Tox'),
              reduction='umap')&theme_map()&NoLegend()
dev.off()



marker.mast = FindMarkers(cart,group.by = 'batch',
                          ident.1 = 'NR',
                          logfc.threshold=0,test.use= 'MAST')

write.csv(marker.mast,'Marker_MAST_CART_NR_VS_R.csv')

marker.wilcoxlimma = FindMarkers(cart,group.by = 'batch',
                          ident.1 = 'NR',
                          logfc.threshold=0,test.use= 'wilcox_limma')

write.csv(marker.mast,'Marker_wilcoxlimma_CART_NR_VS_R.csv')

# devtools::install_github('yanlinlin82/ggvenn')
library(ggvenn)
library(clusterProfiler)
library(org.Mm.eg.db)
marker.mast$gene = rownames(marker.mast)
gseKEGGmarker.mast = mygseKEGG_Cluster(marker.mast,return.obj=T)

write.csv(gseKEGGmarker.mast@result,'GSEKEGG_Marker_MAST_CART_NR_VS_R.csv')


cart = readRDS('./mcherry_cart_from_query.rds')

DimPlot(cart,group.by = 'batch')

cart.cds=readRDS('~/LJX/sc_in/CART/cart.cds.rds')

library(monocle)
plot_cell_trajectory(cart.cds,color_by = 'ident')
plot_complex_cell_trajectory(cart.cds, color_by = "ident")|plot_complex_cell_trajectory(cart.cds, color_by = "State")
cart.cds=orderCells(cart.cds,root_state = 1)
plot_complex_cell_trajectory(cart.cds,root_states = 3, color_by = "ident",)|plot_complex_cell_trajectory(cart.cds, color_by = "State")

selectmarkers= c('Lef1','Ccr7','Klf7','Rorc','Dntt',
                 'Lck','Cd247','mCherry-CAR','Tcf7','Eif1',
                 'Jund','Ccr9','Pdcd4','Cd27',
                 'Nkg7','Mif','Pdcd1','Mki67',
                 'Tacc3','Nelfe','Tcf19','E2f8','Hmmr')

p=plot_genes_in_pseudotime(cart.cds[intersect(selectmarkers,rownames(cart.cds)),])

df = p$data

unique(df$f_id) %>%length()
library(ggrepel)
library(ggpubr)
p=ggplot(df,aes(x=Pseudotime,y=adjusted_expression,color=group))+geom_smooth()+
  facet_grid(f_id~.)+theme_pubr(border = T)+scale_color_manual(values = batchcolor[2:3])
ggsave('Smoothplot_CART_pseudotime_gene.pdf',p,width = 3.5,height = 15)


# > table(cart$batch)
# 
# NR   R 
# 214 463 
# > table(query)
# Error in unique.default(x, nmax = nmax) : unique()只适用于矢量
# > table(query$batch)
# 
# Control      NR       R 
# 9670   10821   12751 
# > 214/10821
# [1] 0.01977636
# > 467/12751
# [1] 0.03662458

DotPlot(cart,features = c('mCherry-CAR','Cd3d'))



#============= NOT RUN =================

TC = subdata$`T Cells`
TC = FindClusters(TC,resolution=2)


DotPlot(TC,features =c('mCherry-CAR','Cd3e','Cd4','Cd8b1','Pdcd1','Il2ra',
                                      'Cd69','Mki67',
                                      'Ccr7','Lef1','Havcr2',
                                      'Tight','Fas','Fasl','Tox',
                                      'Il7r','Cd44','Sell',
                                      'Foxp3','Ctla4'),
        cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  scale_color_gradient2(low = 'steelblue',mid = 'lightgrey',high = 'red')






pdf('./figures/Dotplot_marker_clusters_all.pdf',width =10,height = 6 )
DotPlot(adata,assay = 'RNA',features = markers,cluster.idents = T)+scale_color_gradient2(low = 'blue',mid = 'lightgrey',high = 'red')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+xlab('')+ylab('')
dev.off()

DotPlot(adata,assay = 'RNA',
        # group.by = 'predictions_scanvi',
        features = c('mCherry-CAR','EGFP', 'Acta2','Ptprc',
                     'Cd19','Cd79a','Cd79b','Ighm', 'Ms4a1',
                     'Hba-a1','Hbb-bt',
                     'Dcn','Col1a1','Vwf','Pecam1','Eng','Epcam','Klf9',
                     'S100a8','S100a9','Csf1r','Csf3r',
                     'Adgre1','Il1b','Mrc1','C1qc',
                     'Itgax','Xcr1','Cd63','Cst3','Irf8','H2-Aa',
                     'Cd3e','Cd8a','Cd8b1','Cd4','Il7r','Mki67','Klrk1','Klrb1c','Klra9'),
        cluster.idents = T)+scale_color_gradient2(low = 'blue',mid = 'lightgrey',high = 'red',midpoint = 0)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+xlab('')+ylab('')


cells = names(adata@active.ident)[!adata@active.ident%in%c('22')]
adata =  subset(adata,cells=cells)

immune = subset(adata,subset=(Ptprc>0))

immunecells = names(immune@active.ident)

cells =  names(adata@active.ident)[!adata@active.ident%in%c('4','5','14','19','21')]

adata= subset(adata, cells=cells)

cluster2cell = c(
  '0'='T Cells',
  '1'='T Cells',
  '2'='T Cells',
  '3'='T Cells',
  '4'='EGFP+ Cells',
  '5'='EGFP+ Cells',
  '6'='T Cells',
  '7'='T Cells',
  '8'='T Cells',
  '9'='T Cells',
  '10'='T Cells',
  '11'='T Cells',
  '12'='gdT Cells',
  '13'='NK/NKT Cells',
  '14'='EGFP+ Cells',
  '15'='Myeloid',
  '16'='Fibro/Endo Cells',
  '17'='B Cells',
  '18'='Myeloid',
  '19'='EGFP+ Cells',
  '20'='Epithelium',
  '21'='T Cells'
)


adata = RenameIdents(adata,cluster2cell)

colannotation = data.frame('celltype'= cluster2cell)

colsplit = cluster2cell

df =table(adata$batch,adata$seurat_clusters)
df = apply(df,2,function(x)x/sum(x))
ComplexHeatmap::pheatmap(df,color = hcl.colors(palette = 'viridis',n=20),
                         # column_split= cluster2cell,
                   cluster_rows = F,treeheight_col = 0,
                   # annotation_col = colannotation,
                   cellheight = 20,cellwidth = 20,border_color = 'black')
pdf('./figures/Featureplot_MainCelltype.pdf',width = 10,height =8 )
FeaturePlot(adata,features = c('EGFP','Cd4','Cd8a','Klrb1c','Trdc','Cd19','Lyz2',
                               'Epcam','Col1a1'),reduction = 'umap')&theme_map()&NoLegend()&scale_color_gradientn(colours = c('lightgrey',hcl.colors(palette = 'viridis',n=4)))

dev.off()

immune = subset(adata,subset=(EGFP==0&Ptprc>0))

immunecells = names(immune@active.ident)[immune@active.ident!='EGFP+ Cell']

immune = subset(immune,cells=immunecells)

DimPlot(immune,label = T,repel = T,reduction='umap')

# immune@active.ident =droplevels(immune@active.ident)
#------EGFP- CD45+ ---------
egfp = adata@assays$RNA@counts['EGFP',]
cd45 = adata@assays$RNA@counts['Ptprc',]
CAR=adata@assays$RNA@counts['mCherry-CAR',]
expmeta = data.frame('EGFP'=egfp,'CD45'=cd45,'mCherry'=CAR)
adata = AddMetaData(adata,expmeta)
immune = rownames(adata@meta.data)[adata@meta.data$EGFP==0&adata$CD45>0]
immune =subset(adata,cells=immune)
table(immune$batch)

immune = immune %>%ScaleData()%>%FindVariableFeatures()%>%
  RunPCA()%>%
  RunHarmony(plot_convergence = TRUE, 'batch',
             nclust = 50, max_iter = 10, early_stop = F)%>%
  RunUMAP(reduction='harmony',dims=1:30)%>%FindNeighbors(dims=1:30)%>%FindClusters()


p1= DimPlot(immune,group.by = 'batch',cols = batchcolor,reduction = 'umap')
p2 = DimPlot(immune,group.by = 'seurat_clusters',reduction='umap',label=T,repel = T)
plot_grid(p1+p2)

DotPlot(immune,features = markers,cluster.idents = T)+
  scale_color_gradient2(low = 'steelblue',mid = 'lightgrey',high = 'red')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

DotPlot(adata,features = markers,cluster.idents = T)+
  scale_color_gradient2(low = 'steelblue',mid = 'lightgrey',high = 'red')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


# immune = subset(adata,subset=(('EGFP'==0)&('Ptprc'>0)))

cells = names(immune@active.ident)[!immune@active.ident%in%c('13','19')]

immune = subset(immune,cells = cells)
immune = immune%>%RunHarmony('batch')%>%
  RunUMAP(reduction='harmony',dims=1:30) %>%
  FindNeighbors(dims=1:30,reduction='harmony')%>%FindClusters()

immune = RunUMAP(immune,reduction = 'harmony',dims=1:30)

function(adata){
  adata= ifnb
  ifnb.list <- SplitObject(ifnb, split.by = "stim")
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  features <- SelectIntegrationFeatures(object.list = ifnb.list,nfeatures = 10000)
  features = c('EGFP','mCherry-CAR',features)%>%unqiue()
  
  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
  
  # 生成整合数据
  immune.combined <- IntegrateData(anchorset = immune.anchors)
  
  DefaultAssay(immune.combined) <- "integrated"
  
  
}


#-----final all  --------------
adata =readRDS('./sc_all_29319cells.rds')
adata = h5ad2seurat('./adata_concat_scvi_annotated_2500.h5ad')
adata_concat = h5ad2seurat('./adata_concat1.h5ad')
rename = c('gdT Cell'='T/NK Cells',
           'CD8T Cell'='T/NK Cells',
           'CD4T Cell'='T/NK Cells',
           'Nature Killer Cell'='T/NK Cells',
           'Double-positive T Cell'='T/NK Cells',
           'EGFP+ Cell'='Tumor Cells',
           'B Cell'='B Cells',
           'Myeloid'='Myeloids',
           'Epithelial Cell'='non-immune Cells',
           'Stromal Cell'='non-immune Cells'
           )

adata@active.ident = factor(adata$leiden,levels = unique(adata$leiden)[order(unique(adata$leiden))])

DotPlot(adata_concat,features=c('Cd3e','Cd4','Cd8b1','Trbc1','Trgc2','Klrb1c','Klrk1','Txk','Mki67',
                                'Cd19','Cd79a','Ms4a1',
                                'Lyz2','Csf1r','C1qa','Itgax','Siglech',
                                'Dcn','Col1a1','Eng','Pecam1','Epcam','Krt8','Krt18','EGFP','Acta2','mCherry-CAR'
                                ),group.by = 'leiden',cluster.idents = T)+scale_color_gradient2(low = 'steelblue',mid = 'lightgrey',high = 'red')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

leiden2name =c(
  '0'='Double Positive T Cells',
  '1'='Double Positive T Cells',
  '2'='CD4+ T Cells',
  '3'='CD4+ T Cells',
  '4'='EGFP+ T Cells',
  '5'='CD4+ T Cells',
  '6'='CD8+ T Cells',
  '7'='NK/NKT Cells',
  '8'='Myeloids',
  '9'='gamma-delta T Cells',
  '10'='Double Positive T Cells',
  '11'='gamma-delta T Cells',
  '12'='B/Plasma Cells',
  '13'='EGFP+ Cells',
  '14'='Fibroblasts',
  '15'='gamma-delta T Cells',
  '16' = 'Epithelium',
  '17'='Unknown',
  '18'='Cycling T Cells'
)
adata_concat@active.ident = factor(adata_concat$leiden,levels=unique(adata_concat$leiden)[order(unique(adata_concat$leiden))])

adata_concat = RenameIdents(adata_concat,leiden2name)


library(rcolors)

cellcolors = get_color(col = rcolors$amwg,n = length(levels(adata_concat@active.ident)))
names(cellcolors)= levels(adata_concat@active.ident)
cellcolors['Unknown']<-'lightgrey'

p1=myumapplot(adata_concat,title = 'Single Cell of Mouse with CAR-T therapies',reduction = 'Xumap_',
              colors =get_color(col = rcolors$amwg,n = length(levels(adata_concat@active.ident))),repel = T) 
ggsave('./figures/UMAP_annotated_leiden_celltype.pdf',p1,width = 8,height = 6)

# adata = RenameIdents(adata,rename)
# adata <- RunHarmony(adata, "batch")

#--------remove EGFP+/mCherry+_CD3+ cells---------

obj = subset(adata_concat,subset=(EGFP==0))

obj = subset(obj,subset=(`mCherry-CAR`==0))

obj = obj%>%ScaleData()%>%FindVariableFeatures()

obj = obj%>%FindNeighbors(redcution= 'XscVI_',dims=1:30)%>%RunUMAP(redcution= 'XscVI_',dims=1:30)

DimPlot(obj,reduction='umap',label = T)



#----refine Cell subset-------------
adata.split =  SplitObject(adata_concat)

DPT =adata.split$`Double Positive T Cells`


DPT = DPT %>%ScaleData()%>%FindVariableFeatures()%>%RunPCA()%>%FindNeighbors(dims=1:30)%>%RunUMAP(dims=1:30)

DPT = DPT%>%FindVariableFeatures(nfeatures = 500) %>%RunPCA(features = c('EGFP','mCherry-CAR','Cd4','Cd8b1',VariableFeatures(DPT)))%>%FindNeighbors(dims=1:30)%>%RunUMAP(dims=1:30)

FeatureScatter(DPT,feature1 = 'Cd4',feature2 = 'Cd8b1')





adata =adata%>%NormalizeData%>%ScaleData()%>%FindVariableFeatures(nfeatures = 4000)%>%RunPCA()

myEmbedding = function(obj,res=0.6){
  obj =  adata
  adata <- adata %>%NormalizeData()%>%
    ScaleData()%>%FindVariableFeatures(nfeatures=2500)%>%RunPCA()%>%
    RunHarmony("batch", 
               plot_convergence = TRUE, 
               nclust = 50, max_iter = 10, early_stop = F)
  
  adata = RunUMAP(adata,reduction = 'harmony',dims = 1:30)
  
  # adata = AddMetaData(adata,metadata = adata@active.ident,col.name = 'major Celltype')
  
  
  adata <- adata %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = res)
  return(adata)
}

cells = names(adata@active.ident)[adata@active.ident!=11]
adata = subset(adata,cells = cells)

adata = myEmbedding(adata,res = 1)


adata_concat = subset(adata_concat,cells = rownames(adata@meta.data))

adata_concat@reductions = adata@reductions

adata_concat = AddMetaData(adata_concat,adata@meta.data)

adata_concat = adata_concat%>%NormalizeData()%>%ScaleData()

adata_concat@active.ident = adata@active.ident

adata_concat@active.ident = factor(adata$leiden,
                                   levels=unique(adata$leiden)[order(unique(adata$leiden))])

DimPlot(adata_concat,reduction = 'umap',label = T,repel = T)|
  DimPlot(adata_concat,reduction = 'umap',group.by = 'batch',label = T,
          repel = T,cols = batchcolor)


DotPlot(adata_concat,features = markers,cluster.idents = T)+
  scale_color_gradient2(low = 'steelblue',mid = 'lightgrey',high = 'red')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

p1=DimPlot(adata,group.by = 'leiden',label = T,repel = T)+theme_map()
p2=DimPlot(adata,group.by = 'batch',cols = batchcolor)+theme_map()
p3=DimPlot(adata,group.by = 'majority_voting',label = T,repel = T)+theme_map()
p4=DimPlot(adata,group.by = 'broad_celltype',label = T,repel = T)+theme_map()
plot_grid(p1,p3,p2,p4,ncol = 2)


df =table(adata_concat$batch,adata_concat$leiden)
df = apply(df,2,function(x)x/sum(x))
ComplexHeatmap::pheatmap(df,color = hcl.colors(palette = 'viridis',n=20),
                         # column_split= cluster2cell,
                         cluster_rows = F,treeheight_col = 0,
                         # annotation_col = colannotation,
                         cellheight = 20,cellwidth = 20,border_color = 'black')
adata_concat@active.ident= adata@active.ident

renamelist = list(
  'gdT Cells'=c(2,6,11,12,15,21,23,26),
  'Double-positive T Cells'=c(1,4,7,10,16,19),
  'Nature Killer Cells'=c(14),
  'CD8T Cells'=c(8,9,13),
  'CD4T Cells'=c(0,3,5,25,28),
  'B Cells'=c(18),
  'Myeloid'=c(17,24),
  'Fibroblasts'=c(20,22),
  'Epithelium'=c(27)
)

adata_concat@active.ident = adata@active.ident

rename= reshape2::melt(renamelist)

rename_ident=rename$L1

names(rename_ident)=rename$value

rename_ident=rename_ident[order(rename_ident)]

adata_concat=RenameIdents(adata_concat,rename_ident)

adata_concat=AddMetaData(adata_concat,adata_concat@active.ident,col.name = 'Cell\\ Type')


DimPlot(adata_concat,group.by = 'Cell\\ Type',label = T,repel = T)|DimPlot(adata_concat,label = T,,group.by = 'majority_voting')+ggtitle('Cell Typist')

marker_all = FindAllMarkers(adata_concat,only.pos = T,min.pct =0.3)

marker_all.split = split(marker_all$gene,marker_all$cluster)

library(org.Mm.eg.db)
library(clusterProfiler)

marker.go = lapply(marker_all.split,function(x){
  go= enrichGO(x,OrgDb = org.Mm.eg.db,keyType = 'SYMBOL')
  # return(go@result)
})

# s.features = read.csv('~/LJX/')
adata_concat = CellCycleScoring(adata_concat,s.features = )

#------immune cells -------------


#-----T/NK EGFP-//----------
tc_cells = names(adata@active.ident)[]


#-----myeloid EGFP-//----------
query = PercentageFeatureSet(query,pattern = 'mt-',assay = 'RNA',col.name = 'percent.mt')
DotPlot(query,features=c('Ccr1','Ccr2', grep('Ccl',rownames(query),value = T)))+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

features = c('Ccr1','Ccr2','Ccl2','Ccl6','Ccl7','Ccl8','Ccl9','Ccl12')

DotPlot(query,features=features)+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

macro = c('M1','M2','SPP1+ Mø')

macrocells = names(query@active.ident)[query@active.ident%in%macro]

obj.macro=subset(query,cells=macrocells)
DimPlot(obj.macro,label = T)

obj.macro = obj.macro%>%NormalizeData()%>%ScaleData()

obj.macro <- obj.macro %>% FindVariableFeatures()%>%RunPCA()%>%
  RunHarmony("batch", 
             plot_convergence = TRUE, 
             nclust = 50, max_iter = 10, early_stop = F)

obj.macro = FindNeighbors(obj.macro,reduction = 'harmony',dims=1:10)%>%RunUMAP(dims=1:10,reduction = 'harmony')

DimPlot(obj.macro,reduction = 'umap',label = T)|DimPlot(obj.macro,group.by = 'batch',reduction = 'umap',label = T,cols = batchcolor)

# obj.macro = FindClusters(obj.macro)

obj.macro = AddMetaData(obj.macro,obj.macro@active.ident,col.name = 'predicted_type')

obj.macro = FindClusters(obj.macro)

obj.macro = FindNeighbors(obj.macro,reduction = 'XscVI_',dims=1:20)%>%RunUMAP(reduction = 'XscVI_',dims=1:20)

DimPlot(obj.macro,reduction = 'umap',label = T)|DimPlot(obj.macro,group.by = 'batch',reduction = 'umap',label = T)

markers.macro= FindAllMarkers(obj.macro,only.pos = T,min.pct = 0.3)

markers = c(
            'Lyve1','Cd163','Cd209f','Il4ra','Mrc1','Saa3','Arg1',
            'Ccr1','Ccr2',
            'Ccl2','Ccl6','Ccl7','Ccl8','Ccl9','Ccl12',
            'Ifitm2','Hp','Il1b','Ms4a4c','Ly6c2',#M1
            'Spp1','Fabp4','Fabp5','Ctsd')

DimPlot(obj.macro,reduction = 'umap',label = T,cols = get_color(rcolors$Paired,n = 6))|DimPlot(obj.macro,group.by = 'batch',reduction = 'umap',label = T,cols = batchcolor)
pdf('figures/Dotplot_Macrophage_selectfeatures.pdf',width = 6,height = 3.5)
DotPlot(obj.macro,features = markers,cluster.idents = T)+scale_color_gradient2(low = 'lightgrey',mid = 'white',high = 'red')+theme(legend.position = 'top',axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()
# obj.macro = 
saveRDS(obj.macro,'scMacrophage.RDS')

library(monocle)
# source('~/function.R')
saveRDS(object = obj.macro,'scMacrophage.RDS')

cds= mymonocle(obj = obj.macro,ident = obj.macro@active.ident,group = 'batch',
               cellid = names(obj.macro@active.ident),save.cds ='macrophage.cds' )

cds
cds = orderCells(cds,reverse = F)
saveRDS(cds,'~/LJX/sc_in/rds/monocle_macrophage_cds.rds')

df = data.frame(Pseudotime=cds$Pseudotime,
                State= cds$State,
                cluster = paste0('M_',cds$ident),
                group=obj.macro$batch)
cds$group= obj.macro$batch
pdf('./figures/Monocle_Macrophagesubset.pdf',width =12,height = 4)
plot_cell_trajectory(cds,color_by = 'Pseudotime')+scale_color_viridis_c()|plot_cell_trajectory(cds,color_by = 'ident')+scale_color_manual(values = get_color(rcolors$Paired,n = 6))|plot_cell_trajectory(cds,color_by = 'State')+scale_color_manual(values = get_color(rcolors$Set1,n = 3))|
  plot_cell_trajectory(cds,color_by = 'group')+scale_color_manual(values = batchcolor)
# ggplot(df,aes(x=Pseudotime,fill=cluster,color=cluster))+
#   geom_density()+
#   scale_color_manual(values = get_color(rcolors$Paired,n = 6))+
#   scale_fill_manual(values = get_color(rcolors$Paired,n = 6))+
#   # facet_grid(.~group)+
#   theme_pubclean()+scale_x_continuous(expand = c(0,0))+
#   scale_y_continuous(expand = c(0,0))+facet_grid(.~State)+scale_y_sqrt()
dev.off()


pdf('./figures/DensityMonocle_Macrophagesubset.pdf',width =7,height = 3)
  ggplot(df,aes(x=Pseudotime,fill=cluster,color=cluster))+
  geom_density()+
  scale_color_manual(values = get_color(rcolors$Paired,n = 6))+
  scale_fill_manual(values = get_color(rcolors$Paired,n = 6))+
  # facet_grid(.~group)+
  theme_pubclean()+scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+facet_grid(.~State)+scale_y_sqrt()
dev.off()
#-----Stromal EGFP-//-------
stromal = rownames(adata@meta.data)[adata$Maintype=='Stromal']
stromal = subset(adata,cells=stromal)

stromal = stromal%>%ScaleData()%>%FindVariableFeatures()%>%FindNeighbors(dims=1:20,reduction = 'XscVI_')%>%RunUMAP(dims=1:20,reduction = 'XscVI_')

DimPlot(stromal,reduction = 'umap')

stromal = FindClusters(stromal,resolution = 0.5)

stromalcolor =get_color(col = 'Spectral',n = 5)

p=DimPlot(stromal,label = T,reduction = 'umap',cols = stromalcolor)+ggtitle('Stromal Cell(413)')

df=p$data

p0 = ggplot(df,aes(x=umap_1,y=umap_2,fill=ident))+geom_point(shape=21)+scale_fill_manual(values = stromalcolor)


p0 = LabelClusters(p0,id = 'ident',repel = T,box = T)+theme_cowplot()+ggtitle('Stromal Cell(413)')

pdf('5_UMAP_Stromal.pdf',width = 12,height = 4)
p0|DimPlot(stromal,group.by = 'batch',reduction = 'umap',cols = batchcolor)|
  
  FeaturePlot(stromal,reduction = 'umap',features = c('Pecam1','Vwf','Eng','Epcam',
                                                      'Col1a1','Dcn','Fap','Cxcl10',
                                                      'Cxcl1'),pt.size = 0.1,order = T)&theme_map()&NoLegend()
dev.off()

features = c('Epcam','Krt18','Pecam1','Cxcl2',
             'Vwf','Eng','Fap','Ackr1',
             'Col1a1','Dcn',
             'Cxcl10','Cxcl1')

mat.stromal = AverageExpression(object = stromal,features = features,assays = 'RNA')$RNA

pheatmap::pheatmap(as.matrix(mat.stromal),treeheight_col = 0,
         cellheight = 10,cellwidth = 20,border_color = 'black',
         cluster_rows = F,main = 'Stromal Cell Marker',
         filename = '5_Heatmap_stromal_marker.pdf',
          
         # cluster_cols = F,
         color = colorRampPalette(rev(rcolors::rcolors$RdBu))(50),
         scale = 'row')

df.stromal= table(stromal$batch,stromal@active.ident)%>%as.data.frame()
head(df.stromal)
# install.packages("gg.gap")
library(gg.gap)
library(ggplot2)
library(patchwork)

p1=ggplot(df.stromal,aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat = 'identity', color='black',
           position = position_dodge(),
           show.legend = FALSE) +
  theme_bw() +scale_y_continuous(expand = c(0,0))+theme_pubr(border = T)+
  labs(x = NULL, y = NULL)+scale_fill_brewer(palette = 'Spectral')+ylab('Stromal Cell Number')

p2 =gg.gap(plot = p1,
           segments = c(25, 50),
           tick_width = 40,
           rel_heights = c(0.25, 0, 0.2),
           ylim = c(0, 200)
)

pdf('5_Barplot_stromal_cellnumber.pdf',height = 3,width = 4)
p1
p2
dev.off()
#Example
# #数据
# data <-
#   data.frame(x = c("Alpha", "Bravo", "Charlie", "Delta"),
#              y = c(200, 20, 10, 15))
# #画图
# p1 = ggplot(data, aes(x = x, y = y, fill = x)) +
#   geom_bar(stat = 'identity', position = position_dodge(),show.legend = FALSE) +
#   theme_bw() +
#   labs(x = NULL, y = NULL)
# 
# p1
# p2 =gg.gap(plot = p1,
#            segments = c(25, 190),
#            tick_width = 10,
#            rel_heights = c(0.25, 0, 0.1),
#            ylim = c(0, 200)
# )
# p1+p2


#EC fib
ec = names(adata@active.ident)[adata@active.ident=="Vein EC/Fibro"]
ec=subset(adata,cells=ec)
ec = ec%>%ScaleData()%>%FindVariableFeatures()%>%FindNeighbors(dims=1:10,reduction = 'XscVI_')%>%RunUMAP(dims=1:10,reduction = 'XscVI_')


ec = FindClusters(ec)
DimPlot(ec,label = T,reduction = 'umap')|DimPlot(ec,label = T,group.by = 'batch',reduction = 'umap')|
  
  FeaturePlot(ec,reduction = 'umap',features = c('Pecam1','Vwf','Eng','Epcam',
                                                      'Col1a1','Dcn','Fap','Cxcl10',
                                                      'Cxcl1'),pt.size = 0.1,order = T)&theme_map()&NoLegend()




#-----mapping control to annoated----------
control = h5ad2seurat('adata_concat1.h5ad')
control = SplitObject(control,split.by = 'batch')

control = control$Control

VlnPlot(control,features= c('n_genes_by_counts','pct_counts_mt','doublet_score'),pt.size = 0,group.by = 'batch')

control = subset(control,subset =(pct_counts_mt<10&doublet_score<0.2))

control = control %>%NormalizeData()%>%FindVariableFeatures()

reference = readRDS('~/LJX/sc_in/rds/SC_st_ref_maintype.RDS')

mytransferlabel =function(reference, query){
  
  anchors <- FindTransferAnchors(reference = reference, query = query)
  
  # transfer labels
  predictions <- TransferData(anchorset = anchors, refdata = reference$seurat_annotation)
  query <- AddMetaData(object = query, metadata = predictions)
  return(query)
}


reference$seurat_annotation = reference$celltype
query = mytransferlabel(reference = reference,query = adata)

query = mytransferlabel(reference = reference,query = control)
query$seurat_annotation = query$predicted.id


adata = readRDS('query_33242.RDS')


reference.list = list('CART'=reference,'Control'=query)

pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30) 

integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30) 


RunIntegrated=function(obj){
  pancreas.integrated=obj
  DefaultAssay(pancreas.integrated) <- "integrated" 
  
  # Run the standard workflow for visualization and clustering 
  # 数据标准化 
  pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE) 
  # PCA降维 
  pancreas.integrated <- FindVariableFeatures(pancreas.integrated, npcs = 30, verbose = FALSE)
  pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE) 
  pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
  return(pancreas.integrated)
  
}


integrated = RunIntegrated(integrated)

DimPlot(integrated,group.by='seurat_annotation',label=T,repel=T)

#-------fine immune cell markers-----------
adata = readRDS('query_33242.RDS')
markers =c('EGFP','mCherry-CAR',
           'Cd3d','Cd3e','Cd4',
           'Cd8a','Cd8b1','Il7r','Cd40lg',
           'Ccr7','Sell','Tcf7','Rorc','Ccr9',
           'Foxp3','Il2ra','Clta4',
           'Lef1','Tox','Tox2','Pdcd1','Lag3',
           'Gzma','Gzmb','Klrk1',
           'Mki67','Cd69','Itgae',
           'Trdc','Trgc1','Trgc2',
           'Cd74',
           'Lyz2','Ly6c1','Ly6g','Ly6c2',
           'Spp1','Arg1',
           'S100a8','S100a9','Il1b','Csf3r','Csf1r',
           'Cd14','Adgre1','Mrc1',
           'Itgax','Xcr1','Cst3','Siglech',
           'Cd79a','Cd19','Ms4a1','Ighm','Igha',
           'Dcn','Col1a1',
           'Pecam1','Vegfa','Cd34','Vwf',
           # 'Hba-a2','Hba-a1','Trfc',
           'Epcam','Fxyd3')

p=DotPlot(adata,features = markers,cluster.idents = T)+
  scale_color_gradient2(low = 'steelblue',mid = 'lightgrey',high = 'red')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+ylab('')+xlab('')
ggsave('Dotplot_markers_FineClusters.pdf',p,width = 12,height = 6,device = cairo_pdf)

celltype = levels(adata@active.ident)
celltypelist = list('Myeloid'=c("Neutrophil","M1","M2","SPP1+ Mø","MDSC","Mo-DC","cDC1","cDC2",'pDC'),
                    'B'=c('Plasma','B'),
                    'CD8T'=c("CD8_Tn/cm" ,"CD8_Trm","CD8_Tem","CD8_Tpex"),
                    'CD4T'=c("CD4_Tn/cm","CD4_Tpex","CD4_Th",'Treg'),
                    'γδT'=c("γδ_Teff","γδ_Tex"),
                    'NKT'=c('NKT'),
                    'Prolif T'=c('Tcycling'),
                    'Stromal'=c("Fibro","Vein EC/Fibro","EpC"),
                    'Tumor'='Tumor')

N=intersect(unlist(celltypelist),levels(adata@active.ident)) %>%length()

celltypecolor=get_color(rcolors$amwg,n=27)


names(celltypecolor)=intersect(unlist(celltypelist),levels(adata@active.ident))

adata@active.ident=factor(adata@active.ident,
                          levels=intersect(unlist(celltypelist),levels(adata@active.ident)))



p=DimPlot(adata,cols=celltypecolor,
        # raster = T,
        # label = T,
        # repel = T,
        # raster.dpi = c(256,256)
        )+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  xlab('UMAP_1')+ylab('UMAP_2')+ggtitle('All Cell Clusters (33242 cells)')
p

ggsave('1_Cellcluster_UMAP.pdf',p,device = cairo_pdf,width = 9,height = 6)

p=DimPlot(adata,cols=celltypecolor,
          # raster = T,
          label = T,
          repel = T,
          # raster.dpi = c(256,256)
)+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  xlab('UMAP_1')+ylab('UMAP_2')+ggtitle('All Cell Clusters (33242 cells)')
p

p=DimPlot(adata,cols=celltypecolor,
          # raster = T,
          # label = T,
          # repel = T,
          split.by = 'batch',ncol = 1
          # raster.dpi = c(256,256)
)

p=p&theme_map()+NoLegend()
ggsave('1_Cellcluster_UMAP_split.pdf',p,device = cairo_pdf,width =3 ,height = 10)

p=myumapplot(obj = adata,title = 'CAR-T associated TME',
           reduction = 'Xumap_',colors = celltypecolor)
# dev.off()
ggsave('1_Cellcluster_UMAP_label_cellratio.pdf',p,device=cairo_pdf,width = 10,height = 6)

adata= AddMetaData(adata,metadata = adata@active.ident,col.name = 'subtype')

adata@active.ident = adata$Maintype

p=myumapplot(obj = adata,title = 'CAR-T associated TME',repel = F,
           reduction = 'Xumap_',colors = get_color(rcolors$amwg,n = 9))

ggsave('3_Cellcluster_UMAP_maintype.pdf',p,width = 8,height = 6,device = cairo_pdf)

#immune cellratio for subtype
df =table(adata$subtype,adata$batch)
df.imm = df[grep('Tumor',rownames(df),invert=T),]
df.imm = apply(df.imm,2,function(x)x/sum(x))
df.imm = reshape2::melt(df.imm)
main2sub =reshape2::melt(celltypelist)

colnames(df.imm)=c('subtype','Batch','cellratio')
colnames(main2sub)=c('subtype','main')

df.imm = left_join(df.imm,main2sub,"subtype")

df.imm$main=gsub('CD4','',df.imm$main)
df.imm$main=gsub('CD8','',df.imm$main)
df.imm$main=gsub('NK','',df.imm$main)
df.imm$main=gsub('γδ','',df.imm$main)
df.imm$main=gsub('Prolif ','',df.imm$main)

df.imm$Batch = factor(df.imm$Batch,levels=c('Control','NR','R'),labels = c('Ctr','NR','R'))
df.imm.split = split(df.imm,df.imm$main)

names(celltypecolor)=intersect(unlist(celltypelist),levels(adata@active.ident))
library(prismatic)


plt =lapply(1:length(df.imm.split),function(i){
  df0=df.imm.split[[i]]
  p=ggplot(df0,aes(x=Batch,y=100*cellratio,fill=subtype))+
    geom_bar(stat = 'identity')+ggtitle(names(df.imm.split)[i])+scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values = celltypecolor)+ggpubr::theme_pubr(border = T)+NoLegend()
  if(i==1){
    p=p+xlab('')+ylab('%Immune cells')
    
  }else{
    p=p+xlab('')+ylab('')
      # theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
  }
  return(p)
   
})

p=CombinePlots(plt,ncol = 4)
p
ggsave('Barplot_IMMUNE_Subtype_cellratio.pdf',p,width = 8,height = 4)

plt =lapply(1:length(df.imm.split),function(i){
  df0=df.imm.split[[i]]
  p=ggplot(df0,aes(x=Batch,y=100*cellratio,fill=subtype))+
    geom_bar(stat = 'identity',position = 'fill')+ggtitle(names(df.imm.split)[i])+scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values = celltypecolor)+ggpubr::theme_pubr(border = T)+NoLegend()+ylab(paste0('%Percentage of ',names(df.imm.split)[i],' Cells'))

  return(p)
  
})

p=CombinePlots(plt,ncol = 4)
p
ggsave('Barplot_IMMUNE_Subtype_cellratio2.pdf',p,width = 8,height = 4)


plt =lapply(1:length(df.imm.split),function(i){
  df0=df.imm.split[[i]]
  p=ggplot(df0,aes(x=Batch,y=100*cellratio,fill=subtype))+
    geom_bar(stat = 'identity',position = 'fill')+ggtitle(names(df.imm.split)[i])+scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values = celltypecolor)+ggpubr::theme_pubr(border = T,legend = 'right')+ylab(paste0('%Percentage of ',names(df.imm.split)[i],' Cells'))
  
  return(p)
  
})

p=CombinePlots(plt,ncol = 4)
p
ggsave('Barplot_IMMUNE_Subtype_cellratio3.pdf',p,width = 8,height = 6)






#---2LJX----
mangene= c('Mgat1',
           'Man1a',
           'Man1c1',
           'Srd5a3',
           'B4galt1',
           'Ganab',
           'Fut8',
           'Rpn2'
)


mangene1= c(
  # 'Mgat1',
           'Man1a',
           'Man1c1',
           'Srd5a3',
           'B4galt1',
           # 'Ganab',
           'Fut8',
           'Rpn2'
)


mat = AverageExpression(adata,add.ident='batch',features = mangene)$RNA

group = data.frame(name=colnames(mat),sep=colnames(mat))
group = separate(data = group,col = 'sep',sep = '_',into = c('celltype','batch'))

ComplexHeatmap::pheatmap(t(as.matrix(mat)),row_split= group$celltype,
                         treeheight_row = 0,treeheight_col = 0)

mat0 = mat[,grep('M',colnames(mat))]

mat0 = as.matrix(mat0)
group0 = data.frame(name=colnames(mat0),sep=colnames(mat0))
group0 = separate(data = group0,col = 'sep',sep = '_',into = c('celltype','batch'))

pdf('Pheatmap_N_glycan_gene_in_Myeloid.pdf')

ComplexHeatmap::pheatmap(mat0,scale = 'row',
                         border_color = NA,
                         border=T,
                         cluster_cols = F,
                         cellheight = 10,cellwidth = 20,
                         color = rev(rcolors::rcolors$RdBu),
                         column_split= group0$celltype,
                         treeheight_row = 0,treeheight_col = 0)

dev.off()

mat1 =mat0[,c(1:3)]+cbind(mat0[,c(4:5)],matrix(data=0,ncol = 1,nrow = nrow(mat0)))+mat0[,c(11:13)]

colnames(mat1)=c('Control','NR','RE')


reordergene = c('Man1a','Fut8','Rpn2',
                'Ganab',
                'Mgat1','Srd5a3','Man1c1','B4galt1')

rowsplit = rep('N-Glycan biosynethsis',length(reordergene))


pdf('Pheatmap_N_glycan_gene_in_Myeloid_TAM.pdf')
ComplexHeatmap::pheatmap(mat1[reordergene,],
                         scale='row',
                         border_color = 'black',
                         border=T,main = 'Tumor associated macrophage',
                          # = 'Scaled Avg Exp',
                         cluster_rows = F,
                         cluster_cols = F,
                         row_split=rowsplit,
                         # labels_row ='N-Glycan biosynethsis-related genes', 
                         cellheight = 20,cellwidth = 40,
                         color = rev(rcolors::rcolors$RdBu),
                          
                         # column_split= group0$celltype,
                         treeheight_row = 0,treeheight_col = 0)

dev.off()

# #--------NC_Zhao,2024-----------
# setwd("~/LJX/ValidationInHM/NC_Zhao_2024/")
# # matrix =readr::read_tsv(gzfile('matrix.zip'))
# sample = dir('./matrix/')
#=========CellChat with CAR-T======================
setwd("~/LJX/sc_in/rds/")
cc.nr = readRDS('CellChat_NR.rds')
cc.r = readRDS('CellChat_R.rds')
pathways.show = c('ANNEXIN','TGFb')
pathways.show = c('TGFb')
netVisual_circle(cc.nr)
netVisual_aggregate(cc.nr, signaling = pathways.show,layout = "chord")
netVisual_aggregate(cc.r, signaling = pathways.show, layout = "chord")


pdf('../R_figs/Chord_Cellchat_NR_TGFb.pdf',width=10,height = 8)
netVisual_aggregate(cc.nr, signaling = pathways.show,layout = "chord")
CellChat::plotGeneExpression(cc.nr,signaling='TGFb')
dev.off()

pdf('../R_figs/Chord_Cellchat_R_TGFb.pdf',width=10,height = 8)
CellChat::plotGeneExpression(cc.r,signaling='TGFb')
netVisual_aggregate(cc.r, signaling = pathways.show, layout = "chord")
dev.off()

cc.merge = mergeCellChat(list(cc.nr,cc.r),add.names = c('NR','R'))

rankNet(cc.merge,color.use = batchcolor)

netVisual_bubble(cc.nr,targets.use = 1:4)

netVisual_bubble(cc.r,targets.use = 1:4)
#========myeloid Figure 5==========
setwd('~/LJX/sc_in/rds/')
mye = readRDS('./MYE.RDS')
df= table(mye@active.ident,mye$sample)
df= apply(df,1,function(x)x/sum(x)) %>%reshape2::melt()
p1=ggplot(df,aes(x=Var1,y=value,fill=Var2))+
  geom_bar(stat = 'identity',position = position_dodge(),alpha=0.5)+
  scale_fill_manual(values = c('lightgrey','red'))+theme_pubr()+
  ylab('Relative cell ratio of myeloid cells%')+xlab('')
df2= table(mye@active.ident,mye$sample)%>%as.data.frame()
p2=ggplot(df2,aes(x=Var2,y=Freq,fill=Var1))+geom_bar(stat = 'identity',position ='fill')+
  scale_fill_manual(values = rcolors::get_color(rcolors::rcolors$Set1,n=9))+theme_pubr()+
  ylab('Relative cell ratio of myeloid cells%')+xlab('')+theme(legend.position = 'right')
ggsave('../Barplot_myecelltype_grouped.pdf',p1,width = 8,height = 4)
ggsave('../Barplot_myecelltype_stacked.pdf',p2,width = 4,height = 4)
