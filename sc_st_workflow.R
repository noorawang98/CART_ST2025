devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
n9_path = "~/LJX/st/n9/05.AllheStat/heAuto_level_matrix/subdata"
r5_path ="~/LJX/st/r5/05.AllheStat/heAuto_level_matrix/subdata"

setwd("~/LJX/st/n9/n9")
featurepos= readr::read_tsv("./barcode_pos.tsv",col_names = F)
featurepos[1:3,1:3]
library(ggplot2)
ggplot(featurepos,aes(x=featurepos))

library(dplyr)
library(Seurat)
library(patchwork)
#------------------- my local functions ----------------
myreadmat = function(path){
  setwd(path)
  mat = readMM("matrix.tsv")
  feature = readr::read_delim('features.tsv',col_names = F)
  barcode =readr::read_tsv("barcode.tsv",col_names = F)
  mat@Dimnames = list(feature[[2]],barcode[[1]])
  return(mat)
}

mypreprocess = function(count,meta.data=NULL,min.cells=3,min.features=200){
  #-----------S1-data preprocess -------------------------------------------
  if(!is.null(meta.data=meta.data)){
    pbmc <- CreateSeuratObject(count,meta.data =meta.data ,min.cells = min.cells, min.features = min.features)
  }else{
    pbmc <- CreateSeuratObject(count,min.cells = min.cells, min.features = min.features)
    # pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
    # pbmc <- NormalizeData(pbmc)
  }
   VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  pbmc <-SCTransform(pbmc)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  # FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  pbmc <- FindVariableFeatures(pbmc, x.low.cutoff = 0.0125, y.cutoff = 0.25, do.plot=FALSE)
  pbmc <- RunPCA(pbmc,features = VariableFeatures(pbmc))
  pbmc <- FindNeighbors(pbmc,dim=1:20)
  # pbmc <- FindClusters(pbmc,dim=1:20,resolution = 0.5)
  # pbmc <- RunUMAP(pbmc,dim=1:20)
  
}

bestnpc <- function(data,rep=100){
  data <- JackStraw(data, num.replicate = rep)
  data <- ScoreJackStraw(data, dims = 1:20)
  plot1<-JackStrawPlot(data, dims = 1:15)
  npc = strsplit(gsub(":","",levels(plot1$data$PC.Score)),split = " ") %>%as.data.frame() %>% t()
  npc = apply(npc,2,as.numeric)%>%as.data.frame()
  npc = max(npc$V2[npc$V3<0.05])
  print(paste0("Best N of PC components:",npc))
  return(npc)
} #no sct-transcform


basicFindCluster <- function(data,npc,res=1){
  data = RunPCA(data, verbose = FALSE) 
  data = FindNeighbors(data,dims = 1:npc)
  data =FindClusters(data,resolution = res)
  data =RunTSNE(data,dims = 1:npc)
  data =RunUMAP(data,dims = 1:npc)
  return(data)
}

cal_zoom_rate = function(width, height){
  std_width = 1000
  std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
  if(std_width / std_height > width / height){
    scale = width / std_width
  }
  else{
    scale = height / std_height
  }
  return(scale)
}

mySeurat = function(path,levelpath,hefigpath,
                    projname){
  setwd(path)
  cd = Read10X(levelpath)
  setwd(levelpath)
  pos = readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
  pos= data.frame(pos,row.names = 1)
  colnames(pos)=c("x","y")
  
  hefig =png::readPNG(hefigpath)

  zoom_scale = cal_zoom_rate(ncol(he_fig), nrow(he_fig))
  pos$x=pos$x * zoom_scale
  pos$y=pos$y * zoom_scale
  pos$y=dim(he_fig)[1]-pos$y
  
  setwd("..")
  data = CreateSeuratObject(cd,project = projname,min.cells = 3)
  data= PercentageFeatureSet(data,pattern = "^mt-",col.name = "percent.mt")
  data = SCTransform(data,return.only.var.genes = F,method="glmGamPoi",
                     vars.to.regress = "percent.mt",min_cells=3)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  # data = subset(data, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < 5)
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  # data = RunPCA(data)
  # npc = 20
  # data = basicFindCluster(data,npc=npc)
  data = AddMetaData(data,pos)
  
  return(data)
}


myclstplot = function(he_file,FilePath,
                      Cluster,hetype="png",
                      celltypecol = NULL,pt.size=3.5){
  
  if(hetype=="png"){
    he_fig=png::readPNG(he_file)
  }if(hetype=='tiff'){
    he_fig=tiff::readTIFF(he_file)
  }
  
  w=ncol(he_fig)
  h=nrow(he_fig)
  
  cal_zoom_rate = function(width, height){
    std_width = 1000
    std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
    if(std_width / std_height > width / height){
      scale = width / std_width
    }
    else{
      scale = height / std_height
    }
    return(scale)
  }
  
  zoom_scale = cal_zoom_rate(ncol(he_fig), nrow(he_fig))
  bc_pos_file = gzfile(paste(FilePath,"barcodes_pos.tsv.gz", sep  = "/"),'rt')
  bc_pos = read.table(bc_pos_file, header = FALSE, sep = '\t', quote = '')
  names(bc_pos)=c("Barcode","x","y")
  Cluster = data.frame(cluster=Cluster)
  colnames(Cluster)=c("cluster")
  Cluster$Barcode = rownames(Cluster)
  bc_pos=merge(bc_pos,Cluster,by="Barcode",all=FALSE)
  
  bc_pos$x=bc_pos$x * zoom_scale
  bc_pos$y=bc_pos$y * zoom_scale
  bc_pos$y=dim(he_fig)[1]-bc_pos$y
  
  clstr_plot<-ggplot(data=bc_pos,aes(x ,y)) + 
    ggpubr::background_image(he_fig)+
    # geom_point(aes(colour = cluster),size = opt$point_size, shape=16) + 
    geom_point(aes(fill = cluster,color=cluster),size = pt.size, shape=21) + 
    # scale_fill_brewer(palette = "Set1")+
    # scale_color_manual(values = col)+
    theme_bw() +
    theme(plot.title = element_text(face = 2,size = 50,hjust = 0.5)) +
    theme(axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          panel.border = element_blank(), 
          axis.text.x = element_blank()) +
    xlab('')+ylab('')+
    coord_cartesian(xlim = c(0, dim(he_fig)[2]), ylim = c(0, dim(he_fig)[1]), 
                    expand = FALSE)+
    guides(colour = guide_legend(override.aes = list(size=3.5)))+
    theme(legend.position = "none")
  
  # if(is.null(celltypecol)){
  #   clstr_plot = clstr_plot
  # }else{
  #   clstr_plot = clstr_plot+scale_fill_manual(values= celltypecol)+
  #     scale_color_manual(values= celltypecol)
  # }
  # 
  return(clstr_plot)
} 

st_distance = function(x,y){
  central_x = (max(x)-min(x))/2
  central_y =(max(y)-min(y))/2
  distance = sqrt((x-central_x)^2 + (y-central_y)^2)
}
#-----------end functions ---------------------------

#----1.preprocess single cell data with scanpy---------------
# if 
# For the data preprocessed in Scanpy
setwd("~/LJX/sc_in/")
n9 = Read10X("./N9/")
r5 = Read10X("./R5/")

obj.list = list('Non-response'=n9,'Response'=r5)

obj_pre.list = lapply(1:length(obj.list),function(x){
  mypreprocess(count=obj.list[[x]])
})

names(obj_pre.list)=names(obj_pre.list)
obj.anchors <- FindIntegrationAnchors(object.list = obj_pre.list, dims = 1:30) 
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30) 

obs = readr::read_csv("obs/obs_all_20230815.csv")

n9.cd = obj.list$`Non-response`[,obs$barcodeID[obs$sample=="N9"]]
r5.cd = obj.list
#------read spatial transcriptome data -------------
levelpath = dir(path = n9_path)
levelpath = grep("pdf",levelpath,invert = T,value = T)

#
n9_r2 = mySeurat(path = n9_path,levelpath = levelpath[3],projname = "Non-response")
gc()
n9_r2 = basicFindCluster(data = n9_r2,npc = 20,res = 1)
gc()
r5_r2 = mySeurat(path = r5_path,levelpath = levelpath[3],projname = "Response")
r5_r2 = basicFindCluster(data = r5_r2,npc = 20,res = 1)


n9_r13 = mySeurat(path = n9_path,levelpath = levelpath[2],projname = "Non-response")
gc()
n9_r13 = basicFindCluster(data = n9_r13,npc = 20,res = 1)
gc()
r5_r13 = mySeurat(path = r5_path,levelpath = levelpath[2],projname = "Response")
r5_r13= basicFindCluster(data = r5_r13,npc = 20,res = 1)

objlist = list("N9"=n9_r13,"R5"=r5_r13)


n9_r7 = mySeurat(path = n9_path,levelpath = levelpath[3],projname = "Non-response")
gc()
n9_r7 = basicFindCluster(data = n9_r2,npc = 20,res = 1)
gc()
r5_r7 = mySeurat(path = r5_path,levelpath = levelpath[3],projname = "Response")
r5_r2 = basicFindCluster(data = r5_r2,npc = 20,res = 1)

objlist = list("N9"=n9_r2,"R5"=r5_r2)
rm(n9_r2)
rm(r5_r2)

mm.markers = readxl::read_xlsx("~/LJX/mm/Cell_marker_Mouse_blood_vessels.xlsx",2)
mm.cellmarkers = readxl::read_xlsx("~/LJX/mm/Cell_marker_Mouse.xlsx",1)
features = SelectIntegrationFeatures(object.list = objlist,nfeatures = 5000)
features = c("EGFP","Cd3e","Cd8a","Cd4","Foxp3","Tbx21","Gata1","Gata2","Gata3","Gata4",
             "mCherry-CAR","Il2ra","Gzmb","Gzma","Lag3","Pdcd1","Sell",
             "Siglech","Adgre1","Itgam","Itgax","Cd80","Cd86",'Cd163','Arg1','Mrc1',
             'Cd19','Cd79a','Cd79b','Vpreb3',
             features)
features = unique(features)
features = intersect(features,rownames(objlist$N9@assays$RNA@counts))


n9_r13 = mySeurat(path = n9_path,levelpath = levelpath[2],min.features = 100,max.features = 60000,projname = "Non-response")
gc()
# n9_r13 = basicFindCluster(data = n9_r13,npc = 20,res = 1)
gc()
r5_r13 = mySeurat(path = r5_path,levelpath = levelpath[2],min.features = 100,max.features = 6000,projname = "Response")
# r5_r13= basicFindCluster(data = r5_r13,npc = 20,res = 1)


setwd("L13_cluster")
objlist = list("N9"=n9_r13,"R5"=r5_r13)
immune.markers = c("EGFP","'Col1a1","Col3a1","Cd3d","Cd8b1","Cd4","Foxp3","Tbx21","Gata1","Gata2","Gata3","Gata4",
                   "mCherry-CAR","Il2ra","Gzmb","Gzma","Lag3","Pdcd1","Sell",
                   'Nkg7','Klrb1c',
                   "Ly6g","Fcgr3","Itgam","","Cd14","Siglech","Adgre1","Itgax","Cd80","Cd86",'Cd163','Arg1','Mrc1',
                   'Cd19','Cd79a','Cd79b','Vpreb3',
                   "Kit",'Ly6a',"Cd34","Cd38")
features = SelectIntegrationFeatures(object.list = objlist,nfeatures = 5000)
features = unique(c(features,immune.markers))
features = intersect(features,rownames(objlist$N9@assays$RNA@counts))
immune.anchors <- FindIntegrationAnchors(object.list = objlist, anchor.features = features)
gc()
immune.combined <- IntegrateData(anchorset = immune.anchors) #require min.cell >=3 
# immune.l2=immune.combined
rm(immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE) 
# immune.combined <- RunPCA(immune.combined)
immune.combined = basicFindCluster(immune.combined,npc = 30,res = 1)

group = rep("",nrow(immune.combined@meta.data))

group[grep("1$", rownames(immune.combined@meta.data))] = rep("Non-response",grep("1$", rownames(immune.combined@meta.data),value = T) %>%length())
group[grep("2$", rownames(immune.combined@meta.data))] = rep("Response",grep("2$", rownames(immune.combined@meta.data),value = T) %>%length())
subtypecolor = hcl.colors(levels(immune.combined@active.ident)%>%length(),palette = "Dark 2")
subtypecolor = ArchRPalettes$ironMan
names(subtypecolor)=levels(immune.combined@active.ident) 
immune.combined = AddMetaData(immune.combined,metadata = group,col.name = "Group")
p1 =DimPlot(immune.combined,label = T,cols = subtypecolor)
p2 = DimPlot(immune.combined,group.by = "Group",label = T)+NoLegend()
pdf("./umap_sct_npc30_res1_4556cells.pdf",width=8,height = 4)
p1+p2
dev.off()
pdf("./umap_sct_npc30_res1_4556cells.pdf.pdf",width=8,height = 4)
DimPlot(immune.combined,split.by = "Group",label = T,cols = subtypecolor)
dev.off()

pdf("./featureplot_immunemarkers.pdf",width = 5,height = 4)
FeaturePlot(immune.combined,features = immune.markers,cols = rcolors::rcolors$NCV_jet,
            combine = F,
            label = T,order=T,
            # max.cutoff = "q10"
            )
dev.off()

distance.list = split(immune.combined@meta.data,immune.combined$Group)
dis = lapply(distance.list,function(x){
  df0=x
  dis = st_distance(x = df0$spatial_x,y=df0$spatial_y)
  return(dis)
}) %>%unlist()

immune.combined = AddMetaData(immune.combined,dis,col.name = "distance")

pdf("umap_distance_sct_npc30_res1_4556cells.pdf.pdf",width = 8,height = 7)
FeaturePlot(immune.combined,features = c("distance","mCherry-CAR"),
            cols = rev(viridis_pal()(50)),label = T,split.by = "Group",order = T)
dev.off()

df0 = immune.combined@meta.data

df0$EGFP = as.numeric(immune.combined@assays$RNA@counts["EGFP",])
df0$CAR = as.numeric(immune.combined@assays$RNA@counts["mCherry-CAR",])

# mylineplot =function(df,y,color){ggplot(df,aes(x=distance,y=get(y),color = get(color)))+
#   geom_smooth()+
#   geom_text_repel(aes(label=get(color)))+theme_classic2()+NoLegend()}
# 
# df0.split = split(df0,df0$Group)
# library(ggrepel)
# egfp = lapply(df0.split,function(x){
#   df=x
#   ggplot(df,aes(x=distance,y=EGFP,color = seurat_clusters))+
#     geom_smooth(se = F)+
#     # geom_text_repel(aes(label=seurat_clusters))+
#     theme_classic2()+NoLegend()
# })
# p1 = CombinePlots(egfp)
# car = lapply(df0.split,function(x){
#   df=x
#   ggplot(df,aes(x=distance,y=CAR,color = seurat_clusters))+
#     geom_smooth(se = F)+
#     # geom_text_repel(aes(label=seurat_clusters))+
#     theme_classic2()+NoLegend()
# })
# p2 = CombinePlots(car)
# pdf("./lineplot_seuratclusters.pdf",width = 10,height = 8)
# CombinePlots(plots = list(p1,p2),ncol=1)
# egfp
# car
# dev.off()

#-----subspot annotation -------
immune.markers = list("Tumor"=c("EGFP"),
                      "CAR-T"=c("Cd3d","mCherry-CAR"),
                      "Fibroblast"=c("Col1a1","Col3a1","Dcn"),
                      "CTL"=c("Cd3d","Gzmb","Gzma","Cd8b1","Prf1"),
                      "Te"=c("Cd3d","Cd8b1","Cd44"),
                      "Tex"=c("Cd3d","Cd8b1","Lag3","Pdcd1","Tim3","Tox","Tox2","Tcf7","Tcf1"),
                      "naiveT"=c("Cd3d","Sell"),
                      "Tcm"  =c("Cd3","Il7r","Cd44","Sell"),
                      "activate T"=c("Cd3d","Cd8b1","Cd69"),
                      "Treg"=c("Cd3d","Cd4","Foxp3","Il2ra"),
                      "Th1"=c("Cd3d","Cd4","Tbx21","Ifng"),
                      "Th2"=c("Cd3d","Cd4","Gata3","Il4"),
                      "Th9"=c("Cd3d","Cd4","Il9r","Spi1"),
                      "Th17"=c("Cd3d","Cd4","Il17a"),
                      "Tfh"=c("Cd3d","Cd4","Bcl6","Cxcr5","Il21r"),
                      "Th22"=c("Cd3d","Cd4","Il22"),
                      "NK"=c('CD3d','Nkg7','Klrb1c',"Fcgr3","Ncam1"),
                      "Neutrophil"=c("Ly6g","Itgam","Il1r2"),
                      "Monocyte"=c("Itgam","Fcgr3","Cd14"),
                      "DC"=c("Itgax","Xcr1","Cd74","Clec9a","Siglech"),
                      "M1"=c("Itgam","Adgre1","Cd80","Cd86"),
                      "M2"=c("Itgam","Adgre1",'Cd163','Arg1','Mrc1'),
                      "Bcell"=c("Cd40",'Cd19',"Ighm",'Cd79a','Cd79b','Vpreb3'),
                      "Erthroid"=c("Hbb-bs","Hbb-bt","Tfrc","Hba-a1"),
                      "HSPC"=c("Kit",'Ly6a',"Cd34","Cd38"))



pdf("Featureplot_immune_marker_sct_npc30_res1_4556cells.pdf",width = 5,height = 4)
FeaturePlot(immune.combined,features = unique(unlist(immune.markers)),combine = F,
            order = T,label = T,cols = rcolors::rcolors$NCV_jet)
dev.off()
df= immune.combined@meta.data
df$car =as.numeric(immune.combined@assays$SCT@counts["mCherry-CAR",])
pdf("glmplot_cart_seuratclusters_4556.pdf",width =6,height = 2.5)
ggplot(df,aes(x=distance,y=car,color=seurat_clusters))+geom_smooth(se=F,size=0.5)+
  scale_color_manual(values = subtypecolor)+theme_classic()+facet_grid(.~Group)+ylab("Cmax of CAR-T")
dev.off()

pdf("dotplot_immunemarkers.pdf")
lapply(immune.markers,function(x){
  DotPlot(immune.combined,features =x)+ggtitle(names(x))
})
dev.off()
#-----for celltype annotation ------


cellannotation.markers.l13 = FindAllMarkers(immune.combined,only.pos = T,min.pct = 0.5)

cellannotation.markers.immune = cellannotation.markers.l13[cellannotation.markers.l13$gene%in%immune.markers,]

"0"="Tumor" #
"1"="Tumor" #leukocyte migration
"2"="Tumor/Tcell" #granzyme−mediated programmed cell death signaling pathway
"4"="Fibroblast/Mono/T"
"5"="Fibroblast/Macrophage"
"7"="Tumor/Mono/DC"
"9"="Fibroblast"
"10"="Fibroblast"
"11"="Fibroblast/TAM"
"12"="Fibroblast/HSPC"
"13"="Tumor/Th1"
"14"="Tumor"

plt = lapply(1:length(immune.markers),function(i){
  features = immune.markers[[i]]
  p=DotPlot(immune.combined,features = features)
  df.dot =p$data
  head(df.dot)
  p0=ggplot(df.dot,aes(x=id,y=features.plot,fill=avg.exp.scaled,size = pct.exp))+
    geom_point(shape=21,color="black")+
    scale_fill_gradientn(colours = rev(rcolors$RdBu))+
    ggpubr::theme_pubr(border = T)+ylab(names(immune.markers)[i])+xlab("")
  
  return(p0)
  
})

figheight_length=lapply(plt,function(x)length(x$data[[3]]%>%unique())) %>%unlist

figheight =ifelse(figheight_length>10,6,
                  ifelse(figheight_length>5,4, ifelse(figheight_length>3,3,2)))

pdf("dotplot.expmarker_4556cells.pdf",height = 18)
ggarrange(plotlist =plt,ncol = 1,align = c("hv"),common.legend = T,
          heights =figheight,legend = "right")
dev.off()



#--- for functionEnrichment/GO--------
featuresmarkers.l13 = FindAllMarkers(immune.combined,logfc.threshold = 0)
write.csv(featuresmarkers.l13,"feauresmarkers_sct_npc30_res1_4556.csv")
library(clusterProfiler)
library(org.Mm.eg.db)
featuresmarkers.l13forgo = featuresmarkers.l13[featuresmarkers.l13$p_val_adj<0.05&featuresmarkers.l13$avg_log2FC>0,]
features.split = split(featuresmarkers.l13forgo$gene,featuresmarkers.l13forgo$cluster)
go.split = lapply(features.split,function(x){
  go =enrichGO(x,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "All",pvalueCutoff =1)
})

lapply(go.split,function(x){
  barplot(x,showCategory=15)
})
dev.off()

# Level_2
setwd("~/LJX/st")
dir.create("L2_cluster")
setwd("L2_cluster")
dir.create("L13_cluster")
setwd("L13_cluster")
immune.anchors <- FindIntegrationAnchors(object.list = objlist, anchor.features = features)
gc()
immune.combined <- IntegrateData(anchorset = immune.anchors) #require min.cell >=3 
immune.l2=immune.combined
rm(immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE) 
# immune.combined <- RunPCA(immune.combined)
immune.combined = basicFindCluster(immune.combined,npc = 30,res = 1)

group = rep("",nrow(immune.combined@meta.data))

group[grep("1$", rownames(immune.combined@meta.data))] = rep("Non-response",grep("1$", rownames(immune.combined@meta.data),value = T) %>%length())
group[grep("2$", rownames(immune.combined@meta.data))] = rep("Response",grep("2$", rownames(immune.combined@meta.data),value = T) %>%length())

immune.combined = AddMetaData(immune.combined,metadata = group,col.name = "Group")
#cellcycling
s.gene = read.csv("~/LJX/mm/mm.s.cellcyclegenes.csv")
g2m.gene = read.csv("~/LJX/mm/mm.g2m.cellcyclegenes.csv")

removecc = function(marrow){
  s.genes = s.gene[[2]]
  g2m.genes=g2m.gene[[2]]
  marrow =ScaleData(marrow,genes= rownames(marrow))
  marrow=CellCycleScoring(marrow,s.features = s.genes,g2m.features = g2m.genes)
  
  # Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
  # phase
  marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
  DimPlot(marrow)
  
  marrow <- ScaleData(marrow, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))
  # Now, a PCA on the variable genes no longer returns components associated with cell cycle
  marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)
  # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
  marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
  DimPlot(marrow)
}


n9 <- SCTransform(n9, vars.to.regress = "percent.mt", verbose = FALSE)

n9 <-basicFindCluster(data = n9,npc = 30)
DimPlot(
  n9,label = T
)

FeaturePlot(n9,features = c("EGFP","mCherry-CAR"))
install.packages("ggsignif")
library(ggsignif)
library(ggpubr)
source("~/scripts/singlecell_gene_test.R")
myvlnplot = function(seuratobj,features,group,comps,cols,ncol){
  source("~/scripts/singlecell_gene_test.R")
  A <- singlecell_gene_test(SerautObj = seuratobj, 
                            genes.use = features,
                            group.by = group, 
                            comp = comps)
  anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
  anno_sig <- A$sig
  
  plots_violins <- VlnPlot(seuratobj, 
                           cols = cols,
                           pt.size = 0,
                           group.by =group,
                           features = features, 
                           ncol = ncol, 
                           log = FALSE,
                           combine = FALSE)
  
  for(i in 1:length(plots_violins)) {
    data <- plots_violins[[i]]$data
    colnames(data)[1] <- 'gene'
    plots_violins[[i]] <- plots_violins[[i]] + 
      theme_classic() + 
      theme(axis.text.x = element_text(size = 10,color="black"),
            axis.text.y = element_text(size = 10,color="black"),
            axis.title.y= element_text(size=12,color="black"),
            axis.title.x = element_blank(),
            legend.position='none')+
      # scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
      # scale_x_discrete(labels = comps)+
      geom_signif(annotations = anno_sig[i],
                  y_position = max(data$gene)+0.5,
                  xmin = 1,
                  xmax = 2,
                  tip_length = 0)
  }
  
  cplt=CombinePlots(plots_violins)
  return(cplt)
}

pdf("vlnplot_exhuast_cart.pdf",width = 12,height = 8)

myvlnplot(cart,features = exhuast,group = "sample",comps = levels(cart$sample),cols = ggsci::pal_jco()(2),ncol = 3)

dev.off()

pdf("featureplot_cart_infilter_tumor.pdf",width = 12,height = 6)

FeaturePlot(immune.combined,features = c("EGFP","mCherry-CAR"),
            blend = T,order = T,label = T,split.by="Group",cols=c("lightgrey","green","firebrick"))

dev.off()

allmarkers = FindAllMarkers(immune.combined,only.pos = T)

mm.cellmarkers.lab =readxl::read_xlsx("~/LJX/mm/Cell_marker_Seq.xlsx")

cellmarker = data.frame("gene"=mm.cellmarkers$Symbol,
                        "celltype"=mm.cellmarkers$cell_name,
                        "ref_type"=mm.cellmarkers$tissue_type,
                        "marker_source"= mm.cellmarkers$marker_source,
                        "PMID"=mm.cellmarkers$PMID,"journal"=mm.cellmarkers$journal,"year"=mm.cellmarkers$year)

cellmarker = left_join(allmarkers,cellmarker,"gene")

write.csv(cellmarker,"L2_cluster_allmarkers_Cellmarker.csv")


#------------AUCell---------------
BiocManager::install("mixtools")
BiocManager::install('AUCell')
library(AUCell)
library(GSEABase)
library(mixtools)
setwd("~/LJX/st/L13_cluster/")
geneset = list( "Tumor"=c("EGFP"),
                "CD8_TCell"=c("Cd3e","Cd8a"),
                "CD4_TCell"=c("Cd3e","Cd4"),
                "CAR_T"=c("mCherry-CAR","Cd3e"),
                "FuncT"=c("Foxp3","Tbx21","Gata1","Gata2","Gata3","Gata4","Il2ra","Gzmb","Gzma","Lag3","Pdcd1","Sell"),
                "Mono"=c("Itgam","Ly6g"),
                "Macrophage"=c("Siglech","Adgre1","Itgam","Itgax","Cd80","Cd86",'Cd163','Arg1','Mrc1'),
                "B_Cell"=c('Cd19','Cd79a','Cd79b','Vpreb3'))

mm.geneset = split(mm.cellmarkers$Symbol,mm.cellmarkers$cell_name)

geneSets = c(geneset,mm.geneset)

geneSets = c(immune.markers,mm.geneset)


exprMatrix = as.matrix(GetAssayData(immune.combined,assay = "SCT"))

cells_AUC <- AUCell_run(exprMatrix, geneSets)
save(cells_AUC, file="cells_AUC_l13_sct_4556.RData")

auc.mat = t(cells_AUC@assays@data$AUC) %>% as.data.frame()
immune.combined = AddMetaData(immune.combined,auc.mat)
pdf("AUC_UMAP_seuratlabel.pdf",width = 5,height = 4)
FeaturePlot(immune.combined,label = T,features=colnames(immune.combined@meta.data)[13:37],cols = rcolors::rcolors$GMT_jet,pt.size = 0.1,combine = F)
dev.off()


auc.mat = cbind(auc.mat,immune.combined@meta.data)
egfp = as.numeric(immune.combined@assays$RNA@counts["EGFP",])
car = as.numeric(immune.combined@assays$RNA@counts["mCherry-CAR",])
auc.mat = 


  
save(cells_AUC, file="cells_AUC_l13_sct_4556_immunemarkers.RData")

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank = ceiling(0.3 * nrow(cells_rankings)))
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
save(cells_rankings, file="cells_rankings_l13_sct_4556_immunemarkers.RData")

# exploreThresholds, warning=FALSE, fig.width=7, fig.height=7
set.seed(333)
par(mfrow=c(3,3)) 
pdf('AUC_hist_4556_immunemarkers.pdf')
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]
dev.off()

selectedThresholds = getThresholdSelected(cells_assignment)

pdf("AUCell_UMAP_mmMaker_immunemarkers_l13_sct_4556.pdf",width=8,height = 6)
par(mfrow=c(2,3))
coord=immune.combined@reductions$umap@cell.embeddings
AUCell_plotTSNE(tSNE = coord,exprMat=exprMatrix,
                cellsAUC = cells_AUC,thresholds =selectedThresholds)

dev.off()

pdf("AUCell_UMAP_passThreshold_immunemarkers_l13_sct_4556.pdf",width = 8, height = 6)
cellsTsne=immune.combined@reductions$umap@cell.embeddings
par(mfrow=c(2,3)) # Splits the plot into two rows and three columns
for(geneSetName in names(selectedThresholds)[selectedThresholds>0]){
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0 )
  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    
    # Plot
    plot(cellsTsne, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(cellsTsne)], pch=16) 
  }
}

dev.off()


pdf("AUCell_UMAP_passThreshold_auc0.3_l13_sct_4556.pdf",width = 8, height = 6)
cellsTsne=immune.combined@reductions$umap@cell.embeddings
par(mfrow=c(2,3)) # Splits the plot into two rows and three columns
for(geneSetName in names(selectedThresholds)[selectedThresholds>0.3]){
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0 )
  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    
    # Plot
    plot(cellsTsne, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(cellsTsne)], pch=16) 
  }
}

dev.off()


pdf("AUCell_UMAP_mmMaker_auc0.09_l13_sct_4556.pdf",width=8,height = 6)
par(mfrow=c(2,3))
coord=immune.combined@reductions$umap@cell.embeddings
AUCell_plotTSNE(tSNE = coord,exprMat=exprMatrix,
                cellsAUC = cells_AUC[selectedThresholds>0.09,],thresholds =selectedThresholds[selectedThresholds>0.09])

dev.off()

n=2
path = n9_path
setwd(path)
setwd(levelpath[n])
pos = readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
colnames(pos)=c("id","x","y")
pos$id = paste0(pos$id,"_1")
n9_pos = pos

coord = n9_pos
exp0 = as.data.frame(t(exprMatrix))
exp0$id = rownames(exp0)
exp0 = left_join(n9_pos,exp0,"id")
exp0 =data.frame(exp0,row.names = 1)
exp0 =as.matrix(t(exp0[,3:ncol(exp0)]))
coord = data.frame(coord,row.names = 1) %>% as.matrix()
setwd("~/LJX/st/L13_cluster/")
pdf("AUCell_spatial_n9_mmMaker_auc0.09.pdf",width=8,height = 6)
par(mfrow=c(2,3))
AUCell_plotTSNE(tSNE = coord,exprMat=exp0,
                cellsAUC = cells_AUC[selectedThresholds>0.09,],thresholds =selectedThresholds[selectedThresholds>0.09])
dev.off()


path = r5_path
setwd(path)
setwd(levelpath[n])
pos = readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
colnames(pos)=c("id","x","y")
pos$id = paste0(pos$id,"_2")
r5_pos = pos

coord = r5_pos
exp0 = as.data.frame(t(exprMatrix))
exp0$id = rownames(exp0)
exp0 = left_join(r5_pos,exp0,"id")
exp0 =data.frame(exp0,row.names = 1)
exp0 =as.matrix(t(exp0[,3:ncol(exp0)]))
coord = data.frame(coord,row.names = 1) %>% as.matrix()
setwd("~/LJX/st/L13_cluster/")
pdf("AUCell_spatial_r5_mmMaker_auc0.09.pdf",width=8,height = 6)
par(mfrow=c(2,3))
AUCell_plotTSNE(tSNE = coord,exprMat=exp0,
                cellsAUC = cells_AUC[selectedThresholds>0.09,],thresholds =selectedThresholds[selectedThresholds>0.09])
dev.off()


#annotated by single cell cellmarkers 

load("~/LJX/sc_in/adata_preprocess.RData")

active_ident = adata$Cell_type
names(active_ident)=colnames(adata)
adata@active.ident = active_ident

adata_marker = FindAllMarkers(adata,logfc.threshold = 1,only.pos = T)
geneSets= split(adata_marker$gene,adata_marker$cluster)
geneSets$`CD4T Cell`=c("Cd3e","Cd3g","Cd4",geneSets$`CD4T Cell`)
geneSets$`CD8T Cell` = c("Cd3e","Cd3g","Cd8a","Cd8b1",geneSets$`CD8T Cell`)
geneSets$`T Cell MKI67+`=c("Cd3e","Cd3g","Mki67",geneSets$`T Cell MKI67+`)
geneSets$`T Cell TCF7+` = c("Cd3e","Cd3g","Tcf7",geneSets$`T Cell TCF7+`)
cells_AUC <- AUCell_run(exprMatrix, geneSets)


#------------------------- end AUCell ------------------------------------

#------------------------- Cell type marker annotation--------------

celltype = c("Immune/Tumor",
            "Myeloid/DC",
             "Neutrophil/Bcell",
             "Tcell/Stromal Cell",
             "Fibroblast",
             "TAM/Fibroblast",
             "Neutrophil",
             "Endothelial cell",
             "Myeloid/Endothelial cell",
             "Myeloid/DC",
             "Neutrophil/Epithelial cell",
            " Neutrophoil/Erythoid",
             "M1/Erythoid",
             "Macrophage/MPC",
             "TAM",
             "CD4T/Macrophage/Tumor",
             "Tcell/Stromal Cell",
            "Myeloid/DC")
names(celltype) = levels(immune.combined@active.ident)
immune.combined=RenameIdents(immune.combined,celltype)
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)

lapply(ArchRPalettes,function(x){
  length(x)
})

celltype_colors=ArchRPalettes$ironMan
celltypes = sort(unique(celltype))
celltypes = celltypes[c(5,15,2,8,9,10,11,12,1,7,6,13,14,3,4)]

names(celltype_colors)=celltypes

pdf("umap_Combined_L2_celltype.pdf",width = 8,height = 6)
DimPlot(immune.combined,cols = celltype_colors,repel = T,label = T,label.color = "white")
dev.off()
pdf("umap_Combined_L2_celltype_splited.pdf",width = 12,height = 6)
DimPlot(immune.combined,cols = celltype_colors,split.by = "Group",repel = T,label = T,label.color = "white")
dev.off()


cellmeta = immune.combined@meta.data
cellmeta$celltype_cellmarker = immune.combined@active.ident

cellmeta.split = split(cellmeta,cellmeta$Group)
n=2
path = n9_path
setwd(path)
n9=Read10X(levelpath[n])
egfp = as.numeric(n9["EGFP",])
car = as.numeric(n9["mCherry-CAR",])
setwd(levelpath[n])
pos = readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
colnames(pos)=c("id","x","y")
pos$EGFP = egfp
pos$mCherry = car
st_distance = function(x,y){
  central_x = (max(x)-min(x))/2
  central_y =(max(y)-min(y))/2
  distance = sqrt((x-central_x)^2 + (y-central_y)^2)
}
distance = st_distance(x=pos$x,y=pos$y)
pos$distance = distance
sp = cellmeta.split[[1]]
sp$id =gsub("_1$","",rownames(sp))

# pos = left_join(pos,sp,"id")

setwd("~/LJX/st/L2_cluster")

pdf("spatial_n9_celltype.pdf",width = 8,height = 6)

p_n5=ggplot(pos,aes(x=x,y=y))+geom_point(size=3,color="lightgrey")+
  xlab("Spatial_X")+ylab("Spatial_Y")+
  geom_point(data=sp,aes(x=spatial_x,y=spatial_y,color=seurat_clusters),size=3)+
  # geom_point(data=sp,aes(x=spatial_x,y=spatial_y,color=celltype_cellmarker),size=0.001)+
  # scale_color_manual(values = celltype_colors)+
  cowplot::theme_map()+ggtitle("spatial_L13_celltype_response")

dev.off()
ggsave(p_n9,filename = "spatial_n9_celltype.png",width = 12,height = 10,limitsize = FALSE,dpi=300)

pos = left_join(pos,sp,"id")
table(pos$celltype_cellmarker)
pos$celltype_cellmarker = factor(pos$celltype_cellmarker,levels = names(celltype_colors))
write.csv(pos,"N9_L2_pos_meta.csv")
pdf("smoothplot_N9_celltype_cart.pdf",width=8,height=6)
smplt_np = ggplot(pos,aes(x=distance,y=mCherry,color=celltype_cellmarker))+geom_smooth()+
  scale_color_manual(values = celltype_colors)+theme_classic()
dev.off()

path = r5_path
setwd(path)
r5=Read10X(levelpath[n])
egfp = as.numeric(r5["EGFP",])
car = as.numeric(r5["mCherry-CAR",])

setwd(levelpath[n])
pos = readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
colnames(pos)=c("id","x","y")
pos$EGFP = egfp
pos$mCherry = car
distance = st_distance(x=pos$x,y=pos$y)
pos$distance = distance

sp = cellmeta.split[[2]]
sp$id =gsub("_2$","",rownames(sp))
pos = left_join(pos,sp,"id")
setwd("~/LJX/st/L2_cluster")

p_r5=ggplot(pos,aes(x=x,y=y))+geom_point(size=0.0001,color="lightgrey")+
  xlab("Spatial_X")+ylab("Spatial_Y")+
  geom_point(data=sp,aes(x=spatial_x,y=spatial_y,color=celltype_cellmarker),size=0.0001)+
  scale_color_manual(values = celltype_colors)+
  cowplot::theme_map()+ggtitle("spatial_L2_celltype_response")
ggsave(p_r5,filename = "spatial_r5_celltype_response.png",width = 10,height = 6.5,limitsize = FALSE,dpi=300)
write.csv(pos,"R5_L2_pos_meta.csv")

dreference <- Reference(counts, cell_types, nUMI)

puck <- SpatialRNA(coords, counts, nUMI)

myRCTD <- create.RCTD(puck, reference, max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
str(myRCTD)

#https://ccsm.uth.edu/SPASCER/download.html database


#annotated with STdeconvolve
remotes:: install_github ('JEFworks-Lab/STdeconvolve')
library(STdeconvolve)
#example
# data(mOB)
# pos <- mOB$pos
# cd <- mOB$counts
# counts <- cleanCounts(cd, min.lib.size = 100)
# corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
# ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
# optLDA <- optimalModel(models = ldas, opt = 3)
# results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
# deconProp <- results$theta
# corMtx <- getCorrMtx(m1 = as.matrix(deconProp), m2 = as.matrix(deconProp), type = "t")
# rownames(corMtx) <- paste0("X", seq(nrow(corMtx)))
# colnames(corMtx) <- paste0("X", seq(ncol(corMtx)))
# correlationPlot(mat = corMtx, title = "Proportional correlation", annotation = TRUE) +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))

#annotated with RCTD
devtools::install_github("dmcable/spacexr",build_vignettes = False)
devtools::install_github("dmcable/RCTD", build_vignettes = False)


# 
# #ggfig
# n9_l2=Read10X("./L2_heAuto/")
# n9_pos=readr::read_tsv(gzfile("./L2_heAuto/barcodes_pos.tsv.gz"),col_names = F)
# 
# pdf("N9_L2_TUMOR_BURDERN.pdf")
# df$mCherry =as.numeric(n9_l2['mCherry-CAR',])
# ggplot(df,aes(x=X2,y=X3,color=EGFP))+geom_point(size=0.1)+
#        geom_polygon(data=df,aes(x=X2,y=X3,fill=mCherry),stat="density_2d", alpha = .3, color = NA)+
#   scale_fill_viridis_c(option = "A")+scale_color_viridis_c()+cowplot::theme_map()\
# 
# ggplot(df,aes(x=X2,y=X3,fill=EGFP))+
#        geom_polygon(stat="density_2d", alpha = .3, color = NA)+scale_fill_viridis_c()+cowplot::theme_map()
# 
# 
# dev.off()


sc_nUMI = colSums(sc_counts)

reference = Reference(sc_counts, cell_types, sc_nUMI)

#读取空间数据，这里选择的是level 13的数据

#位置信息
path = ""
coords = read.table(gzfile(paste0(path,'/05.AllheStat/heAuto_level_matrix/subdata/L13_heAuto/barcodes_pos.tsv.gz')),    header = F) %>% dplyr::rename(barcodes = V1, xcoord = V2, ycoord = V3)

rownames(coords) <- coords$barcodes; coords$barcodes <- NULL

#表达量矩阵

expr <- Read10X('/xxx/05.AllheStat/BSTViewer_project/subdata/L13_heAuto/', cell.column = 1)

sp_data <- CreateSeuratObject(counts = expr,assay = "Spatial")

sp_counts <- as_matrix(sp_data[['Spatial']]@counts)

#nUMI

sp_nUMI <- colSums(sp_counts)

#构建空间实验集

puck <- SpatialRNA(coords, sp_counts, sp_nUMI)


#Annoataed by spolight/spodeconvolution

library(STdeconvolve)
pos = readr::read_tsv(gzfile("L13_heAuto/barcodes_pos.tsv.gz"),col_names = F)
barcode = pos$X1
pos = data.frame(x=pos$X2,y=pos$X3,row.names = pos$X1)
cd =Read10X("./L13_heAuto/")
cd =as.matrix(cd)
car= cd["mCherry-CAR",]
EGFP =cd["EGFP",]
celltype = ifelse(EGFP>0,"Tumor","Non-Tumor")
celltype = ifelse(car>0,"CAR-T",celltype)
annot =factor(celltype,levels=unique(celltype))
central_x = (max(pos$x)-min(pos$x))/2
central_y = (max(pos$y)-min(pos$x))/2
pos$distance = sqrt((pos$x-central_x)^2 + (pos$y-central_y)^2) 
pos$CAR=car
pos$group = rep("Non-Response",nrow(pos))
pos_merge=rbind(pos,r5_pos)
ggplot(pos_merge,aes(x=distance,y=CAR,color=group))+geom_smooth()+scale_color_brewer(palette = "Set1")+
  ggpubr::theme_pubr(border =T )+
  ylab("Cmax of CAR-T")+ggtitle("Distribution of CAR-T infiltering in two groups")

counts <- cleanCounts(counts=cd,
                      min.lib.size = 100,
                      min.reads = 1,
                      min.detected = 1,verbose = T
                      )

corpus <- restrictCorpus(counts,
                         removeAbove = 1.0,
                         removeBelow = 0.05,
                         alpha = 0.05,
                         plot = TRUE,
                         verbose = TRUE)

ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)

optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA,
                        perc.filt = 0.05,
                        betaScale = 1000)


deconProp <- results$theta
deconGexp <- results$beta


# proxy theta for the annotated layers
mobProxyTheta <- model.matrix(~ 0 + annot)
rownames(mobProxyTheta) <- names(annot)
# fix names
colnames(mobProxyTheta) <- unlist(lapply(colnames(mobProxyTheta), function(x) {
  unlist(strsplit(x, "annot"))[2]
}))

mobProxyGexp <- counts %*% mobProxyTheta

corMtx_beta <- getCorrMtx(# the deconvolved cell-type `beta` (celltypes x genes)
  m1 = as.matrix(deconGexp),
  # the reference `beta` (celltypes x genes)
  m2 = t(as.matrix(mobProxyGexp)),
  # "b" = comparing beta matrices, "t" for thetas
  type = "b")

rownames(corMtx_beta) <- paste0("decon_", seq(nrow(corMtx_beta)))

correlationPlot(mat = corMtx_beta,
                # colLabs (aka x-axis, and rows of matrix)
                colLabs = "Deconvolved cell-types",
                # rowLabs (aka y-axis, and columns of matrix)
                rowLabs = "Ground truth cell-types",
                title = "Transcriptional correlation", annotation = TRUE) +
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))


corMtx_theta <- getCorrMtx(# deconvolved cell-type `theta` (pixels x celltypes)
  m1 = as.matrix(deconProp),
  # the reference `theta` (pixels x celltypes)
  m2 = as.matrix(mobProxyTheta),
  # "b" = comparing beta matrices, "t" for thetas
  type = "t")

rownames(corMtx_theta) <- paste0("decon_", seq(nrow(corMtx_theta)))

correlationPlot(mat = corMtx_theta,
                # colLabs (aka x-axis, and rows of matrix)
                colLabs = "Deconvolved cell-types",
                # rowLabs (aka y-axis, and columns of matrix)
                rowLabs = "Ground truth cell-types",
                title = "Proportional correlation", annotation = TRUE) +
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))

pairs <- lsatPairs(t(corMtx_theta))
m <- t(corMtx_theta)[pairs$rowix, pairs$colsix]

correlationPlot(mat = t(m), # transpose back
                # colLabs (aka x-axis, and rows of matrix)
                colLabs = "Deconvolved cell-types",
                # rowLabs (aka y-axis, and columns of matrix)
                rowLabs = "Ground truth cell-types",
                title = "Transcriptional correlation", annotation = TRUE) +
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))

dev.off()

pdf("r5_deconvolve.pdf",width =6,height = 4)
vizAllTopics(deconProp, pos, 
             groups = annot, 
             group_cols = rainbow(length(levels(annot))),
             r=5)
dev.off()


distance


#------annotated with ssGSEA -------------------------
mobProxyLayerMarkers <- list()

## make the tissue layers the rows and genes the columns
gexp <- t(as.matrix(mobProxyGexp))

for (i in seq(length(rownames(gexp)))){
  celltype <- i
  ## log2FC relative to other cell-types
  ## highly expressed in cell-type of interest
  highgexp <- names(which(gexp[celltype,] > 10))
  ## high log2(fold-change) compared to other deconvolved cell-types and limit to top 200
  log2fc <- sort(log2(gexp[celltype,highgexp]/colMeans(gexp[-celltype,highgexp])), decreasing=TRUE)[1:200]
  
  ## for gene set of the ground truth cell-type, get the genes
  ## with log2FC > 1 (so FC > 2 over the mean exp of the other cell-types)
  markers <- names(log2fc[log2fc > 1])
  mobProxyLayerMarkers[[ rownames(gexp)[celltype] ]] <- markers
}


#------------match cluster point on png ----------------------------

myclstplot = function(he_file,FilePath,
                      Cluster,
                      celltypecol = NULL,pt.size=3.5){
  
  he_fig=png::readPNG(he_file)
  w=ncol(he_fig)
  h=nrow(he_fig)
  
  cal_zoom_rate = function(width, height){
    std_width = 1000
    std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
    if(std_width / std_height > width / height){
      scale = width / std_width
    }
    else{
      scale = height / std_height
    }
    return(scale)
  }
  
  zoom_scale = cal_zoom_rate(ncol(he_fig), nrow(he_fig))
  bc_pos_file = gzfile(paste(FilePath,"barcodes_pos.tsv.gz", sep  = "/"),'rt')
  bc_pos = read.table(bc_pos_file, header = FALSE, sep = '\t', quote = '')
  names(bc_pos)=c("Barcode","x","y")
  Cluster = data.frame(cluster=Cluster)
  Cluster$Barcode = rownames(Cluster)
  bc_pos=merge(bc_pos,Cluster,by="Barcode",all=FALSE)
  
  bc_pos$x=bc_pos$x * zoom_scale
  bc_pos$y=bc_pos$y * zoom_scale
  
  clstr_plot<-ggplot(data=bc_pos,aes(x ,dim(he_fig)[1]-y)) + 
    background_image(he_fig)+
    # geom_point(aes(colour = cluster),size = opt$point_size, shape=16) + 
    geom_point(aes(fill = cluster,color=cluster),size = pt.size, shape=21) + 
  # scale_fill_brewer(palette = "Set1")+
  # scale_color_manual(values = col)+
  theme_bw() +
    theme(plot.title = element_text(face = 2,size = 50,hjust = 0.5)) +
    theme(axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          panel.border = element_blank(), 
          axis.text.x = element_blank()) +
    xlab('')+ylab('')+
    coord_cartesian(xlim = c(0, dim(he_fig)[2]), ylim = c(0, dim(he_fig)[1]), 
                    expand = FALSE)+
    guides(colour = guide_legend(override.aes = list(size=3.5)))+
    theme(legend.position = "none")
  
  if(is.null(celltypecol)){
    clstr_plot = clstr_plot
  }else{
      clstr_plot = clstr_plot+scale_fill_manual(values= celltypecol)+scale_color_manual(values= celltypecol)
      }
  
  return(clstr_plot)
} 

#level_13 non-responese
clusterlist = split(immune.combined@active.ident,immune.combined$Group)
he_file = "~/LJX/st/n9/05.AllheStat/allhe/he_roi_small.png"
FilePath = paste(n9_path,levelpath[2],sep="/")
cluster = clusterlist$`Non-response`
names(cluster)=gsub("_1$","",names(cluster))

he_fig=png::readPNG(he_file)
w=ncol(he_fig)
h=nrow(he_fig)
setwd("~/LJX/st/L13_cluster/")
pdf("he_cluster_n9_sct_4556.pdf",width = w/100,height =h/100 )
clstr_plot = myclstplot(he_file =he_file,Cluster = cluster,FilePath = FilePath,celltypecol = subtypecolor )
clstr_plot
dev.off()

#level_13 response
he_file = "~/LJX/st/r5/05.AllheStat/allhe/he_roi_small.png"
FilePath = paste(r5_path,levelpath[2],sep="/")
cluster = clusterlist$Response
names(cluster)=gsub("_2$","",names(cluster))

he_fig=png::readPNG(he_file)
w=ncol(he_fig)
h=nrow(he_fig)
clstr_plot = myclstplot(he_file =he_file,Cluster = cluster,FilePath = FilePath,celltypecol = subtypecolor )
setwd("~/LJX/st/L13_cluster/")
pdf("he_cluster_r5_4556.pdf",width = w/100,height =h/100 )
clstr_plot
dev.off()

#level_2 non-response

clusterlist = split(immune.l2@active.ident,immune.l2$Group)
he_file = "~/LJX/st/n9/05.AllheStat/allhe/he_roi_small.png"
FilePath = paste(n9_path,levelpath[3],sep="/")
cluster = clusterlist$`Non-response`
names(cluster)=gsub("_1$","",names(cluster))


setwd("~/LJX/st/L2_cluster/")
clstr_plot = myclstplot(he_file =he_file,Cluster = cluster,
                        FilePath = FilePath,
                        celltypecol = celltype_colors,pt.size = 0.4)
he_fig=png::readPNG(he_file)
w=ncol(he_fig)
h=nrow(he_fig)
png("he_cluster_n9_L2.png",width = w*5,height =h*5,res = 300)
clstr_plot
dev.off()

#level_2 response

he_file = "~/LJX/st/r5/05.AllheStat/allhe/he_roi_small.png"
FilePath = paste(r5_path,levelpath[3],sep="/")
cluster = clusterlist$Response
names(cluster)=gsub("_2$","",names(cluster))

clstr_plot = myclstplot(he_file =he_file,Cluster = cluster,
                        FilePath = FilePath,
                        celltypecol = celltype_colors,pt.size = 0.4)
setwd("~/LJX/st/L2_cluster/")
he_fig=png::readPNG(he_file)
w=ncol(he_fig)
h=nrow(he_fig)
png("he_cluster_r5_L2.png",width = w*5,height =h*5,res = 300)
clstr_plot
dev.off()


#--------------CAR-T subtypes reculstering ---------------
#level_2
subspot =immune.l2@assays$RNA@scale.data["mCherry-CAR",]
subspot = names(subspot)[subspot>0]
immune.l2.sub = subset(immune.l2,cells=subspot)
immune.l2.sub = FindVariableFeatures(immune.l2.sub,nfeatures = 10000)
immune.l2.sub = SCTransform(immune.l2.sub,variable.features.n = 10000)

immune.l2.sub =ScaleData(immune.l2.sub)
immune.l2.sub = basicFindCluster(immune.l2.sub,npc = 10)

markers = FindAllMarkers(immune.l2.sub,logfc.threshold = 0)

write.csv(markers,"cart_fc.csv")

myclstplot()

#level_2 and surrounding 7-point


#---Fig1 reorder generel EGFP and CART--------
setwd('~/LJX/st/L13_cluster/')
setwd(n9_path)

coord.n9 = gzfile("../../level_matrix/level_13/barcodes_pos.tsv.gz",'rt')
coord.n9 = read.table(coord.n9,header = FALSE,sep="\t",quote = "")
L13_st.n9 = Read10X("./L13_heAuto/")

df=L13_st.n9[c("mCherry-CAR",'EGFP','Cdh5','Vegfa'),] %>% t() %>%as.data.frame()

df$mV =ifelse(df$Cdh5*df$Vegfa>0,"mV","other")
df$label = ifelse(df$EGFP*df$`mCherry-CAR`>0,"CART-Tumor",ifelse(df$`mCherry-CAR`>0,"CAR-T",
                                                                 ifelse(df$EGFP>0,"Tumor",df$mV)))

table(df$label)
table(df$mV)

he_fig = png::readPNG("../../allhe/he_roi_small.png")
myclstplot = function(he_fig,coord,
                      metadata,fill,color,filename=NULL,
                      celltypecol = NULL,pt.size=3.5){
  # 
  # if(hetype=="png"){
  #   he_fig=png::readPNG(he_file)
  # }if(hetype=='tiff'){
  #   he_fig=tiff::readTIFF(he_file)
  # }
  # 
  w=ncol(he_fig)
  h=nrow(he_fig)
  
  cal_zoom_rate = function(width, height){
    std_width = 1000
    std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
    if(std_width / std_height > width / height){
      scale = width / std_width
    }
    else{
      scale = height / std_height
    }
    return(scale)
  }
  
  zoom_scale = cal_zoom_rate(ncol(he_fig), nrow(he_fig))
  # bc_pos_file = gzfile(paste(FilePath,"barcodes_pos.tsv.gz", sep  = "/"),'rt')
  # bc_pos = read.table(bc_pos_file, header = FALSE, sep = '\t', quote = '')
  colnames(coord)=c("Barcode","x","y")
  # Cluster = data.frame(cluster=Cluster)
  # colnames(Cluster)=c("cluster")
  metadata$Barcode = rownames(metadata)
  bc_pos=merge(coord,metadata,by="Barcode",all=FALSE)
  
  bc_pos$x=bc_pos$x * zoom_scale
  bc_pos$y=bc_pos$y * zoom_scale
  bc_pos$y =dim(he_fig)[1]-bc_pos$y
  
  clstr_plot<-ggplot(data=bc_pos,aes(x ,y)) + 
    ggpubr::background_image(he_fig)+
    # geom_point(aes(colour = cluster),size = opt$point_size, shape=16) + 
    geom_point(aes(fill = fill,color=color),size = pt.size, shape=21) + 
    # scale_fill_brewer(palette = "Set1")+
    # scale_color_manual(values = col)+
    theme_bw() +
    theme(plot.title = element_text(face = 2,size = 50,hjust = 0.5)) +
    theme(axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          panel.border = element_blank(), 
          axis.text.x = element_blank()) +
    xlab('')+ylab('')+
    coord_cartesian(xlim = c(0, dim(he_fig)[2]), ylim = c(0, dim(he_fig)[1]), 
                    expand = FALSE)+
    guides(colour = guide_legend(override.aes = list(size=3.5)))+
    theme(legend.position = "none")
  
  # if(is.null(celltypecol)){
  #   clstr_plot = clstr_plot
  # }else{
  #   clstr_plot = clstr_plot+scale_fill_manual(values= celltypecol)+
  #     scale_color_manual(values= celltypecol)
  # }
  # 
  
  if(!is.null(filename)){
    write.csv(bc_pos,filename)
    return(clstr_plot)
  }else{
    print(clstr_plot)
    return(bc_pos)
  }
  
} 

n9_bc_pos = myclstplot(he_fig = he_fig,coord = coord.n9,metadata = df,fill =df$label,color=df$mV)
df0=df
df=df[df$mV=="mV",]
df=df0
df$label = factor(df$label,levels=c("CAR-T","Tumor","CART-Tumor","mV","other"))
p = myclstplot(he_fig = he_fig,coord = coord.n9,metadata =df,fill =df$label,pt.size = 2,color=df$mV,filename ="he_coord_anno_n9_mV.csv")

mVcolors = c("lightgrey","firebrick")
names(mVcolors)=unique(df$mV)
pdf("N9_general_he.pdf",width = 6,height = 6)
p=p+scale_fill_brewer(palette = "Paired")+scale_color_manual(values = mVcolors)+theme(legend.position = c(200,300))
LabelClusters(p,id="label",box = T,position = "median",repel = T)+scale_fill_brewer(palette = "Set1")
dev.off()

st_distance = function(x,y){
  central_x = (max(x)-min(x))/2
  central_y =(max(y)-min(y))/2
  distance = sqrt((x-central_x)^2 + (y-central_y)^2)
}

df=df0[df0$EGFP>0,]
df=df0
p = myclstplot(he_fig = he_fig,coord = coord.n9,metadata =df,fill =as.numeric(df$EGFP),pt.size = 2,color=NULL,filename ="he_coord_anno_n9tumor.csv")
pdf("N9_general_he_tumor.pdf",width = 6,height = 6)
p+scale_fill_viridis_c()+scale_color_manual(values = NA)
# LabelClusters(p,id="label",box = T,position = "median",repel = T)+scale_fill_brewer(palette = "Set1")

dev.off()

pdf("N9_general_he_cart.pdf",width = 6,height = 6)
p = myclstplot(he_fig = he_fig,coord = coord.n9,metadata =df,fill =as.numeric(df$`mCherry-CAR`),pt.size = 2,color=NULL,filename ="he_coord_anno_n9cart.csv")
p+scale_fill_viridis_c(option = "A")+scale_color_manual(values = NA)
# p
dev.off()
