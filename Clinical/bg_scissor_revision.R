library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(homologene)
library(Scissor)
source('~/scripts/Scissor5.R')
source('~/LJX/code/schard-main/R/functions.R')
source('~/LJX/code/schard-main/R/h5ad_util.R')

setwd('/home/Data/ST_20241110/spatial_input/h5ad/recluster_obj/')
output='Scissor_rdata'
dir.create(output)
workpath = '/home/Data/ST_20241110/spatial_input/h5ad/recluster_obj'
filelist = dir('./') 
filelist =grep('_sc2st_DestVI_destvi_recluster.h5ad',filelist,value = T)
filelist=  grep('RData',filelist,invert = T,value = T)


nr2 = h5ad2seurat('./NR2_sc2st_DestVI_destvi_recluster.h5ad')
nr2$barcode =paste0(nr2$orig.ident,nr2$sample)
genes = rownames(nr2)
mouse2human = homologene(genes, inTax = 10090, outTax = 9606)

mm2hm=function(obj){
  meta.features = obj@assays$RNA@meta.features
  meta.features$`10090`=meta.features$gene_names
  meta.features.homo = left_join(meta.features,mouse2human,"10090")
  meta.features.homo =meta.features.homo[!duplicated(meta.features.homo$`_index`),]
  rownames(meta.features.homo)=meta.features.homo$gene_names
  obj@assays$RNA@meta.features=meta.features.homo
  count=FetchData(obj,vars = rownames(meta.features.homo)[!duplicated(meta.features.homo$`9606`)])
  colnames(count)=meta.features.homo$`9606`[!duplicated(meta.features.homo$`9606`)]
  count=t(count)
  count = count[na.omit(rownames(count)),]
  st_dataset= Seurat_preprocessing(count, verbose = F,resolution = 1)
  return(st_dataset)
}



add_scissorlab=function(sc_dataset,infos5,cols="scissor"){
  Scissor_select <- rep(0, ncol(sc_dataset))
  names(Scissor_select) <- colnames(sc_dataset)
  Scissor_select[infos5$Scissor_pos] <- 'NDR'
  Scissor_select[infos5$Scissor_neg] <- 'DR'
  sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name =cols )
  return(sc_dataset)
}


# dir.create("./Scissor_rdata"
library(data.table)

bulk_dir <- "/home/Data/ST_20241110/spatial_input/h5ad/scissor_no_relapse_exprllist/bulk"

bulk_files <- dir(
  bulk_dir,
  pattern = "^bulk_.*\\.txt$",
  full.names = TRUE
)

names(bulk_files) <- gsub("^bulk_|\\.txt$", "", basename(bulk_files))

# dir.create("./Scissor_rdata", showWarnings = FALSE, recursive = TRUE)

meta_merge.sub = read.csv('/home/Data/ST_20241110/spatial_input/h5ad/scissor_no_relapse_exprllist/pheno.csv')

meta_merge.sub$label=factor(meta_merge.sub$label,levels=c('DR','NDR'))

phenotypeall <- as.numeric(meta_merge.sub$label) - 1

tag <- c("DR", "NDR")
names(phenotypeall)=meta_merge.sub$sample_id

# infos2_hm2mm_list = lapply(1:length(filelist),function(i){
infos2_hm2mm_list <- lapply(names(bulk_files), function(cohort_name) {
  
  message("========== Running cohort: ", cohort_name, " ==========")
  
  bulk_file <- bulk_files[[cohort_name]]
  bulk_df <- fread(bulk_file, data.table = FALSE)
  
  gene_col <- colnames(bulk_df)[1]
  
  bulk_df[[gene_col]] <- as.character(bulk_df[[gene_col]])
  bulk_df <- bulk_df[!is.na(bulk_df[[gene_col]]) & bulk_df[[gene_col]] != "", ]
  
  expr_df <- bulk_df[, -1, drop = FALSE]
  
  for (j in seq_along(expr_df)) {
    expr_df[[j]] <- as.numeric(expr_df[[j]])
  }
  
  expr_df[[gene_col]] <- bulk_df[[gene_col]]
  
  expr_df <- expr_df %>%
    dplyr::group_by(.data[[gene_col]]) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  gene_name <- expr_df[[gene_col]]
  
  expr_mat <- expr_df[, -1, drop = FALSE]
  expr_mat <- as.matrix(expr_mat)
  rownames(expr_mat) <- gene_name
  
  expr_mat <- expr_mat[rowSums(expr_mat, na.rm = TRUE) > 0, , drop = FALSE]
  
  common_samples = intersect(names(phenotypeall),colnames(expr_mat))
  expr_mat <- expr_mat[, common_samples, drop = FALSE]
  
  phenotype <- phenotypeall[colnames(expr_mat)]
  
  if (any(is.na(phenotype))) {
    stop(
      "Missing phenotype labels for samples: ",
      paste(colnames(expr_mat)[is.na(phenotype)], collapse = ", ")
    )
  }
  
  infos2_hm2mm_list2 <- lapply(seq_along(filelist), function(i) {
    
    file <- filelist[i]
    message("Running file: ", file)
    
    sample_name <- gsub(
      "_sc2st_DestVI_destvi_recluster.h5ad$",
      "",
      basename(file)
    )
    
    st_obj <- h5ad2seurat(file)
    
    st_obj@meta.data$barcode <- paste0(
      rownames(st_obj@meta.data),
      "-",
      sample_name
    )
    
    st_obj <- mm2hm(obj = st_obj)
    
    save_file <- paste0(
      "./Scissor_rdata/",
      cohort_name, "_",
      sample_name,
      ".Scissor_260525_cutoff0.2_alpha0.005.RData"
    )
    
    meta_file <- paste0(
      "./Scissor_rdata/meta_",
      cohort_name, "_",
      sample_name,
      ".Scissor_260520_alpha0.005.csv"
    )
    
    infos2 <- Scissor5(
      bulk_dataset = expr_mat,
      sc_dataset = st_obj,
      phenotype = phenotype,
      tag = c("DR", "NDR"),
      family = "binomial",
      alpha=0.005,
      cutoff = 0.2,
      Save_file = save_file
    )
    
    st_obj <- add_scissorlab(
      st_obj,
      infos2,
      cols = paste0(cohort_name, "_scissor")
    )
    
    write.csv(st_obj@meta.data, meta_file)
    
    return(infos2)
  })
  
  names(infos2_hm2mm_list2) <- gsub(
    "_sc2st_DestVI_destvi_recluster.h5ad$",
    "",
    basename(filelist)
  )
  
  return(infos2_hm2mm_list2)
 
})

names(infos2_hm2mm_list) <- names(bulk_files)

save(
  infos2_hm2mm_list,
  file = "./Scissor_rdata/Scissor_5bulk_all_samples.RData"
)

length(infos2_hm2mm_list)

#---------TEST FOR Scissor labels -----------







#--------scissor_region----------
sample = c('Vehicle','NR1','NR2','NR3','R1','R2','R3')
scanpy_tab20 <- c(
  "#1f77b4", "#aec7e8",
  "#ff7f0e", "#ffbb78",
  "#2ca02c", "#98df8a",
  "#d62728", "#ff9896",
  "#9467bd", "#c5b0d5",
  "#8c564b", "#c49c94",
  "#e377c2", "#f7b6d2",
  "#7f7f7f", "#c7c7c7",
  "#bcbd22", "#dbdb8d",
  "#17becf", "#9edae5"
)

nice_cols <- c(
  "#5DA5DA", "#FAA43A", "#60BD68", "#F17CB0", "#B2912F",
  "#B276B2", "#DECF3F", "#F15854", "#4D4D4D", "#9C755F",
  "#59A14F", "#EDC948", "#AF7AA1", "#FF9DA7", "#76B7B2",
  "#E15759", "#BAB0AC", "#8CD17D", "#B6992D", "#499894",
  "#D37295", "#FABFD2", "#86BCB6", "#F1CE63", "#D4A6C8"
)

scissor_count = lapply(sample,function(n){
  count = data.frame(readxl::read_xlsx('./results.xlsx',n),row.names = 1)
  count = count[,5:16]
  colnames(count)=c('seurat_cluster',"CC2025_scissor_label" ,"GSE153437_scissor_label","GSE153438_scissor_label", "GSE197977_scissor_label",
                    "GSE248835_scissor_label" ,"CC2025_scissor_score","GSE153437_scissor_score","GSE153438_scissor_score","GSE197977_scissor_score",
                    "GSE248835_scissor_score","scissor_score")
  return(count) 
}
)
names(scissor_count)=sample
workpath = '/home/Data/ST_20241110/spatial_input/h5ad/recluster_obj/'
stobjlist = lapply(sample,function(n){
  obj=h5ad2seurat(paste0(workpath,n,'_sc2st_DestVI_destvi_recluster.h5ad'))
  obj = AddMetaData(obj,scissor_count[[n]])
  return(obj)
}
 
  )

dir.create('./figures')
lapply(1:length(stobjlist),function(i){
  p=DimPlot(stobjlist[[i]],reduction = 'Xspatial_',
            # cols=nice_cols,
            cols= scanpy_tab20,
            group.by = 'leiden_0.8',label = T)+coord_fixed()+scale_y_reverse()|DimPlot(stobjlist[[i]],reduction = 'Xspatial_',
                                                                                           group.by = 'scissor_score',cols = c('steelblue','skyblue','lightgrey','pink','red'))+coord_fixed()+scale_y_reverse()
  ggsave(paste0('figures/',sample[i],'_Spatial_scissor_leiden.pdf'),p,width =8,height =6)
  
})

lapply(1:length(stobjlist),function(i){
  p=DimPlot(stobjlist[[i]],reduction = 'Xspatial_',
            # cols=nice_cols,
            cols= scanpy_tab20,
            group.by = 'seurat_cluster',label = T)+coord_fixed()+scale_y_reverse()|DimPlot(stobjlist[[i]],reduction = 'Xspatial_',
                                                                                       group.by = 'scissor_score',cols = c('steelblue','skyblue','lightgrey','pink','red'))+coord_fixed()+scale_y_reverse()
  ggsave(paste0('figures/',sample[i],'_Spatial_scissor_seurat.pdf'),p,width =8,height =6)
  
})

DimPlot(stobjlist[[2]],reduction = 'Xspatial_',
        group.by = 'leiden_0.8',label = T,
        cols=rcolors::get_color(rev(rcolors::rcolors$Set1),n=12))+coord_fixed()+scale_y_reverse()|DimPlot(stobjlist[[2]],reduction = 'Xspatial_',
                                                                                   group.by = 'scissor_score',cols = c('blue','skyblue','lightgrey','pink','red'))+coord_fixed()+scale_y_reverse()

DimPlot(stobjlist[[3]],reduction = 'Xspatial_',
        group.by = 'leiden_0.8',label = T)+coord_fixed()+scale_y_reverse()|DimPlot(stobjlist[[3]],reduction = 'Xspatial_',
                                                                                   group.by = 'scissor_score',cols = c('blue','skyblue','lightgrey','pink','red'))+coord_fixed()+scale_y_reverse()

DimPlot(stobjlist[[4]],reduction = 'Xspatial_',
        group.by = 'leiden_0.8',label = T)+coord_fixed()+scale_y_reverse()|DimPlot(stobjlist[[4]],reduction = 'Xspatial_',
                                                                                   group.by = 'scissor_score',cols = c('blue','skyblue','lightgrey','pink','red'))+coord_fixed()+scale_y_reverse()

DimPlot(stobjlist[[5]],reduction = 'Xspatial_',
        group.by = 'leiden_0.8',label = T)+coord_fixed()+scale_y_reverse()|DimPlot(stobjlist[[5]],reduction = 'Xspatial_',
                                                                                   group.by = 'scissor_score',cols = c('blue','skyblue','lightgrey','pink','red'))+coord_fixed()+scale_y_reverse()
DimPlot(stobjlist[[6]],reduction = 'Xspatial_',
        group.by = 'leiden_0.8',label = T)+coord_fixed()+scale_y_reverse()|DimPlot(stobjlist[[6]],reduction = 'Xspatial_',
                                                                                   group.by = 'scissor_score',cols = c('blue','skyblue','lightgrey','pink','red'))+coord_fixed()+scale_y_reverse()
DimPlot(stobjlist[[7]],reduction = 'Xspatial_',
        group.by = 'leiden_0.8',label = T)+coord_fixed()+scale_y_reverse()|DimPlot(stobjlist[[7]],reduction = 'Xspatial_',
                                                                                   group.by = 'scissor_score',cols = c('blue','skyblue','lightgrey','pink','red'))+coord_fixed()+scale_y_reverse()


lapply(stobjlist,function(obj){
  obj <- SetIdent(obj, value = "leiden_0.8")
  # obj = obj%>%ScaleDatsa()%>%FindVariableFeatures()%>%RunPCA()%>%FindNeighbors()
  FindAllMarkers(obj,
                 features=c('Ccr1','Ccr2','Ccr5','Man1a1','Slc3a2'),only.pos = T)
}
 )


obj = stobjlist[[1]]

FindMarkers(obj,ident.1 = '2',only.pos=T)


FeaturePlot(stobjlist[[1]],
            features = c('Ccr1','Ccr2','Ccr5','Man1a','Slc3a2','Fap','Vim','Dcn','Pecam1','Cdh5','mCherry-CAR'),
            reduction='Xspatial_')&scale_y_reverse()&coord_fixed()&scale_color_viridis_c(limits = c(0, 2))


lapply(1:length(stobjlist),function(i){
  p=FeaturePlot(stobjlist[[i]],
              features = c('Ccr1','Ccr2','Ccr5','Man1a','Spp1','Cd14','Cd68',
                           'Mrc1','Slc3a2','Fap','Vim','Dcn',
                           'Col1a1','Cdh5','mCherry-CAR','Havcr2'),
              reduction='Xspatial_')&scale_y_reverse()&coord_fixed()&scale_color_gradientn(colours = c('lightgrey',
                                                                                                       rcolors::rcolors$MPL_viridis),
                                                                                           limits = c(0, 2))
  ggsave(paste0('figures/',sample[i],'_Spatial_markers_distribution.pdf'),p,width =15,height =12)
  
})








