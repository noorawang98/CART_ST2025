#=======test Create Spatial seurat object===================
setwd('~/LJX/st/')
library(Seurat)

expr_sample_cor = Read10X('./Control_L13_heAuto/')
fig = png::readPNG('./Control_L13_heAuto/he_roi_small.png')
pos_sample_df_cor = read_tsv(gzfile('./Control_L13_heAuto/barcodes_pos.tsv.gz'),col_names = F)
colnames(pos_sample_df_cor)=c('spot','xcoord','ycoord')
colnames(pos_sample_df_cor)=c('spot','ycoord','xcoord')
pos_sample_df_cor = data.frame('spot'=pos_sample_df_cor$spot,
                               'xcoord'=pos_sample_df_cor$ycoord,
                               'ycoord'=pos_sample_df_cor$xcoord)

ctr=myCreateSpatialObject(expr_sample_cor = expr_sample_cor,
                          sample = 'control',
                          # cfig=fig,
                          pos_sample_df_cor = pos_sample_df_cor)

ctr = ctr%>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures()%>%RunPCA()%>%FindNeighbors()

ctr =FindClusters(ctr)

SpatialDimPlot(ctr,pt.size.factor = 3)
SpatialFeaturePlot(ctr,features = c('mCherry-CAR','Spp1','Ccr1'),pt.size.factor = 3)


# pos_sample_df_cor$ycoord =dim(fig)[2]-pos_sample_df_cor$ycoord

# pos_sample_df_cor$xcoord =dim(fig)[1]-pos_sample_df_cor$xcoord

# pos_sample_df_cor= data.frame(pos_sample_df_cor,row.names = 1)


#Prepare istar enhanced data
setwd('/home/Data/')
expr_sample_cor = Read10X('./NR1_L7/L7_heAuto/')
fig = png::readPNG('./Control_L13_heAuto/he_roi_small.png')
pos_sample_df_cor = read_tsv(gzfile('./NR1_L7/L7_heAuto/barcodes_pos.tsv.gz'),col_names = F)
colnames(pos_sample_df_cor)=c('spot','xcoord','ycoord')

pos_sample_df_cor = data.frame('spot'=pos_sample_df_cor$spot,
                               'xcoord'=pos_sample_df_cor$ycoord,
                               'ycoord'=pos_sample_df_cor$xcoord)


features.list = list('T Cell'=c('Cd3e'),
             'CD8 T cell'=c('Cd8a','Cd8b1'),
             'CD4 T cell'=c('Cd4'),
             'B Cell'=c('Cd19','Cd79a','Cd79b','Ms4a1'),
             'Plasma'=c('Sdc1','Ms4a1'),
             'NK'=c('Klrk1','Klrb1c'),
             'Myeloid'=c('Lyz2'),
             'Fibroblast'=c('Fap','Col1a1','Col1a2'),
             'M1'=c('Cd14','Csf1r','Adgre1','Il1b','Tnf','Ifng'),
             'M2'=c('Mrc1','Cd163','Vegfa','Arg1'),
             'Neu'=c('S100a8','S100a9','Il1b','Csf3r'),
             'cDC'=c('Itgax','Itgae','Cd74','Xcr1','Clec10a'),
             'pDC'=c('Siglech','Cd74'),
             'Tumor'=c('EGFP'),
             'EpC'=c('Epcam','Krt18'),
             'VEC'=c('Pecam1','Vwf'),
             'LEC'=c('Pecam1','Lyve1','Prox1','Thy1'),
             'CART'=c('mCherry-CAR'),
             'Cytotoxic'=c('Gzma','Gzmb','Gzmk','Prf1'),
             'Memory'=c('Tcf7','Sell','Il7r'),
             'Proliferation'=c('Mki67','Il2ra','Cd69'),
             'Exhaustion'=c('Pdcd1','Tigit','Havcr2','Tox'),
             'Treg'=c('Foxp3','Il2ra','Ctla4'),
             'CCR1 MEB'=c('Ccr1','Spp1','Col1a1','Lyz2','Mrc1')
             )
marker_to_celltype = reshape2::melt(features.list)

colnames(marker_to_celltype)=c('gene','label')

CCR1_MEB=data.frame('gene'=c('Lyz2','Mrc1','Ccr1','Spp1','Col1a1','Fap'))

readr::write_tsv(CCR1_MEB,'./data/markers/signature-score-template.txt',col_names = F)


obj = myCreateSpatialObject(expr_sample_cor = expr_sample_cor,
                                sample = 'NR1',
                                # cfig=fig,
                                pos_sample_df_cor = pos_sample_df_cor)

obj = obj%>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures()%>%RunPCA()%>%FindNeighbors()

obj =FindClusters(obj)
cellcolors= get_color(rcolors::rcolors$GMT_paired,n = 12)
names(cellcolors)=levels(obj@active.ident)

SpatialDimPlot(obj,cols =cellcolors,label.size = 1.2)
SpatialFeaturePlot(obj,features = c('Ccr1','Spp1','Col1a1','Lyz2','Mrc1','Fap'))
SpatialFeaturePlot(obj,features = c('Pecam1','Vwf','Lyve1',
                                    'EGFP','mCherry-CAR','Cd79b'))

cnts = GetAssayData(obj)

cnts = cnts[rowSums(cnts)>100,]
genenames =rownames(cnts)
# cnts=readr::read_tsv('../istar-master/data/demo/cnts.tsv')

intersect(marker_to_celltype$gene,rownames(cnts))
cnts  =data.frame(t(cnts))
colnames(cnts) = genenames
cnts$spot= rownames(cnts)

cnts =cnts[,c(ncol(cnts),1:(ncol(cnts)-1))]
readr::write_tsv(cnts,'./data/cnts.tsv')

pos = pos_sample_df_cor

colnames(pos)=c('spot','x','y')

readr::write_tsv(pos,'data/locus-raw.tsv')

# Define Tumor border
library(SPIAT)

# createSpatialObject = 
create_spatial_experiment <- function(count_matrix, spatial_coords, image_data = NULL) {
  # Check if required packages are installed and loaded
  required_packages <- c("SpatialExperiment", "SummarizedExperiment")
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste(package, "package is not installed. Please install it before using this function."))
    }
    library(package, character.only = TRUE)
  }
  
  # Create a SummarizedExperiment object from the count matrix
  se <- SummarizedExperiment(assays = list(counts = count_matrix))
  
  # Create the SpatialExperiment object
  spe <- SpatialExperiment(
    assay = se,
    spatialCoords = spatial_coords
  )
  
  # Add image data if provided
  if (!is.null(image_data)) {
    imgData(spe) <- image_data
  }
  
  return(spe)   
}

