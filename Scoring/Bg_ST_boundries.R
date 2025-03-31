#-----Function and library--------
library(Cottrazm)
library(Seurat)
library(infercnv)
library(dplyr)

print('STPreProcess')

mySTCNV=function (TumorST = TumorST, assay = c("Morph", "Spatial"), OutDir = NULL, gene_order_file=NULL,
          Sample = Sample, num_threads = 30) 
{
  if (is.null(OutDir) == TRUE) {
    OutDir <- paste(getwd(), "/", Sample, "/", sep = "")
    dir.create(OutDir)
  }
  matrix <- Seurat::GetAssayData(TumorST, slot = "counts", 
                                 assay = assay) %>% as.matrix()
  annotation_file <- paste(OutDir, "InferCNV/CellAnnotation.txt", 
                           sep = "")
  NormalCluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(split(TumorST@meta.data[, 
                                                                                               c("seurat_clusters", "NormalScore")], TumorST@meta.data[, 
                                                                                                                                                       c("seurat_clusters", "NormalScore")]$seurat_clusters), 
                                                                       function(test) mean(test$NormalScore))), decreasing = T)[1]]
  ref_cluster <- NormalCluster
  if(is.null(gene_order_file)==TRUE){
    gene_order_file <- system.file("extdata/gencode_v38_gene_pos.txt", 
                                   package = "Cottrazm")
  }
  infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = matrix, 
                                                 annotations_file = annotation_file, delim = "\t", gene_order_file = gene_order_file, 
                                                 ref_group_names = ref_cluster)
  infercnv_obj <- infercnv::run(infercnv_obj, cutoff = 0.1, 
                                out_dir = paste(OutDir, "InferCNV/output_", assay, sep = ""), 
                                cluster_by_groups = F, analysis_mode = "subclusters", 
                                write_phylo=T,write_expr_matrix = T,
                                denoise = T, HMM = T, tumor_subcluster_partition_method = "random_trees", 
                                HMM_type = "i6", BayesMaxPNormal = 0, num_threads = num_threads)
  return(infercnv_obj)
}

# mySTCNVScore=function (TumorST = TumorST, assay = c("Saptial", "Morph"), OutDir = NULL, 
#           Sample = Sample) 
# {
#   if (is.null(OutDir) == TRUE) {
#     OutDir <- paste(getwd(), "/", Sample, "/", sep = "")
#     dir.create(OutDir)
#   }
#   cnv_outdir = paste(OutDir, "InferCNV/output_", assay, sep = "")
#   cell_groupings <- read.tree(file = paste(cnv_outdir, "/infercnv.observations_dendrogram.txt", 
#                                            sep = ""))
#   infercnv.label <- dendextend::cutree(cell_groupings, k = 8)
#   infercnv.label <- as.data.frame(infercnv.label)
#   infercnv.label <- rbind(infercnv.label, data.frame(row.names = rownames(TumorST@meta.data)[!rownames(TumorST@meta.data) %in% 
#                                                                                                infercnv.label$row.names], infercnv.label = rep("Normal", 
#                                                                                                                                                length(rownames(TumorST@meta.data)[!rownames(TumorST@meta.data) %in% 
#                                                                                                                                                                                     infercnv.label$row.names]))))
#   TumorST@meta.data$CNVLabel <- infercnv.label$infercnv.label[match(rownames(TumorST@meta.data), 
#                                                                     rownames(infercnv.label))]
#   .cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", 
#                      "#882E72", "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", 
#                      "#E78AC3", "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", 
#                      "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", "#666666", 
#                      "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", 
#                      "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", 
#                      "#ffff00")
#   pdf(paste(OutDir, Sample, "_cnv_label.pdf", sep = ""), width = 7, 
#       height = 7)
#   p <- SpatialDimPlot(TumorST, group.by = "CNVLabel", cols = .cluster_cols, 
#                       pt.size.factor = 1, alpha = 0.6) + scale_fill_manual(values = .cluster_cols)
#   print(p)
#   dev.off()
#   #pdf(paste(OutDir, Sample, "_reduction_cnvlabel.pdf", sep = ""), 
#   #    width = 7, height = 7)
#   #p <- DimPlot(TumorST, group.by = "CNVLabel", cols = .cluster_cols) + 
#   #    scale_fill_manual(values = .cluster_cols)
#   #print(p)
#   #dev.off()
#   cnv_table <- read.table(paste(cnv_outdir, "/infercnv.observations.txt", 
#                                 sep = ""), header = T)
#   cnv_score_table <- as.matrix(cnv_table)
#   cnv_score_tableA <- abs(cnv_score_table - 3)
#   cell_scores_CNV <- as.data.frame(colSums(cnv_score_tableA)) %>% 
#     set_colnames(., c("cnv_score"))
#   rownames(cell_scores_CNV) <- gsub("\\.", "-", rownames(cell_scores_CNV))
#   TumorST@meta.data$cnv_score <- cell_scores_CNV$cnv_score[match(rownames(TumorST@meta.data), 
#                                                                  rownames(cell_scores_CNV))]
#   TumorST@meta.data$cnv_score <- ifelse(TumorST@meta.data$CNVLabel == 
#                                           "Normal", 0, TumorST@meta.data$cnv_score)
#   cell_scores_CNVA <- TumorST@meta.data[, c("CNVLabel", "cnv_score")]
#   pdf(paste(OutDir, Sample, "_cnv_observation_vlnplot.pdf", 
#             sep = ""), width = 6, height = 4)
#   p <- ggplot(cell_scores_CNVA, aes(x = CNVLabel, y = cnv_score, 
#                                     fill = CNVLabel)) + geom_violin(alpha = 0.5) + geom_boxplot(stat = "boxplot", 
#                                                                                                 alpha = 1, width = 0.5, outlier.size = 0.5) + labs(y = "CNV_scores") + 
#     ggpubr::stat_compare_means() + scale_fill_manual(values = .cluster_cols) + 
#     theme(axis.text = element_text(colour = "black"), panel.background = element_blank(), 
#           panel.grid = element_blank(), legend.title = element_text(face = "bold"), 
#           axis.line = element_line(colour = "black"), legend.text = element_text(size = 10), 
#           title = element_text(face = "bold")) + labs(title = "CNV Scores") + 
#     NoLegend()
#   print(p)
#   dev.off()
#   return(TumorST)
# }

#-----Run-----
setwd('/home/Data/')
TumorST = readRDS('ST_boundries_inDir/TumorST.rds')
NormalFeatures <-c("Ccr1","Lyz2","Spp1","Col1a2","Mrc1",'H2-Aa',"H2-K1", 
                   'Gzmf',"Gzmg","Tcf4","Csf1r",'Actb',
                   "Arg1" ,"C1qc","C1qa","C1qb",'Cxcl2',
                   'S100a8','Ighm','Clec10a',"Cd74")
TumorST@meta.data$NormalScore <- apply(TumorST@assays$Spatial$data[rownames(TumorST@assays$Spatial$data) %in% 
                                                                     NormalFeatures, ], 2, mean)
nr1.infercnv=mySTCNV(TumorST = TumorST,
                     gene_order_file = 'ST_boundries_inDir/mm_gene_orderfile.txt',
                     assay = 'Spatial',Sample = 'NR1',
                     OutDir = 'ST_boundries_outDir/')
Sample='NR1'
OutDir='ST_boundries_outDir/'

TumorST <-
  STCNVScore(
    TumorST = TumorST,
    assay = "Spatial",
    Sample = Sample,
    OutDir = OutDir
  )

saveRDS(TumorST,'ST_boundries_outDir/TumorST_step2.rds')

TumorSTn=BoundaryDefine(TumorST,
                        # MalLabel =,
                        OutDir = OutDir,
                        Sample=Sample)



