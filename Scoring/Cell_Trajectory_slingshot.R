#Slingshot
# library(slingshot)
.libPaths()
# BiocManager::install('slingshot')
library('slingshot')
library(uwot)
#--function----

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}


#-----create object---------
counts = cd8t@assays$SCT@data
counts = tc@assays$RNA@data
sce <- SingleCellExperiment(assays = List(counts = counts))

# rd <- slingshotExample$rd
# cl <- slingshotExample$cl

#----upstream analysis-----------

geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
length(geneFilter)
sce <- sce[geneFilter, ]
#normalise
assays(sce)$norm <- FQnorm(assays(sce)$counts)
#dimension reduction
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)


#RunUMAP
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
reducedDims(sce)<- SimpleList(PCA = rd1, UMAP = rd2,UMAP2=tc@reductions$umap@cell.embeddings)

# reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

#cluster cells
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1
colData(sce)$cluster = tc$seurat_clusters


library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP')

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

par(mfrow=c(1,2))
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')


plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[sce$cluster], pch=5, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')


#--------downstream analysis--------
library(tradeSeq)

# fit negative binomial GAM
sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)

# topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:2000]
cd8markers

cd8tmarkers = FindAllMarkers(cd8t,only.pos =T)

cd8tmarkers[cd8tmarkers$pct.1>2*cd8tmarkers$pct.2&cd8tmarkers$pct.1>0.4,]

topgenes=cd8tmarkers[cd8tmarkers$pct.1>2.5*cd8tmarkers$pct.2&cd8tmarkers$pct.1>0.3,]$gene

topgenes = grep('^Gm|^Rps|Rpl|AW',topgenes,invert = T,value = T)


DoHeatmap(cd8t,topgenes)+scale_fill_gradientn(colours = rev(rcolors::rcolors$RdBu),limits=c(-2.5,2.5))|DimPlot(cd8t,label = T,repel=T)

pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
# cell.ord = colnames(assays(sce))
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$cluster[pst.ord]


heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

