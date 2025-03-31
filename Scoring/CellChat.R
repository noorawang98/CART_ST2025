#CellChat for spatial transcriptome
library('CellChat')
mcherry = readRDS('../st/L2_cluster/st.l2.mcherry.nr.RDS')
# data.input = Seurat::GetAssayData(visium.brain, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = Idents(visium.brain), row.names = names(Idents(visium.brain))) 
# spatial.locs = Seurat::GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol")) 
data.input =mcherry@assays$SCT@data
spatial.locs = mcherry@reductions$spatial@cell.embeddings

cart = subset(mcherry,subset=((Cd3e>0|Cd3d>0|Cd3g>0)&`mCherry-CAR`>0))

cart =ScaleData(cart)%>%NormalizeData()%>%RunPCA(npcs = 10)%>%FindNeighbors(dims=1:10)%>%RunUMAP(dims=1:10)

table(cart$celltype,cart$exhaustion)

levels(cart$celltype)

st.integrated = readRDS("./")

carttype = ifelse(cart$exhaustion=='High Exhaustion','CART_PR','CART_GR')

st.integrated =AddMetaData(st.integrated,st.integrated@active.ident,col.name = 'ST_region')


st.integrated = SCTransform(st.integrated)

st.integrated = FindNeighbors(st.integrated)
DefaultAssay(st.integrated)<-'integrated'

st.integrated = FindClusters(st.integrated,resolution = 2)

st.markers = FindAllMarkers(st.integrated)

DotPlot(st.integrated,features = marker,cluster.idents = T)

names(carttype)=names(cart@active.ident)

carttype = data.frame(celltype=carttype)
carttype$cells = rownames(carttype)

celltype = data.frame(mcherry@active.ident)
celltype$cells = rownames(celltype)

celltype = left_join(celltype,carttype,"cells")
celltype$celltype =ifelse(is.na(celltype$celltype),
                          as.character(celltype$mcherry.active.ident),
                          celltype$celltype)
meta = data.frame(labels = celltype$celltype,row.names = celltype$cells)

meta$labels=gsub('Red cell','Erythroid',meta$labels)
unique(meta$labels)

scale.factors = list(spot.diameter = 7, spot = 1 # these two information are required
                     # fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)


cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signalin

cellchat@DB <- CellChatDB.use

#preprocession data
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# execution.time = Sys.time() - ptm
# print(as.numeric(execution.time, units = "secs"))
    
ptm = Sys.time()

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.01, 
                              interaction.length = 10,
                              distance.use = TRUE, 
                              # interaction.range = 250,
                              scale.distance = 1)

cellchat <- computeCommunProbPathway(cellchat)


cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
par(mfrow=c(1,1))
pathways=cellchat@netP$pathways
par(mfrow=c(2,2))
plts = lapply(pathways,function(x){
  pathway.show=x
  netVisual_aggregate(cellchat, signaling = pathway.show, layout = "circle",                    
                      edge.width.max = 2, 
                      vertex.size.max = 1, alpha.image = 0.2, 
                      vertex.label.cex = 1
                      )
})

length(pathways)
names(plts) = pathways


plts[c('CCL','CSF','CXCL','IL10','TNF','IL6','IL4','IL1','CSF3','IL2')]

plts$SPP1
pathways.show = 'IL10'
pathways.show =c('CCL','CSF','CXCL','IL10','TNF','IL6','IL4','IL1','CSF3','IL2')
# netAnalysis_computeCentrality()
cellchat=netAnalysis_computeCentrality(cellchat)
pathways.show = "MIF"
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "spatial", 
                    edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, 
                    vertex.label.cex = 1)
                                          
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "circle", 
                    edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, 
                    vertex.label.cex = 1.2)

netVisual_aggregate(cellchat, signaling = pathways.show,layout = "hierarchy", 
                    edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, 
                    vertex.label.cex =1.2)
pathways.show='CCL'
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 10, height = 2.5, font.size = 10)

netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = LR.show, 
                     vertex.receiver = vertex.receiver)


