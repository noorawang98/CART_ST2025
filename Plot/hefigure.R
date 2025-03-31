#Figure HE FIGURES with annotated 
he.n9 ='~/LJX/st/n9/05.AllheStat/allhe/he_roi_small.png'
he.r5 = '~/LJX/st/r5/05.AllheStat/allhe/he_roi_small.png'

df =DimPlot(st.integrated,cols = stcolors,reduction = 'spatial')
df = df$data
colnames(df)=c('x','y','group')
df$cluster=st.integrated$seurat_clusters
df$cluster = factor(df$cluster,levels=0:12)
df.split = split(df,df$group)


he_fig =png::readPNG(he.n9)

heplot(data = df.split$NR,color = stcolors,he_fig = he_fig)+NoLegend()
pdf('HE_plot_100um_seuratclusters_NR.pdf',width = 10.29/2,height = 10.24/2)
heplot(data = df.split$NR,color = stcolors,he_fig = he_fig)+NoLegend()
dev.off()


he_fig =png::readPNG(he.r5)

heplot(data = df.split$R,color = stcolors,he_fig = he_fig)+NoLegend()
pdf('HE_plot_100um_seuratclusters_R.pdf',width = 10.29/2,height = 10.24/2)
heplot(data = df.split$R,color = stcolors,he_fig = he_fig)+NoLegend()
dev.off()

nmcolor = rep(NA,length(stcolors))
nmcolor[c(3,5)]<-stcolors[c(3,5)]

he_fig = png::readPNG(he.n9)
pdf('HE_plot_100um_nm2_4_NR.pdf',width = 10.29/2,height = 10.24/2)
heplot(data = df.split$NR,color = nmcolor,he_fig = he_fig)+NoLegend()
dev.off()

he_fig = png::readPNG(he.r5)
pdf('HE_plot_100um_nm2_4_R.pdf',width = 10.29/2,height = 10.24/2)
heplot(data = df.split$R,color = nmcolor,he_fig = he_fig)+NoLegend()
dev.off()


# df=FeaturePlot(st.integrated,features = c('Pseudotime',
#                                           'Ccr2','Ccr1','Spp1','mCherry-CAR'
#                                           ),
#                reduction = 'spatial',split.by = 'sample')
# head(df$data)
# data =df$data
# data$sample =st.integrated$sample

intersect(c('Ccr2','Ccr1','Spp1','mCherry-CAR')
          ,rownames(st.integrated@assays$RNA@data))

mat = st.integrated@assays$RNA@data[c('Ccr2','Ccr1','Spp1','mCherry-CAR'),]
mat = t(mat)
mat = as.data.frame(mat)


df=st.integrated@meta.data[,c('x','y','sample',"Pseudotime")]

df=cbind(df,mat)

he_fig = png::readPNG(he.n9)

df.split = split(df,df$sample)
 
mat.split = split(mat,df$sample)

for (i in 1:ncol(mat)){
  
  fignames = paste0('FeatureHeplot_',colnames(mat)[i],'_NR.pdf')
  df.in = df.split$NR
  df.in$gene = mat.split$NR[[i]]
  
  he_fig=png::readPNG(he.n9)
  p=featureheplot(data = df.in,he_fig = he_fig)+NoLegend()
  ggsave(fignames,p,width = 10.29/2,height = 10.24/2)
  
  fignames = paste0('FeatureHeplot_',colnames(mat)[i],'_R.pdf')
  df.in = df.split$R
  df.in$gene =mat.split$R[[i]]
  he_fig=png::readPNG(he.r5)
  p=featureheplot(data = df.in,he_fig = he_fig)+NoLegend()
  ggsave(fignames,p,width = 10.29/2,height = 10.24/2)
}
# featureheplot()


featureheplot = function(data,he_fig,pt.size=1,color='viridis'){
    p<-ggplot(data=data,aes(x ,y)) + 
      ggpubr::background_image(he_fig)+
      # geom_point(aes(colour = cluster),size = opt$point_size, shape=16) + 
      geom_point(aes(fill = gene,color=gene),size = pt.size, shape=21) + 
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
      guides(colour = guide_legend(override.aes = list(size=3.5)))
    # theme(legend.position = "none")
    
    # if(is.null(celltypecol)){
    #   clstr_plot = clstr_plot
    # }else{
    #   clstr_plot = clstr_plot+scale_fill_manual(values= celltypecol)+
    #     scale_color_manual(values= celltypecol)
    # }
    # 
    if(is.null(color)){
      p=p
    }else{
      p+scale_fill_gradientn(colours = c('lightgrey',hcl.colors(palette = color,n = 100)))+
        scale_color_gradientn(colours = c('lightgrey',hcl.colors(palette = color,n = 100)))
    }
    # return(clstr_plot)
  }

heplot = function(data,he_fig,pt.size=1,color=NULL){
  p<-ggplot(data=data,aes(x ,y)) + 
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
    guides(colour = guide_legend(override.aes = list(size=3.5)))
  # theme(legend.position = "none")
  
  # if(is.null(celltypecol)){
  #   clstr_plot = clstr_plot
  # }else{
  #   clstr_plot = clstr_plot+scale_fill_manual(values= celltypecol)+
  #     scale_color_manual(values= celltypecol)
  # }
  # 
  if(is.null(color)){
    p=p
  }else{
    p+scale_fill_manual(values = color)+scale_color_manual(values = color)
  }
  # return(clstr_plot)
}


