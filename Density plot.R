#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
#DENSITY Plot Fig 1E
getwd()

seuraObj <- seurnodblct #change according the name of the seurat object subject of analysis

DefaultAssay(seuraObj) <- 'RNA'
Idents(seuraObj) <- "annot"

subset_H <- subset(seuraObj, subset = sample == "EPJ3")
DimPlot(subset_H)
subset_LLC <- subset(seuraObj, subset = sample == "EPJ4")
DimPlot(subset_LLC)
dim(seuraObj)
dim(subset_LLC)
dim(subset_H)
createDensityPlot <- function(sobj){
  
  tmp.all<-as.data.frame(Embeddings(object = sobj, reduction = "umap"))
  p <- ggplot(tmp.all, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(colour="#00000000") +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon", bins=50) +
    scale_fill_gradientn(colors = c("#4169E100","royalblue", "darkolivegreen3","goldenrod1","red")) +
    theme_classic() + 
    theme(legend.position="right") +
    coord_fixed(ratio=1) +
    xlim(-30,30) + 
    ylim(-20,20) +
    NoAxes() 
  
  return(p)
}
createDensityPlot(subset_LLC)
dev.off()

createDensityPlot(subset_H)
dev.off()