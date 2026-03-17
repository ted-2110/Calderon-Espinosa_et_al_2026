
#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
#Supervised annotation#

res<-0.6
dims.use<-30
DEgenes
seurnodbl

features.use=unlist(lapply(DEgenes, function(x) { head(x[x$avg_log2FC>0,]$gene,3)})) #select the top 3 DEgenes from each cluster and unlist this to one vector
length(features.use) # or 3x21
names(features.use)=NULL
features.use=features.use[!duplicated(features.use)] # make sure there are no overlapping genes between clusters
DotPlot(seurnodbl, features = features.use)+RotatedAxis()
dev.off()
# Log transform the data 
seurnodbl$nCount_RNA_log=log2(seurnodbl$nCount_RNA)
seurnodbl$nFeature_RNA_log=log2(seurnodbl$nFeature_RNA)

features.use <- c((features.use), "perc.mito", "doublet.score", "nCount_RNA_log","nFeature_RNA_log")

DotPlot(seurnodbl, features = features.use,group.by = "RNA_snn_PC30_res.0.6")+RotatedAxis()+theme(axis.title = element_blank(),legend.position = "bottom")
dev.off()

#DoHeatmap is a Seurat function to make heatmaps. 
#first scale all genes 
seurnodbl$annot <- Idents(seurnodbl)
seurnodbl_scaled <- ScaleData(seurnodbl, features = rownames(seurnodbl))
DoHeatmap(seurnodbl_scaled, features = features.use, group.by="annot",  size = 3.5,  angle = 75)+theme(legend.position = "bottom")
dev.off() 



#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# If needed to rename; first rename clusters from 1:21;
Idents(seurnodbl)=  paste0("RNA_snn_PC",dims.use,"_res.",res)
Idents(seurnodbl)=  factor(Idents(seurnodbl),levels = 1:(length(unique(Idents(seurnodbl)))))
#Reorder the levels
levels(x = seurnodbl) <-c(21,11,2,4,20,13,19,1,10,8,12,3,6,9,7,14,5,18,17,15,16)
Idents(seurnodbl) <- plyr::mapvalues(x = Idents(seurnodbl), from = c(21,11,2,4,20,13,19,1,10,8,12,3,6,9,7,14,5,18,17,15,16),
                                     to = 1:21)
Idents(seurnodbl) <- plyr::mapvalues(x = Idents(seurnodbl), from = 1:21,
                                     to=c("Early pro-B cell","Late pro-B cell","Immature B cell","Immature B cell_2","B cell",
                                          "NK/T cell","Basophil/eosinophil prog", "Neutrophil_1", "Neutrophil_2", "Neutrophil_3","Neutrophil_4", "GP",
                                          "GMP","HPC","Promonocyte","DC","Monocyte","Macrophage_1","Macrophage_2",
                                          "MEP","Proerythroblast"))
DimPlot(seurnodbl, cols=cols.use, label= TRUE, label.size = 3)
dev.off()
DimPlot(seurnodbl, cols=cols.use)
dev.off()


seurnodbl$annot <- Idents(seurnodbl)

seurnodbl_scaled <- ScaleData(seurnodbl, features = rownames(seurnodbl))
DoHeatmap(seurnodbl_scaled, features = features.use, group.by="annot",  size = 3.5,  angle = 75)+theme(legend.position = "bottom")
dev.off() 

DimPlot(seurnodbl, group.by ="sample", cols =c('cadetblue3',"lightpink1"),shuffle=TRUE)+ggtitle("Dimension reduction plot - grouped by sample")
dev.off()

DimPlot(seurnodbl, cols = cols.use) +ggtitle("Dimension reduction plot - Annotation") #Extended Data Fig 1B
dev.off()


# DimPlot(seur_clean,repel =T,label=T,cols=cols.use) 
DimPlot(seurnodbl,repel =T,label=T, group.by = "CMO_call", split.by = "sample")
DimPlot(seurnodbl,repel =T,label=T, group.by = "Assignment", split.by = "sample")

#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#

#Annotate Clusters according the cell type they belong: HP-related, Neu-related, MoDC-related, Ery-related, BaEo-related, B cell-related, NK/T cell-related

Idents(seurnodbl)=  paste0("RNA_snn_PC",dims.use,"_res.",res)
Idents(seurnodbl)=  factor(Idents(seurnodbl),levels = 1:(length(unique(Idents(seurnodbl)))))
#Reorder the levels
levels(x = seurnodbl) <-c(21,11,2,4,20,13,19,1,10,8,12,3,6,9,7,14,5,18,17,15,16)
Idents(seurnodbl) <- plyr::mapvalues(x = Idents(seurnodbl), from = c(21,11,2,4,20,13,19,1,10,8,12,3,6,9,7,14,5,18,17,15,16),
                                     to = 1:21)
Idents(seurnodbl) <- plyr::mapvalues(x = Idents(seurnodbl), from = 1:21,
                                     to=c("B cell","B cell","B cell","B cell","B cell",
                                          "NK/T cell","Basophil/eosinophil prog", "Neutrophil", "Neutrophil", "Neutrophil","Neutrophil", "Neutrophil",
                                          "Neutrophil","Progenitor","Monocyte/DC","Monocyte/DC","Monocyte/DC","Monocyte/DC","Monocyte/DC",
                                          "Erythroid","Erythroid"))


#create a new Seurat Object containing the annotation from the cell type
seurnodblct <- seurnodbl
seurnodblct$annot <- Idents(seurnodblct)
DimPlot(seurnodblct, cols=cols.use, label= TRUE, label.size = 3) #Fig 1A

#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
#### Feature plot caracterising each cell type# Fig 1C
#HP-related
FeaturePlot(seurnodblct,"Ms4a3", repel = TRUE, cols = c("grey80","darkred"))+ NoLegend()+NoAxes()
dev.off()

#Neu-related 
FeaturePlot(seurnodblct,features=c("Cd177"),repel = TRUE, cols = c("grey80","darkred"))+ NoLegend()+NoAxes()
dev.off()

#MoDC-related
FeaturePlot(seurnodblct,"Ms4a6c", repel = TRUE, cols = c("grey80","darkred"))+ NoLegend()+NoAxes()
dev.off()

#Ery-related
FeaturePlot(seurnodblct,"Car2", repel = TRUE, cols = c("grey80","darkred"))+ NoLegend()+NoAxes()
dev.off()

#BaEo-related
FeaturePlot(seurnodblct,"Ms4a2", repel = TRUE, cols = c("grey80","darkred"))+ NoLegend()+NoAxes()
dev.off()

#B cell-related
FeaturePlot(seurnodblct,"Cd79a", repel = TRUE, cols = c("grey80","darkred"))+ NoLegend()+NoAxes()
dev.off()

#NK/T cell-related
FeaturePlot(seurnodblct,"Cd3e", repel = TRUE, cols = c("grey80","darkred"))+ NoLegend()+NoAxes()
dev.off()

#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
#Heatmap DEG - Fig 1D

DEgenes=list()
for (i in levels(Idents(seurnodblct))){
  DEgenes[[i]]<-FindMarkers(seurnodblct, ident.1 = i,min.cells.group=2,pseudocount.use = 0.01, max.cells.per.ident = 1000,verbose = F) 
  #The max.cells.per.ident is set lower to speed up this calculation, but should be Inf to do it 'completely'
  #logfc.threshold=0.2 and min.diff.pct=0.2
  #If you want only marker genes that are positive (Log2FC>0 or upregulated); only.pos=TRUE
  DEgenes[[i]]$cluster=i
  DEgenes[[i]]$pseudocount=0.01
  DEgenes[[i]]$max.cells.per.ident=1000
  DEgenes[[i]]$gene=row.names(DEgenes[[i]])
}


features.use=unlist(lapply(DEgenes, function(x) { head(x[x$avg_log2FC>0,]$gene,3)})) #select the top 3 DEgenes from each cluster and unlist this to one vector
length(features.use) # or 3x21
names(features.use)=NULL
features.use=features.use[!duplicated(features.use)] # make sure there are no overlapping genes between clusters
DotPlot(seurnodblct, features = features.use)+RotatedAxis()
dev.off()
# Log transform the data 
seurnodblct$nCount_RNA_log=log2(seurnodblct$nCount_RNA)
seurnodblct$nFeature_RNA_log=log2(seurnodblct$nFeature_RNA)

features.use <- c((features.use), "perc.mito", "doublet.score", "nCount_RNA_log","nFeature_RNA_log")

#DoHeatmap is a Seurat function to make heatmaps. 
#first scale all genes 
seurnodblct$annot <- Idents(seurnodblct)
seurnodblct_scaled <- ScaleData(seurnodblct, features = rownames(seurnodblct))

#Do Heatmap annotated clusters 
DoHeatmap(seurnodblct_scaled, features = features.use, group.by="annot",  size = 3.5,  angle = 75)+theme(legend.position = "bottom")
dev.off() 


#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
#Volcano Plot - Fig 4A
DefaultAssay(seursub) <- "RNA"
# DE gene expression analysis on the whole subset group 
DE <- FindMarkers(seurnodblct, group.by = "sample", ident.1 = "EPJ4",logfc.threshold=0) 
#pseudocount.use = 0.1 )

DE[["avg_log2FC"]] <- as.numeric(DE[["avg_log2FC"]])
keyvals <- ifelse(
  DE$avg_log2FC < -0.8, 'dodgerblue3',
  ifelse(DE$avg_log2FC > 0.8, 'firebrick3',
         'slategray'))
keyvals[is.na(keyvals)] <- 'slategray'
names(keyvals)[keyvals == 'firebrick3'] <- 'upregulated'
names(keyvals)[keyvals == 'slategray'] <- 'ns'
names(keyvals)[keyvals == 'dodgerblue3'] <- 'downregulated'

library(EnhancedVolcano)
EnhancedVolcano(DE,
                lab = row.names(DE),
                x = "avg_log2FC",
                y = "p_val_adj",
                # xlim = c(min(resPlot$log2FoldChange), max(resPlot$log2FoldChange)),
                # ylim = c(min(-log(resPlot$log2FoldChange)), max(-log(resPlot$log2FoldChange))),
                title = 'Orthotopic LLC-bearing vs Healthy C57BL/6',
                subtitle = paste("Differential gene expression -","Erythroid"),
                caption = paste("fold change cutoff, 0.8; p-value cutoff, 0.1"),
                pCutoff = 0.1,
                FCcutoff = 0.8,
                pointSize = 3.0,
                labSize = 6,
                colCustom = keyvals,
                max.overlaps=Inf,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                arrowheads = FALSE)
dev.off()

DimPlot(seurnodblct, group.by ="sample", cols =c('cadetblue3',"lightpink1"),shuffle=TRUE)+ggtitle("Dimension reduction plot - grouped by sample")
dev.off()


#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# CACOA - expressionshiftmagnitudes ####
# https://github.com/kharchenkolab/cacoa
# http://pklab.med.harvard.edu/viktor/cacoa/walkthrough_short.html

BiocManager::install(c("clusterProfiler", "DESeq2", "DOSE", "EnhancedVolcano", "enrichplot", "fabia", "GOfuncR", "Rgraphviz"))
install.packages("remotes")
remotes::install_github("kharchenkolab/cacoa")
library(devtools)
library(cacoa)
sample.per.cell <- as.vector(seurnodblct$Assignment)
names(sample.per.cell) <- rownames(seurnodblct@meta.data)


cell.groups <- as.vector(seurnodblct@meta.data[["annot"]])
names(cell.groups) <- rownames(seurnodblct@meta.data)

sample.groups <- as.factor(c("Healthy","Healthy","Healthy","Healthy","LLC","LLC","LLC","LLC"))
names(sample.groups) <- c("Healthy_1", "Healthy_2", "Healthy_3", "Healthy_4", "LLC_5","LLC_6", "LLC_7", "LLC_14")
ref.level <- "Healthy"
target.level <- "LLC"

cao <- Cacoa$new(seurnodblct, sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=sample.per.cell, 
                 ref.level=ref.level, target.level=target.level,embedding = seurnodblct@reductions$umap@cell.embeddings, n.cores=1,graph.name = "RNA_snn_PC30")

cao$plot.params <- list(size=0.1, alpha=0.1, font.size=c(2, 3))
cao$plot.theme <- cao$plot.theme + theme(legend.background=element_blank())
#install.packages(c('coda.base', 'psych'))
cao$estimateCellLoadings()
cao$estimateExpressionShiftMagnitudes(min.cells.per.sample = 1)
cao$estimateExpressionShiftMagnitudes()

#Cell abundance - CACOA#
cao$plotCellLoadings(show.pvals=FALSE, ordering = "loading",jitter.size=0)+
  geom_boxplot(color= 'gray20', fill='gray74')+
  theme(axis.text=element_text(size=11))+
  scale_x_continuous(limits=c(-0.8, 0.8))+
  scale_x_continuous(guide = guide_axis(angle = 45))+
  scale_y_discrete(guide = guide_axis(angle = 0), 
                   limits = c("Erythroid","Monocyte/DC","Progenitor","Neutrophil","Basophil/eosinophil prog","NK/T cell","B cell"))
dev.off()

#Expression shift magnitudes- CACOA#
cao$plotExpressionShiftMagnitudes()
dev.off()

cao$plotExpressionShiftMagnitudes(type='box',notch=FALSE,show.pvalues=c("adjusted")) +
  geom_boxplot(color="gray20", fill='gray74')+
  coord_flip() +scale_y_continuous(limits=c(-0.08, 0.08))+ 
  scale_y_continuous(guide = guide_axis(angle = 45))+
  scale_x_discrete(guide = guide_axis(angle = 0), 
                   limits = c("Erythroid","Monocyte/DC","Progenitor","Neutrophil","Basophil/eosinophil prog","NK/T cell","B cell"))

dev.off()



#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#### Abundance of cells - Hashing - Fig 1F ####
seurnodblct$Assignment <- factor(seurnodblct$Assignment, levels=c("Healthy_1","Healthy_2","Healthy_3","Healthy_4","LLC_5","LLC_6", "LLC_7","LLC_14"))

#proportions of cells per sample
library(dplyr)
meta.data<-seurnodblct@meta.data %>%
  dplyr::group_by(annot,sample,Assignment) %>%
  dplyr::summarise(count=n())%>% suppressMessages()%>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(perc.per.sample = (count / sum(count))*100)    %>% 
  dplyr::group_by(annot) %>% 
  dplyr::mutate(perc.per.clust = (count / sum(count))*100)
meta.data$annot=factor(meta.data$annot, levels=levels(seurnodblct$annot))
meta.data$Assignment=factor(meta.data$Assignment, levels=levels(seurnodblct$Assignment))


# Redo the analysis of above, now leaving out the Multiplets, Blanks and Unassigned
Idents(seurnodblct) <- "Assignment"
object = subset(seurnodblct, idents = c("Healthy_1","Healthy_2","Healthy_3","Healthy_4","LLC_5","LLC_6","LLC_7","LLC_14"))

meta.data2<-object@meta.data %>%
  dplyr::group_by(annot,sample,Assignment) %>%
  dplyr::summarise(count=n())%>% suppressMessages()%>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(perc.per.sample = (count / sum(count))*100)    %>% 
  dplyr::group_by(annot) %>% 
  dplyr::mutate(perc.per.clust = (count / sum(count))*100)
meta.data2$annot=factor(meta.data2$annot, levels=levels(object$annot))
meta.data2$Assignment=factor(meta.data2$Assignment, levels=levels(object$Assignment))

meta.data2$sample <- factor(meta.data2$sample)
levels(meta.data2$sample) <- list(Healthy = "EPJ3", LLC = "EPJ4")


All<-ggplot(meta.data2, aes(x=perc.per.clust, y= forcats::fct_rev(annot),fill=forcats::fct_rev(Assignment)))+geom_bar(stat = "identity",width=0.75)+ 
  guides(fill = guide_legend(reverse = TRUE))+
  theme_bw()+scale_fill_manual(values=c("lightpink","lightpink1","lightpink2","lightpink3",'cadetblue1','cadetblue2','cadetblue3','cadetblue'))+
  theme(axis.title.x=element_blank(),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title =element_text(hjust = 0.5),
  )+ ylab(paste0("% cells per cluster"))+RotatedAxis()+ggtitle ("Proportion of cells per sample")


All
dev.off()

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
####Volcano plots - Heatmaps up- and down-regulated genes per cell type - Fig 4D

# Differential Gene expression nodb ####
clustersDEG <- levels(Idents(seurnodblct))  ### Specify here which clusters to use
# clustersDEG <- clustersDEG[clustersDEG!="Fibroblasts"]
# clustersDEG <- "Neutrophil Camp Lcn2"
DEgenes = list()
for ( i in clustersDEG) {
  DEgenes[[i]]=FindMarkers(seurnodblct,subset.ident=i, pseudocount.use = 0.01, group.by = "sample", ident.1 = "EPJ4",logfc.threshold=0.8) 
  DEgenes[[i]]$cluster=i
  DEgenes[[i]]$group.1="EPJ4"
  DEgenes[[i]]$group.2="EPJ3"
}
DEnames <- c()
for ( i in clustersDEG) {
  sig <- DEgenes[[i]]$p_val_adj < 0.1
  DEnames <- c(DEnames, rownames(DEgenes[[i]][sig,]))
}
DEnames <- unique(DEnames)
library(EnhancedVolcano)
n=0
for ( i in names(DEgenes)) {
  n=n+1
  DEnames <- unique(DEnames)
  DEgenes[[i]][["avg_log2FC"]] <- as.numeric(DEgenes[[i]][["avg_log2FC"]])
  keyvals <- ifelse(
    DEgenes[[i]][["avg_log2FC"]] < -0.8, 'dodgerblue3',
    ifelse(DEgenes[[i]][["avg_log2FC"]] > 0.8, 'firebrick3',
           'slategray'))
  keyvals[is.na(keyvals)] <- 'slategray'
  names(keyvals)[keyvals == 'firebrick3'] <- 'upregulated'
  names(keyvals)[keyvals == 'slategray'] <- 'ns'
  names(keyvals)[keyvals == 'dodgerblue3'] <- 'downregulated'
  print(EnhancedVolcano(DEgenes[[i]],
                        lab = row.names(DEgenes[[i]]),
                        x = "avg_log2FC",
                        y = "p_val_adj",
                        title = 'Orthotopic LLC-bearing vs Healthy C57BL/6',
                        subtitle = paste("Differential gene expression -",i),
                        caption = paste("fold change cutoff, 0.8; p-value cutoff, 0.1"),
                        pCutoff = 0.1,
                        FCcutoff = 0.8,
                        pointSize = 6.0,#normally 3
                        labSize = 9, #normally 6
                        colCustom = keyvals,
                        max.overlaps=20,
                        drawConnectors = TRUE,
                        widthConnectors = 1.0,
                        colConnectors = 'black',
                        arrowheads = FALSE))
  dev.off()
}

DEnames2 <- c()
for ( i in clustersDEG) {
  sig <- abs(DEgenes[[i]]$"avg_log2FC") > 0.8 & DEgenes[[i]]$p_val_adj < 0.1
  DEnames2 <- c(DEnames2, rownames(DEgenes[[i]][sig,]))
}
DEnames2 <- unique(DEnames2)

####pheatmap- UPregulated and DOWNregulated genes####
#Create a data frame with all the genes per cluster that are upregulated and downregulated

#DEgenes was calculated for volcano plots
Bcell <- row.names(DEgenes[["B cell"]])
NKT <- row.names(DEgenes[["NK/T cell"]])
Baso<- row.names(DEgenes[["Basophil/eosinophil prog"]])
Neu<- row.names(DEgenes[["Neutrophil"]])
Prog<- row.names(DEgenes[["Progenitor"]])
Mono<-row.names(DEgenes[["Monocyte/DC"]])
Ery<-row.names(DEgenes[["Erythroid"]])

df = list('Progenitor'= Prog,'Neutrophil'= Neu,'Monocyte/DC'= Mono,'Erythroid'=Ery,'Basophil/eosinophil prog'=Baso,'B cell'= Bcell, 'NK/T cell' = NKT   
)
attributes(df) = list(names = names(df),
                      row.names=1:max(length(Prog),length(Neu),length(Mono),length(Ery), length(Baso),length(Bcell), length(NKT)
                      ), class='data.frame')
#create a binary data frame in which includes the 134 genes with the 8 clusters
#fill it with 1 and 0, if the name of the gene is present in that cluster or not
bin<- data.frame(matrix(ncol=length(DEnames),nrow=ncol(df)))

rownames(bin) <- colnames(df)
colnames(bin) <- DEnames
for (i in clustersDEG){
  for( j in 1:length(DEnames)){ 
    bin[i,j]<- ifelse(DEnames[j] %in% df[,i] ,1,0)}}
bin
library(pheatmap)
#binary pheatmap:red and blue 
pheatmap(bin, cluster_rows = F, cluster_cols = F, cellwidth = 8, cellheight = 8)
dev.off()
#Replacement of bin (previous data frame) with log fold change values
for (i in clustersDEG){
  for (x in 1:length(DEnames)){
    bin[i,x]<- ifelse(DEnames[x] %in% rownames(DEgenes[[i]]) ,DEgenes[[i]][which(rownames(DEgenes[[i]]) %in% DEnames[x]),2],0)}}

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length 
pallettelengthsmyBreaks <- c(seq(min(-2), 0, length.out=ceiling(paletteLength/2) + 1), 
                             seq(max(4.5)/paletteLength, max(4.5), length.out=floor(paletteLength/2)))
pheatmap(bin,cluster_rows = F,cluster_cols = F,cellwidth = 8,cellheight = 8,color=myColor,breaks = pallettelengthsmyBreaks)
dev.off()

####Heatmap up and downregulated genes - ranked from fewer to greater appearance (expression)####
gene_count <- sapply(DEnames, function(gene) sum(sapply(df, function(cluster_genes) gene %in% cluster_genes)))

# Sort DEnames based on gene_count (ascending order)
sorted_genes <- DEnames[order(gene_count)]

# Create a binary data frame (bin) with the sorted gene order
bin <- data.frame(matrix(ncol=length(sorted_genes),nrow=length(df)))

rownames(bin) <- names(df)
colnames(bin) <- sorted_genes

# Fill it with 1 and 0, if the name of the gene is present in that cluster or not
for (i in names(df)){
  for (j in 1:length(sorted_genes)){
    bin[i,j] <- ifelse(sorted_genes[j] %in% df[[i]], 1, 0)
  }
}
# Replace the binary matrix with log fold change (LFC) values, using the sorted genes
for (i in names(df)){
  for (x in 1:length(sorted_genes)){
    gene <- sorted_genes[x]
    bin[i,x] <- ifelse(gene %in% rownames(DEgenes[[i]]), 
                       DEgenes[[i]][which(rownames(DEgenes[[i]]) == gene), 2], 
                       0)
  }
}
# Define the color palette and breaks
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(-2), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(4.5)/paletteLength, max(4.5), length.out=floor(paletteLength/2)))

# Plot the heatmap
pheatmap(bin, cluster_rows = F, cluster_cols = F, cellwidth = 8, cellheight = 8, color=myColor, breaks = myBreaks)
dev.off()

####Choose only the genes expressed in at least two clusters, ranked from the fewest to the greatest
# Count how many clusters each gene appears in
gene_count <- sapply(DEnames, function(gene) sum(sapply(df, function(cluster_genes) gene %in% cluster_genes)))

# Filter genes that are present in at least 2 clusters
filtered_genes <- DEnames[gene_count >= 2]

# Sort the filtered genes by occurrence (ascending order)
sorted_genes <- filtered_genes[order(gene_count[gene_count >= 2])]

# Create a binary data frame (bin) with the sorted gene order
bin <- data.frame(matrix(ncol=length(sorted_genes), nrow=length(df)))

rownames(bin) <- names(df)
colnames(bin) <- sorted_genes

# Fill it with 1 and 0, indicating if the gene is present in the respective cluster
for (i in names(df)){
  for (j in 1:length(sorted_genes)){
    bin[i, j] <- ifelse(sorted_genes[j] %in% df[[i]], 1, 0)
  }
}

# Replace the binary matrix with log fold change (LFC) values, using the sorted genes
for (i in names(df)){
  for (x in 1:length(sorted_genes)){
    gene <- sorted_genes[x]
    bin[i, x] <- ifelse(gene %in% rownames(DEgenes[[i]]), 
                        DEgenes[[i]][which(rownames(DEgenes[[i]]) == gene), 2], 
                        0)
  }
}

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#Heatmap of gene of interest - Fig 3B 
# ---- Select the gene of interest ----
gene_of_interest <- "Mki67"

# Ensure the gene exists
if (!gene_of_interest %in% DEnames) {
  stop("Mki67 not found in DEnames")
}

# Create matrix
bin <- data.frame(matrix(ncol = 1, nrow = length(df)))
rownames(bin) <- names(df)
colnames(bin) <- gene_of_interest

# ---- Fill with log fold change values FIRST ----
for (i in names(df)) {
  bin[i, 1] <- ifelse(
    gene_of_interest %in% rownames(DEgenes[[i]]),
    DEgenes[[i]][gene_of_interest, 2],
    0
  )
}


rownames(bin) <- cluster_map[rownames(bin)]

# ---- Define colors and breaks ----
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

myBreaks <- c(
  seq(-2, 0, length.out = ceiling(paletteLength/2) + 1),
  seq(4.5/paletteLength, 4.5, length.out = floor(paletteLength/2))
)

# ---- Plot heatmap ----
tiff(
  file = "Results ECE/Integration-nodoublets/20230613/Graphs/wholedataset-Mki67_HEatmap.tiff",
  units = "in",
  width = 6,
  height = 4,
  res = 400
)

pheatmap(
  bin,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 40,
  cellheight = 15,
  color = myColor,
  breaks = myBreaks
)

dev.off()

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#

#Feature plots#### Select the gene of interest - Fig  3A, 3C, 4G, Extended Data Fig 3; Extended Data Fig 4 
FeaturePlot(seurnodblct,features = "Lcn2",split.by = "sample")
dev.off()
