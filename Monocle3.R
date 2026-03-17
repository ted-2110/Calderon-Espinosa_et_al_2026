
####Monocle-seurnodblct####
set.seed(1234)
devtools::install_github('cole-trapnell-lab/monocle3', lib="/path/to/your/personal/R-libraries")
remotes::install_github('satijalab/seurat-wrappers')
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)


# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3
DefaultAssay(seurnodblct) <- "RNA"
DimPlot(seurnodblct)
# subset cells ####

subset_H <- subset(seurnodblct, subset = sample == "EPJ3")
DimPlot(subset_H)
subset_LLC <- subset(seurnodblct, subset = sample == "EPJ4")
DimPlot(subset_LLC)
dim(seurnodblct)
#12242 13541
dim(subset_LLC)
# 12242  6867
dim(subset_H)
#12242  6674
#transform subset_LLC to cdsLLCLLC
cdsLLC <- as.cell_data_set(subset_LLC, assay ="RNA")
cdsLLC <- estimate_size_factors(cdsLLC)


cdsLLC

# to get cell metadata
colData(cdsLLC)
# to gene metdata
fData(cdsLLC)
rownames(fData(cdsLLC))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cdsLLC)$gene_short_name <- rownames(fData(cdsLLC))

# to get counts
counts(cdsLLC)


# see the partitions
plot_cells(cdsLLC, reduction_method = "UMAP", color_cells_by = 'partition', label_cell_groups=FALSE)
# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have
#cdsLLC <- partitionCells(cdsLLC)
# assign partitions
reacreate.partition <- c(rep(1,length(cdsLLC@colData@rownames))) #partition can only be 1, if not errors below
names(reacreate.partition) <- cdsLLC@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cdsLLC@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- subset_LLC@active.ident
cdsLLC@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cdsLLC@int_colData@listData$reducedDims$UMAP <- subset_LLC@reductions$umap@cell.embeddings

# ...3. Learn trajectory graph ------------------------
cdsLLC <- learn_graph(cdsLLC, use_partition = TRUE)
#cdsLLC <- learn_graph(cdsLLC, RGE_method = 'SimplePPT')
# learn trajectory graph
# visualise the learned trajectory

plot_cells(cdsLLC,
           color_cells_by = 'annot', #cell_type
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 3)
dev.off()

# expression of marker genes across the sample
plot_cells(cdsLLC, genes=c('Kit','S100a8','Flt3','Cd79b','Cd3d','Gata2','Car2'), reduction_method = "UMAP")
dev.off()


# specifying root cells: `root_pr_nodes` argument - check the principal points
plot_cells(cdsLLC,
           color_cells_by = "cluster", #cell_type
           label_cell_groups=FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_principal_points = TRUE,       # set this to TRUE
           graph_label_size=3)

#You can see now the principal points and their labels in the form Y_number. 
#Pick the principal point in the cluster that you expect to be the beginning of the trajectory
#and type its name in the root_pr_nodes argument when calling order_cells() function

# specifying root cells: `root_pr_nodes` argument - use the relevant principal point
cdsLLC <- order_cells(cdsLLC, root_pr_nodes='Y_146') #always check with previos graph
#other method is below
#cdsLLC <- order_cells(cdsLLC, reduction_method = 'UMAP', root_cells = colnames(cdsLLC[,clusters(cdsLLC) == "G0"]))

# plot cells in pseudotime

plot_cells(cdsLLC,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()
# cells ordered by monocle3 pseudotime


# access pseudotime calculated for each cell and store it alongside cell metadata
pseudotime <- pseudotime(cdsLLC) 	
cdsLLC@colData$pseudotime <- pseudotime 		 
data.pseudo <- as.data.frame(colData(cdsLLC))


# pseudotime(cdsLLC)
# cdsLLC$monocle3_pseudotime <- pseudotime(cdsLLC)
# data.pseudo <- as.data.frame(colData(cdsLLC))

#ggplot(data.pseudo, aes(monocle3_pseudotime,cdsLLC$annot))+geom_boxplot()
ggplot(data.pseudo, aes(pseudotime, reorder(cdsLLC$annot, pseudotime, median)))+ #instead of annot was cell_type
  geom_boxplot()+ ylab(paste0("Neutrophil cluster"))+RotatedAxis()+ggtitle ("Pseudotime")
#dev.off()


# make the subset cdsLLC to analyse only DE genes obtained from the analysis prior volcano plot
test_genes <- c("S100a9")
cdsLLC_subset <- cdsLLC[rowData(cdsLLC)$gene_short_name %in% test_genes,]

# produce violin plots
#tiff(file = "Results ECE/Integration-nodoublets/20230613/Graphs/Subcluster/Neutrophils/Violinplot-genes.tiff", units="in", width=7, height=3.5, res=400)
plot_genes_violin(cdsLLC_subset, group_cells_by="annot", ncol=2) #instead of annot is cell_type
#dev.off()
# test_genes2= DEnames[101:134]
# cdsLLC_subset2 <- cdsLLC[rowData(cdsLLC)$gene_short_name %in% test_genes,]
# 
# # produce violin plots
# plot_genes_violin(cdsLLC_subset2, group_cells_by="cell_type", ncol=2)
plot_genes_in_pseudotime(cdsLLC_subset, color_cells_by="annot", min_expr=0.5)
dev.off()

a <- plot_genes_in_pseudotime(cdsLLC_subset, color_cells_by="annot", min_expr=0.5)
a
# ...5. Finding genes that change as a function of pseudotime --------------------
DEg <- graph_test(cdsLLC, neighbor_graph = 'principal_graph', cores = 4)

DEg %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()
DEg

# Mrpl15      OK       0              70.56813 0.2480182          Mrpl15       0
# Rrs1        OK       0              66.25256 0.2327379            Rrs1       0
# Snhg6       OK       0              58.87542 0.2068614           Snhg6       0
# Arfgef1     OK       0              38.46631 0.1351336         Arfgef1       0
# Tram1       OK       0              38.75084 0.1361385           Tram1       0
# Rpl7        OK       0             132.14356 0.4644797            Rpl7       0
# visualizing pseudotime in seurat

subset_LLC$pseudotime <- pseudotime(cdsLLC)
#Idents(subset_LLC) <- subset_LLC$seurat_clusters

FeaturePlot(subset_LLC, features = "pseudotime", label = T)
dev.off()

###########################################################################################################
##########################################################################################################
#transform subset_H to cdsH
cdsH <- as.cell_data_set(subset_H, assay ="RNA")
cdsH <- estimate_size_factors(cdsH)


cdsH

# to get cell metadata
colData(cdsH)
# to gene metdata
fData(cdsH)
rownames(fData(cdsH))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cdsH)$gene_short_name <- rownames(fData(cdsH))

# to get counts
counts(cdsH)
# # clustering
# cdsH_auto_cluster <- cluster_cells(cdsH, reduction_method = "UMAP")
# # see the clusters
# plot_cells(cdsH_auto_cluster, reduction_method = "UMAP", color_cells_by = 'cluster', group_label_size = 5)     
# # clustering with changed resolution value
# cdsH_clustered <- cluster_cells(cdsH, reduction_method = "UMAP", resolution = 0.0008)
# # see the new clusters
# plot_cells(cdsH_clustered, reduction_method = "UMAP", color_cells_by = 'cluster', group_label_size = 5) 	
# 

#NEW*****************************************************************************************************************
# create a new column ‘cell_type’ and initialise it with clusters values
#colData(cdsH)$cell_type <- as.character(clusters(cdsH)) 	
# check the annotation
# tiff(file = "Results ECE/Integration-nodoublets/20230613/Graphs/Subcluster/Neutrophils/Cell-type_Monocle.tiff", units="in", width=7, height=3.5, res=400)
# plot_cells(cdsH, color_cells_by="cell_type", label_cell_groups=FALSE)    
# dev.off()
# 
# see the partitions
plot_cells(cdsH, reduction_method = "UMAP", color_cells_by = 'partition', label_cell_groups=FALSE)
# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have
#cdsH <- partitionCells(cdsH)
# assign partitions
reacreate.partition <- c(rep(1,length(cdsH@colData@rownames))) #partition can only be 1, if not errors below
names(reacreate.partition) <- cdsH@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cdsH@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- subset_H@active.ident
cdsH@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cdsH@int_colData@listData$reducedDims$UMAP <- subset_H@reductions$umap@cell.embeddings

# ...3. Learn trajectory graph ------------------------
cdsH <- learn_graph(cdsH, use_partition = TRUE)
#cdsH <- learn_graph(cdsH, RGE_method = 'SimplePPT')
# learn trajectory graph
# visualise the learned trajectory

plot_cells(cdsH,
           color_cells_by = 'annot', #cell_type
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 3)
dev.off()

# expression of marker genes across the sample
plot_cells(cdsH, genes=c('Kit','S100a8','Flt3','Cd79b','Cd3d','Gata2','Car2'), reduction_method = "UMAP")
dev.off()

# plot_cell_trajectory(cdsH,
#                      color_by = "cell_type2") +
#   scale_color_manual(values = cell_type_color)

# ...4. Order the cells in pseudotime -------------------

# specifying root cells: `root_pr_nodes` argument - check the principal points
plot_cells(cdsH,
           color_cells_by = "cluster", #cell_type
           label_cell_groups=FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_principal_points = TRUE,       # set this to TRUE
           graph_label_size=3)

#You can see now the principal points and their labels in the form Y_number. 
#Pick the principal point in the cluster that you expect to be the beginning of the trajectory
#and type its name in the root_pr_nodes argument when calling order_cells() function

# specifying root cells: `root_pr_nodes` argument - use the relevant principal point
cdsH <- order_cells(cdsH, root_pr_nodes='Y_39')
#other method is below
#cdsH <- order_cells(cdsH, reduction_method = 'UMAP', root_cells = colnames(cdsH[,clusters(cdsH) == "G0"]))

# plot cells in pseudotime

plot_cells(cdsH,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()
# cells ordered by monocle3 pseudotime


# access pseudotime calculated for each cell and store it alongside cell metadata
pseudotime <- pseudotime(cdsH) 	
cdsH@colData$pseudotime <- pseudotime 		 
data.pseudo <- as.data.frame(colData(cdsH))


# pseudotime(cdsH)
# cdsH$monocle3_pseudotime <- pseudotime(cdsH)
# data.pseudo <- as.data.frame(colData(cdsH))

#ggplot(data.pseudo, aes(monocle3_pseudotime,cdsH$annot))+geom_boxplot()
ggplot(data.pseudo, aes(pseudotime, reorder(cdsH$annot, pseudotime, median)))+ #instead of annot was cell_type
  geom_boxplot()+ ylab(paste0("Neutrophil cluster"))+RotatedAxis()+ggtitle ("Pseudotime")
dev.off()


# make the subset cdsH to analyse only DE genes obtained from the analysis prior volcano plot
#test_genes <- c("S100a9")

cdsH_subset <- cdsH[rowData(cdsH)$gene_short_name %in% test_genes,]

# produce violin plots
plot_genes_violin(cdsH_subset, group_cells_by="annot", ncol=2) #instead of annot is cell_type
dev.off()
# test_genes2= DEnames[101:134]
# cdsH_subset2 <- cdsH[rowData(cdsH)$gene_short_name %in% test_genes,]
# 
# # produce violin plots
# plot_genes_violin(cdsH_subset2, group_cells_by="cell_type", ncol=2)
plot_genes_in_pseudotime(cdsH_subset, color_cells_by="annot", min_expr=0.5)
dev.off()

b <- plot_genes_in_pseudotime(cdsH_subset, color_cells_by="annot", min_expr=0.5)



pseudodfb <- data.frame(expression=b[["plot_env"]][["cds_exprs"]][["expression"]],pseudotime=b[["data"]][["pseudotime"]], annot=  b[["data"]][["annot"]],condition="Healthy")
pseudodfa <- data.frame(expression=a[["plot_env"]][["cds_exprs"]][["expression"]],pseudotime=a[["data"]][["pseudotime"]], annot =  a[["data"]][["annot"]], condition="LLC")
pseudodf <- rbind(pseudodfb,pseudodfa)
ggplot(pseudodf, aes(x=pseudotime, y=expression, color=condition)) + geom_point()+ geom_smooth( aes(fill=condition))


ggplot(pseudodf) +geom_point(aes(x=pseudotime, y=expression, color=annot))+geom_smooth( aes(x=pseudotime, y=expression, color=condition))

ggplot(pseudodf) +geom_smooth( aes(x=pseudotime, y=expression, color=condition))
dev.off()


