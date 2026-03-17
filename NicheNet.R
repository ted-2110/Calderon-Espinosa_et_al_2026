#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#

#####NicheNet analyis#### Fig 5E and Fig 5F
#https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md

getwd()
install.packages("remotes")
remotes::install_github("saeyslab/nichenetr")
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
seuratObj <- seurnodblct
seuratObj@meta.data %>% head()
seuratObj <- UpdateSeuratObject(seuratObj)


seuratObj@meta.data$annot %>% table() # note that the number of cells of some cell types is very low and should preferably be higher for a real application

DimPlot(seuratObj, reduction = "umap")
seuratObj@meta.data$sample %>% table()
## Healthy  LLC
## EPJ3     EPJ4 
## 6674     6867
DimPlot(seuratObj, reduction = "umap", group.by = "sample")


organism = "mouse"

if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

lr_network = lr_network %>% distinct(from, to)
head(lr_network)
## # A tibble: 6 × 2
##   from          to   
##   <chr>         <chr>
## 1 2300002M23Rik Ddr1 
## 2 2610528A11Rik Gpr15
## 3 9530003J23Rik Itgal
## 4 a             Atrn 
## 5 a             F11r 
## 6 a             Mc1r
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
#
#NicheNet analysis - Receiver: Progenitor
nichenet_output <- nichenet_seuratobj_aggregate(
  assay_oi = "RNA",
  seurat_obj = seuratObj,
  receiver = "Progenitor",
  condition_colname = "sample",
  condition_oi = "EPJ4",
  condition_reference = "EPJ3",
  geneset = "up",
  sender ="all",
  ligand_target_matrix = ligand_target_matrix,
  lr_network =lr_network,
  weighted_networks = weighted_networks,
  #organism ="mouse", #human and mouse are allowed
  expression_pct = 0.05) #gene expressed at least in 5% of the cell type

# Warning: Scaling data with a low number of groups may produce misleading results
# Warning message:
#   In nichenet_seuratobj_aggregate(seurat_obj = seuratObj, receiver = "Neutrophil",  :
#        Seurat object is result from the Seurat integration 
#        workflow. Make sure that the way of defining expressed 
#        and differentially expressed genes in this wrapper is appropriate for your 
#        integrated data.

nichenet_output$ligand_activity_target_heatmap
nichenet_output$ligand_receptor_heatmap
