
# Hashed scRNAseq analysis - murine BM from healthy and Lewis lung carcinoma(LLC)-bearing mice  ####
# EPJ3 corresponds to healthy mice
# EPJ4 corresponds to LLC bearing mice 


# the folder containing the unzipped datafolder
# setwd("")

# memory.limit(size=24000)
#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# Install packages####
install.packages("Seurat")
BiocManager::install(pkgs ="scater")
install.packages("ggvenn")
install.packages("ggpubr")
install.packages("dplyr")
BiocManager::install(pkgs ="DropletUtils") 
install.packages("VennDiagram")
BiocManager::install(pkgs ="scDblFinder")
install.packages("SingleR")
BiocManager::install(pkgs = "SingleR")
BiocManager::install(pkgs ="celldex")
install.packages("HGNChelper")
install.packages("clustree")
install.packages("cowplot")
install.packages("Matrix")
library(Matrix)

# Load packages ####
library(Seurat)
library(scater)
library(ggvenn)
library(ggpubr)
library(dplyr)
library(DropletUtils)
library(VennDiagram)
library(scDblFinder)
library(SingleR)
library(celldex)
library(HGNChelper)
library(clustree)
library(cowplot)

dataset.name="EPJ3_4"


#Select colors you want to use in the dimension plot
cols.use=c("brown1","goldenrod3","goldenrod1","darkorchid","mediumpurple",
           "darkolivegreen1","magenta","dodgerblue","turquoise3","lightblue",
           "gold4","coral","violetred","grey","green3","goldenrod1","palevioletred",
           "yellow", "brown3", "pink", "bisque", "yellowgreen", "grey30","purple","lightgreen","royalblue","pink")


#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# Loading input data from Cellranger output ####
## We work wih data from Cellranger v3
path_data="data/"
sample.names<-c("EPJ3","EPJ4")


#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# We use the Read10X function to read in the raw unfiltered expression matrix from the CellRanger output
#Adjust the Read10X function a bit, in case you struggle with Csparsematrix ####
Read10X <- function(
    data.dir,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
) {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    gene.loc <- file.path(run, 'genes.tsv')
    features.loc <- file.path(run, 'features.tsv.gz')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc) ) {
      stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    }
    data <- Matrix::readMM(file = matrix.loc)
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    } else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, , drop = FALSE])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "CsparseMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}
#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# Filtering empty droplets ####
# First create an empty list
seur=list()
for (i in 1:length(sample.names)){
  seur[[i]]<-CreateSeuratObject(counts = Read10X(paste0(path_data,sample.names[i],"_cDNA/outs/raw_feature_bc_matrix")))
  seur[[i]]$sample=sample.names[i]
  print(sample.names[i])
  print(dim(seur[[i]]))
}

# Save this list for easier work in future analysis

# Distinguish between droplets with cells and ambient RNA
# Set seed to ensure  ensures that any results that rely on randomness, e.g. subsampling or permutations, are reproducible.
set.seed(2022)

for (i in 1:length(sample.names)){
  results =DropletUtils:: emptyDrops(GetAssayData(seur[[i]], slot="counts"), niters = 5000)
  seur[[i]]$ED.cell <- results$FDR <= 0.01
  print(sample.names[i])
  print(table(Sig=seur[[i]]$ED.cell, Limited=results$Limited))
}

# Create the knee plot
# A useful diagnostic for droplet-based data is the barcode rank plot, 
# which shows the (log-)total UMI count for each barcode on the y-axis and the (log-)rank on the x-axis. 
# This is effectively a transposed empirical cumulative density plot with log-transformed axes. 
# It is useful as it allows users to examine the distribution of total counts across barcodes, focusing on those with the largest counts.

for (i in 1:length(sample.names)){
  br.out <- DropletUtils::barcodeRanks(GetAssayData(seur[[i]], slot="counts"))
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total", main=sample.names[i])
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=min(br.out$total[seur[[i]]$ED.cell], na.rm=T), col="red", lty=2)
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), 
         legend=c("knee", "inflection","ED-min"))
  # dev.off()
}

#remove variables that are unnecessary in the further analysis; to save memory space
rm(o)
rm(br.out)
rm(results)

#Filter out the empty droplets
for (i in 1:length(sample.names)){
  print(sample.names[i])
  seur[[i]]<-seur[[i]][, seur[[i]]$ED.cell &!is.na(seur[[i]]$ED.cell)] # here you select all columns(droplets) which have value TRUE and discard all values with FALSE and NAa
  print(dim(seur[[i]]))
}


#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# Integration of both samples (EPJ3 and EPJ4) ####
#https://satijalab.org/seurat/articles/integration_introduction.html
# normalize and identify variable features for each dataset independently
seur <- lapply(X = seur, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seur)
# We then identify anchors using the FindIntegrationAnchors() function, 
# which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
immune.anchors <- FindIntegrationAnchors(object.list = seur, anchor.features = features)
# this command creates an 'integrated' data assay
# creates a combined large seurat object
seur.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seur.combined) <- "RNA"


# rename all values in the ED.cell from TRUE to NULL
seur.combined$ED.cell=NULL
# show the amount of droplets (columns) in sample EPJ3 and EPJ4
table(seur.combined$sample)
#rename the seuratobject from seur.combined to seur
#rename the list of both samples to seur.list 
seur.list <- seur
seur <- seur.combined

#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# QC cells ####
# Find outliers for total UMI counts, number of genes and percent mitochondrial genes per cell
# Outliers are determined based on the median absolute deviation (MAD).
# The goal is to remove the long tail before/after the peak in the ditribution 
# of the QC metric. It is possible to modify the nmads parameter 
# (minimum number of MADs away from median, required for a value to be called an outlier), 
# or to set the threshold manually (e.g. remove all cells with percent mitochondrial genes above 40).


#Calculate UMI count per cell 
# UMI counts are stored in the nCount_RNA column in the Seurat object
# outlier treshold can be adjusted via the nmads variable

outliers=c()
for ( i in unique(seur$sample)){
  outliers=c(outliers, 
             scater::isOutlier(seur$nCount_RNA[seur$sample==i], nmads=3, type="lower", log=TRUE)
  )
}
seur$nUMI.outlier.low <- outliers[colnames(seur)]
cat("Outliers:",sum(seur$nUMI.outlier.low))

# Create histograms of UMIs per cell 
for ( i in unique(seur$sample)){
  hist(seur$nCount_RNA[seur$sample==i],
       breaks = 100,xlab="nCount_RNA",
       main=paste0("Total UMI counts per cell: ",i))
  if(sum(seur$sample==i & seur$nUMI.outlier.low)!=0)
    abline(v = max(seur$nCount_RNA[seur$sample==i & seur$nUMI.outlier.low]), col = "red")
  # dev.off()
}

# Create violin plot
for ( i in unique(seur$sample)){
  print(ggplot(as.data.frame(seur[[]])[seur$sample==i,], aes(1, nCount_RNA)) + 
          geom_violin(fill="gray80") +theme_classic()+ theme(axis.title.x = element_blank())+
          geom_jitter(height = 0, width = 0.3, aes(col=nUMI.outlier.low)) +
          scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle(paste0("Total UMI counts per cell: ",i)))
  #dev.off()
}

# Number of genes per cell
outliers=c()
for ( i in unique(seur$sample)){
  outliers=c(outliers, 
             scater::isOutlier(seur$nFeature_RNA[seur$sample==i], nmads=3, type="lower", log=TRUE) # number of genes is saved in the nFeature_RNA column
  )
}
seur$nGene.outlier.low <- outliers[colnames(seur)]
cat("Outliers:",sum(seur$nGene.outlier.low))

# Create a histogram 
for ( i in unique(seur$sample)){
  hist(seur$nFeature_RNA[seur$sample==i],
       breaks = 100,xlab="nCount_RNA",
       main=paste0("Number of genes per cell: ",i))
  if(sum(seur$sample==i & seur$nGene.outlier.low)!=0)
    abline(v = max(seur$nFeature_RNA[seur$sample==i & seur$nGene.outlier.low]), col = "red")
  #  dev.off()
}

# Create violin plots
for ( i in unique(seur$sample)){
  print(ggplot(as.data.frame(seur[[]])[seur$sample==i,], aes(1, nFeature_RNA)) + 
          geom_violin(fill="gray80") +theme_classic()+ theme(axis.title.x = element_blank())+
          geom_jitter(height = 0, width = 0.3, aes(col=nGene.outlier.low)) +
          scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle(paste0("Number of genes per cell: ",i)))
  #  dev.off()
}


# Proportion of mitochondrial genes per cell
# Calculate % mitochondrial genes per cell
seur[["perc.mito"]] <- PercentageFeatureSet(seur, pattern = "^mt-",assay = 'RNA')  #PercentageFeatureSet uses the raw counts matrix, which is not present in the integrated assay. If you pass assay = "RNA", then PercentageFeatureSet will run correctly
# We include all genes starting with "mt-" for calculating this QC satistic. Depending on the organism, the search pattern might need to be modified, e.g. to "MT-" for human
summary(seur$perc.mito)

outliers=c()
for ( i in unique(seur$sample)){
  outliers=c(outliers, 
             scater::isOutlier(seur$perc.mito[seur$sample==i], nmads=3, type="higher", log=TRUE)
  )
}
seur$mito.outlier.high <- outliers[colnames(seur)]
cat("Outliers:",sum(seur$mito.outlier.high))

# Create Histograms
for ( i in unique(seur$sample)){
  hist(seur$perc.mito[seur$sample==i],
       breaks = 100,xlab="perc.mito",
       main=paste0("% mito genes per cell: ",i))
  if(sum(seur$sample==i & seur$mito.outlier.high)!=0)
    abline(v = min(seur$perc.mito[seur$sample==i & seur$mito.outlier.high]), col = "red")
  #  dev.off()
}

#Create violin plots
for ( i in unique(seur$sample)){
  print(ggplot(as.data.frame(seur[[]])[seur$sample==i,], aes(1, perc.mito)) + 
          geom_violin(fill="gray80") +theme_classic()+ theme(axis.title.x = element_blank())+
          geom_jitter(height = 0, width = 0.3, aes(col=mito.outlier.high)) +
          scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle(paste0("% mito genes per cell: ",i)))
  #  dev.off()
}


# Overlap of cells, outliers for UMI counts, number of genes and % mitochodrial genes per cell
for ( i in unique(seur$sample)){
  v <-venn.diagram(
    list (nUMI=colnames(seur)[seur$nUMI.outlier.low & seur$sample==i],
          nGene=colnames(seur)[seur$nGene.outlier.low & seur$sample==i],
          perc.mito=colnames(seur)[seur$mito.outlier.high & seur$sample==i]),
    filename=NULL,main=i,
    alpha = c( 0.5,0.5,0.5),
    fill = c("darkgreen","orangered","mediumpurple")
  )
  grid.newpage()
  grid.draw(v)
  dev.off()
  #rm(v)
}

#For now we only remove the perc.mito outliers
seur_clean<-seur[, !seur$mito.outlier.high]
print(dim(seur))
print(dim(seur_clean))
table(seur_clean$sample)

#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# QC genes ####

# To define a cutoff of lowly-abundant genes, we plot the distribution of log-means across all genes.
# The cutoff is placed in middle of the rectangular component of the graph before the peak.
ave.counts=list()
for ( i in 1: length(unique(seur_clean$sample))) {
  ave.counts[[i]] <- rowMeans(as.matrix(GetAssayData(seur_clean)[,seur_clean$sample==unique(seur$sample)[i]]))
}

thresholds<-c(0.005,0.005)
for ( i in 1: length(unique(seur_clean$sample))) {
  # tiff(paste("ResultsKI/integclean/Graphs/UMIperGene-Histogram_sample",i,".tiff",sep = ""), units="in", width=6, height=5, res=400)
  hist(log10(ave.counts[[i]]), breaks=100, main=paste0("Histogram of mean UMI counts per gene: ",unique(seur$sample)[i]), col="grey80",
       xlab=expression(Log[10]~"mean count per gene"))
  abline(v=log10(thresholds[i]), col="blue", lwd=2, lty=2)
  dev.off()
}

#number of genes to keep
usegenes=list()
for ( i in 1: length(unique(seur_clean$sample))) {
  print(unique(seur$sample)[i])
  usegenes[[i]]<-ave.counts[[i]]>thresholds[i]
  print(table(usegenes[[i]]))
}

#Filter out the lowly-abundant genes that overlap between all samples
i=1
genes.filter=names(usegenes[[i]][! usegenes[[i]]])
for ( i in 2: length(unique(seur_clean$sample))) {
  genes.filter=intersect(genes.filter, names(usegenes[[i]][! usegenes[[i]]])) #find the intersect between lowly abundant genes
}
usegenes.final=!rownames(seur_clean) %in% genes.filter
table(usegenes.final)

seur_clean<-seur_clean[usegenes.final, ]
print(dim(seur))
print(dim(seur_clean))
table(seur_clean$sample)

# Doublet finder
doublet.score<- scDblFinder::scDblFinder(as.SingleCellExperiment(seur_clean,assay="RNA"), samples="sample",returnType="table") ###  BPPARAM=BiocParallel::MulticoreParam(3) Using 3 cores

seur_clean <- AddMetaData(object = seur_clean, metadata = doublet.score$scDblFinder.score, col.name = "doublet.score")
seur_clean <- AddMetaData(object = seur_clean, metadata = doublet.score$scDblFinder.class, col.name = "doublet.class")

table(seur_clean$doublet.class)
summary(seur_clean$doublet.score)


#---#---#---#---#---#---#---#---#---#---#---#---#
#---#---#---#---#---#---#---#---#---#---#---#---#
# Demultiplexing ####

CMO_tags <-readxl::read_excel(paste0(getwd(),"CellPlex_CMO_tags")
                              CMO_tags$`CMO-tag`=paste0("CMO",CMO_tags$`CMO-tag`)
                              
                              # For each of both groups put assignment confidence table in a list
                              data=list()
                              for (i in 1:length(sample.names)){
                                data[[i]]<-read.csv(paste0(getwd(),sample.names[i],"_multi/outs/multi/multiplexing_analysis/assignment_confidence_table.csv"))
                                data[[i]]$sample=sample.names[i]
                              }
                              #put both list in one big dataframe, pasting data as rows
                              data <- dplyr::bind_rows(data)
                              
                              head(data$Barcodes) #one column of data
                              head(colnames(seur_clean)) #colnames of seurat object
                              
                              meta <- seur_clean[[]] #take all metadata from the seurat object
                              meta$barcode <-sapply(strsplit(rownames(meta),"_"),"[[",1) #save the rownames in a column, but cut of the _1 part
                              meta$barcode <- paste0(meta$sample,"_", meta$barcode) #paste the samplenames (=groupname) before the barcode
                              data$barcode <- paste0(data$sample,"_", data$Barcodes) # paste the samplenames before the barcode in the assignment confidence table 
                              length(intersect(data$barcode,meta$barcode)) # how many "cells" can be matched between both datasets
                              nrow(meta)
                              nrow(data) # there is a difference here due to filtering in the seurat object
                              
                              
                              meta$cell <- row.names(meta)
                              # paste some of the information of the assignment confidence table in the metadata of the cells 
                              meta <- merge(meta,data[,c("barcode","Assignment","Assignment_Probability")], by="barcode",all.x = T ) 
                              row.names(meta) <- meta$cell
                              # add the gathered metadata to the metadata in seurat object
                              seur_clean <- AddMetaData(seur_clean, metadata = meta)
                              table(seur_clean$Assignment)
                              
                              data=list()
                              for (i in 1:length(sample.names)){
                                data[[i]]<-read.csv(paste0(getwd(),sample.names[i],"_multi/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv"))
                                data[[i]]$sample=sample.names[i]
                              }
                              data <- dplyr::bind_rows(data)
                              
                              data$barcode <- sapply(strsplit(data$cell_barcode,"_"),"[[",1)
                              data$barcode <- paste0(data$sample,"_", data$cell_barcode)
                              length(intersect(data$barcode,meta$barcode))
                              colnames(data) <- c("cell_barcode","num_features","CMO_call","nCount_CMO","sample","barcode") #change of colnames needed, otherwise not possible to match the data and meta file
                              nrow(meta)
                              nrow(data)
                              
                              meta <- merge(meta,data[,c('cell_barcode',"barcode","CMO_call", "nCount_CMO" )], by="barcode",all.x = T )
                              row.names(meta) <- meta$cell
                              
                              seur_clean <- AddMetaData(seur_clean, metadata = meta)
                              table(seur_clean$Assignment)
                              
                              
                              seur_clean$Assignment <- plyr::mapvalues(seur_clean$Assignment, from=CMO_tags$`CMO-tag`, to=CMO_tags$subsamples)
                              seur_clean$Assignment[is.na(seur_clean$Assignment)] <- "Unassigned"
                              seur_clean$CMO_call[is.na(seur_clean$CMO_call)] <-"Unassigned"
                              seur_clean$nCount_CMO[is.na(seur_clean$nCount_CMO)] <-0
                              table(seur_clean$Assignment)
                              
                              
                              
                              
                              seur_clean$nCount_CMO_call <- seur_clean$nCount_CMO
                              temp <- seur_clean$nCount_CMO_call[grep("[|]",seur_clean$nCount_CMO_call)]
                              temp <- sapply(strsplit(temp,"[|]"),as.numeric)
                              temp <- sapply(temp,sum)
                              identical(seur_clean@meta.data[grep("[|]",seur_clean$nCount_CMO),"cell"],names(temp))
                              seur_clean$nCount_CMO[grep("[|]",seur_clean$nCount_CMO)] <- temp
                              
                              seur_clean$nCount_CMO_log <- log1p(as.numeric(seur_clean$nCount_CMO))
                              FeaturePlot(seur_clean,c("nCount_CMO_log"), min.cutoff = "q10")
                              meta <- seur_clean@meta.data[,c("Assignment" ,"Assignment_Probability","CMO_call"  ,"nCount_CMO_call" ,"nCount_CMO"  ,  "nCount_CMO_log" ,  "nCount_CMO"    )]
                              
                              # DimPlot(seur_clean,repel =T,label=T, group.by = "Assignment", split.by = "sample", cols=cols.use) 
                              # 
                              # DimPlot(seurnodbl,cells.highlight = WhichCells(seurnodbl, expression = Assignment %in% c("Multiplet", "Unassigned", "Blanks"), invert=T), sizes.highlight = 0.2)+ggtitle("Singlets") +NoLegend()
                              # DimPlot(seurnodbl,cells.highlight = WhichCells(seurnodbl, expression = Assignment=="Multiplet"), sizes.highlight = 0.2)+ggtitle("Multiplet") +NoLegend()
                              # FeaturePlot(seurnodbl,"doublet.score")
                              # DimPlot(seurnodbl,cells.highlight = WhichCells(seurnodbl, expression = Assignment=="Blanks"), sizes.highlight = 0.2)+ggtitle("Blanks") +NoLegend()
                              # DimPlot(seurnodbl,cells.highlight = WhichCells(seurnodbl, expression = Assignment=="Unassigned"), sizes.highlight = 0.2)+ggtitle("Unassigned") +NoLegend()
                              
                              
                              
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              DefaultAssay(seur_clean) <- "integrated"
                              # Standard workflow for visualization and clustering ####
                              seur_clean <- ScaleData(seur_clean, verbose = FALSE)
                              seur_clean <- RunPCA(seur_clean, npcs = 60,  verbose = FALSE)
                              
                              
                              #PC selection
                              ElbowPlot(object = seur_clean,ndims =50)
                              DimHeatmap(seur_clean, dims = 25:42, cells = 5000, balanced = TRUE)
                              
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              # UMAP DIMENSIONALITY REDUCTION ####
                              ### PC selection for downstream analysis
                              dims.use<-30
                              ### Run UMAP
                              seur_clean <- RunUMAP(seur_clean, dims = 1:dims.use, verbose=F)
                              
                              # Visualize QC metrics
                              # tiff("ResultsKI/integclean/Graphs/QCUMAPS_PC30.tiff", units="in", width=8, height=7, res=400)
                              FeaturePlot(object = seur_clean,
                                          features = c("nCount_RNA", "nFeature_RNA", "perc.mito","doublet.score"))
                              dev.off()
                              
                              #example UMAP plot with one gene
                              FeaturePlot(seurnodbl,features=c("S100a9"),split.by = "sample")
                              
                              #split UMAP based on sample origin
                              # tiff("ResultsKI/integclean/Graphs/Dimplot_splitsample.tiff", units="in", width=10, height=8, res=400)
                              DimPlot(seur_clean, group.by = "sample")
                              dev.off()
                              
                              DimPlot(seur_clean,repel =T,label=T, group.by = "Assignment", split.by = "sample", cols=cols.use)
                              DimPlot(seur_clean,cells.highlight = WhichCells(seur_clean, expression = Assignment %in% c("Multiplet", "Unassigned", "Blanks"), invert=T), sizes.highlight = 0.2)+ggtitle("Singlets") +NoLegend()
                              DimPlot(seur_clean,cells.highlight = WhichCells(seur_clean, expression = Assignment=="Multiplet"), sizes.highlight = 0.2)+ggtitle("Multiplet") +NoLegend()
                              FeaturePlot(seur_clean,"doublet.score")
                              DimPlot(seur_clean,cells.highlight = WhichCells(seur_clean, expression = Assignment=="Blanks"), sizes.highlight = 0.2)+ggtitle("Blanks") +NoLegend()
                              DimPlot(seur_clean,cells.highlight = WhichCells(seur_clean, expression = Assignment=="Unassigned"), sizes.highlight = 0.2)+ggtitle("Unassigned") +NoLegend()
                              
                              DimPlot(seur_clean,repel =T,label=T, group.by = "sample")
                              # DimPlot(seur_clean,repel =T,label=T,cols=cols.use) 
                              DimPlot(seur_clean,repel =T,label=T, group.by = "CMO_call", split.by = "sample")
                              DimPlot(seur_clean,repel =T,label=T, group.by = "Assignment", split.by = "sample")
                              
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              # remove the Multiplets/Blanks/Unassigned ####
                              DefaultAssay(seur_clean) <- "RNA"
                              seur_clean <- SetIdent(seur_clean, value = seur_clean@meta.data$Assignment)
                              Idents(seur_clean)
                              clusters <- c("Healthy_1", "Healthy_2", "Healthy_3", "Healthy_4", "LLC_5","LLC_6", "LLC_7", "LLC_14")
                              clusters[!clusters %in% Idents(seur_clean)] #check if some of the cluster names is not correct
                              #Subset the cells
                              cells=WhichCells(seur_clean,idents = clusters)
                              seur<-seur_clean[,cells]
                              dim(seur_clean)
                              dim(seur)
                              
                              
                              seur <- SplitObject(seur, split.by = "sample")
                              
                              # #https://satijalab.org/seurat/articles/integration_introduction.html
                              # # normalize and identify variable features for each dataset independently
                              # seurnodbl_l <- lapply(X = seurnodbl_l, FUN = function(x) {
                              #   x <- NormalizeData(x)
                              #   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                              # })
                              # select features that are repeatedly variable across datasets for integration
                              features <- SelectIntegrationFeatures(object.list = seur)
                              
                              # We then identify anchors using the FindIntegrationAnchors() function, 
                              # which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
                              immune.anchors <- FindIntegrationAnchors(object.list = seur, anchor.features = features)
                              # this command creates an 'integrated' data assay
                              # creates a combined large seurat object
                              seur.combined <- IntegrateData(anchorset = immune.anchors)
                              
                              # specify that we will perform downstream analysis on the corrected data note that the
                              # original unmodified data still resides in the 'RNA' assay
                              DefaultAssay(seur.combined) <- "integrated"
                              
                              seurnodbl <- seur.combined
                              # Run the standard workflow for visualization and clustering
                              
                              seurnodbl <- ScaleData(seurnodbl, verbose = FALSE)
                              # seur_clean <- FindVariableFeatures(seur_clean,verbose=F)
                              # seur_clean <- RunPCA(seur_clean, npcs = 60, features =VariableFeatures(seur_clean), verbose = FALSE)
                              seurnodbl <- RunPCA(seurnodbl, npcs = 60,  verbose = FALSE)
                              
                              
                              #PC selection
                              ElbowPlot(object = seurnodbl,ndims =50)
                              # DimHeatmap(seurnodbl, dims = 25:42, cells = 5000, balanced = TRUE)
                              
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              # UMAP DIMENSIONALITY REDUCTION ####
                              
                              ### PC selection for downstream analysis
                              dims.use<-30
                              
                              seurnodbl <- RunUMAP(seurnodbl, dims = 1:dims.use, verbose=F)
                              
                              # CLUSTERING nodbl####
                              # We run Leiden clustering with varying the clustering resolution
                              seurnodbl <- FindNeighbors(seurnodbl, dims = 1:dims.use, 
                                                         graph.name = paste0("RNA_snn_PC",dims.use), verbose=F)
                              
                              for ( i in seq(0.1,0.8, 0.1)) #Will take a while to run this piece of code
                                seurnodbl <- FindClusters(seurnodbl, resolution = i, algorithm = 4, graph.name=paste0("RNA_snn_PC",dims.use), verbose=F,method="igraph") # algorithm=4 -Leiden clustering
                              
                              clustree::clustree(seurnodbl, prefix = paste0("RNA_snn_PC",dims.use,"_res."))+
                                ggtitle(paste("PC =",dims.use))
                              dev.off()
                              
                              
                              plot<-list()
                              for ( res in c(0.5,0.6,0.7,0.8))
                                plot[[as.character(res)]]<-DimPlot(seurnodbl,label=T,repel=T, group.by = paste0("RNA_snn_PC",dims.use,"_res.",res))+
                                ggtitle(paste("PC =",dims.use,"res=",res))
                              plot_grid(plotlist=plot)
                              dev.off()
                              
                              #check where there is still highest density of possible doublets
                              FeaturePlot(seurnodbl,"doublet.score")
                              dev.off()
                              
                              
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              #---#---#---#---#---#---#---#---#---#---#---#---#
                              # chose PC=30 res=0.6 before we go to FindMarkers ####
                              #We choose 0.6 as basophils are splitted from MEP at this resolution
                              res=0.6
                              Idents(seurnodbl)=  paste0("RNA_snn_PC",dims.use,"_res.",res)
                              Idents(seurnodbl)=  factor(Idents(seurnodbl),levels = 1:(length(unique(Idents(seurnodbl)))))
                              
                              DimPlot(seurnodbl, group.by =paste0("RNA_snn_PC",dims.use,"_res.",res),label=T,repel=T, cols=cols.use)
                              dev.off()
                              # For performing differential expression after integration, we switch back to the original data
                              DefaultAssay(seurnodbl) <- "RNA"
                              # Find differentially expressed genes per cluster nodb ####
                              #FindMarkers function is the one that takes a lot of RAM and time
                              DEgenes=list()
                              for (i in levels(Idents(seurnodbl))){
                                DEgenes[[i]]<-FindMarkers(seurnodbl, ident.1 = i,min.cells.group=2,pseudocount.use = 0.01, max.cells.per.ident = 1000,verbose = F) 
                                #The max.cells.per.ident is set lower to speed up this calculation, but should be Inf to do it 'completely'
                                #logfc.threshold=0.2 and min.diff.pct=0.2
                                #If you want only marker genes that are positive (Log2FC>0 or upregulated); only.pos=TRUE
                                DEgenes[[i]]$cluster=i
                                DEgenes[[i]]$pseudocount=0.01
                                DEgenes[[i]]$max.cells.per.ident=1000
                                DEgenes[[i]]$gene=row.names(DEgenes[[i]])
                              }
                              
                              
                              