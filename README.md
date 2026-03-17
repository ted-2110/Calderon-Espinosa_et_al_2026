# Calderon-Espinosa_et_al_2026
This repository contains the R scripts used to analyze the 10x Genomics hashed single‑cell RNA‑sequencing data supporting the findings of Calderon‑Espinosa et al., 2026. The dataset includes bone marrow samples from healthy (EPJ3) and Lewis lung carcinoma (LLC)–bearing mice (EPJ4).
**Data structure and Input Requirements**
The script expects the following directory structure for loading 10X Genomics data:
EPJ3_cDNA/outs/raw_feature_bc_matrix/
EPJ4_cDNA/outs/raw_feature_bc_matrix/
Each raw_feature_bc_matrix folder must contain: 
-barcodes.tsv
-features.tsv
-matrix.mtx
**Repository Contents**
**1. Loading Data_QC_integration**
This script loads the raw scRNA‑seq data into R, performs quality control, and integrates the datasets. It prepares the Seurat object used throughout the downstream analyses.
**2. Annotation_Cell type**
This script reproduces the clustering and reclustering of the Seurat object. It includes:
Marker‑based cell type annotation, volcano plots, heatmaps, cell abundance and transcriptomic shift analysis using the Cacoa package, and feature plots for gene expression visualization
**3. Density plot**
The script generates density plots according to sample condition. 
**4. Monocle3**
This script performs pseudotime trajectory analysis using Monocle3, including visualization of gene expression dynamics along differentiation trajectories.
**5. NicheNet**
This script runs NicheNet to infer ligand–receptor interactions and potential cell–cell communication networks within the bone marrow microenvironment.
