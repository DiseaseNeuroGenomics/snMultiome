# -----------------------------------------------------------------------------
# Project: Single-cell Multi-Modal Analysis - filtering and clustering
# Date: August 2024
# Script name: filter_cluster.R
# Description: This script loads a pre-processed Seurat object containing 
#             multi-modal data (RNA and ATAC-seq), performs filtering, 
#             normalization, dimensional reduction, and clustering on the data,
#             and saves the processed object for downstream analysis.
# 
# outputs: Potential plots for Fig. 1b,e 
# -----------------------------------------------------------------------------

###################### LIBRARIES ###################### 

# Load necessary libraries
library(Seurat)  # Comprehensive toolkit for single-cell RNA-seq analysis
library(Signac)  # Extends Seurat for single-cell chromatin accessibility analysis (e.g., ATAC-seq)
library(Matrix)  # Sparse and dense matrix classes and methods
library(future)  # Parallel and distributed processing

# Set global options
options(future.globals.maxSize = 50000 * 1024^2)  # Increase the maximum allowed size for future objects to 50 GB

###################### LOAD DATA ###################### 

# Load the preprocessed Seurat object that contains multi-modal data (RNA and ATAC-seq)
obj_all <- readRDS('./data_processed/ALL_snMultiome.rds')

###################### FILTERING & CLUSTERING ###################### 

##---------------------------------- Filtering Cells ----------------------------------##
# Function to filter cells based on quality control metrics
# The function filters cells using RNA-seq and ATAC-seq specific metrics

filter_cells <- function(
    object, 
    min_nCount_RNA = 200, 
    max_nCount_RNA = 50000, 
    min_nCount_ATAC = 200, 
    max_nCount_ATAC = 100000, 
    max_percent_mt = 5, 
    max_nucleosome_signal = 3, 
    min_TSS_enrichment = 1) {
  filtered_object <- subset(
    x = object,
    subset = nCount_RNA > min_nCount_RNA &    # Keep cells with more than min_nCount_RNA RNA counts
      nCount_RNA < max_nCount_RNA &           # Keep cells with fewer than max_nCount_RNA RNA counts
      nCount_ATAC > min_nCount_ATAC &         # Keep cells with more than min_nCount_ATAC ATAC counts
      nCount_ATAC < max_nCount_ATAC &         # Keep cells with fewer than max_nCount_ATAC ATAC counts
      percent.mt < max_percent_mt &           # Remove cells with more than max_percent_mt% mitochondrial gene expression
      nucleosome_signal < max_nucleosome_signal & # Remove cells with nucleosome signal > max_nucleosome_signal
      TSS.enrichment > min_TSS_enrichment     # Keep cells with TSS enrichment > min_TSS_enrichment
  )
  return(filtered_object)
}

# Filter the cells in the Seurat object
obj_filt <- filter_cells(obj_all)

##---------------------------------- Processing and Clustering ----------------------------------##
# Function to perform data processing and clustering
# This function performs normalization, dimensional reduction, and clustering separately on RNA-seq and ATAC-seq data
process_and_cluster <- function(object, n.pc = 30, n.lsi = 10) {
  
  ##----- Filter features (genes/peaks) detected in fewer than 10 cells -----
  filter_features <- function(object, assay_name) {
    feature_counts <- Matrix::rowSums(object[[assay_name]]@counts > 0)  # Calculate the number of cells in which each feature is detected
    filtered_features <- names(which(feature_counts >= 10))  # Keep only features detected in at least 10 cells
    object[[assay_name]] <- subset(object[[assay_name]], features = filtered_features)  # Subset the assay by filtered features
    return(object)
  }
  
  object <- filter_features(object, "RNA")  # Filter RNA features
  object <- filter_features(object, "ATAC") # Filter ATAC features
  
  ##----- RNA-seq Processing and Clustering -----
  DefaultAssay(object) <- 'RNA'  # Set the default assay to RNA for RNA-seq analysis
  object <- SCTransform(object)  # Normalize RNA data using SCTransform, which adjusts for sequencing depth and other confounders
  object <- RunPCA(object)  # Perform Principal Component Analysis (PCA) for dimensionality reduction
  object <- RunUMAP(object, reduction = 'pca', dims = 1:n.pc, assay = 'SCT',
                    reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')  # Run UMAP for visualization based on PCA
  object <- FindNeighbors(object, reduction = 'pca', dims = 1:n.pc, assay = 'SCT')  # Identify nearest neighbors based on PCA
  object <- FindClusters(object, graph.name = 'SCT_snn', algorithm = 3, resolution = 0.2)  # Cluster cells based on nearest neighbors
  
  ##----- ATAC-seq Processing and Clustering -----
  DefaultAssay(object) <- 'ATAC'  # Set the default assay to ATAC for ATAC-seq analysis
  object <- RunTFIDF(object, method = 3)  # Apply Term Frequency-Inverse Document Frequency (TF-IDF) normalization to ATAC data
  object <- FindTopFeatures(object, min.cutoff = 'q75')  # Select the top features based on a quantile cutoff (top 25%)
  object <- RunSVD(object)  # Perform Singular Value Decomposition (SVD) for dimensionality reduction
  object <- RunUMAP(object, reduction = 'lsi', dims = 2:n.lsi, assay = 'ATAC',
                    reduction.name = "umap.atac", reduction.key = "atacUMAP_")  # Run UMAP for visualization based on LSI
  object <- FindNeighbors(object, reduction = 'lsi', dims = 2:n.lsi, assay = 'ATAC')  # Identify nearest neighbors based on LSI
  object <- FindClusters(object, graph.name = 'ATAC_snn', algorithm = 3, resolution = 0.2)  # Cluster cells based on nearest neighbors
  
  ##----- Multi-Modal Integration and Clustering -----
  object <- FindMultiModalNeighbors(object,
                                    reduction.list = list("pca", "lsi"),
                                    dims.list = list(1:n.pc, 2:n.lsi),
                                    modality.weight.name = 'RNA.weight')  # Find nearest neighbors across both RNA and ATAC modalities
  object <- RunUMAP(object, nn.name = "weighted.nn", 
                    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")  # Run UMAP using the weighted nearest neighbor graph
  object <- FindClusters(object, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.2)  # Cluster cells using the WNN graph
  
  ##----- Optional: Multi-Modal Integration with Harmony (Batch Correction) -----
  # Harmony can be applied to integrate data from different batches or conditions
  object <- FindMultiModalNeighbors(object,
                                    reduction.list = list("harmony.pca", "harmony.lsi"),
                                    dims.list = list(1:n.pc, 2:n.lsi),
                                    modality.weight.name = 'RNA.weight.harmony',
                                    knn.graph.name = 'wknn.harmony',
                                    snn.graph.name = 'wsnn.harmony',
                                    weighted.nn.name = 'weighted.nn.harmony')  # Repeat WNN analysis with Harmony-corrected data
  
  DefaultAssay(object) <- 'SCT'  # Set default assay back to RNA for further RNA-based analysis
  
  return(object)  # Return the processed and clustered Seurat object
}

# Process and cluster the filtered Seurat object
obj_filt <- process_and_cluster(obj_filt, n.pc = 30, n.lsi = 10)

# Save the processed and clustered Seurat object for downstream analysis
saveRDS(obj_filt, file = './data_processed/OBJ_FILT_clustered.rds')


# -----------------------------------------------------------------------------
# END OF SCRIPT
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Summary of the Code
# -----------------------------------------------------------------------------
# 1. **Load Necessary Libraries**:
#    - Load R libraries essential for multi-modal single-cell analysis, including Seurat for data processing and clustering,
#      Signac for ATAC-seq data, Matrix for handling matrix operations, and future for parallel processing.

# 2. **Load Data**:
#    - Load a preprocessed Seurat object containing multi-modal data (both RNA-seq and ATAC-seq) from a specified file path.

# 3. **Filtering Cells**:
#    - **filter_cells**: Define and apply a function to filter cells based on specified quality control metrics:
#      - RNA and ATAC counts within specified ranges.
#      - Maximum percentage of mitochondrial gene expression.
#      - Maximum nucleosome signal and minimum TSS enrichment.
#    - Filter the Seurat object to retain high-quality cells.

# 4. **Processing and Clustering**:
#    - **process_and_cluster**: Define and apply a function to perform the following steps:
#      - **RNA-seq Processing**:
#        - Normalize RNA data using SCTransform.
#        - Perform PCA and UMAP for dimensionality reduction and visualization.
#        - Identify nearest neighbors and cluster cells.
#      - **ATAC-seq Processing**:
#        - Normalize ATAC data with TF-IDF.
#        - Select top features and perform SVD.
#        - Apply UMAP and clustering.
#      - **Multi-Modal Integration**:
#        - Integrate RNA and ATAC data using weighted nearest neighbor (WNN) analysis.
#        - Optionally, apply Harmony for batch correction and re-run WNN analysis.
#    - Process and cluster the filtered Seurat object based on the defined parameters.

# 5. **Save Processed Object**:
#    - Save the final processed and clustered Seurat object to a specified file path for further downstream analysis.

# -----------------------------------------------------------------------------
