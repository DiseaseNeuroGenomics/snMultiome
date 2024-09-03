# -----------------------------------------------------------------------------
# Project: Single-cell Multi-Modal Analysis
# Date: August 2024
# Script name: preprocess_merge.R
# Description: This script pre-processes single-cell RNA-seq and ATAC-seq 
#              data from multiple samples using Seurat and Signac, merges the 
#          pre-processed data into a single Seurat object, and saves the output.
# -----------------------------------------------------------------------------

#---------------------------------------------------------------------
# LIBRARIES
#---------------------------------------------------------------------

# Load necessary libraries for multi-modal single-cell analysis
library(Seurat)              # Comprehensive toolkit for single-cell RNA-seq analysis
library(Signac)              # Extends Seurat for single-cell chromatin accessibility analysis (e.g., ATAC-seq)
library(ggplot2)             # Visualization library for creating plots
library(harmony)             # Tool for batch correction and data integration across different conditions
library(EnsDb.Hsapiens.v86)  # Ensembl-based annotation package for the human genome (version 86)
library(dplyr)               # Grammar of data manipulation
library(tidyr)               # Tools for tidying up data

#---------------------------------------------------------------------
# HELPER FUNCTIONS
#---------------------------------------------------------------------
### Function: create_seurat_object
# Description: Creates a Seurat object from 10x Genomics data for both RNA and ATAC modalities.
create_seurat_object <- function(sample_id, data_dir, genome_version = 'hg38') {
  setwd(file.path(data_dir, sample_id))
  
  inputdata_10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
  metadata <- read.csv('per_barcode_metrics.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  
  seurat_obj <- CreateSeuratObject(counts = inputdata_10x$`Gene Expression`, 
                                   assay = 'RNA', 
                                   project = sample_id, 
                                   meta.data = metadata)
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  atac_counts <- filter_standard_chromosomes(inputdata_10x$Peaks)
  seurat_obj[['ATAC']] <- create_atac_assay(atac_counts, genome_version)
  seurat_obj <- calculate_atac_qc_metrics(seurat_obj, genome_version)
  
  return(seurat_obj)
}

### Function: filter_standard_chromosomes
# Description: Filters ATAC-seq peaks to include only those from standard chromosomes.
filter_standard_chromosomes <- function(atac_counts) {
  grange_counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  
  return(atac_counts[as.vector(grange_use), ])
}

### Function: create_atac_assay
# Description: Creates an ATAC-seq ChromatinAssay for integration into the Seurat object.
create_atac_assay <- function(atac_counts, genome_version) {
  annotations <- prepare_annotations(genome_version)
  frag_file <- "atac_fragments.tsv.gz"
  
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = genome_version,
    fragments = frag_file,
    min.cells = 10,
    annotation = annotations
  )
  
  return(chrom_assay)
}

### Function: prepare_annotations
# Description: Prepares gene annotations for the specified genome version.
prepare_annotations <- function(genome_version) {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- genome_version
  
  return(annotations)
}

### Function: calculate_atac_qc_metrics
# Description: Calculates quality control metrics for ATAC-seq data.
calculate_atac_qc_metrics <- function(seurat_obj, genome_version) {
  DefaultAssay(seurat_obj) <- 'ATAC'
  
  seurat_obj <- NucleosomeSignal(seurat_obj)
  seurat_obj <- TSSEnrichment(seurat_obj, fast = FALSE)
  seurat_obj$pct_reads_in_peaks <- seurat_obj$atac_peak_region_fragments / seurat_obj$atac_fragments * 100
  seurat_obj$blacklist_fraction <- FractionCountsInRegion(seurat_obj, assay = 'ATAC', regions = blacklist_hg38_unified)
  
  peaks <- call_and_filter_peaks(seurat_obj, genome_version)
  macs_count <- FeatureMatrix(fragments = Fragments(seurat_obj),
                              features = peaks,
                              cells = colnames(seurat_obj))
  
  seurat_obj[['ATAC']] <- CreateChromatinAssay(
    counts = macs_count,
    sep = c(":", "-"),
    genome = genome_version,
    fragments = "atac_fragments.tsv.gz",
    min.cells = 10,
    annotation = prepare_annotations(genome_version)
  )
  
  return(seurat_obj)
}

### Function: call_and_filter_peaks
# Description: Calls ATAC-seq peaks using MACS2 and filters out blacklist regions.
call_and_filter_peaks <- function(seurat_obj, genome_version) {
  peaks <- CallPeaks(seurat_obj, assay = 'ATAC', macs2.path = '/home/xxx/anaconda3/bin/macs2')
  peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  return(peaks)
}

#---------------------------------------------------------------------
# MAIN SCRIPT
#---------------------------------------------------------------------

# Define the root directory where the data for all samples is stored
data_dir <- "../data"

# Create and preprocess Seurat objects for each sample
obj_example_1 <- create_seurat_object(sample_id = "example_1", data_dir = data_dir)
save(obj_example_1, file = file.path(data_dir, "obj_example_1_processed.rda"))

obj_example_2 <- create_seurat_object(sample_id = "example_2", data_dir = data_dir)
obj_example_3 <- create_seurat_object(sample_id = "example_3", data_dir = data_dir)
obj_example_4 <- create_seurat_object(sample_id = "example_4", data_dir = data_dir)

# Merge the Seurat objects from all samples into a single Seurat object
obj_all <- merge(obj_example_1, y = c(obj_example_2, obj_example_3, obj_example_4), 
                 project = "postnatal_brain")

# Save the merged Seurat object for downstream analysis
save(obj_all, file = file.path(data_dir, "data_preprocessed/ALL_snMultiome.rda"))
saveRDS(obj_all, file.path(data_dir, "data_processed/ALL_snMultiome.rds"))

# -----------------------------------------------------------------------------
# END OF SCRIPT
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Summary of the Code
# -----------------------------------------------------------------------------
# 1. **Load Necessary Libraries**:
#    - Load R libraries required for single-cell multi-modal analysis, including Seurat for single-cell analysis,
#      Signac for ATAC-seq data, ggplot2 for visualization, harmony for batch correction, EnsDb.Hsapiens.v86 for
#      genome annotations, dplyr and tidyr for data manipulation.

# 2. **Define Helper Functions**:
#    - **create_seurat_object**: Creates a Seurat object from 10x Genomics data for both RNA and ATAC modalities,
#      including reading data, creating assays, and calculating quality control metrics.
#    - **filter_standard_chromosomes**: Filters ATAC-seq peaks to include only those from standard chromosomes.
#    - **create_atac_assay**: Creates an ATAC-seq ChromatinAssay for integration into the Seurat object.
#    - **prepare_annotations**: Prepares gene annotations for the specified genome version.
#    - **calculate_atac_qc_metrics**: Calculates quality control metrics for ATAC-seq data, including nucleosome signal,
#      TSS enrichment, and fraction of reads in peaks.
#    - **call_and_filter_peaks**: Calls ATAC-seq peaks using MACS2 and filters out blacklist regions.

# 3. **Main Script**:
#    - Define the root directory where the data for all samples is stored.
#    - Create and preprocess Seurat objects for each sample using the `create_seurat_object` function.
#    - Save each preprocessed Seurat object to file for later use.
#    - Merge all Seurat objects from different samples into a single Seurat object to combine data from multiple samples.
#    - Save the merged Seurat object for downstream analysis in both `.rda` and `.rds` formats.

# -----------------------------------------------------------------------------
