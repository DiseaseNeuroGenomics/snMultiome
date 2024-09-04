# Project Overview - snMultiome
This repository contains scripts for the snMultiome postnatal project, focused on the analysis and integration of multi-modal single-cell data, specifically RNA-seq and ATAC-seq. Below is an overview of the scripts used in the project. Further information, including relevant input data files, can be found in synapse under SynID: **syn62738189**.


## SnMultiome Analysis - Preprocessing
- **Date:** August 2024
- **Script Name:** `preprocess_merge.R`
- **Description:** Pre-processes single-cell RNA-seq and ATAC-seq data from multiple samples using Seurat and Signac, merges data into a single Seurat object, and saves the output.

## SnMultiome Analysis - Filtering and Clustering
- **Date:** August 2024
- **Script Name:** `filter_cluster.R`
- **Description:** Loads a pre-processed Seurat object with multi-modal data (RNA and ATAC-seq), performs filtering, normalization, dimensional reduction, and clustering, then saves the processed object.
- **Outputs:** Potential plots for Fig. 1b, e.

## Variance Partition and crumblr Analysis
- **Date:** August 2024
- **Script Name:** `VarPart_crumblr.R`
- **Description:** Performs VariancePartition analysis and crumblr.
- **Outputs:** Potential plots for Fig. 1e, Fig. 5b.


## MAGMA and LDsc Analysis
- **Date:** August 2024
- **Script Name:** `magma_ldsc.R`
- **Description:** Analyzes enrichment (MAGMA / LDsc) and TOBIAS results.
- **Outputs:** Plots for Fig. 2c, 2e, 3a, 6c, and S6.

## ABC-Max Analysis
- **Date:** August 2024
- **Script Name:** `abcmax.R`
- **Description:** Processes genomic data for various GWAS traits and linkage disequilibrium (LD) buddies, focusing on ABC-MAX analysis (peak overlap with GWAS SNPs and causal genes).
- **Outputs:** Plots for Fig. 3-4, Fig. S8-S17, Table S8-S10.

## Pseudotime analysis using Monocle3
- **Date:** August 2024
- **Script Name:** `pseudotime.R`
- **Description:** Processes a Seurat object, converts it to a Monocle3 `cell_data_set`, performs trajectory inference, and identifies temporal genes using Moran's I test.
- **Outputs:** Potential plots for Fig. 5a, c-d.

## Hotspot assignment and SCENIC+ correlation with pseudotime
- **Date:** August 2024
- **Script Name:** `OPCOligo_Hotspot_SCENIC.R`
- **Description:** Performs module scoring and downstream analysis for Hotspot results and pseudotime, including correlation analysis between pseudotime and SCENIC results.
- **Outputs:** Potential plots for Fig. 5e-f, Fig. 6 a-c.
