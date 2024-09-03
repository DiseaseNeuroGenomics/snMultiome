# ------------------------------------------------------------------------------
# Project: Single-cell Multi-Modal Analysis - filtering and clustering
# Date: August 2024
# Script name: VarPart_crumblr.R
# Description: VariancePartion analysis and crumblr
#
# outputs: Potential plots for Fig. 1e, Fig.5 b
# ------------------------------------------------------------------------------

# -------------------------------
# Setup and Library Loading
# -------------------------------

# Set library paths (optional)
.libPaths()

# Load Required Libraries

## Dreamlet Universe: For multi-omics analysis
library(dreamlet)      # Dreamlet package for mixed models
library(crumblr)       # Crumblr for correlation modeling
library(zenith)        # Zenith for data integration

## Data Input/Output: For handling single-cell data
library(SingleCellExperiment)  # SingleCellExperiment object
library(zellkonverter)         # Conversion between SingleCellExperiment and other formats
library(tidyr)                 # Data tidying

## Plotting: For visualization
library(ggplot2)       # Data visualization
library(dplyr)         # Data manipulation
library(RColorBrewer)  # Color palettes
library(cowplot)       # Plot composition
library(ggtree)        # Tree visualization
library(aplot)         # Advanced plotting

## Metadata Analysis: For meta-analysis
library(muscat)        # Multi-level analysis of single-cell data
library(metafor)       # Meta-analysis
library(broom)         # Convert statistical analysis objects into tidy tibbles
library(tidyverse)     # Collection of R packages for data science

## Gene Set Enrichment Analysis (GSEA): For pathway analysis
library(topGO)         # Gene ontology enrichment
library(enrichR)       # Enrichr API

# Initialize Enrichr settings
listEnrichrSites()             # List available Enrichr sites
setEnrichrSite("Enrichr")      # Set the Enrichr site to Human genes
websiteLive <- TRUE             # Flag to check if Enrichr website is accessible
dbs <- listEnrichrDbs()         # List available databases in Enrichr

# Check if Enrichr website is accessible
if (is.null(dbs)) {
  websiteLive <- FALSE
}

# Display available databases if website is live
if (websiteLive) {
  head(dbs)
}

# Additional libraries for enrichment analysis
library(clusterProfiler) # ClusterProfiler for enrichment analysis
library(enrichplot)      # Enrichment plot visualization

# -------------------------------
# Function Definitions
# -------------------------------

#' Perform Meta-Analysis on a List of Data Tables
#'
#' This function performs meta-analysis on a list of data tables, computing effect sizes
#' and adjusting p-values for multiple testing.
#'
#' @param tabList A list of data frames, each containing gene-level statistics.
#' @return A tibble with meta-analysis results, including effect sizes and FDR-adjusted p-values.
meta_analysis <- function(tabList) {
  
  # Assign default names if none are provided
  if (is.null(names(tabList))) {
    names(tabList) <- as.character(seq(length(tabList)))
  }
  
  # Add a 'Dataset' column to each data frame in the list
  for (key in names(tabList)) {
    tabList[[key]]$Dataset <- key
  }
  
  # Combine all data frames into one
  df <- do.call(rbind, tabList) 
  
  # Perform meta-analysis for each assay
  results <- df %>%
    as_tibble() %>%
    group_by(assay) %>%
    do(tidy(rma(yi = logFC, sei = logFC / t, data = ., method = "FE"))) %>%
    select(-term, -type) %>%
    ungroup() %>%
    mutate(FDR = p.adjust(p.value, method = "fdr")) %>% 
    mutate(log10FDR = -log10(FDR))
  
  return(results)
}

#' Plot Phylogenetic Tree
#'
#' This function plots a phylogenetic tree using ggtree with customizable color scales.
#'
#' @param tree A phylogenetic tree object.
#' @param low Color for low values (default: "grey90").
#' @param mid Color for mid values (default: "red").
#' @param high Color for high values (default: "darkred").
#' @param xmax.scale Scale factor to increase the x-axis limit (default: 1.5).
#' @return A ggplot object representing the tree.
plotTree <- function(tree, low = "grey90", mid = "red", high = "darkred", xmax.scale = 1.5) {
  
  # Create base tree plot
  fig <- ggtree(tree, branch.length = "none") + 
    geom_tiplab(color = "black", size = 4, hjust = 0, offset = 0.2) +
    theme(legend.position = "top left", plot.title = element_text(hjust = 0.5))
  
  # Get current maximum x-axis value
  xmax <- layer_scales(fig)$x$range$range[2]
  
  # Adjust x-axis limits
  fig <- fig + xlim(0, xmax * xmax.scale) 
  
  return(fig)
}

#' Plot Effect Sizes
#'
#' This function creates a plot of effect sizes with error bars and significance indicators.
#'
#' @param tab Data frame containing the results (must include 'estimate', 'std.error', 'FDR', and 'assay' columns).
#' @param coef The coefficient name used in the analysis.
#' @param fig.tree A tree object used for ordering the cell types.
#' @param low Color for low FDR (default: "grey90").
#' @param mid Color for mid FDR (default: "red").
#' @param high Color for high FDR (default: "darkred").
#' @param ylab Label for the y-axis.
#' @param tick_size Size of the axis tick labels.
#' @return A ggplot object representing the effect sizes.
plotCoef2 <- function(tab, coef, fig.tree, low = "grey90", mid = "red", high = "darkred", ylab, tick_size) {
  
  # Prepare data
  tab$logFC <- tab$estimate
  tab$celltype <- factor(tab$assay, levels = rev(ggtree::get_taxa_name(fig.tree)))
  tab$se <- tab$std.error
  
  # Create effect size plot
  fig.es <- ggplot(tab, aes(x = celltype, y = logFC)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +
    geom_errorbar(aes(ymin = logFC - 1.96 * se, ymax = logFC + 1.96 * se), width = 0) +
    geom_point2(aes(color = pmin(4, -log10(FDR)), size = pmin(4, -log10(FDR)))) + 
    scale_color_gradient2(
      name = expression(-log[10]~FDR),
      limits = c(0, 4),
      low = low,
      mid = mid,
      high = high,
      midpoint = -log10(0.01)
    ) +
    scale_size_area(name = expression(-log[10]~FDR), limits = c(0, 4)) +
    geom_text2(aes(label = '+', subset = FDR < 0.05), color = "white", size = 6, vjust = 0.3, hjust = 0.5) +
    theme_classic() +
    coord_flip() +
    xlab('') + 
    ylab(ylab) +
    theme(
      axis.text.y = element_blank(),
      axis.text = element_text(size = 12),
      axis.ticks.y = element_blank(),
      text = element_text(size = tick_size)
    ) +
    scale_y_continuous(breaks = scales::breaks_pretty(3))
  
  return(fig.es)
}

#' Select Top Genes for Enrichment Analysis
#'
#' This function selects the top N genes for each covariate based on the specified criteria.
#'
#' @param VP.LST List containing variance partition data.
#' @param cova Vector of covariate names.
#' @param topN_genes Number of top genes to select for each covariate.
#' @param analysis Description of the analysis being performed.
#' @return A data frame with selected genes and metadata.
SelectGenes <- function(VP.LST, cova, topN_genes, analysis) {
  
  DF <- data.frame()
  
  for (i in cova) {
    print(i)
    
    # Select top N genes based on the covariate
    top_genes <- VP.LST %>%
      as_tibble() %>%
      arrange(desc(!!sym(i))) %>%
      slice_head(n = topN_genes)
    
    Gene <- top_genes$gene
    
    # Create a data frame for the selected genes
    df <- data.frame(
      Gene = Gene,
      covariate = i,
      analysis = analysis,
      topN_genes = topN_genes,
      stringsAsFactors = FALSE
    )
    
    DF <- rbind(DF, df)
  }
  
  # Note: The following line seems to overwrite 'topN_genes' with 'gene', which might be unintended.
  # Verify if this is desired or if it should be removed/commented out.
  # DF$topN_genes <- DF$gene
  
  return(DF)
}

#' Retrieve Expression Values for a Gene
#'
#' This function extracts expression values for a specified gene from a data frame.
#'
#' @param dataframe Data frame containing gene expression data.
#' @param Gene_name Name of the gene to retrieve.
#' @return A list of expression values for the specified gene.
ListExpression <- function(dataframe, Gene_name) {
  
  # Subset the data frame for the specified gene
  test <- subset(dataframe, subset = Gene == Gene_name)
  
  # Remove the last column if it's non-expression data (adjust as needed)
  test.df <- test[head(seq_len(ncol(test)), -1)]
  
  # Transpose and convert to a list
  test.list <- as.list(as.data.frame(t(test.df)))
  
  return(test.list)
}

#' Prepare Data for Dotplot Visualization
#'
#' This function filters and selects genes for dotplot visualization based on specified criteria.
#'
#' @param Df Data frame containing enrichment results.
#' @param Database Name of the enrichment database to filter.
#' @param numGenes Number of genes to select.
#' @param numGenesCutoff Maximum number of genes to include after cutoff.
#' @return A filtered data frame ready for dotplot visualization.
PrepforDotplot <- function(Df, Database, numGenes, numGenesCutoff) {
  
  # Filter data for the specified database
  Df_bio_proc <- Df[which(Df$database == Database), ]
  
  # Count the number of genes per row (assuming 'Genes' are separated by ';')
  CNT <- sapply(Df_bio_proc$Genes, function(gens) length(strsplit(gens, split = ";")[[1]]))
  Df_bio_proc$Gene_count <- CNT
  
  D0 <- data.frame()
  
  # Select top genes for each 'Top_DE'
  for (i in unique(Df_bio_proc$Top_DE)) {
    d <- Df_bio_proc[which(Df_bio_proc$Top_DE == i), ]
    d0 <- d %>% top_n(n = -numGenes, wt = Adjusted.P.value)
    print(paste("Number of genes after initial filtering for", i, ":", nrow(d0)))
    
    if (nrow(d0) > numGenesCutoff) {
      d0 <- head(d0, numGenesCutoff)
    }
    
    D0 <- rbind(D0, d0)
  }
  
  return(D0)
}

#' Run Gene Ontology (GO) Enrichment Analysis
#'
#' This function performs GO enrichment analysis for each covariate and visualizes the results.
#'
#' @param df_corr Data frame containing correlated genes and covariates.
#' @return A combined data frame with enrichment results.
runGOEnrichment <- function(df_corr) {
  
  Df <- data.frame()
  
  # Iterate over each covariate
  for (i in unique(df_corr$covariate)) {
    temp <- df_corr[df_corr$covariate == i, ]
    Cluster_DE <- i
    print(paste("Performing enrichment for:", Cluster_DE))
    GENES <- temp$Gene
    
    # Define enrichment databases
    dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021",
             "WikiPathway_2021_Human", "Elsevier_Pathway_Collection", "Reactome_2022")
    
    # Perform enrichment if the website is accessible
    if (websiteLive) {
      enriched <- enrichr(GENES, dbs)
    } else {
      warning("Enrichr website is not accessible. Skipping enrichment analysis.")
      next
    }
    
    # Process each enrichment result
    for (j in seq_along(enriched)) {
      enriched[[j]]$Top_DE <- Cluster_DE
      enriched[[j]]$database <- dbs[j]
      
      # Extract GO terms and numbers
      TEXT <- sapply(enriched[[j]]$Term, function(term) strsplit(term, split = "\\(")[[1]][1])
      GOT <- sapply(enriched[[j]]$Term, function(term) {
        parts <- strsplit(term, split = "\\(")[[1]]
        if (length(parts) > 1) {
          strsplit(parts[2], split = "\\)")[[1]][1]
        } else {
          NA
        }
      })
      
      enriched[[j]]$GO_text <- TEXT
      enriched[[j]]$GO_number <- GOT
    }
    
    # Combine all enriched data
    Df <- rbind(Df, do.call(rbind, enriched))
    
    # Generate enrichment plots
    if (websiteLive) {
      p0 <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, 
                       y = "Count", orderBy = "Adjusted.P.value")
      p1 <- plotEnrich(enriched[[2]], showTerms = 20, numChar = 50, 
                       y = "Count", orderBy = "Adjusted.P.value")
      p3 <- plotEnrich(enriched[[3]], showTerms = 20, numChar = 50, 
                       y = "Count", orderBy = "Adjusted.P.value")
      p4 <- plotEnrich(enriched[[4]], showTerms = 20, numChar = 50, 
                       y = "Count", orderBy = "Adjusted.P.value")
      p5 <- plotEnrich(enriched[[5]], showTerms = 20, numChar = 50, 
                       y = "Count", orderBy = "Adjusted.P.value")
      
      # Set plot dimensions
      options(repr.plot.width = 20, repr.plot.height = 8)
      
      # Combine and print plots
      print(p0 + p1)
      print(p3 + p4)
      print(p1 + p5)
    }
  }
  
  return(Df)
}



# --------------------------------------------
# Main Analysis Workflow for VariancePartition
# --------------------------------------------

# --- 1) Prepare data and convert to pseudobulk

# Load data; from Seurat object - DietSeurat() to single assay and then save as h5ad


# Create pseudobulk data by specifying cluster_id and sample_id
# Count data for each cell type is then stored in the `assay` field
# assay: entry in assayNames(sce) storing raw counts
# cluster_id: variable in colData(sce) indicating cell clusters
# sample_id: variable in colData(sce) indicating sample id for aggregating cells
pb <- aggregateToPseudoBulk(sce,
                            assay = "counts",
                            cluster_id = "cell",
                            sample_id = "id",
                            verbose = FALSE
)

# one 'assay' per cell type
assayNames(pb)

# --- 2) Perform preprocessing

# Define the model formula with fixed and random effects
form <- ~ (Age) + (1 | Brain.region) + (1 | Sex) + (1 | Brain.bank) + (1 | Individual)

# Normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, form, min.count = 5)

# the resulting object of class dreamletProcessedData stores
# normalized data and other information
res.proc

# view details of dropping samples
details(res.proc)

# show voom plot for each cell clusters
plotVoom(res.proc)

#----- 3) Variance Partition analysis

# run variance partitioning analysis
vp.lst <- fitVarPart(res.proc, form)


# --- 4) Generate plots from Variance Partition

# Set plot dimensions for variance partition plot
options(repr.plot.width = 10, repr.plot.height = 10)

# Plot variance partition results
plotVarPart(sortCols(vp.lst), label.angle = 60, ncol = 4) 






# --------------------------------------------
# Main Analysis Workflow for crumblr
# --------------------------------------------

# --- 1) Prepare data for pseudobulk

pb_celltype <- pb #use pseudobulk from VariancePartition analysis 

# --- 2) Crumblr on Pseudobulk Data

# Create crumblr object from cell counts
cobj <- crumblr(cellCounts(pb_celltype))

# Define the model formula with fixed and random effects
form <- ~ (Age) + (1 | Brain.region) + (1 | Sex) + (1 | Brain.bank) + (1 | Individual)

# Fit variance partition model and extract results
vp.c <- fitExtractVarPartModel(cobj, form, colData(pb_celltype))

# Set plot dimensions for variance partition plot
options(repr.plot.width = 6, repr.plot.height = 5)

# Plot variance partition results for crumblr
plotVarPart(sortCols(vp.c), label.angle = 60, ncol = 4) 

# Create percent bars plot for variance partition
fig.vp <- plotPercentBars(sortCols(vp.c))
fig.vp


# --- 3) Analysis with dream()

# Define the model formula (redundant; already defined earlier)
# form <- ~ (Age) + (1 | Brain.region) + (1 | Sex) + (1 | Brain.bank) + (1 | Individual)

# Fit dream model using the crumblr object and the formula
fit <- dream(cobj, form, colData(pb_celltype))

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# --- 4) Meta-Analysis

# Define prefix for output files
prefix <- 'Age_analysis_'

# Extract top table for the 'Age' coefficient
res.age <- topTable(fit, coef = 'Age', number = Inf) %>%
  rownames_to_column('assay')

# Perform meta-analysis on the 'Age' results
res.age <- meta_analysis(list(res.age))
head(res.age)

# View unique coefficients (for verification)
unique(fit$coefficients)

# Build cluster tree from pseudobulk data
hc <- buildClusterTreeFromPB(pb_celltype)

# Plot the phylogenetic tree with increased x-axis scale and adjust legend position
fig.tree <- plotTree(ape::as.phylo(hc), xmax.scale = 2.2) + theme(legend.position = "bottom")

# --- 4) Effect Size Plotting

# Create effect size plot for 'Age' meta-analysis results
fig.es1 <- plotCoef2(res.age, coef = 'Age', fig.tree, ylab = 'Age', tick_size = 12)

# --- 5) Combine Plots

# Set plot dimensions
options(repr.plot.width = 8, repr.plot.height = 5)

# Combine effect size plot and tree plot
combined_plot <- fig.es1 %>% insert_left(fig.tree, width = 1.3) 

# Save the combined plot as a PDF
ggsave(filename = paste0(prefix, '_controls_age_sex.pdf'), plot = combined_plot, width = 13, height = 5)




# --------------------------------------------------------------------
# Main Analysis Workflow for Gene Set Enrichment Analysis (GSEA)
# --------------------------------------------------------------------

# The following steps for GSEA are outlined via function definitions above.
# To execute GSEA, you would typically:

# 1. Select top genes using `SelectGenes()`
# 2. Retrieve their expression values using `ListExpression()`
# 3. Prepare data for visualization using `PrepforDotplot()`
# 4. Run enrichment analysis using `runGOEnrichment()`

# Example (assuming you have appropriate data):
# selected_genes <- SelectGenes(VP.LST = vp.c, cova = c("Age"), topN_genes = 100, analysis = "Age Analysis")
# expression_values <- ListExpression(dataframe = your_expression_data, Gene_name = "GeneX")
# dotplot_data <- PrepforDotplot(Df = your_enrichment_results, Database = "GO_Molecular_Function_2021", numGenes = 20, numGenesCutoff = 50)
# enrichment_results <- runGOEnrichment(df_corr = your_correlated_genes)

# Adjust and execute these steps based on your specific data and analysis needs.


# GSEA analysis functions

SelectGenes <- function(VP.LST, cova, topN_genes, analysis) {
  DF <- data.frame()
  for (i in cova) {
    print(i)
    VP.LST %>%
      as_tibble() %>%
      arrange(desc(!!sym(i))) %>%
      slice_head(n = topN_genes) -> top100
    
    Gene <- top100$gene
    df <- as.data.frame(Gene)
    df$covariate <- i
    df$analysis <- analysis
    df$topN_genes <- topN_genes
    DF <- rbind(DF, df)
  }
  
  DF$topN_genes <- DF$gene
  return(DF)
}

ListExpression <- function(dataframe, Gene_name) {
  test <- subset(dataframe, subset = Gene == Gene_name)
  test.df <- test[head(seq_len(ncol(test)), -1)]
  test.list <- as.list(as.data.frame(t(test.df)))
  return(test.list)
}

PrepforDotplot <- function(Df, Database, numGenes, numGenesCutoff) {
  Df_bio_proc <- Df[which(Df$database == Database),]
  CNT <- c()
  for (i in 1:nrow(Df_bio_proc)) {
    gens <- Df_bio_proc[i,]$Genes
    CNT <- append(CNT, count.fields(textConnection(gens), sep = ";"))
  }
  
  Df_bio_proc$Gene_count <- CNT
  D <- Df_bio_proc
  D0 <- data.frame()
  for (i in unique(D$Top_DE)) {
    d <- D[which(D$Top_DE == i),]
    d %>% top_n(n = -numGenes, wt = Adjusted.P.value) -> d0
    if (nrow(d0) > numGenes) {
      d0 <- head(d0, numGenesCutoff)
    }
    D0 <- rbind(D0, d0)
  }
  return(D0)
}

runGOEnrichment <- function(df_corr) {
  Df <- data.frame()
  
  for (i in unique(df_corr$covariate)) {
    temp <- df_corr[which(df_corr$covariate == i),]
    Cluster_DE <- i
    print(Cluster_DE)
    GENES <- temp$Gene
    
    # EnrichR and GO
    dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021",
             "WikiPathway_2021_Human", "Elsevier_Pathway_Collection", "Reactome_2022")
    if (websiteLive) {
      enriched <- enrichr(GENES, dbs)
    }
    
    enriched[[1]]$Top_DE <- Cluster_DE
    enriched[[1]]$database <- 'GO_Molecular_Function_2021'
    enriched[[1]]$GO_text <- sapply(enriched[[1]]$Term, function(x) strsplit(x, split = "\\(")[[1]][1])
    enriched[[1]]$GO_number <- sapply(enriched[[1]]$Term, function(x) strsplit(strsplit(x, split = "\\(")[[1]][2], split = "\\)")[[1]][1])
    
    # Similar processing for other databases
    for (db_index in 2:5) {
      enriched[[db_index]]$Top_DE <- Cluster_DE
      enriched[[db_index]]$database <- names(dbs)[db_index]
      enriched[[db_index]]$GO_text <- sapply(enriched[[db_index]]$Term, function(x) strsplit(x, split = "\\(")[[1]][1])
      enriched[[db_index]]$GO_number <- sapply(enriched[[db_index]]$Term, function(x) strsplit(strsplit(x, split = "\\(")[[1]][2], split = "\\)")[[1]][1])
    }
    
    Df <- rbind(Df, do.call(rbind, enriched))
    
    if (websiteLive) {
      p0 <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Adjusted.P.value")
      p1 <- plotEnrich(enriched[[2]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Adjusted.P.value")
      p3 <- plotEnrich(enriched[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Adjusted.P.value")
      p4 <- plotEnrich(enriched[[4]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Adjusted.P.value")
      p5 <- plotEnrich(enriched[[5]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Adjusted.P.value")
      options(repr.plot.width = 20, repr.plot.height = 8)
      
      print(p0 + p1)
      print(p3 + p4)
      print(p1 + p5)
    }
  }
  
  return(Df)
}


# -------------------------------
# End of Script
# -------------------------------




# -----------------------------------------------------------------------------
# Summary of the Code
# -----------------------------------------------------------------------------

# **Libraries**:
# - Various R libraries for data processing, meta-analysis, and plotting are loaded,
#   including 'dreamlet', 'crumblr', and 'zenith' for specialized analysis, and 'ggplot2', 'dplyr', 'ggtree' for visualization.

# **Functions**:
# - `meta_analysis()`: Performs meta-analysis on provided data, adjusting p-values and calculating FDR.
# - `plotTree()`: Visualizes a tree structure with customizable color gradients and x-axis scaling.
# - `plotCoef2()`: Creates a plot of effect sizes with error bars and color-coded significance.

# **Analysis Workflow**:
# Variance Partition Analysis
# 1. Prepare data and convert to pseudobulk format
# 2. Normalize data and apply voom or voomWithDreamWeights
# 3. Fit the variance partition model and visualize results

# Crumblr Analysis
# 1. Perform crumblr analysis on data
# 2. Interpret and visualize results
# 3. *Dream Analysis*: Fit linear models and perform eBayes analysis.
# 4. *Meta-Analysis*: Conduct meta-analysis on results and plot effect sizes.
# 5. *GSEA Analysis*: Perform gene set enrichment analysis using various databases, and plot results if Enrichr is available.

# **Output**:
# - Various plots and saved results from meta-analysis and GSEA analysis.

# -----------------------------------------------------------------------------








