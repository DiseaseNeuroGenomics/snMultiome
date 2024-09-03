# -----------------------------------------------------------------------------
# Title: Downstream analysis from oligodendrogenesis pseudotime
#         correlated with hotspot and SCENIC+ results
# Date: August 2024
# Script name: OPCOligo_Hotspot_SCENIC.R
# Description: This script performs module scoring and downstream analysis for
#             Hotspot results and pseudotime. The script also includes 
#             correlation analysis between pseudotime and eRegulon activity
#             using SCENIC results.
# -----------------------------------------------------------------------------

# Load necessary libraries
library(Seurat)    # For single-cell RNA-seq data analysis
library(Signac)    # For single-cell ATAC-seq data analysis
library(ggplot2)   # For data visualization
library(cowplot)   # For combining ggplot objects
library(reshape2)  # For data reshaping

# Load preprocessed Seurat object
opc_olig3 <- readRDS('OPCOligObject.rds')  # Replace with the correct path if needed

# Load external cell annotations
annot_file <- '3DG_OPCOlig_annotation_Sep2023.csv'
annot.use <- read.csv(annot_file, row.names = 1)
rownames(annot.use) <- annot.use$Colnames

# Add annotations to the Seurat object
opc_olig3 <- AddMetaData(opc_olig3, metadata = annot.use[Cells(opc_olig3), ])

# Load gene modules from Hotspot analysis
hs_file <- 'hs_results_modules.csv'
hs.modules <- read.csv(hs_file, row.names = 1)

# Convert modules into a list excluding non-informative module '-1'
hs.list <- lapply(setdiff(sort(unique(hs.modules$Module)), -1), function(x) {
  rownames(hs.modules)[hs.modules$Module == x]
})

# Add module scores to Seurat object
opc_olig3 <- AddModuleScore(opc_olig3, features = hs.list)

# Rename module score columns for clarity
colnames(opc_olig3@meta.data)[grep('^Cluster', colnames(opc_olig3@meta.data))] <- paste0('HS_', 1:length(hs.list))

# Merge gene modules with visually similar patterns into combined modules
hs.new <- list(c(1, 3, 8), c(2, 6), c(4, 5, 7, 11, 12), c(9, 10))
hs.list.2 <- lapply(hs.new, function(x) {
  rownames(hs.modules)[hs.modules$Module %in% x]
})

# Recalculate module scores for the newly combined modules
opc_olig3 <- AddModuleScore(opc_olig3, features = hs.list.2)

# Rename the new module score columns for clarity
colnames(opc_olig3@meta.data)[grep('^Cluster', colnames(opc_olig3@meta.data))] <- paste0('HS2_', 1:length(hs.list.2))

# Set cell identities based on predefined annotations
opc_olig3$fgroup <- factor(opc_olig3$anno_clus_dreamorigBB_bioarx_v2,
                           levels = c('TypeI_OPC', 'TypeII_OPC', 'COP', 
                                      'TypeI_Oligo_OPALIN', 'TypeII_Oligo_RBFOX1'))
Idents(opc_olig3) <- 'fgroup'

# Plot scatter plots for each combined gene module score vs. pseudotime
plot.list <- lapply(1:length(hs.list.2), function(x) {
  FeatureScatter(opc_olig3, feature1 = 'pseudotime_sc_OPC1', 
                 feature2 = paste0('HS2_', x), pt.size = 0.6)
})

# Display scatter plots in a grid layout
cowplot::plot_grid(plotlist = plot.list)

# Predict trend lines for hotspot modules over pseudotime using LOESS regression
pred.use <- c()
for (i in 1:length(hs.list)) {
  data.use <- FetchData(opc_olig3, vars = c('pseudotime_sc_OPC1', paste0('HS_', i)))
  colnames(data.use) <- c('pt', 'feature')
  
  # Fit a LOESS model for the trend line
  model.use <- loess(feature ~ pt, data = data.use)
  pred.use <- cbind(pred.use, model.use$fitted)
}

# Convert predictions to a data frame for plotting
pred.use <- as.data.frame(pred.use)
colnames(pred.use) <- paste0('HS_', 1:ncol(pred.use))
pred.use$pt <- data.use$pt
rownames(pred.use) <- Cells(opc_olig3)

# Plot trend lines for each merged gene module
plot.list <- lapply(1:length(hs.new), function(x) {
  plot.data <- reshape2::melt(pred.use[, paste0('HS_', hs.new[[x]])])
  colnames(plot.data) <- c('Module', 'Value')
  plot.data$pt <- rep(pred.use$pt, times = length(hs.new[[x]]))
  
  ggplot(plot.data, aes(x = pt, y = Value, color = Module)) +
    geom_line(size = 1.2) + theme_bw() + ggtitle(paste0('Trend_', x))
})

# Display trend line plots in a grid layout
cowplot::plot_grid(plotlist = plot.list)

# -------- Pseudotime vs eRegulon Activity Correlation Analysis --------

# Load the processed Seurat object and SCENIC results
opc_olig3 <- readRDS('OPCOligObject.rds')
ereg.use <- readRDS('scenic_sc.rds')

# Ensure only shared cells between Seurat and SCENIC results are used
cells.use <- rownames(ereg.use)
opc_olig3 <- subset(opc_olig3, cells = cells.use)

# Add SCENIC results as a new assay
opc_olig3[['scenic']] <- CreateAssayObject(data = t(ereg.use))
DefaultAssay(opc_olig3) <- 'scenic'

# Scale SCENIC assay data
opc_olig3 <- ScaleData(opc_olig3, features = rownames(opc_olig3))

# Randomly sample cells for plotting, stratified by pseudotime bins
set.seed(42)
plot.cells <- lapply(levels(opc_olig3$pt_bin), function(x) {
  sample(Cells(opc_olig3)[opc_olig3$pt_bin == x], size = min(1000, sum(opc_olig3$pt_bin == x)))
})
plot.cells <- unlist(plot.cells)
plot.cells <- plot.cells[order(opc_olig3@meta.data[plot.cells, 'pseudotime_sc_OPC1'])]

# Correlation between pseudotime and eRegulon AUC scores
cor.ereg <- sapply(1:ncol(ereg.use), function(x) 
  cor(opc_olig3$pseudotime_sc_OPC1, ereg.use[cells.use, x])
)
names(cor.ereg) <- colnames(ereg.use)
cor.ereg <- sort(cor.ereg)

# Identify top 10 positively and negatively correlated eRegulons
plot.features <- c(head(names(cor.ereg), 10), tail(names(cor.ereg), 10))
plot.features <- gsub('_', '-', plot.features)

# Plot heatmap of eRegulon activity
DoHeatmap(opc_olig3, features = plot.features, group.by = 'is_cell', 
          cells = plot.cells, label = FALSE)

# Switch to RNA assay and plot TF expression for top correlated eRegulons
DefaultAssay(opc_olig3) <- 'SCT'
plot.features <- unname(sapply(plot.features, function(x) 
  unlist(strsplit(x, '-'))[1]
))
opc_olig3 <- GetResidual(opc_olig3, plot.features)

# Plot heatmap of TF expression
DoHeatmap(opc_olig3, features = plot.features, group.by = 'is_cell', 
          cells = plot.cells, label = FALSE)

# Save the final Seurat object
saveRDS(opc_olig3, file = 'OPC_Olig_final.rds')

# -----------------------------------------------------------------------------
# END OF SCRIPT
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Summary of the Code
# -----------------------------------------------------------------------------
# 1. **Load Necessary Libraries**: 
#    - Load required R libraries for single-cell RNA-seq and ATAC-seq data analysis, visualization, and data manipulation.

# 2. **Load Data**:
#    - Load the preprocessed Seurat object (`OPCOligObject.rds`) containing the OPC and oligodendrocyte data.
#    - Load and integrate external cell annotations (`3DG_OPCOlig_annotation_Sep2023.csv`) into the Seurat object.

# 3. **Process Gene Modules from Hotspot **:
#    - Load gene module information from Hotspot analysis (`hs_results_modules.csv`).
#    - Convert gene modules into a list, excluding non-informative modules.
#    - Calculate module scores for each cell based on these gene modules and add them to the Seurat object.

# 4. **Combine Gene Modules**:
#    - Merge visually similar gene modules into new combined modules.
#    - Recalculate module scores for these newly combined modules and update the Seurat object.

# 5. **Set Cell Identities**:
#    - Set cell identities based on predefined annotations and update the Seurat object's identity classes.

# 6. **Plot Gene Module Scores vs. Pseudotime**:
#    - Create scatter plots to visualize the relationship between pseudotime and each gene module score.
#    - Display these plots in a grid layout.

# 7. **Predict Trend Lines**:
#    - Fit LOESS regression models to predict trend lines for each gene module score over pseudotime.
#    - Plot these trend lines to visualize trends.

# 8. **Pseudotime vs. eRegulon Activity Correlation Analysis**:
#    - Load SCENIC results (`scenic_sc.rds`) and ensure compatibility with the Seurat object.
#    - Add SCENIC results as a new assay in the Seurat object and scale the data.
#    - Randomly sample cells for plotting and calculate correlations between pseudotime and eRegulon AUC scores.
#    - Identify the top 10 positively and negatively correlated eRegulons and plot heatmaps of their activity and transcription factor expression.

# 9. **Save Final Object**:
#    - Save the updated Seurat object with processed data and results to an RDS file (`OPC_Olig_final.rds`).

# -----------------------------------------------------------------------------
