# -----------------------------------------------------------------------------
# Project: Analyzing Single-Cell RNA-seq Data using Seurat and Monocle3
# Date: August 2024
# Script name: pseudotime.R
# Description: This script processes a Seurat object, converts it to a Monocle3
#             cell_data_set, performs trajectory inference, and identifies temporal
#             genes using Moran's I test.
# -----------------------------------------------------------------------------


# Load necessary libraries
library(Seurat) # for single-cell RNA-seq analysis
library(SeuratWrappers) # for interfacing Seurat with other tools like Monocle3
library(monocle3) # for trajectory inference in single-cell data

#---------------------------------------------------------------------
# 1. Load and Prepare Seurat Object
#---------------------------------------------------------------------

# Load the processed Seurat object from a .rds file
opc_olig3 <- readRDS('OPCOlig_allprocessed.rds')

# Create a new Seurat object 'o0' specifically for Monocle3 processing, 
# using the counts data and metadata from the loaded Seurat object
o0 <- CreateSeuratObject(counts = opc_olig3[['SCT']]@counts, meta.data = opc_olig3@meta.data)

# Set the 'data' slot of 'o0' with the normalized data from the original Seurat object
o0 <- SetAssayData(o0, slot = 'data', new.data = opc_olig3[['SCT']]@data)

# Set the 'scale.data' slot of 'o0' with the scaled data from the original Seurat object
o0 <- SetAssayData(o0, slot = 'scale.data', new.data = opc_olig3[['SCT']]@scale.data)

# Add UMAP embeddings to the 'o0' object, using the embeddings from the original Seurat object
# and assigning a key prefix 'UMAP_' to the dimension reduction object
o0[['umap']] <- CreateDimReducObject(embeddings = Embeddings(opc_olig3, 'wnn.dream2BB.umap'), key = 'UMAP_')

# Convert the Seurat object 'o0' to a Monocle3 cell_data_set object 'cds' for trajectory analysis
cds <- as.cell_data_set(o0)

#---------------------------------------------------------------------
# 2. Trajectory Inference with Monocle3
#---------------------------------------------------------------------

# Preprocess the 'cds' object by normalizing and reducing dimensionality
cds <- preprocess_cds(cds)

# Cluster the cells in the 'cds' object based on the reduced dimensions
cds <- cluster_cells(cds)

# Learn the trajectory graph in the 'cds' object without partitioning cells 
# and without closing loops in the graph
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# Plot the inferred trajectory, coloring the cells by their original cluster annotations
# Branch points, leaves, and roots are not labeled; group labels are also omitted
plot_cells(cds, color_cells_by = 'anno_clus_dreamorigBB',
           label_branch_points = FALSE, label_leaves = FALSE, label_roots = FALSE,
           group_label_size = 4, label_groups_by_cluster = FALSE)

#---------------------------------------------------------------------
# 3. Pseudotime Calculation
#---------------------------------------------------------------------

# Define a function to select the earliest principal node in the trajectory graph
# based on the identity of a specific cell cluster (default is "OPC 1")
get_earliest_principal_node <- function(cds, ident="OPC 1") {
  
  # Identify cells belonging to the specified cluster
  cell_ids <- which(colData(cds)[, "anno_clus_dreamorigBB"] == ident)
  
  # Retrieve the closest vertices (nodes) in the trajectory graph for each cell
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  # Determine the root principal node as the most frequent vertex among the specified cells
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  # Return the root principal node
  root_pr_nodes
}

# Calculate pseudotime by ordering the cells starting from the selected root node
root_pr_nodes <- get_earliest_principal_node(cds, ident="OPC 1")
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

# Store the pseudotime values in the 'pseudotime' slot of the 'cds' object
cds$pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

# Plot the cells in the trajectory, coloring them by their pseudotime values
plot_cells(cds, color_cells_by = 'pseudotime',
           label_branch_points = FALSE, label_leaves = FALSE, label_roots = FALSE,
           group_label_size = 4, label_groups_by_cluster = FALSE)

#---------------------------------------------------------------------
# 4. Save Pseudotime and Graph Information
#---------------------------------------------------------------------

# Save the calculated pseudotime values to an RData file
pt.monocle3.sc <- pseudotime(cds)
save(pt.monocle3.sc, file = 'pt.RData')

# Save the entire Monocle3 cell_data_set object for later use
saveRDS(cds, file = 'cds.rds')

# Extract and save the graph information (UMAP coordinates) from the 'cds' object
y_to_cells <- principal_graph_aux(cds)$UMAP
saveRDS(y_to_cells, file = 'cds_graphInfo.rds')

#---------------------------------------------------------------------
# 5. Temporal Gene Identification using Moran's I Test
#---------------------------------------------------------------------

# Perform Moran's I test on the principal graph to identify 'temporal genes'
# that are correlated with pseudotime. Use 4 cores for parallel processing.
# Note: Make sure to correct 'rBind' to 'rbind' in the graph_test function if needed.
cds_moran_test <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

# Remove rows with missing values (NA) from the Moran's I test results
cds_moran_test <- cds_moran_test[complete.cases(cds_moran_test),]

# Save the Moran's I test results to an .rds file
saveRDS(cds_moran_test, file = 'cds_moranTest.rds')

#---------------------------------------------------------------------
# End of Script
#---------------------------------------------------------------------




# Summary of Steps:
# 1. Load necessary libraries for single-cell RNA-seq analysis and trajectory inference.
# 2. Load a pre-processed Seurat object and prepare it for Monocle3 by creating a new Seurat object
#    with necessary data (counts, normalized data, scaled data, and UMAP embeddings).
# 3. Convert the Seurat object into a Monocle3 cell_data_set, perform trajectory inference by 
#    pre-processing the data, clustering the cells, and learning a trajectory graph.
# 4. Calculate pseudotime by determining the root node of the trajectory graph and ordering the 
#    cells based on their position in the graph.
# 5. Save the pseudotime values, the entire Monocle3 object, and graph information for later analysis.
# 6. Identify temporal genes using Moran's I test on the principal graph, filtering out any incomplete
#    data, and save the results for further study.
