library(Seurat)
library(SeuratWrappers)
library(monocle3)

opc_olig3 <- readRDS('OPCOlig_allprocessed.rds')

# create cds object with SCT + WNN UMAP from Seurat
o0 <- CreateSeuratObject(counts = opc_olig3[['SCT']]@counts,meta.data = opc_olig3@meta.data)
o0 <- SetAssayData(o0,slot = 'data',new.data = opc_olig3[['SCT']]@data)
o0 <- SetAssayData(o0,slot = 'scale.data',new.data = opc_olig3[['SCT']]@scale.data)
o0[['umap']] <- CreateDimReducObject(embeddings = Embeddings(opc_olig3,'wnn.dream2BB.umap'),key = 'UMAP_')
cds <- as.cell_data_set(o0)

cds <- preprocess_cds(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds,use_partition = F,close_loop = F)

# plot inferred trajectory
plot_cells(cds,color_cells_by = 'anno_clus_dreamorigBB',
           label_branch_points = F,label_leaves = F,label_roots = F,
           group_label_size = 4,label_groups_by_cluster = F)

# select root node
get_earliest_principal_node <- function(cds, ident="OPC 1"){
  cell_ids <- which(colData(cds)[, "anno_clus_dreamorigBB"] == ident)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
# calculate pseudotime
root_pr_nodes <- get_earliest_principal_node(cds,ident="OPC 1")
cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)
cds$pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
# plot pseudotime
plot_cells(cds,color_cells_by = 'pseudotime',
           label_branch_points = F,label_leaves = F,label_roots = F,
           group_label_size = 4,label_groups_by_cluster = F)

pt.monocle3.sc <- pseudotime(cds)
save(pt.monocle3.sc,file = 'pt.RData')
saveRDS(cds,file = 'cds.rds')
y_to_cells <-  principal_graph_aux(cds)$UMAP
saveRDS(y_to_cells,file = 'cds_graphInfo.rds')

# identify 'temporal genes' in monocle3 & nominate candidates
# need to change 'rBind' to 'rbind' in graph_test
cds_moran_test <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
cds_moran_test <- cds_moran_test[complete.cases(cds_moran_test),]
saveRDS(cds_moran_test,file = 'cds_moranTest.rds')

