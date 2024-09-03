# ------------------------------------------------------------------------------
# Project: MAGMA and LDsc analysis
# Date: August 2024
# Script name: magma_ldsc.R
# Description: This script describes the enrichment (MAGMA / LDsc) and TOBIAS analysis

# outputs: Plots for Fig. 2C, 2E, 3A, 6C and S6

# ------------------------------------------------------------------------------


# -------------------------------
# Setup and Library Loading
# -------------------------------

library(data.table)   # Fast data manipulation
library(ggplot2)      # Create graphics with a flexible system
library(plyr)         # Data splitting, applying, and combinin
library(ggsci)        # Scientific color palettes for ggplot2
library(dplyr)        # Data manipulation with a consistent syntax
library(tidyverse)    # Collection of data science packages
library(ggrepel)      # Prevent text label overlap in plots


########################################################################################
##### CONFIG ###########################################################################

{
  ROOT = "~/Desktop/3dg_code/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!!
  INPUT_MARKERS_GENES = file.path(ROOT, "inputs", "markers_genes_celltype.csv") # Marker genes from CellRanger; # copied from /sc/arion/projects/CommonMind/3D_Genome/single_cell/snMultiome_NG_revisions2024/clarence/All_CELLTYPE_markers_log2FC0.20_pct0.15_SCT_2024_04_23.csv
  INPUT_MARKERS_PEAKS = file.path(ROOT, "inputs", "marker_peaks_celltype.csv")  # Marker peaks from CellRanger; # copied from /sc/arion/projects/CommonMind/3D_Genome/single_cell/snMultiome_NG_revisions2024/clarence/All_CELLTYPE_markers_log2FC0.20_pct0.15_ATAC_2024_04_23.csv
  INPUT_CELLSPEC_EREGULONS = file.path(ROOT, "inputs", "02.eRegulon_activor_CTonly_metadata_filtered_44eRegs.csv") # All eRegulons # copied from /sc/arion/projects/CommonMind/3D_Genome/single_cell/snMultiome_NG_revisions2024/Xuan/0.20240520_NG_revision/02.eRegulon_activor_CTonly_metadata_filtered_44eRegs.csv
  INPUT_ALL_EREGULONS = file.path(ROOT, "inputs", "02.eRegulon_activor_146_metadata_filtered.csv")                 # Cell type specific eRegulons # copied from /sc/arion/projects/CommonMind/3D_Genome/single_cell/snMultiome_NG_revisions2024/Xuan/0.20240520_NG_revision/02.eRegulon_activor_146_metadata_filtered.csv
  INPUT_FIG_2C_MAGMA = file.path(ROOT, "inputs", "Fig_2c_MAGMA.csv") # Pre-calculated MAGMA results for cell-specific eRegulons # copied from /sc/arion/projects/CommonMind/3D_Genome/single_cell/analysisJaro/magma_eRegulons_fig2_cellspec/meta-files/geneSetResults.tsv.gz
  INPUT_FIG_2C_LDSC = file.path(ROOT, "inputs", "Fig_2c_LDSC.csv")   # Pre-calculated LD-sc results for cell-specific eRegulons # copied from /sc/arion/projects/CommonMind/3D_Genome/single_cell/analysisJaro/ldsc_eRegulons_fig2_cellspec/meta-files/aggregatedPartInfo.tsv
  INPUT_FIG_3A_MAGMA = file.path(ROOT, "inputs", "Fig_3a_MAGMA.csv") # Pre-calculated MAGMA results for cell type marker genes # copied from /sc/arion/projects/CommonMind/3D_Genome/single_cell/results/cellspec_genes_rnaseq_sctransformed_default_ensemblProtCodGenes35kb10kbAutosomesNoBmhc/meta-files/geneSetResults.tsv.gz 
  INPUT_FIG_3B_LDSC = file.path(ROOT, "inputs", "Fig_3b_LDsc.csv")   # Pre-calculated LD-sc results for cell type marker peaks # copied from /sc/arion/projects/CommonMind/3D_Genome/single_cell/results/ldsc_marker_cell/meta-files/aggregatedPartInfo.tsv
  INPUT_FIG_6C = "/sc/arion/projects/CommonMind/3D_Genome/single_cell/analysisJaro/magma_eRegulons_fig6/meta-files/geneSetResults.tsv.gz" # Pre-calculated MAGMA results for all eRegulons # cp from /sc/arion/projects/CommonMind/3D_Genome/single_cell/analysisJaro/magma_eRegulons_fig6/meta-files/geneSetResults.tsv.gz
  
  TOBIAS_RESULTS_PATHS = list( # Pre-calculated TOBIAS results with footprinting scores
    "Mg_all_IKZF1" = file.path(ROOT, "inputs", "tobias", "bindetect_all", "IKZF1F4_M3450_1.02_overview.txt"),
    "Mg_all_IRF8" = file.path(ROOT, "inputs", "tobias", "bindetect_all", "IRF8_M5578_1.02_overview.txt"),
    "Mg_all_NFATC2" = file.path(ROOT, "inputs", "tobias", "bindetect_all", "NFATC1234_M1520_1.02_overview.txt"),
    "Mg_all_RUNX1" = file.path(ROOT, "inputs", "tobias", "bindetect_all", "RUNX123_M5794_1.02_overview.txt"),
    "Mg_all_SPI1" = file.path(ROOT, "inputs", "tobias", "bindetect_all", "SPI1SPICCTD-2545M3.6_M6484_1.02_overview.txt"),
    
    "Mg_isect_IKZF1" = file.path(ROOT, "inputs", "tobias", "bindetect_isect", "IKZF1F4_M3450_1.02_overview.txt"),
    "Mg_isect_IRF8" = file.path(ROOT, "inputs", "tobias", "bindetect_isect", "IRF8_M5578_1.02_overview.txt"),
    "Mg_isect_NFATC2" = file.path(ROOT, "inputs", "tobias", "bindetect_isect", "NFATC1234_M1520_1.02_overview.txt"),
    "Mg_isect_RUNX1" = file.path(ROOT, "inputs", "tobias", "bindetect_isect", "RUNX123_M5794_1.02_overview.txt"),
    "Mg_isect_SPI1" = file.path(ROOT, "inputs", "tobias", "bindetect_isect", "SPI1SPICCTD-2545M3.6_M6484_1.02_overview.txt"))
  
  # Conversion between GWAS trait abbreviation and desired final GWAS tratit name
  convTraitNames = list("adhd_ipsych"="ADHD","alzBellenguez"="AD", "asd"="ASD", "bip2" = "BD", "ms" = "MS", "pd_without_23andMe" = "PD", "als"="ALS",   "anorexia"="AN", "mdd_ipsych" = "MDD", "sz3"="SCZ") #alzKunkleNoapoe
}

########################################################################################
##### HELPER FUNCTIONS / CONSTANTS #####################################################

{
  # For diagonally split tiles in heatmap
  make_triangles = function(x, y, point = "up") {
    x <- as.integer(as.factor((x)))
    y <- as.integer(as.factor((y)))
    
    if (point == "up") {
      newx <- sapply(x, function(x) {
        c(x - 0.5, x - 0.5, x + 0.5)
      }, simplify = FALSE)
      newy <- sapply(y, function(y) {
        c(y - 0.5, y + 0.5, y + 0.5)
      }, simplify = FALSE)
    } else if (point == "down") {
      newx <- sapply(x, function(x) {
        c(x - 0.5, x + 0.5, x + 0.5)
      }, simplify = FALSE)
      newy <- sapply(y, function(y) {
        c(y - 0.5, y - 0.5, y + 0.5)
      }, simplify = FALSE)
    }
    data.frame(x = unlist(newx), y = unlist(newy))
  }
  
  myPalette = colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"), space = "Lab")
}

########################################################################################
##### FIG. 2C :: ENRICHMENT FOR CELL TYPE EREGULONS ####################################

{
  # Load cell type eRegulons & all eRegulons
  ctypeRegulonsDf = read.csv(INPUT_CELLSPEC_EREGULONS)
  allRegulonsDf = read.csv(INPUT_ALL_EREGULONS)
  background = unique(allRegulonsDf$Gene_ID)
  
  # Convert data frame to list, each item is named after eRegulon TF and contains regulatory targets (Ensembl Gene IDs)
  sets = lapply(unique(ctypeRegulonsDf$TF), function(tf_parent) {
    unique(ctypeRegulonsDf[(ctypeRegulonsDf$TF == tf_parent),"motif_target_EnsemblID"])
  })
  names(sets) = unique(ctypeRegulonsDf$TF)
  
  print(paste0("***** FIG. 2C :: INPUT FOR MAGMA / LD-sc *****"))
  print(paste0("+ Set of genes for cell specific eRegulons: ", paste0(sort(sapply(names(sets), function(setName) paste0(setName, " = ", length(sets[[setName]]), " genes"))), collapse="; ")))
  print(paste0("+ Background of all genes (present in at least one eRegulon): ", length(background), " genes"))
  
  # Comment on pre-calculation of MAGMA/LD-sc
  print("We skip this step here as it requires having MAGMA / LDSC installed and having processed GWAS summary stats available which cannot be done within GitHub repository.") 
  
  # Load LD-sc results
  ldsc = read.csv(INPUT_FIG_2C_LDSC)
  ldsc$adj.P.value = p.adjust(ldsc$p_regression, method="BH")
  ldsc$minusLog = -log10(ldsc$p_regression)
  ldsc$P = ldsc$p_regression
  ldsc$combo = paste0(ldsc$annoID, "_", ldsc$gwasAcronym)
  
  # Load MAGMA results
  magmaConditioning = read.csv(INPUT_FIG_2C_MAGMA)
  magmaConditioning$combo = paste0(magmaConditioning$annoID, "_", magmaConditioning$gwasAcronym)
  magmaConditioning$adj.P.value = p.adjust(magmaConditioning$P, method="BH")
  magmaConditioning$minusLog = -log10(magmaConditioning$P)
  magmaConditioning = magmaConditioning[order(magmaConditioning$combo),]
  
  # Merging LD-sc and MAGMA results
  isectCols = intersect(colnames(ldsc), colnames(magmaConditioning))
  df = cbind.data.frame(ldsc[order(ldsc$combo),setdiff(isectCols, "combo")], magmaConditioning[order(magmaConditioning$combo),setdiff(isectCols, "combo")])
  colnames(df) = c(paste0("ldsc_", setdiff(isectCols, "combo")), paste0("magma_", setdiff(isectCols, "combo")))
  df$magma_label = ifelse(df$magma_adj.P.value < 0.05, "#", ifelse(df$magma_P < 0.05, "·", ""))
  df$ldsc_label = ifelse(df$ldsc_adj.P.value < 0.05, "#", ifelse(df$ldsc_P < 0.05, "·", ""))
  
  # Calculation of coordinates for diagonal tiles
  newcoord_up = make_triangles(df$ldsc_annoID, df$ldsc_gwasAcronym)
  newcoord_down = make_triangles(df$ldsc_annoID, df$ldsc_gwasAcronym, point = "down")
  newcoord_down = newcoord_down %>% select(xdown = x, ydown = y) # just a dirty trick for renaming
  repdata = map_df(1:nrow(df), function(i) df[rep(i, 3), ])
  newdata = bind_cols(repdata, newcoord_up, newcoord_down)
  newdata$magma_label = sapply(1:nrow(newdata), function(i) {
    ifelse(i %% 3 == 1, newdata$magma_label[i], "")
  })
  newdata$ldsc_label = sapply(1:nrow(newdata), function(i) {
    ifelse(i %% 3 == 1, newdata$ldsc_label[i], "")
  })
  df$ldsc_gwasAcronym = ordered(df$ldsc_gwasAcronym, levels=unique(df$ldsc_gwasAcronym))
  
  # Make final plot Fig. 2c
  pdf(file.path(ROOT, "outputs", "Fig_2c.pdf"), width = 8 + 0.6 * length(unique(df$ldsc_gwasAcronym)), height = 3 + 0.2 * length(unique(df$ldsc_annoID)))
  print(ggplot(newdata) +
    geom_polygon(aes(x = x, y = y, fill = magma_minusLog, group = interaction(ldsc_annoID, ldsc_gwasAcronym)), color = "white") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 10)) + geom_text(aes(x = x, y = y, label=magma_label),alpha=0.7, nudge_x=0.25, nudge_y=0.70) + 
    ggnewscale::new_scale_fill() +
    geom_polygon(aes(x = xdown, y = ydown, fill = ldsc_minusLog, group = interaction(ldsc_annoID, ldsc_gwasAcronym)), color = "white") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 10)) + geom_text(aes(x = x, y = y, label=ldsc_label),alpha=0.7, nudge_x=0.7, nudge_y=0.25) + 
    scale_x_continuous(breaks = seq_along(levels(df$ldsc_annoID)), labels = levels(df$ldsc_annoID)) +
    scale_y_continuous(breaks = seq_along(unique(df$ldsc_gwasAcronym)), labels = unique(df$ldsc_gwasAcronym)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    coord_equal())
  dev.off()
}

########################################################################################
##### FIG. 2E :: FOOTPRINTING TOBIAS ANALYSIS ##########################################

{
  tobiasDetail = lapply(names(TOBIAS_RESULTS_PATHS), function(x) {
    df = read.csv(TOBIAS_RESULTS_PATHS[[x]], sep="\t");
    colnames(df)[11] = "score"
    colnames(df)[12] = "bound"
    rownames(df) = paste0(df$peak_id, "_", df$TFBS_start)
    df$type = x;
    df$combo = rownames(df)
    df
  })
  names(tobiasDetail) = names(TOBIAS_RESULTS_PATHS)
  
  tobiasScore = do.call("rbind.data.frame", lapply(tobiasDetail, function(x) { x[,c("score", "type", "combo", "bound")] }))
  tobiasScore$type = as.factor(tobiasScore$type)
  tobiasScore$typeGeneral = sapply(as.character(tobiasScore$type), function(x) {
    ifelse(startsWith(x, prefix="Mg_all_"), "all", "isect")
  })
  tobiasScore = tobiasScore[tobiasScore$typeGeneral %in% c("all", "isect"),]
  
  violin = ggplot(tobiasScore, aes(x=typeGeneral, y=score)) + theme_bw() + geom_violin(width = 1) + geom_boxplot(width=0.1, color="grey", alpha=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + ylim(c(0, 65))
  boxplot = ggplot(tobiasScore, aes(x=typeGeneral, y=score)) + theme_bw() + geom_boxplot(width=0.1) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + ylim(c(0, 65)) 
  
  print(paste0("> KS test: ", ks.test(tobiasScore[tobiasScore$typeGeneral=="all", "score"], tobiasScore[tobiasScore$typeGeneral=="isect", "score"], alternative="two.sided")$p.value))
  
  pdf(file.path(ROOT, "outputs", "Fig_2e_violin.pdf"), width=3, height=3); print(violin); dev.off()
  pdf(file.path(ROOT, "outputs", "Fig_2e_boxplot.pdf"), width=3, height=3); print(boxplot); dev.off()
}

########################################################################################
##### FIG. 3A :: ENRICHMENT FOR CELL TYPE MARKER GENES & PEAKS #########################

{
  # Load marker genes for MAGMA analysis
  markerGenes = read.csv(INPUT_MARKERS_GENES)
  magmaSet = lapply(unique(markerGenes$cluster), function(celltype) {
    unique(markerGenes[(markerGenes$cluster == celltype),"gene"])
  })
  names(magmaSet) = unique(markerGenes$cluster)
  
  # Load marker peaks for LDsc analysis
  markerPeaks = read.csv(INPUT_MARKERS_PEAKS)
  ldscSet = lapply(unique(markerPeaks$cluster), function(celltype) {
    unique(markerPeaks[(markerPeaks$cluster == celltype),"gene"])
  })
  names(ldscSet) = unique(markerPeaks$cluster)
  
  print(paste0("***** INPUT FOR MAGMA / LD-sc *****"))
  print(paste0("+ Set of genes for MAGMA: ", paste0(sort(sapply(names(magmaSet), function(setName) paste0(setName, " = ", length(magmaSet[[setName]]), " genes"))), collapse="; ")))
  print(paste0("+ Set of peaks for LD-sc: ", paste0(sort(sapply(names(ldscSet), function(setName) paste0(setName, " = ", length(ldscSet[[setName]]), " peaks"))), collapse="; ")))
  
  # Comment on pre-calculation of MAGMA/LD-sc
  print("We skip this step here as it requires having MAGMA / LDSC installed and having processed GWAS summary stats available which cannot be done within GitHub repository.") 
  
  ####
  # Load & plot MAGMA results
  magma = read.csv(INPUT_FIG_3A_MAGMA)
  magma = magma[magma$VARIABLE %in% c("GABAergic", "Glutamatergic", "OPC", "Astrocyte", "Oligodendrocyte", "Microglia", "Endothelial"),]
  magma$adj.P.value = p.adjust(magma$P, method="BH")
  magma$minusLogP = -log10(magma$P)
  magma$name_full = ordered(magma$name_full, levels=unique(magma$name_full))
  magma$size = ifelse(magma$adj.P.value < 0.05, 2, 1)
  
  magma$name_full = as.character(magma$name_full)
  magmaPlot = ggplot(magma, aes(x = gwasAcronym, y = minusLogP, color=name_full)) + geom_point(aes(fill = magma$minusLogP, size=magma$size)) + coord_flip() + theme_classic() + 
    theme(strip.text=element_text(face="bold"), axis.text.x = element_text(angle = 45)) + scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0.5)) + 
    geom_hline(yintercept=-log10(0.05), color="gray", linetype="dashed") + ylim(c(0, round(max(magma$minusLogP) + 1))) + xlab("")
  pdf(file.path(ROOT, "outputs", "Fig_3a_left.pdf"), width = 6, height = 6); print(magmaPlot); dev.off()
  
  ####
  # Load & plot LD-sc results
  ldsc = read.csv(INPUT_FIG_3B_LDSC)
  ldsc$adj.P.value = p.adjust(ldsc$p_regression, method="BH")
  ldsc$minusLogP = -log10(ldsc$p_regression)
  ldsc$P = ldsc$p_regression
  ldsc$size = ifelse(ldsc$adj.P.value < 0.05, 2, 1)
  
  ldsc$annoName = as.character(ldsc$annoName)
  ldscPlot = ggplot(ldsc, aes(x = gwasAcronym, y = minusLogP, color=annoName)) + geom_point(aes(fill = ldsc$minusLogP, size=ldsc$size)) + coord_flip() + theme_classic() + 
    theme(strip.text=element_text(face="bold"), axis.text.x = element_text(angle = 45)) + scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0.5)) + 
    geom_hline(yintercept=-log10(0.05), color="gray", linetype="dashed") + ylim(c(0, round(max(ldsc$minusLogP) + 1))) + xlab("")
  
  # Plotting final Fig. 3b
  pdf(file.path(ROOT, "outputs", "Fig_3b_right.pdf"), width = 6, height = 6); print(ldscPlot); dev.off()
}

########################################################################################
##### FIG. 6C :: ENRICHMENT FOR OLIG EREGULONS #########################################

{
  allRegulonsDf = read.csv(INPUT_ALL_EREGULONS)
  background = unique(allRegulonsDf$Gene_ID)
  
  # Get olig eRegulons
  eOLa = c("ASCL1", "ETV5", "OLIG2", "PRRX1", "SMAD3", "ZNF385D", "KLF12", "NFIB", "TCF7L1", "TCF7L2")
  lOLa = c("MITF", "FOXN2", "CPEB1", "ZNF189", "MYRF", "ATF7", "NFIX", "TFEB", "CREB5", "ZNF536")
  
  sets = lapply(c(eOLa, lOLa), function(tf_parent) { na.omit(unique(allRegulonsDf[(allRegulonsDf$TF == tf_parent),"Gene_ID"])) })
  names(sets) = c(eOLa, lOLa)
  
  print(paste0("***** INPUT FOR MAGMA *****"))
  print(paste0("+ Set of genes for early-stage and late-stage OLa eRegulons: ", paste0(sort(sapply(names(sets), function(setName) paste0(setName, " = ", length(sets[[setName]]), " genes"))), collapse="; ")))
  print(paste0("+ Background of all genes (present in at least one eRegulon): ", length(background), " genes"))
  
  # Comment on pre-calculation of MAGMA
  print("We skip this step here as it requires having MAGMA installed and having processed GWAS summary stats available which cannot be done within GitHub repository.") 
  
  # Load MAGMA results
  magmaConditioning = read.csv(INPUT_FIG_6C)
  magmaConditioning$annoID = ordered(gsub(".500bp.all", "", gsub("NKX6\\.", "NKX6-", magmaConditioning$VARIABLE)), levels=rev(c("ASCL1", "ETV5", "OLIG2", "PRRX1", "SMAD3", "ZNF385D", "KLF12", "NFIB", "TCF7L1", "TCF7L2", "MITF", "FOXN2", "CPEB1", "ZNF189", "MYRF", "ATF7", "NFIX", "TFEB", "CREB5", "ZNF536")))
  magmaConditioning$adj.P.value = p.adjust(magmaConditioning$P, method="BH")
  magmaConditioning$minusLog = -log10(magmaConditioning$P)
  
  # Plotting final Fig. 6c
  df = magmaConditioning
  df$plotLabel = ""
  df$plotLabel[df$P < 0.05] = "·"
  df$plotLabel[df$adj.P.value < 0.05] = "#"
  df$tempScoreCol = -log10(df$P)
  zz = ggplot(df, aes(gwasTrait, annoID, fill = tempScoreCol)) + geom_tile() + scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + ylab("Trait") + 
    xlab("Annotation") + theme_classic(base_size = 11) + theme(axis.text = element_text(colour = "black")) + 
    coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.title = element_text(size = 10, face = "bold"))
  pdf(file.path(ROOT, "outputs", "Fig_6c.pdf"), width = 6, height = 8)
  print(zz + geom_tile(aes(fill = tempScoreCol)) + scale_fill_gradientn(colours = myPalette(100), name = "-logP") + geom_text(aes(label = plotLabel), size = 11 * 0.55))
  dev.off()
}

########################################################################################
##### FIG. S6 :: ENRICHMENT FOR ALL EREGULONS ##########################################

{
  allRegulonsDf = read.csv(INPUT_ALL_EREGULONS)
  background = unique(allRegulonsDf$Gene_ID)
  
  sets = lapply(unique(allRegulonsDf$TF), function(tf_parent) { na.omit(unique(allRegulonsDf[(allRegulonsDf$TF == tf_parent),"Gene_ID"])) })
  names(sets) = unique(allRegulonsDf$TF)
  
  print(paste0("***** INPUT FOR LD-sc *****"))
  print(paste0("+ Set of genes for early-stage and late-stage OLa eRegulons: ", paste0(sort(sapply(names(sets), function(setName) paste0(setName, " = ", length(sets[[setName]]), " genes"))), collapse="; ")))
  
  LDSC_RESULTS = list("Astrocyte" = file.path(ROOT, "inputs", "Fig_S6_Astrocyte.tsv"), 
                       "Endothelial" = file.path(ROOT, "inputs", "Fig_S6_Endothelial.tsv"), 
                       "GABA" = file.path(ROOT, "inputs", "Fig_S6_GABA.tsv"), 
                       "GLU" = file.path(ROOT, "inputs", "Fig_S6_GLU.tsv"), 
                       "Microglia" = file.path(ROOT, "inputs", "Fig_S6_Microglia.tsv"), 
                       "Oligo" = file.path(ROOT, "inputs", "Fig_S6_Oligo.tsv"),
                      "OPC" = file.path(ROOT, "inputs", "Fig_S6_OPC.tsv"))
                      
  # Comment on pre-calculation of LD-sc
  print("We skip this step here as it requires having LD-sc installed and having processed GWAS summary stats available which cannot be done within GitHub repository.") 
  
  # Load pre-calculated LD-sc results 
  ldscAll = list("Astrocyte" = read.csv(LDSC_RESULTS[["Astrocyte"]], sep="\t"),
                 "Endothelial" = read.csv(LDSC_RESULTS[["Endothelial"]], sep="\t"),
                 "GABA" = read.csv(LDSC_RESULTS[["GABA"]], sep="\t"),
                 "GLU" = read.csv(LDSC_RESULTS[["GLU"]], sep="\t"),
                 "Microglia" = read.csv(LDSC_RESULTS[["Microglia"]], sep="\t"),
                 "Oligo" = read.csv(LDSC_RESULTS[["Oligo"]], sep="\t"),
                 "OPC" = read.csv(LDSC_RESULTS[["OPC"]], sep="\t"))
  
  df = data.frame(do.call("rbind.data.frame", ldscAll))
  df$celltype = sapply(rownames(df), function(x) strsplit(x, "\\.")[[1]][1])
  df$minusLog = -log10(df$p_regression)
  df$P = df$p_regression
  df = df[order(df$p_regression),]
  df = df[df$gwasAcronym %in% names(convTraitNames),]
  df$adj.P.value = p.adjust(df$p_regression, method="BH")
  df$size = ifelse(df$adj.P.value < 0.05, 2, 1)
  df$color = ifelse(df$adj.P.value < 0.05, "gray", "green")
  df$tf = df$origPeakSets
  df$tf_label = ifelse(df$adj.P.value < 0.05, df$tf, NA)
  df$tf_label = sapply(1:nrow(df), function(i) {
    ifelse((which(df[(df$celltype == df[i,"celltype"]) & (df$gwasAcronym == df[i,"gwasAcronym"]),"tf"] == df[i,"tf"]) <= 5), df[i,"tf_label"], "")
  })
  
  # Plotting Fig. S6
  df$size = ifelse(df$p_regression > 0.05, "#CCCCCC", ifelse(df$adj.P.value < 0.05, "#0072b2", "#d55e00"))
  ldscPlot = ggplot(df, aes(x = gwasAcronym, y = minusLog, color=size, label=tf_label)) + geom_point(aes(fill = minusLog, size=minusLog, label=tf_label)) + coord_flip() + theme_classic() + 
    theme(strip.text=element_text(face="bold"), axis.text.x = element_text(angle = 45), plot.margin = unit(c(1, 1, 1, 1), "cm")) + scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0.5)) + 
    geom_hline(yintercept=-log10(0.05), color="gray", linetype="dashed") + ylim(c(0, round(max(df$minusLog) + 1))) + 
    geom_text_repel(aes(label = tf_label), max.overlaps = Inf) + facet_wrap(~celltype, nrow=2)
  pdf(file.path(ROOT, "outputs", "Fig_S6.pdf"), width = 14, height = 8); print(ldscPlot); dev.off()
}

# -------------------------------
# End of Script
# -------------------------------
