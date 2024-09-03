# -----------------------------------------------------------------------------
# Project: snMultiome
# Date: August 2024
# Script name: abcmax.R
# Description: This script is designed for processing genomic data in the context of 
#              different GWAS traits and analyzing their linkage disequilibrium (LD) 
#              buddies. It loads various data sources, processes genomic information, 
#              and generates useful output for further analysis.
#              This script describes ABC-MAX analysis (connection between peaks 
#              overlapping GWAS SNPs and causal genes)
#
# outputs: Plots for Fig. 3-4, Fig. S8-S17, Table S8-S10
#
# -----------------------------------------------------------------------------

# -------------------------------
# Setup and Library Loading
# -------------------------------

library(plyr)  # Data manipulation and aggregation
library(data.table)  # Fast data manipulation for large datasets
library(GenomicRanges)  # Handling genomic intervals
library(bedr)  # Working with BED files and genomic ranges
library(liftOver)  # Converting genomic coordinates between assemblies
library(ggplot2)  # Creating customizable plots
library(edgeR)  # Differential expression analysis for count data
library(scales)  # Formatting and scaling for plots
library(ggsci)  # Color palettes for ggplot2


########################################################################################
##### CONFIG ###########################################################################

{
  ROOT = "~/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!!
  ABC_THRESHOLD = 0.015     # Minimal ABC score to consider it peak-gene as valid EP 
  MAX_LINKS_PER_PEAK = 1    # Maximum number of E-P links detected by ABC per peak in GWAS loci (per each cell type separately) 
}

########################################################################################
##### PROLOGUE: PATHS AND SETTING ANALYSIS TYPE ########################################

{
  # List of paths to table with index SNPs for analyzed GWAS traits 
  GWAS_LDBUDDIES_ROOT = file.path(ROOT, "inputs", "abcmax_gwas")
  GWAS = list("SCZ" = list(trait="SCZ", path=file.path(GWAS_LDBUDDIES_ROOT, "pgc3_scz_table.csv"), chr="Chromosome", left="range.left", right="range.right", index="top.pos", index_rs="top.index", p="top.P"),
              "BD" = list(trait="BD", path=file.path(GWAS_LDBUDDIES_ROOT, "BD_GWAS_mullins_2020.tsv"), chr="CHR", left="range.left", right="range.right", index="BP", index_rs="SNP", p="P"),
              "MDD" = list(trait="MDD", path=file.path(GWAS_LDBUDDIES_ROOT, "mdd_GWAS_2023.csv"), chr="CHR", left="NONE", right="NONE", index="BP.index.SNP", index_rs="leadSNP", p="P.x.index.SNP"),
              "AN" = list(trait="AN", path=file.path(GWAS_LDBUDDIES_ROOT, "anorexia_GWAS_2019.csv"), chr="CHR", left="left", right="right", index="BP", index_rs="SNP", p="P"),
              "ASD" = list(trait="ASD", path=file.path(GWAS_LDBUDDIES_ROOT, "ASD_GWAS_grove_2019.csv"), chr="CHR", left="left", right="right", index="BP", index_rs="SNP", p="P"),
              "ADHD" = list(trait="ADHD", path=file.path(GWAS_LDBUDDIES_ROOT, "adhd_GWAS_2022.csv"), chr="Chr", left="NONE", right="NONE", index="Position_bp", index_rs="rsID", p="P_value"),
              "AD" = list(trait="AD", path=file.path(GWAS_LDBUDDIES_ROOT, "ad_GWAS_meta_shea_2023_bellenguezSubset.csv"), chr="CHR", left="left", right="right", index="BP", index_rs="SNP", p="P"),
              "ALS" = list(trait="ALS", path=file.path(GWAS_LDBUDDIES_ROOT, "als_GWAS_2021_vanRheenen.csv"), chr="CHR", left="NONE", right="NONE", index="BP", index_rs="SNP", p="P"), # https://www.nature.com/articles/s41588-021-00973-1
              "PD" = list(trait="PD", path=file.path(GWAS_LDBUDDIES_ROOT, "pd_GWAS_2023_kim.csv"), chr="CHR", left="NONE", right="NONE", index="BP", index_rs="rsID", p="P_RE"), # https://www.medrxiv.org/content/10.1101/2022.08.04.22278432v1.full-text
              "MS" = list(trait="MS", path=file.path(GWAS_LDBUDDIES_ROOT, "ms_GWAS_2019.csv"), chr="CHR", left="NONE", right="NONE", index="POS", index_rs="SNP", p="EUR") # https://www-science-org.eresources.mssm.edu/doi/10.1126/science.aav7188#supplementary-materials
  )
  
  # Gene expression matrices from cell type psedobulks
  COUNT_MATRIX_SNRNASEQ_TO_BULK = list(
    "Astrocyte" = file.path(ROOT, "inputs", "abcmax_materials", "Astrocyte_all.exon.geneID.txt"),
    "Endothelial" = file.path(ROOT, "inputs", "abcmax_materials", "Endothelial_all.exon.geneID.txt"),
    "GABA" = file.path(ROOT, "inputs", "abcmax_materials", "GABA_all.exon.geneID.txt"),
    "GLU" = file.path(ROOT, "inputs", "abcmax_materials", "GLU_all.exon.geneID.txt"),
    "Microglia" = file.path(ROOT, "inputs", "abcmax_materials", "Microglia_all.exon.geneID.txt"),
    "OPC" = file.path(ROOT, "inputs", "abcmax_materials", "OPC_all.exon.geneID.txt"),
    "Oligo" = file.path(ROOT, "inputs", "abcmax_materials", "Oligo_all.exon.geneID.txt")
  )
  
  LD_BUDDIES = file.path(ROOT, "inputs", "abcmax_materials", "EUR_geno.txt.gz")  # LD buddies # from https://zenodo.org/record/3404275#.YWL0eVPMJdB
  TRUBETSKOY_PRIORITIZED_GENES = file.path(ROOT, "inputs", "Trubetskoy_Supplementary_Table_12.csv") # Table S12 from PMID:35396580
  TRUBETSKOY_PRIORITIZED_GENES2 = file.path(ROOT, "inputs", "Trubetskoy_Supplementary_Table_12b.csv") # Table S12 from PMID:35396580
  CELLTYPE_MARKERS_GENERAL = file.path(ROOT, "inputs", "abcmax_materials", "TopALL_3dg_ALL_Celltype_SCT_logfc0.25_minpct0.25.csv") # Cell type gene markers 
  COUNTMATRIX_PEAKS_GENERAL = file.path(ROOT, "inputs", "abcmax_materials", "readCountMatrices.RData")                             # Cell type peak markers 
  COUNTMATRIX_GENES_GENERAL = file.path(ROOT, "inputs", "abcmax_materials", "ps_3dg_RNA_celltype_brain.region_gene.csv")           # Count matrix on cell type pseudobulks 
  CYTOGENIC_LOCATIONS = file.path(ROOT, "inputs", "abcmax_materials", "gen.location.RDS")          # Convertor of genes to cytogenic locations  
  ENSEMBL_INFO_hg19 = file.path(ROOT, "inputs", "abcmax_materials", "muchEnsemblInfo_hg19.tsv.gz") # Downloaded Ensembl hg19 
  ENSEMBL_INFO_hg38 = file.path(ROOT, "inputs", "abcmax_materials", "muchEnsemblInfo_hg38.tsv.gz") # Downloaded Ensembl hg38 
  ENSEMBL_INFO = file.path(ROOT, "inputs", "abcmax_materials", "muchEnsemblInfo_hg38.tsv.gz")      # muchEnsemblInfo_hg38_v101.gz
  GENCODE = file.path(ROOT, "inputs", "abcmax_materials", "gencode.v30.annotation.gtf.gz")         # REFERENCE_HG38_GENCODE30/gencode.v30.annotation.gtf.gz
  
  allCells_env = new.env(); load(file.path(ROOT, "inputs", "abcmax_materials", "abc_and_peaks.RData"), envir=allCells_env) # Pre-calculated ABC E-P links 
  
  outDir = file.path(ROOT, "outputs", "abcmax")  # Set output dir
  dir.create(outDir, recursive=T)
  
  TX = list()
  CELL_TYPES = c("GABA", "GLU", "Oligo", "OPC", "Astrocyte", "Endothelial", "Microglia") # Vector of cell analyzed cell types
  TRAITS = c("AD", "ADHD", "ALS", "AN", "ASD", "BD", "MDD", "MS", "PD", "SCZ")           # Vector of analyzed GWAS traits
}

########################################################################################
##### HELPER FUNCTIONS #################################################################

{
  # Get cytogenic location
  get_location <- function(gene.name) {
    cyt.location <- readRDS(CYTOGENIC_LOCATIONS)
    test <- lapply(cyt.location, FUN = function(x){as.logical(which(x==gene.name))})
    ifelse (
      length(as.character(na.omit(names(test)[test==T]))) !=0 ,
      stringr::str_replace_all(
        as.character(na.omit(names(test)[test==T])),
        "\\d+",
        function(m) sprintf("%02d", as.integer(m)) ), 
      NA )
  }
  
  # Convert Ensembl gene name to Ensembl gene ID 
  findEnsemblNames = function(inputNames, keep_NA=F) {
    ensemblIds = sapply(inputNames, function(x) 
      ifelse(!is.na(which(ensemblInfo$Associated.Gene.Name==x)[1]), ensemblInfo[which(ensemblInfo$Associated.Gene.Name==x)[1], "Ensembl.Gene.ID"], 
             ifelse((which(ensemblInfo$HGNC.symbol==x)[1]), ensemblInfo[which(ensemblInfo$HGNC.symbol==x)[1], "Ensembl.Gene.ID"], 
                    ifelse((which(ensemblInfo$Associated.Transcript.Name==x)[1]), ensemblInfo[which(ensemblInfo$Associated.Transcript.Name==x)[1], "Ensembl.Gene.ID"], NA)))
    )
    if(keep_NA) {
      return(ensemblIds)
    } else {
      return(ensemblIds[!is.na(ensemblIds)])
    }
  }
  
  untransform = function(x) { 
    return(2^as.numeric(x))
  }
  
  # Plotting/CSV-exporting functions
  mpdf = function(x, width=7,height=7, outDir=outDir, onefile=T)eval.parent(substitute({ pdf(paste0(outDir,"/plot_",make.names(x),".pdf"),useDingbats=F,width=width,height=height,onefile=onefile) })) #outDir must be defined as a global var
  qplotly = function(x,command,outDir=outDir, width=700,height=700)eval.parent(substitute({
    if(!getOption("disablePlotly",default=F)){
      plotlyDir=paste0(outDir,"/plotly/") 
      dir.create(plotlyDir,showWarnings=F,recursive=T)
      htmlwidgets::saveWidget(as_widget(ggplotly(command, width=width, height=height)), paste0(plotlyDir, x, ".html"), selfcontained=F)
    }
  }))
}

########################################################################################
##### PROLOGUE :: PREPARE ENVIRONMENT ##################################################

{
  # Load Ensembl hg38
  ensemblInfo = read.delim(ENSEMBL_INFO_hg38, stringsAsFactors=F)
  ensemblGeneInfo = unique(ensemblInfo[,c("Ensembl.Gene.ID", "Strand", "Transcript.start..bp.", "Transcript.end..bp.", "Associated.Gene.Name", "Gene.type", "Status..gene.", "Source..gene.", "Transcript.count", "Description", "Chromosome.Name","HGNC.symbol", "Transcription.Start.Site..TSS."),])
  z = ddply(ensemblGeneInfo,c("Ensembl.Gene.ID"),summarize,HGNC.symbol=paste(HGNC.symbol,collapse=";"))
  ensemblGeneInfo$HGNC.symbol = z$HGNC.symbol[match(ensemblGeneInfo$Ensembl.Gene.ID,z$Ensembl.Gene.ID)]
  ensemblGeneInfo = unique(ensemblGeneInfo)
  rownames(ensemblGeneInfo) = ensemblGeneInfo$Ensembl.Gene.ID
  
  #####
  # Load Ensembl hg19
  ensemblInfo_hg19 = read.delim(ENSEMBL_INFO_hg19, stringsAsFactors=F)
  ensemblGeneInfo_hg19 = unique(ensemblInfo_hg19[,c("Ensembl.Gene.ID", "Strand", "Transcript.Start..bp.", "Transcript.End..bp.", "Associated.Gene.Name", "Gene.type", "Status..gene.", "Source..gene.", "Transcript.count", "Description", "Chromosome.Name","HGNC.symbol"),])
  z = ddply(ensemblGeneInfo_hg19,c("Ensembl.Gene.ID"),summarize,HGNC.symbol=paste(HGNC.symbol,collapse=";"))
  ensemblGeneInfo_hg19$HGNC.symbol = z$HGNC.symbol[match(ensemblGeneInfo_hg19$Ensembl.Gene.ID,z$Ensembl.Gene.ID)]
  ensemblGeneInfo_hg19 = unique(ensemblGeneInfo_hg19)
  
  ensemblGeneInfoShort = ensemblGeneInfo[,c("Chromosome.Name", "Transcript.start..bp.", "Transcript.end..bp.", "Transcription.Start.Site..TSS.", "Ensembl.Gene.ID", "Strand")]
  colnames(ensemblGeneInfoShort) = c("seqnames", "start", "end", "tss", "name", "strand")
  ensemblGeneInfoShort$seqnames = paste0("chr", ensemblGeneInfoShort$seqnames)
  ensemblGeneInfoShort = ensemblGeneInfoShort[ensemblGeneInfoShort$seqnames %in% paste0("chr", c(1:22, "X")),]
  ensemblGeneInfoShort$strand = ifelse(ensemblGeneInfoShort$strand == "1", "+", "-")
  rownames(ensemblGeneInfoShort) = ensemblGeneInfoShort$name
  
  #####
  # Load cell type marker genes
  ctypeDiff = read.csv(CELLTYPE_MARKERS_GENERAL)
  ctypeDiff$gene_id = findEnsemblNames(ctypeDiff$gene, keep_NA = T)

  # Load gene count matrix (pseudobulk)
  df = data.frame(t(read.csv(COUNTMATRIX_GENES_GENERAL)))
  celltypes = c("Astrocyte", "Endothelial", "GABAergic", "Glutamatergic", "Microglia", "OPC", "Oligodendrocyte")
  colnames(df) = df[1,]
  df = df[-1,]  
  geneNames = rownames(df)
  df = data.frame(apply(df, 2, untransform))
  rownames(df) = geneNames
  geneMatrix = do.call("cbind.data.frame", lapply(celltypes, function(ctype) { rowMeans(df[,startsWith(colnames(df), prefix=ctype)]) }))

  best_cellType = sapply(1:nrow(geneMatrix), function(i) { names(which.max(geneMatrix[i,])) })
  best_cellType_ratio = sapply(1:nrow(geneMatrix), function(i) { mx = which.max(geneMatrix[i,]); geneMatrix[i,mx] / sum(geneMatrix[i,]) })
  
  qcGeneAnno = cbind.data.frame(rownames(geneMatrix), best_cellType, best_cellType_ratio)
  colnames(qcGeneAnno) = c("GeneID", "best_cellType", "best_cellType_ratio")
}

########################################################################################
##### LOAD PEAK COUNT MATRIX (PSEUDOBULK) ##############################################

{
  tmpEnv = new.env(); load(COUNTMATRIX_PEAKS_GENERAL, envir=tmpEnv)
  peakMatrix = data.frame(apply(tmpEnv$initialVoomObj$E, 2, untransform))
  qcPeakAnno = allCells_env$mergedPeaksDf
  colnames(peakMatrix) = gsub("_all$", "", colnames(peakMatrix))
  tmp = allCells_env$barcodes[!duplicated(allCells_env$barcodes$sampleID),c("sampleID", "age.group")]
  rownames(tmp) = tmp$sampleID
     
  for(ctype in names(allCells_env$abcResults_cell)) {
    allCells_env$abcResults_cell[[ctype]]$CellType_general = ctype
    allCells_env$abcResults_cell[[ctype]]$SampleID = sapply(allCells_env$abcResults_cell[[ctype]]$CellType, function(x) strsplit(x, "_")[[1]][1])
    allCells_env$abcResults_cell[[ctype]]$AgePeriod = tmp[allCells_env$abcResults_cell[[ctype]]$SampleID, "age.group"]
  }
  rownames(peakMatrix) = rownames(tmpEnv$initialVoomObj)
  
  qcPeakAnno$best_cellType = sapply(1:nrow(qcPeakAnno), function(i) { names(which.max(peakMatrix[i,])) }) # For each peak, find cell type in which it it the most expressed (~accessible)
  qcPeakAnno$best_cellType_ratio = sapply(1:nrow(qcPeakAnno), function(i) { mx = which.max(peakMatrix[i,]); peakMatrix[i,mx] / sum(peakMatrix[i,]) }) # ... and proportion of expression for that most expressed cell type (compared to the other cell types)
  cellTypeRatiosDf = do.call("cbind.data.frame", lapply(colnames(peakMatrix), function(ctype) {
    sapply(1:nrow(peakMatrix), function(i) {
      peakMatrix[i,ctype] / sum(peakMatrix[i,])
    })
  }))
  colnames(cellTypeRatiosDf) = colnames(peakMatrix)
  qcPeakAnno = cbind.data.frame(qcPeakAnno, cellTypeRatiosDf)
}

########################################################################################
##### LOAD GWAS INDEX SNPS AND LD BUDDIES ##############################################

{
  ld_buddies_list = list()
  lociTable_list = list()
  ld_main_and_buddies_list = list()
  
  for(gwas in GWAS) {
    print(paste0("*** PROCESSING ", gwas$trait, " GWAS"))
    
    # Load GWAS info
    lociTable = read.csv(gwas$path, sep="\t")
    lociTable_list[[gwas$trait]] = lociTable[,c(gwas$index_rs, gwas$chr, gwas$index, gwas$index, gwas$p)]
    colnames(lociTable_list[[gwas$trait]]) = c("rs", "chr", "start", "end", "p")
    lociTable_list[[gwas$trait]]$start = lociTable_list[[gwas$trait]]$start - 1
    lociTable_list[[gwas$trait]]$trait = gwas$trait
    
    lociTableTmp = lociTable_list[[gwas$trait]][,c("chr", "start", "end", "rs", "chr", "start", "end", "rs", "rs", "p")]
    colnames(lociTableTmp) = c("chr", "start", "end", "rs", "buddy_chr", "buddy_start", "buddy_end", "buddy_rs", "r2", "p")
    lociTableTmp$r2 = NA
    lociTableTmp$trait = gwas$trait
    lociTableTmp$is_index_snp = T
    
    # Find LD buddies for GWAS index SNP locations
    ld_buddies = tabix(region=c(paste0(lociTable[,gwas$chr], ":", lociTable[,gwas$index]-1, "-", lociTable[,gwas$index])), LD_BUDDIES)
    if(!is.null(ld_buddies) == 0) {
      next
    }
    colnames(ld_buddies) = c("chr", "start", "end", "rs", "buddy_chr", "buddy_start", "buddy_end", "buddy_rs", "r2")
    
    TX[[paste0("ld_buddy_search_", gwas$trait, "_1")]] = paste0("GWAS ", gwas$trait, ": We found at least 1 LD buddy for ", sum(unique(lociTable[,gwas$index_rs]) %in% unique(ld_buddies$rs)), " out of ", length(unique(lociTable[,gwas$index])), " GWAS index SNPs.")
    TX[[paste0("ld_buddy_search_", gwas$trait, "_2")]] = paste0("GWAS ", gwas$trait, ": We found ", min(table(ld_buddies$rs)), "-", max(table(ld_buddies$rs)), " LD buddies (avg/median=", round(mean(nrow(ld_buddies)) / nrow(lociTable)), "/", median(table(ld_buddies$rs)), ") for GWAS index SNPs; total buddies = ", length(unique(ld_buddies$buddy_rs)))
    print(TX[[paste0("ld_buddy_search_", gwas$trait, "_1")]])
    print(TX[[paste0("ld_buddy_search_", gwas$trait, "_2")]])
    
    ld_buddies$trait = gwas$trait
    ld_buddies = ld_buddies[(ld_buddies$rs %in% lociTableTmp$rs),]
    ld_buddies$p = lociTableTmp[match(ld_buddies$rs, lociTableTmp$rs),"p"]
    ld_buddies$is_index_snp = F
    ld_buddies_list[[gwas$trait]] = ld_buddies
    
    ld_main_and_buddies_list[[gwas$trait]] = rbind.data.frame(ld_buddies, lociTableTmp)
  }
  
  ld_buddies_df = data.frame(do.call("rbind", ld_buddies_list))
  lociTable_df = data.frame(do.call("rbind", lociTable_list))
  ld_main_and_buddies_df = data.frame(do.call("rbind", ld_main_and_buddies_list))
  ld_main_and_buddies_df$combo = paste0(ld_main_and_buddies_df$buddy_rs, "_", ld_main_and_buddies_df$chr, "_", ld_main_and_buddies_df$start, "_", ld_main_and_buddies_df$end)
  
  ld_main_and_buddies_df$rsAlt = paste0(ld_main_and_buddies_df$chr, "_", ld_main_and_buddies_df$start) # that's because some buddies don't have rsID so they wouldn't be counted
  ld_main_and_buddies_df$buddy_rsAlt = paste0(ld_main_and_buddies_df$buddy_chr, "_", ld_main_and_buddies_df$buddy_start) # that's because some buddies don't have rsID so they wouldn't be counted
  
  ld_main_and_buddies_df$rs = sapply(1:nrow(ld_main_and_buddies_df), function(i) {
    ifelse(startsWith(ld_main_and_buddies_df[i,"rs"], "rs"), ld_main_and_buddies_df[i,"rs"], paste0(ld_main_and_buddies_df[i,"chr"], "_", ld_main_and_buddies_df[i,"end"]))
  })
  ld_main_and_buddies_df$buddy_rs = sapply(1:nrow(ld_main_and_buddies_df), function(i) {
    ifelse(startsWith(ld_main_and_buddies_df[i,"buddy_rs"], "rs"), ld_main_and_buddies_df[i,"buddy_rs"], paste0(ld_main_and_buddies_df[i,"chr"], "_", ld_main_and_buddies_df[i,"buddy_end"]))
  })
  uniqIndexSnps = length(unique(ld_main_and_buddies_df$rsAlt))
  uniqIndexBuddies = length(unique(ld_main_and_buddies_df$buddy_rsAlt))
  write.table(ld_main_and_buddies_df, file="./ld_buddies.tsv", quote=F, row.names=F, sep="\t")
  
  TX[["buddy_search_global"]] = paste0("Initially, we assembled a set of ", uniqIndexSnps, " genome-wide significant variants linked to select neuropsychiatric and neurodegenerative traits (with a significance threshold of P<5×10-8) and expanded this set to include ", uniqIndexBuddies, " variants based on their high linkage disequilibrium (LD; R2≥0.8). [[total of ", length(c(unique(ld_main_and_buddies_df$rsAlt), unique(ld_main_and_buddies_df$buddy_rsAlt))), "]]")
  print(TX[["buddy_search_global"]])
}

########################################################################################
##### FUNCTIONS FOR CONVERTING PEAKS (FROM ABC-LINKS) TO HG19 ##########################

{
  # Loading liftover and chains for conversion
  chPath = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  ch = import.chain(chPath)
  
  chPath2 = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
  ch2 = import.chain(chPath2)
  
  # Helper function for converting hg38 cell-type to hg19 (lifted peak cannot have size shorter or larger than 10% compared to original one)
  toHg19 = function(x) {
    x$name = paste0(x@seqnames, "_", x@ranges@start, "_", x@ranges@width, "_", x$TargetGene)
    x@ranges@NAMES = x$name
    hg19 = liftOver(x, ch)@unlistData
    hg19 = unlist(range(split(hg19, ~name)))
    
    idsInBothChrVersions = as.character(intersect(names(hg19), x$name))
    differences = (hg19[idsInBothChrVersions,]@ranges@width - x[idsInBothChrVersions,]@ranges@width) / x[idsInBothChrVersions,]@ranges@width
    idsTooDifferent = idsInBothChrVersions[which(abs(differences) > 0.10)]
    
    hg19 = hg19[!hg19@ranges@NAMES %in% idsTooDifferent,]
    hg19$gene = x[hg19@ranges@NAMES,]@elementMetadata$TargetGene
    hg19$gene_name = x[hg19@ranges@NAMES,]@elementMetadata$TargetGene_name
    hg19$peak = x[hg19@ranges@NAMES,]@elementMetadata$PeakID
    hg19$ABC.Score = x[hg19@ranges@NAMES,]@elementMetadata$ABC.Score
    hg19$AgePeriod = x[hg19@ranges@NAMES,]@elementMetadata$AgePeriod
    hg19$SampleID = x[hg19@ranges@NAMES,]@elementMetadata$SampleID
    hg19$CellType_general = x[hg19@ranges@NAMES,]@elementMetadata$CellType_general 
    hg19$TargetGeneTSS = x[hg19@ranges@NAMES,]@elementMetadata$TargetGeneTSS
    names(hg19) = NULL
    
    return(hg19)
  }
  
  # Helper function for converting hg38 cell-type to hg19 (lifted peak cannot have size shorter or larger than 10% compared to original one)
  toHg38 = function(x) {
    x@ranges@NAMES = x$name
    hg38 = liftOver(x, ch2)@unlistData
    return(hg38)
  }
}

########################################################################################
##### LOAD PEAK-TO-GENES LINKS #########################################################

{
  # Convert ABC E-P links to GenomicRanges
  epLinks = lapply(allCells_env$abcResults_cell, function(x) { x$start = x$PeakID_start; x$end = x$PeakID_end; gr = makeGRangesFromDataFrame(x, keep.extra.columns=T) })
  names(epLinks) = names(allCells_env$abcResults_cell)

  # Remove VLMC_Pericyte links
  epLinks$VLMC_Pericyte = NULL
  epLinks$VLMC = NULL
  allCells_env$abcResults_cell_age_df = allCells_env$abcResults_cell_age_df[!(allCells_env$abcResults_cell_age_df$cell_type %in% c("VLMC", "VLMC_Pericyte")),]
  allCells_env$abcResults_cell_df = allCells_env$abcResults_cell_df[!(allCells_env$abcResults_cell_df$cell_type %in% c("VLMC", "VLMC_Pericyte")),]
  
  # Convert ABC E-P links from hg38 to hg19
  epLinks_hg19 = lapply(epLinks, function(x) { toHg19(x) })
  names(epLinks_hg19) = names(epLinks)
  
  # Get summary counts of EP links
  tx_linksUnique = unique(unlist(lapply(epLinks, function(x) unique(x$combo))))
  tx_linksAll = sum(unlist(lapply(epLinks, function(x) length(unique(x$combo)))))
  tx_genesLinkedUnique = length(unique(unlist(lapply(epLinks, function(x) unique(x$TargetGene)))))
  tx_peaksLinkedUnique = length(unique(unlist(lapply(epLinks, function(x) unique(x$PeakID)))))
  tx_totalEnhancers = nrow(allCells_env$mergedPeaksDf)
  TX$maps_1 = paste0("Collectively, we identified ", length(tx_linksUnique), " enhancer–gene connections for ", tx_genesLinkedUnique, " expressed genes and ", tx_peaksLinkedUnique, " unique enhancers (note - total enhancers: ", tx_totalEnhancers , ").")
  print(TX$maps_1)
  
  # Find out how many ABC E-P links are present in at least two developmental periods
  atLeastTwoDevelPeriods = sapply(CELL_TYPES, function(ctype) {
    epLinksList = sapply(allCells_env$abcResults_cell_age[startsWith(names(allCells_env$abcResults_cell_age), prefix=ctype)], function(df) {
      df$combo
    })
    freqTable = table(table(unlist(epLinksList)))
    sum(freqTable[2:length(freqTable)]) / sum(freqTable)
  })
  tx_linksFreq = table(unlist(lapply(epLinks, function(x) unique(x$combo))))
  tx_linksFreqSum = table(tx_linksFreq)
  atLeastTwoCellTypes = sum(tx_linksFreqSum[2:length(tx_linksFreqSum)]) / sum(tx_linksFreqSum)
  TX$maps_2 = paste0("On average, ", round(atLeastTwoCellTypes * 100), "% of the E–P links were shared in at least two cell types, whereas ", round(min(atLeastTwoDevelPeriods) * 100), "–", round(max(atLeastTwoDevelPeriods) * 100), "% of E–P links were shared across regions, when comparing within the same cell type, with a high correlation of ABC score (demanding-to-calculate).")
  print(TX$maps_2)

  peaksMelted = allCells_env$abcResults_cell_age_df
  peaks_count_celltype_age = ggplot(data=peaksMelted, aes(x=category, y=abcLinks_count)) + geom_bar(stat="identity") + xlab("Cell type & Developmental period") + ylab("Number of E-P links") +
    theme_classic() + theme(axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1))
  mpdf(paste0("abc_celltype_age"), width=12, height=7); print(peaks_count_celltype_age); dev.off()
}

########################################################################################
##### FORMAT OUTPUT DATAFRAME (WITH PEAK-GENE LINKS) ###################################

{
  # Get unique rsIDs (skip ".")
  shared_rs = names(which(table(ld_buddies_df$buddy_rs) > 0)) 
  shared_rs = shared_rs[shared_rs != "."]
  
  # List of associations between "rsIDs - associated disease GWAS" 
  allcomb = list()
  allcomb[unique(ld_main_and_buddies_df$rs)] = ld_main_and_buddies_df[match(unique(ld_main_and_buddies_df$rs), ld_main_and_buddies_df$rs), "trait"]
  for(rsAlt in shared_rs) {
    allcomb[[rsAlt]] = paste0(sort(ld_main_and_buddies_df[which(ld_main_and_buddies_df$buddy_rsAlt == rsAlt),"trait"]), collapse=",")
    buddies_rsSource = ld_main_and_buddies_df[which(ld_main_and_buddies_df$buddy_rsAlt==rsAlt), "rsAlt"]
    allcomb[buddies_rsSource] = paste0(sort(ld_main_and_buddies_df[which(ld_main_and_buddies_df$buddy_rs==rsAlt),"trait"]), collapse=",")
  }
  lociTable = ld_main_and_buddies_df[, c("buddy_chr", "buddy_start", "buddy_end", "buddy_rs", "trait", "rs")]
  colnames(lociTable) = c("chr", "start", "end", "rs", "trait", "rsx")
  lociTable$end = as.numeric(lociTable$start)
  lociTable$start = as.numeric(lociTable$start) - 1
  lociTable$chr = gsub("chr23", "chrX", lociTable$chr)
  lociTableGr = makeGRangesFromDataFrame(lociTable, keep.extra.columns=T)
  
  lociTable$id = ifelse(lociTable$rs != ".", lociTable$rs, paste0(lociTable$chr, "_", lociTable$start))
  lociTable$rsAlt = paste0(lociTable$chr, "_", lociTable$start)
  lociTable = lociTable[!duplicated(lociTable$id),]
  rownames(lociTable) = lociTable$id
  
  lociTable_hg38 = data.frame(toHg38(makeGRangesFromDataFrame(lociTable, keep.extra.columns=T)))
  lociTable_hg38$desc = paste0(lociTable_hg38$seqnames, "_", lociTable_hg38$start, "_", lociTable_hg38$end)
  rownames(lociTable_hg38) = lociTable_hg38$id
  lociTable$pos_hg38 = lociTable_hg38[lociTable$id, "end"]
  
  # Export locations of GWAS rsIDs and their buddies
  exportLociTable = lociTable_hg38[,c("seqnames", "start", "end", "id","trait")]
  exportLociTable$score = ifelse(exportLociTable$id %in% unique(ld_main_and_buddies_df$rs[(ld_main_and_buddies_df$rs==ld_main_and_buddies_df$buddy_rs)]), 1000, 500)
  write.table(exportLociTable, row.names=F, sep="\t", quote=F, file=file.path(outDir, "rs_locations.bed"))
  
  for(trait in unique(exportLociTable$trait)) {
    df = exportLociTable[(exportLociTable$trait == trait) & (exportLociTable$score == 1000),]
    write.table(df[,!colnames(df) %in% c("trait")], row.names=F, sep="\t", quote=F, file=file.path(outDir, paste0("rs_locations_", trait, "_index.bed")))
    df = exportLociTable[(exportLociTable$trait == trait) & (exportLociTable$score == 500),]
    write.table(df[,!colnames(df) %in% c("trait")], row.names=F, sep="\t", quote=F, file=file.path(outDir, paste0("rs_locations_", trait, "_buddy.bed")))
  }
  
  allResults = list()
  for(ctype in names(epLinks_hg19)) {
    print(ctype)
    hits = findOverlaps(lociTableGr, epLinks_hg19[[ctype]])
    results = cbind.data.frame(ld_main_and_buddies_df[hits@from,], epLinks_hg19[[ctype]][hits@to,])
    colnames(results) = c("chr", "start", "end", "rs", "trait_buddy_chr", "trait_buddy_start", "trait_buddy_end", "trait_buddy_rs", "r2", "trait", "p", "is_index_snp", "combo", "rsAlt", "buddy_rsAlt", "linked_chr", "linked_start", "linked_end", "linked_width", "linked_strand", "gene", "gene_name", "peak", "ABC.Score", "AgePeriod", "SampleID", "CellType_general", "TargetGeneTSS_hg38")
    results$linked_chr = qcPeakAnno[results$peak, "chr"]
    results$linked_start = qcPeakAnno[results$peak, "start"]
    results$linked_end = qcPeakAnno[results$peak, "end"]
    
    results$p = as.numeric(results$p)
    results$trait_gene = paste0(results$gene, "_", results$trait)
    results$trait_gene_ldbuddy = paste0(results$gene, "_", results$trait, "_", results$trait_buddy_rs)
    results$trait_gene_rs = paste0(results$gene, "_", results$trait, "_", results$rs)
    results[which("AC018685.2" == results$gene), "cyt_loc"] = "chr02p25" # get_location("EIPR1"); get_location("MYT1L") #fixme
    results[which("AC096570.1" == results$gene), "cyt_loc"] = "chr02p24"
    results[which("AC156455.1" == results$gene), "cyt_loc"] = "chr12q24"
    results[which("U91319.1" == results$gene), "cyt_loc"] = "chr16p13" 
    results[which("AC138647.1" == results$gene), "cyt_loc"] = "chr08q24"
    results[which("HIST1H2AK" == results$gene), "cyt_loc"] = "chr08q24"
    results = results[!duplicated(results$trait_gene_ldbuddy),]
    dim(results[which(results$peak=="Peak_505316"),])
    
    results$ctypeMarker = sapply(results$gene, function(geneName) {
      subset = ctypeDiff[which(geneName==ctypeDiff$gene_id),]
      if(nrow(subset) > 0) {
        paste0(subset[order(subset$p_val, -subset$avg_log2FC),"cluster"], collapse=",")
      } else {
        NA
      }
    })
    results = results[results$peak %in% rownames(qcPeakAnno),]
    results$is_marker_gene = ifelse(is.na(results$ctypeMarker), F, T)
    results$peak_best_cellType = qcPeakAnno[results$peak, "best_cellType"]
    results$peak_best_cellType_ratio = qcPeakAnno[results$peak, "best_cellType_ratio"]
    
    results$best_cellType = qcGeneAnno[results$gene_name, "best_cellType"]
    results$best_cellType_ratio = qcGeneAnno[results$gene_name, "best_cellType_ratio"]
    
    allResults[[ctype]] = results
  }
  
  final = data.frame(do.call("rbind.data.frame", allResults))
  final$chr = gsub("chr23", "chrX", final$chr)
  final = final[final$linked_chr == final$chr,] # this looks as a weird operation but hg38->hg19 can even change chromosome, those cases are rare so let's drop them
  final$gene_biotype = ensemblGeneInfo[match(final$gene, ensemblGeneInfo$Ensembl.Gene.ID),"Gene.type"]
  final$gene_name = sapply(final$gene, function(x) { ensemblInfo[match(x, ensemblInfo$Ensembl.Gene.ID),"Associated.Gene.Name"] })
  final$cyt_loc = sapply(final$gene_name, function(x) get_location(x))
}

########################################################################################
##### FILTER OUT E-P LINKS THAT ARE NOT EXPRESSED IN THE CELL TYPE #####################

{
  countMatrix = do.call("cbind.data.frame", lapply(names(COUNT_MATRIX_SNRNASEQ_TO_BULK), function(ctype) {
    print(ctype)
    cReadCount = read.csv(COUNT_MATRIX_SNRNASEQ_TO_BULK[[ctype]], sep="\t", header=T, skip=1, stringsAsFactors=F)[,7]
    cReadCount
  }))
  colnames(countMatrix) = names(COUNT_MATRIX_SNRNASEQ_TO_BULK)
  head(countMatrix)
  df = read.csv(COUNT_MATRIX_SNRNASEQ_TO_BULK[[ctype]], sep="\t", header=T, skip=1, stringsAsFactors=F)
  countMatrix$Geneid = sapply(df$Geneid, function(x) { strsplit(x, "\\.")[[1]][1] })
  countMatrix = countMatrix[!duplicated(countMatrix$Geneid),]
  rownames(countMatrix) = countMatrix$Geneid
  countMatrix = countMatrix[,1:(length(countMatrix)-1)]
  
  ensemblInfo = read.csv(ENSEMBL_INFO, sep="\t", header=T, stringsAsFactors=F)
  ensemblInfo = ensemblInfo[,c("Ensembl.Gene.ID", "Transcript.length..including.UTRs.and.CDS.", "Gene.type", "Gene...GC.content")]
  colnames(ensemblInfo) = c("Ensembl.Gene.ID", "Length", "transcript_biotype", "gcContent")
  gtf = rtracklayer::import(GENCODE)
  gtf_df = as.data.frame(gtf)
  gtf_df = data.frame(gtf_df)
  gtf_df_new = gtf_df[which(gtf_df$type=="gene"),]
  gtf_df_new = gtf_df_new[,c("gene_id","seqnames","start","end","strand","gene_name","transcript_name","width","gene_type")]
  gtf_df_new$exonLength = ensemblInfo[match(sapply(strsplit(as.character(gtf_df_new[,"gene_id"]), "\\."), "[[", 1), ensemblInfo[,"Ensembl.Gene.ID"]),]$Length
  gtf_df_new$biotype = ensemblInfo[match(sapply(strsplit(as.character(gtf_df_new[,"gene_id"]), "\\."), "[[", 1), ensemblInfo[,"Ensembl.Gene.ID"]),]$transcript_biotype
  gtf_df_new$gcContent = ensemblInfo[match(sapply(strsplit(as.character(gtf_df_new[,"gene_id"]), "\\."), "[[", 1), ensemblInfo[,"Ensembl.Gene.ID"]),]$gcContent
  gtf_df_new = gtf_df_new[!is.na(match(sapply(strsplit(as.character(gtf_df_new[,"gene_id"]), "\\."), "[[", 1), ensemblInfo[,"Ensembl.Gene.ID"])),]
  qcPeakAnnoGene = gtf_df_new
  qcPeakAnnoGene$gene_id = sapply(qcPeakAnnoGene$gene_id, function(x) { strsplit(x, "\\.")[[1]][1] })
  qcPeakAnnoGene = qcPeakAnnoGene[!duplicated(qcPeakAnnoGene$gene_id),]
  rownames(qcPeakAnnoGene) = qcPeakAnnoGene$gene_id
  
  datExpr = data.frame(countMatrix)
  expObjAll = DGEList(counts=datExpr, genes=qcPeakAnnoGene[match(rownames(datExpr),qcPeakAnnoGene$PeakID),])
  
  MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM = 0.5
  expressedGenes = lapply(names(COUNT_MATRIX_SNRNASEQ_TO_BULK), function(ctype) {
    fracSamplesWithMinCPM = rowMeans(cpm(expObjAll[,ctype]) >= 1)
    isNonLowExpr = fracSamplesWithMinCPM >= MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM
    expObjNonLow = expObjAll[isNonLowExpr, , keep.lib.sizes=F]
    rownames(expObjNonLow$genes)
  })
  names(expressedGenes) = names(COUNT_MATRIX_SNRNASEQ_TO_BULK)

  # Not necessary as these links were already filtered when ABC was calculated
  flags = sapply(1:nrow(final), function(i) { 
    final[i,"gene"] %in% expressedGenes[[final$CellType_general[i]]]
  })
  final = final[flags,]
}

###################################################################################################################################################
##### ADD INDICATION WHETHER GENE (FROM PEAK-GENE LINKE) IS CELL TYPE MARKER (FROM EXTERNAL STUDIES) OR WHETHER PEAK IS OVERLAPPING FINEMAP SNP ###

{
  #####
  # Prioritized genes from PGC3 study (Supplementary Table 12)
  pgcGenes = read.csv(TRUBETSKOY_PRIORITIZED_GENES, sep="\t")
  pgcGenes2 = read.csv(TRUBETSKOY_PRIORITIZED_GENES2, sep="\t")
  pgcGenes$rsGene = paste0(pgcGenes$Index.SNP, "_", pgcGenes$Ensembl.ID)
  pgcGenes2$rsGene = paste0(pgcGenes2$Index.SNP, "_", pgcGenes2$Ensembl.ID)
   
  final$pgc3_prioritized1 = sapply(1:nrow(final), function(i) {
    ifelse(final$rs[i] %in% pgcGenes$Index.SNP, T, F)
  })
  table(final$pgc3_prioritized1)

  final$pgc3_prioritized2 = sapply(1:nrow(final), function(i) {
    ifelse(final$rs[i] %in% pgcGenes2$Index.SNP, T, F)
  })
  table(final$pgc3_prioritized2)

  final$pgc3_prioritized1_genes = sapply(1:nrow(final), function(i) {
    if(final$pgc3_prioritized1[i]) {
      paste0(sort(unique(pgcGenes[(final$rs[i] == pgcGenes$Index.SNP),"Ensembl.ID"])), collapse=",")
    } else {
      NA
    }
  })
  final$pgc3_prioritized2_genes = sapply(1:nrow(final), function(i) {
    if(final$pgc3_prioritized2[i]) {
      paste0(sort(unique(pgcGenes[(final$rs[i] == pgcGenes$Index.SNP),"Ensembl.ID"])), collapse=",")
    } else {
      NA
    }
  })
  final$peak_gene = paste0(final$peak, "_", final$gene)

  #####
  # Distance between index SNPs or LD buddies and the distance from peak edge (it should be 0 if PEAK_PADDING == 0)
  final$trait_gene_rs_ctype = paste0(final$trait_gene_rs, "_", final$CellType_general)
  final$trait_gene_rs_ctype_agePeriod = paste0(final$trait_gene_rs, "_", final$CellType_general, "_", final$AgePeriod)
  final$trait_gene_rs_ctype_sampleId = paste0(final$trait_gene_rs, "_", final$CellType_general, "_", final$SampleID)
  final$trait_gene_rs_ctype_agePeriod_sampleId = paste0(final$trait_gene_rs, "_", final$CellType_general, "_", final$AgePeriod, "_", final$SampleID)
  final$replications = sapply(final$trait_gene_rs_ctype_agePeriod, function(id) { sum(final$trait_gene_rs_ctype_agePeriod==id) })
  final$trait_buddy_end_hg38 = lociTable_hg38[final$trait_buddy_rs,"end"]
  final = final[order(final$rs == final$trait_buddy_rs, decreasing=T),]
  
  #####
  # Calculated distances between rs SNP and target gene & number of skipped genes
  peakToGeneRegions = final[,c("chr", "trait_buddy_end_hg38", "TargetGeneTSS_hg38", "trait_gene_rs_ctype_agePeriod", "gene")]
  colnames(peakToGeneRegions) = c("chr", "start", "end", "name", "geneID")
  
  for(i in 1:nrow(peakToGeneRegions)) {
    if(i %% 10000 == 0){
      print(i)
    }
    if(peakToGeneRegions[i, "start"] > peakToGeneRegions[i, "end"]) { 
      x = peakToGeneRegions[i, "end"]
      peakToGeneRegions[i, "end"] = peakToGeneRegions[i, "start"]; 
      peakToGeneRegions[i, "start"] = x
    }
    if(peakToGeneRegions[i,"geneID"] %in% ensemblGeneInfo$Ensembl.Gene.ID) {
      peakToGeneRegions[i, "end"] = max(c(peakToGeneRegions[i, "end"], ensemblGeneInfo[match(peakToGeneRegions[i, "geneID"], ensemblGeneInfo$Ensembl.Gene.ID),"Transcript.start..bp."], 
                                          ensemblGeneInfo[match(peakToGeneRegions[i, "geneID"], ensemblGeneInfo$Ensembl.Gene.ID),"Transcript.start..bp."] ))
    }
  }
  
  ensemblGeneInfoOnlyExpressedGr = makeGRangesFromDataFrame(ensemblGeneInfoShort[rownames(ensemblGeneInfoShort) %in% unique(unlist(expressedGenes)),], keep.extra.columns=T)
  peakToGeneRegionsGr = makeGRangesFromDataFrame(peakToGeneRegions, keep.extra.columns=T)
  hits = data.frame(findOverlaps(peakToGeneRegionsGr, ensemblGeneInfoOnlyExpressedGr))
  final$genesInBetween = sapply(1:nrow(final), function(i) { max(0, sum(hits$queryHits == i) - 1) })
  
  final$distance_rsToGene = sapply(1:nrow(final), function(i) { # fixme now
    snpPos = lociTable_hg38[final[i,"trait_buddy_rs"], "end"] # lociTable_hg38[final[i,"rs"], "end"]
    tss = final$TargetGeneTSS_hg38[i] 
    abs(snpPos-tss)
  })
  
  TX$abcSummary = paste0("ABC filtering: ", nrow(final), " E-P links for ", length(unique(final$gene)), " expressed genes and ", length(unique(final$peak)), " unique enhancers; ", length(unique(final$trait_gene_rs_ctype)), " trait-gene-loci-cell combinations; ", length(unique(final$trait_gene_rs)), " trait-gene-loci combinations; ", length(unique(final$trait_buddy_rs)), " loci")
  print(TX$abcSummary)
}

###################################################################################################################################################
##### FILTER OUT E-P LINKS USING ADDITINAL CRITERIA (ACCESSIBLE PEAKS, EXISTENCE OF E-P IN MIN. NUMBER OF SAMPLES, ETC)  ##########################

{
  # Some stats
  print(paste0("> unique cytogenic locations = ", length(unique(final$trait_buddy_rs)) ))
  print(paste0("> rsSNP -to- genes = ", nrow(final)))
  print(paste0("> rsSNP -to- genes = ", length(unique(final$rs))))
  print(paste0("> unique peaks -to- genes = ", length(unique(final$peak_gene))))
  print(paste0("> unique mapped genes = ", length(table(final$gene))))
  
  print(paste0("> All of EP links overlaping GWAS loci: ", nrow(final)))
  final = final[!duplicated(final$trait_gene_rs_ctype_agePeriod_sampleId),]
  print(paste0("> All of EP links with non-duplicated link for LD buddy (i.e. separately for index SNP and its buddy): ", nrow(final)))
  
  final = final[sapply(1:nrow(final), function(i) {final[i,"gene"]  %in% expressedGenes[[final[i,"CellType_general"]]] }),]
  print(paste0("> All of EP links with gene that is sufficiently expressed in cell type in which EP link was identified: ", nrow(final)))
  
  final = final[(final$peak %in% rownames(tmpEnv$initialVoomObj_filtered)),]
  print(paste0("> All of EP links (dtto) with sufficiently accessible peak (existence of peak in filtered matrix): ", nrow(final)))
  
  final = final[sapply(1:nrow(final), function(i) { qcPeakAnno[final$peak[i],final$CellType_general[i]] > 0.05 }),]
  print(paste0("> All of EP links (dtto) with sufficiently accessible peak in the given cell type (needs to be at least 5% from the bulk): ", nrow(final)))
}

###################################################################################################################################################
##### ABC-MAX FILTER (Table S8, Fig. 4c, Fig. S8-S17) #############################################################################################

{
  finalAbc = do.call("rbind.data.frame", lapply(setdiff(unique(final$CellType_general), "VLMC_Pericyte"), function(cellType) {
    print(cellType)
    df = final[which(cellType==final$CellType_general),]
    
    dfx = do.call("rbind.data.frame", lapply(unique(df$peak), function(peakId) {
      df3 = df[(df$peak == peakId),]
      df3$gene_biotype = ordered(df3$gene_biotype, levels=c("protein_coding", setdiff(unique(df3$gene_biotype), "protein_coding")))
      df3 = df3[order(df3$gene_biotype, -df3$ABC.Score),]
      df3 = df3[1:(min(MAX_LINKS_PER_PEAK, nrow(df3))),] # make change here if you want to change the allowed amount of E-P links per peak
      df3 = df3[!duplicated(df3$rs),]
      df3$order = 1:nrow(df3)
      df3
    }))
    dfx$combo2 = paste0(dfx$trait, "_", dfx$gene)
    
    df = do.call("rbind.data.frame", lapply(unique(dfx$trait), function(trait) {
      dfz = dfx[dfx$trait == trait,]
      
      do.call("rbind.data.frame", lapply(unique(dfz$combo2), function(combo2) {
        dfy = dfz[dfz$combo2 == combo2,]
        xx = dfy[order(dfy$ABC.Score, decreasing=T),][1,]
        xx$rs_all = paste0(unique(dfy$rs), collapse=",")
        xx$trait_buddy_rs_all = paste0(unique(dfy$trait_buddy_rs), collapse=",")
        xx$peak_all = paste0(unique(dfy$peak), collapse=",")
        xx$AgePeriod_all = paste0(unique(dfy$trait_buddy_rs), collapse=",")
        xx
      }))
    }))
    df = df[sapply(1:nrow(df), function(i) { qcPeakAnno[df$peak[i],df$CellType_general[i]] > 0.05 }),]
    df
  }))
  finalAbc$gene_name[finalAbc$gene_name==""] = unlist(sapply(finalAbc$gene[finalAbc$gene_name == ""], function(x) { ensemblInfo[match(x, ensemblInfo$Ensembl.Gene.ID),"Associated.Gene.Name"] }))
  finalAbc$cyt_loc = sapply(finalAbc$gene_name, function(x) get_location(x))
  
  finalAbc$start = as.numeric(finalAbc$start)
  finalAbc$end = as.numeric(finalAbc$end)
  finalAbc$trait_buddy_start = as.numeric(finalAbc$trait_buddy_start)
  finalAbc$trait_buddy_end = as.numeric(finalAbc$trait_buddy_end)
  finalAbc$plotRange = sapply(1:nrow(finalAbc), function(i) {
    paste0(finalAbc[i,"chr"], ":", min(finalAbc[i,"linked_start"], ensemblGeneInfo[finalAbc[i,"gene"],"Transcript.start..bp."]), "-", max(finalAbc[i,"linked_end"], ensemblGeneInfo[finalAbc[i,"gene"],"Transcript.end..bp."]))
  })
  
  # Printing Table S8: Enhancer-promoter links connected to GWAS loci as predicted by the ABC-MAX approach.
  print(paste0("ABC-MAX filtering: ", nrow(finalAbc), " E-P links for ", length(unique(finalAbc$gene)), " expressed genes and ", length(unique(finalAbc$peak)), " unique enhancers; ", length(unique(finalAbc$trait_gene_rs_ctype)), " trait-gene-loci-cell combinations; ", length(unique(finalAbc$trait_gene_rs)), " trait-gene-loci combinations; ", length(unique(finalAbc$trait_buddy_rs)), " loci"))
  exportAbcMax = finalAbc[,c("cyt_loc", "rs", "trait_buddy_rs", "r2", "trait", "CellType_general", "ABC.Score", "gene", "gene_name", "peak", "linked_chr", "linked_start", "linked_end")]
  colnames(exportAbcMax) = c("Cytogenic location", "rsID", "rsID (LD buddy)", "r2 (LD buddy)", "GWAS abbreviation", "Cell type", "ABC score", "Ensembl Gene ID", "Gene name", "Peak name", "Peak_chr", "Peak_start", "Peak_end")
  write.csv(exportAbcMax, row.names=F, file=file.path(ROOT, "outputs", "Table_S8.csv"))
  
  # Export all valid ABC EP links (not ABC-MAX!)
  exportAbc = final[,c("cyt_loc", "rs", "trait_buddy_rs", "r2", "trait", "CellType_general", "ABC.Score", "gene", "gene_name", "peak", "linked_chr", "linked_start", "linked_end")]
  colnames(exportAbc) = c("Cytogenic location", "rsID", "rsID (LD buddy)", "r2 (LD buddy)", "GWAS abbreviation", "Cell type", "ABC score", "Ensembl Gene ID", "Gene name", "Peak name", "Peak_chr", "Peak_start", "Peak_end")
  write.csv(exportAbc, row.names=F, file=file.path(ROOT, "outputs", "abcmax", "abc_filtered.csv"))
  abcResult = final
  dir.create(file.path(ROOT, "outputs", "abcmax", "genes"))
  
  for(trait in unique(abcResult$trait)) {
    traitDf = abcResult[(abcResult$trait == trait) & (abcResult$gene_name != ""),]
    genesDf = ensemblGeneInfo[match(traitDf$gene, rownames(ensemblGeneInfo)),]
    genesDf = genesDf[(genesDf$Gene.type == "protein_coding"),]
    traitDf = traitDf[traitDf$gene %in% genesDf$Ensembl.Gene.ID,]
    
    # Export genes to the file & save it to the list for further exploration
    for(ctype in unique(traitDf$CellType_general)) {
      genes = unique(traitDf[traitDf$CellType_general == ctype, "gene"])
      geneListForGsea[["ABC_MAX"]][[paste0(trait, "_", ctype)]] = genes
      write.table(genes, file.path(ROOT, "outputs", "abcmax", "genes", paste0(trait, "_", ctype , ".csv")), quote=F, row.names=F)
    }
    write.table(unique(traitDf$gene), file.path(outDir, "genes", paste0(trait, "_ALL.cvs")), quote=F, row.names=F)
  
    df = do.call("rbind.data.frame", lapply(unique(traitDf$gene_name), function(gene) {
      df = (traitDf[traitDf$gene_name == gene,])
      df = df[order(df$ABC.Score, decreasing=T),]
      df = df[!duplicated(df$CellType_general),]
      rownames(df) = df$CellType_general
      xx = sapply(unique(traitDf$CellType_general), function(ctype) {
        df[ctype,"ABC.Score"]
      })
      xx
    }))
    dim(df)
    rownames(df) = unique(traitDf$gene_name)
    colnames(df) = unique(traitDf$CellType_general)
    df$gene_name = unique(traitDf$gene_name)
    df[is.na(df)] = 0
    
    df$dist = sapply(rownames(df), function(gene_name) {
      mean(na.omit(traitDf[(traitDf$gene_name == gene_name),"distance_rsToGene"]))
    })
    df$closestGene = sapply(rownames(df), function(gene_name) {
      ifelse(min(traitDf[(traitDf$gene_name == gene_name),"genesInBetween"]) == 0, T, F)
    })
    
    # Perform hierarchical clustering
    dist_mat = dist(df[,1:length(unique(traitDf$CellType_general))])
    hc = hclust(dist_mat)
    ordered_df = df[order.dendrogram(as.dendrogram(hc)), ]
    ordered_df$gene_name = ordered(ordered_df$gene_name, levels=ordered_df$gene_name)
    ordered_df_melt = melt(ordered_df[,1:(length(unique(traitDf$CellType_general))+1)], id.vars = c("gene_name"))
    
    ordered_df_melt$gene_name = ordered(ordered_df_melt$gene_name, levels=ordered_df$gene_name)
    ordered_df_melt$variable = ordered(ordered_df_melt$variable, levels=CELL_TYPES)
    ordered_df_melt$value = log(ordered_df_melt$value+1)
    
    colors = c("white", colorRampPalette(c("lightsalmon", "salmon", "darkred"))(100))
    hmap = ggplot(data = ordered_df_melt, aes(y=variable, x=gene_name, fill=value)) + geom_tile() + scale_x_discrete(position = "top", guide = guide_axis(angle = 90)) + 
      scale_fill_gradientn(colors = colors, values = rescale(c(0, 0.02, 0.2, 1))) + theme_bw() + labs(x="", y="", fill="ABC score", title=trait) + 
      theme(axis.ticks.x=element_blank(), axis.title.y=element_blank(), title=element_blank(), axis.ticks.y=element_blank()) # scale_fill_material(colors = colors, values = rescale(c(0, 0.02, 0.2, 1)))
    hmap
    mpdf(paste0("fig_S8-S17_", trait), outDir=file.path(ROOT, "outputs", "abcmax"), hmap, height=2, width=floor(0.2 * length(unique(traitDf$gene_name))+1.5)); print(hmap); dev.off()
    
    colors = c("white", colorRampPalette(c("white", "salmon", "darkred"))(1000))
    ordered_df$temp = "Distance"
    ordered_df$distLog = log10(ordered_df$dist)
    hmap2 = ggplot(data = ordered_df, aes(x=gene_name, y=temp, fill=distLog)) + geom_tile() + scale_x_discrete(guide = guide_axis(angle = 90)) + 
      theme_bw() + labs(x="", y="", fill="Distance", title=trait) + scale_fill_material("light-blue") + 
      theme(axis.ticks.x=element_blank(), axis.title.y=element_blank(), title=element_blank(), axis.ticks.y=element_blank())
    mpdf(paste0("fig_S8-S17_", trait, "_distance"), outDir=file.path(ROOT, "outputs", "abcmax"), hmap, height=2, width=floor(0.2 * length(unique(traitDf$gene_name))+1.5)); print(hmap2); dev.off()
    
    colors = c("white", "black")
    ordered_df$temp = "ClosestGene"
    hmap3 = ggplot(data = ordered_df, aes(x=gene_name, y=temp, fill=closestGene)) + geom_tile() + scale_x_discrete(guide = guide_axis(angle = 90)) + 
      scale_fill_manual(values=c("white", "black")) + theme_bw() + labs(x="", y="", fill="Distance", title=trait) +
      theme(axis.ticks.x=element_blank(), axis.title.y=element_blank(), title=element_blank(), axis.ticks.y=element_blank())
    mpdf(paste0("fig_S8-S17_", trait, "_closestGeneIndicator"), outDir=file.path(ROOT, "outputs", "abcmax"), hmap, height=2, width=floor(0.2 * length(unique(traitDf$gene_name))+1.5)); print(hmap3); dev.off()
  }
  
  TX[["abc_abcMax"]] = paste0("Then, to narrow down the focus on GWAS traits, we decided to keep only E-P links overlapping GWAS index SNPs or their LD buddies (R2≥0.8), resulting in ", nrow(final), "E-P links. To limit the number of causal genes per GWAS locus, we applied the ABC-MAX approach [33828297], allowing only one E-P link per peak, resulting in ", nrow(finalAbc), " E-Ps.")
  print(TX[["abc_abcMax"]])
  
  TX[["abcmax_results"]] = paste0("Across all cell types, we nominated a total of ", length(unique(finalAbc$gene)), " unique genes that were putatively associated with ", round(length(unique(finalAbc$rs)) / uniqIndexSnps * 100, 0), "% of analyzed loci.")
  print(TX[["abcmax_results"]])
  
  TX[["abcmax_ms"]] = paste0("Our approach yielded predictions for ", length(unique(finalAbc[(finalAbc$trait == "MS"),"rs"])), " loci associated with multiple sclerosis (MS), nominating a total of ", length(unique(finalAbc[(finalAbc$trait == "MS"),"gene"])), " unique genes (Table traitsByCtype).")
  print(TX[["abcmax_ms"]])
}

###################################################################################################################################################
##### SUMMARIZE RESULTS (E-P OVERLAPPING GWAS) BY GENES ###########################################################################################

{
  finalx = finalAbc[order(finalAbc$cyt_loc),]
  
  dir.create(file.path(ROOT, "outputs", "abcmax", "indiv_gene_tables"))
  finalx = finalAbc[order(finalAbc$cyt_loc),]
  
  finalx = do.call("rbind", lapply(unique(finalx$gene), function(gene) {
    old_row = finalx[(finalx$gene==gene),]
    old_row = old_row[order(old_row$ABC.Score, decreasing=T),]
    write.csv(old_row, file=file.path(ROOT, "outputs", "abcmax", "indiv_gene_tables", paste0(gene, ".csv")), row.names=F)
    
    row = old_row[1,c("cyt_loc", "trait_buddy_rs", "trait", "gene", "gene_name")]
    row$trait_buddy_rs = paste0(unique(old_row$trait_buddy_rs),collapse=",")
    row$rs = paste0(unique(old_row$rs),collapse=",")
    row$no_buddies = length(which(finalx$gene == gene))
    row$no_traits = length(unique(old_row$trait))
    row$no_peaks = length(unique(old_row$peak))
    row$trait = paste0(sort(unique(old_row$trait)),collapse=",")
    row$cyt_loc = paste0(unique(old_row$cyt_loc),collapse=",")
    row$peak = paste0(unique(old_row$peak),collapse=",")
    row$cellType_note = paste0(sapply(unique(old_row$CellType_general)[order(table(unique(old_row$CellType_general)), decreasing=T)], function(ctype) paste0(ctype, ":", length(old_row[old_row$CellType_general==ctype,"trait"]))), collapse="; ")
    row$cellType_note2 = paste0(sapply(unique(old_row$trait), function(trt) {
      paste0(trt, ": ", paste0(old_row[old_row$trait == trt,"CellType_general"], collapse=" > "))
    }), collapse=" ;; ")
    row
  }))
  
  finalByGenes = finalx
  
  ld_main_and_buddies_df$combo = paste0(ld_main_and_buddies_df$buddy_rs, "_", ld_main_and_buddies_df$chr, "_", ld_main_and_buddies_df$start, "_", ld_main_and_buddies_df$end)
  comboLength = length(unique(ld_main_and_buddies_df$combo))
  
  exportByGenes = finalByGenes[,c("cyt_loc", "rs", "trait_buddy_rs", "trait", "gene", "gene_name", "no_traits", "no_peaks", "no_buddies", "peak")]
  colnames(exportByGenes) = c("Cytogenic location", "rsIDs", "rsID (LD buddies)", "GWAS abbreviations", "Ensembl Gene ID", "Gene name", "Number of overlapped GWAS traits", "Number of overlapped peaks", "Number of overlapped rsIDs", "Overlapped peaks")
  write.csv(exportByGenes, row.names=F, file=file.path(ROOT, "outputs", "abcmax", "table_exportByGenes.csv"))
  
  TX$multiDiseases = paste0("Finally, our approach identified a set of ", sum(finalByGenes$no_traits>=2), " genes associated with multiple diseases, ...")
  print(TX$multiDiseases)
}

###################################################################################################################################################
##### GSEA ANALYSIS FOR MS (Table S10) ############################################################################################################

{
  gseaResultsDf = read.csv(file.path(ROOT, "inputs", "Table_S10.csv"), sep="\t")
  gseaResultsDf$minusLogP = -log10(gseaResultsDf$pval)
  gseaResultsDfOrig = gseaResultsDf
  
  #####
  # Generate Table S10, i.e. scatterplot for MS and all cell types (only those with at least one BH-significant in top 10 most enriched are shown)
  dfx = gseaResultsDfOrig[(startsWith(gseaResultsDfOrig$Set, prefix="MS_")),]
  dfx$BH_AdjP = p.adjust(dfx$pval, method="BH")
  dfx = dfx[(dfx$name_full %in% unique(dfx$name_full)[1:10]),]
  dfx$name_full = ordered(dfx$name_full, levels=unique(dfx$name_full))
  dfx$size = ifelse(dfx$BH_AdjP < 0.05, 2, 1)
  dfx$Set = gsub("^MS_", "", dfx$Set)
  gseaMsDf = dfx[,c("name_full", "Set", "pval", "BH_AdjP", "odds.ratio", "Jaccard", "intersection", "union")]
  colnames(gseaMsDf) = c("Pathway name", "Cell type", "P-value", "BH-corrected P-value", "Odds ratio", "Jaccard", "Intersection", "Union")
  write.csv(gseaMsDf, row.names=F, file=file.path(ROOT, "outputs", "Table_S10.csv"))
}

###################################################################################################################################################
##### ALS-RELATED ANALYSIS ########################################################################################################################

{
  alsNominated = unique(finalAbc[finalAbc$trait == "ALS","rs"])
  alsTotal = unique(lociTable_df[lociTable_df$trait == "ALS","rs"])
  TX$als = paste0("Our approach was able to nominate causal genes for ", length(alsNominated), " out of ", length(alsTotal), " loci (", round(length(alsNominated) / length(alsTotal) * 100, 1), "%)")
  print(TX$als)
  
  df = finalAbc[(finalAbc$trait == "ALS"),]
  ALS_LIU_OMIM = file.path(ROOT, "inputs", "liu_als_omim.csv")
  ALS_LIU_GWAS = file.path(ROOT, "inputs", "liu_als_gwas.csv")
  ALS_VAN_RHEENEN = file.path(ROOT, "inputs", "als_van_Rheenen.csv")
  als_omim = read.csv(ALS_LIU_OMIM, sep="\t")
  als_gwas = read.csv(ALS_LIU_GWAS, sep="\t")
  als_gwas_vanRheenen = read.csv(ALS_VAN_RHEENEN, sep="\t")
  als_gwas_vanRheenenList = unique(unlist(lapply(als_gwas_vanRheenen$Prioritized_gene, function(genes) { as.character(sapply(unlist(strsplit(genes, ",")[[1]]), trimws)) })))
  df$existSupport = df$gene_name %in% c(als_omim$Gene.symbol, als_gwas$Gene.symbol)
  
  isect_omim = intersect(unique(df$gene_name), als_omim$Gene.symbol)
  isect_gwas = intersect(unique(df$gene_name), als_gwas$Gene.symbol)
  isect_singleGwas = intersect(unique(df$gene_name), als_gwas_vanRheenenList)
  
  TX$als_omim = paste0("> Overlap between ALS OMIM (", nrow(als_omim) , " genes) and our predictions (", length(unique(df$gene)), " genes) is in ", length(isect_omim), " genes (", length(unique(df[df$gene_name %in% isect_omim,"rs"])), " loci), i.e. ", paste0(sort(isect_omim), collapse=","))
  print(TX$als_omim)
  TX$als_gwas = paste0("> Overlap between ALS GWAS (", nrow(als_gwas) , " genes) and our predictions (", length(unique(df$gene)), " genes) is in ", length(isect_gwas), " genes (", length(unique(df[df$gene_name %in% isect_gwas,"rs"])), " loci), i.e. ", paste0(sort(isect_gwas), collapse=","))
  print(TX$als_gwas)
  TX$als_gwas2 = paste0("> Overlap between ALS GWAS (", length(als_gwas_vanRheenen) , " genes) and our predictions (", length(unique(df$gene)), " genes) is in ", length(isect_singleGwas), " genes (", length(unique(df[df$gene_name %in% isect_singleGwas,"rs"])), " loci), i.e. ", paste0(sort(isect_singleGwas), collapse=","))
  print(TX$als_gwas2)
}

###################################################################################################################################################
##### SCZ-RELATED ANALYSIS ########################################################################################################################

{
  finalAbc$pgc3_prioritized1_matchOur = sapply(1:nrow(finalAbc), function(i) {
    if(finalAbc$pgc3_prioritized1[i]) {
      #ifelse(length(intersect(unique(pgcGenes[(final$rs[i] == pgcGenes$Index.SNP),"Ensembl.ID"]), unique(finalAbc[(finalAbc$trait == "SCZ") & (finalAbc$rs == final$rs[i]),"gene"]))) > 0, T, F)
      ifelse(length(intersect(unique(pgcGenes[(finalAbc$rs[i] == pgcGenes$Index.SNP),"Ensembl.ID"]), unique(finalAbc[(finalAbc$trait == "SCZ") & (finalAbc$rs == finalAbc$rs[i]),"gene"]))) > 0, T, F)
    } else {
      NA
    }
  })
  table(finalAbc$pgc3_prioritized1_matchOur)
  
  finalAbc$pgc3_prioritized2_matchOur = sapply(1:nrow(finalAbc), function(i) {
    if(finalAbc$pgc3_prioritized2[i]) {
      ifelse(length(intersect(unique(pgcGenes2[(finalAbc$rs[i] == pgcGenes2$Index.SNP),"Ensembl.ID"]), unique(finalAbc[(finalAbc$trait == "SCZ") & (finalAbc$rs == finalAbc$rs[i]),"gene"]))) > 0, T, F)
    } else {
      NA
    }
  })
  table(finalAbc$pgc3_prioritized2_matchOur)
  
  #### Comparison with prioritized genes from SCZ PGC3 GWAS (there are two versions of prioritized genes from PGC and we consider using only protein coding-vs-all)
  isectPgc3prioritized_oursSczPrioritized_pc = intersect(unique(finalAbc[(finalAbc$trait=="SCZ") & (finalAbc$gene_biotype == "protein_coding"),"rs"]), pgcGenes[pgcGenes$gene_biotype=="protein_coding", "Index.SNP"])
  isectPgc3prioritized_oursSczPrioritized_sameOrNot_pc = sapply(isectPgc3prioritized_oursSczPrioritized_pc, function(rs) { any(finalAbc[finalAbc$rs == rs,"pgc3_prioritized1_matchOur"]) })
  TX$pgc3_overlap1_pc = paste0("> Out of ", length(unique(ld_buddies_df[ld_buddies_df$trait=="SCZ","rs"])), " SCZ loci, PGC3(v1) made predictions of protein coding target gene for ", length(unique(pgcGenes[(pgcGenes$gene_biotype == "protein_coding"), "Index.SNP"])), " loci while our study made a prediction for ", length(unique(finalAbc[(finalAbc$trait == "SCZ") & (finalAbc$gene_biotype == "protein_coding"), "rs"])), " loci. For ", length(isectPgc3prioritized_oursSczPrioritized_pc), " loci with available predictions in both studies, the same gene is prioritized in ", sum(isectPgc3prioritized_oursSczPrioritized_sameOrNot_pc), " cases.")
  print(TX$pgc3_overlap1_pc)
  
  isectPgc3prioritized_oursSczPrioritized = intersect(unique(finalAbc[(finalAbc$trait=="SCZ"),"rs"]), pgcGenes$Index.SNP)
  isectPgc3prioritized_oursSczPrioritized_sameOrNot = sapply(isectPgc3prioritized_oursSczPrioritized, function(rs) { any(finalAbc[finalAbc$rs == rs,"pgc3_prioritized1_matchOur"]) })
  TX$pgc3_overlap1 = paste0("> Out of ", length(unique(ld_buddies_df[ld_buddies_df$trait=="SCZ","rs"])), " SCZ loci, PGC3(v1) made predictions of all (not just protein coding) target gene for ", length(unique(pgcGenes$Index.SNP)), " loci while our study made a prediction for ", length(unique(finalAbc[(finalAbc$trait == "SCZ"), "rs"])), " loci. For ", length(isectPgc3prioritized_oursSczPrioritized), " loci with available predictions in both studies, the same gene is prioritized in ", sum(isectPgc3prioritized_oursSczPrioritized_sameOrNot), " cases.")
  print(TX$pgc3_overlap1)
  
  isectPgc3prioritized_oursSczPrioritized_pc = intersect(unique(finalAbc[(finalAbc$trait=="SCZ") & (finalAbc$gene_biotype == "protein_coding"),"rs"]), pgcGenes2[pgcGenes2$gene_biotype=="protein_coding", "Index.SNP"])
  isectPgc3prioritized_oursSczPrioritized_sameOrNot_pc = sapply(isectPgc3prioritized_oursSczPrioritized_pc, function(rs) { any(finalAbc[finalAbc$rs == rs,"pgc3_prioritized2_matchOur"]) })
  TX$pgc3_overlap2_pc = paste0("> Out of ", length(unique(ld_buddies_df[ld_buddies_df$trait=="SCZ","rs"])), " SCZ loci, PGC3(v2) made predictions of protein coding target gene for ", length(unique(pgcGenes2[(pgcGenes2$gene_biotype == "protein_coding"), "Index.SNP"])), " loci while our study made a prediction for ", length(unique(finalAbc[(finalAbc$trait == "SCZ") & (finalAbc$gene_biotype == "protein_coding"), "rs"])), " loci. For ", length(isectPgc3prioritized_oursSczPrioritized_pc), " loci with available predictions in both studies, the same gene is prioritized in ", sum(isectPgc3prioritized_oursSczPrioritized_sameOrNot_pc), " cases.")
  print(TX$pgc3_overlap2_pc)
  
  isectPgc3prioritized_oursSczPrioritized2 = intersect(unique(finalAbc[(finalAbc$trait=="SCZ"),"rs"]), pgcGenes2$Index.SNP)
  isectPgc3prioritized_oursSczPrioritized_sameOrNot2 = sapply(isectPgc3prioritized_oursSczPrioritized2, function(rs) { any(finalAbc[finalAbc$rs == rs,"pgc3_prioritized2_matchOur"]) })
  TX$pgc3_overlap2 = paste0("> Out of ", length(unique(ld_buddies_df[ld_buddies_df$trait=="SCZ","rs"])), " SCZ loci, PGC3(v2) made predictions of all (not just protein coding) target gene for ", length(unique(pgcGenes2$Index.SNP)), " loci while our study made a prediction for ", length(unique(finalAbc[(finalAbc$trait == "SCZ"), "rs"])), " loci. For ", length(isectPgc3prioritized_oursSczPrioritized2), " loci with available predictions in both studies, the same gene is prioritized in ", sum(isectPgc3prioritized_oursSczPrioritized_sameOrNot2), " cases.")
  print(TX$pgc3_overlap2)
  
  TX$abcmax_scz_epDistance = paste0("For example, ABC-MAX made predictions for ", length(unique(finalAbc[(finalAbc$trait=="SCZ"),"rs"])), " SCZ loci, nominating ", length(unique(finalAbc[(finalAbc$trait=="SCZ"),"gene"])), " unique genes.")
  print(TX$abcmax_scz_epDistance)
  TX$abcmax_ms_epDistance = paste0("For example, ABC-MAX made predictions for ", length(unique(finalAbc[(finalAbc$trait=="MS"),"rs"])), " MS loci, nominating ", length(unique(finalAbc[(finalAbc$trait=="MS"),"gene"])), " unique genes.")
  print(TX$abcmax_ms_epDistance)
  
  nonConflictingRs = names(isectPgc3prioritized_oursSczPrioritized_sameOrNot)[isectPgc3prioritized_oursSczPrioritized_sameOrNot == T] # rsID with shared EP
  conflictingRs = names(isectPgc3prioritized_oursSczPrioritized_sameOrNot)[isectPgc3prioritized_oursSczPrioritized_sameOrNot == F]    # rsID with non-shared EP
  conflictingGenesPgc = pgcGenes2[(pgcGenes2$Index.SNP %in% conflictingRs),]
  conflictingGenesPgc$gene_rs = paste0(conflictingGenesPgc$Ensembl.ID, "_SCZ_", conflictingGenesPgc$Index.SNP)
  existingInAbc = final[final$trait_gene_rs %in% conflictingGenesPgc$gene_rs,"trait_gene_rs"]
  print(paste0("To further demonstrate the accuracy of our approach, we investigated ", length(unique(finalAbc[(finalAbc$trait=="SCZ"),"rs"])), " SCZ loci for which we have predicted E-P(ABC-MAX) links, covering ", round(length(unique(finalAbc[(finalAbc$trait=="SCZ"),"rs"])) / length(unique((ld_buddies_df[(ld_buddies_df$trait=="SCZ"),"rs"]))) * 100, 1), "% of all SCZ loci"))
  print(paste0("When comparing our predicted causal genes with those prioritized by SMR/FINEMAP in the latest SCZ GWAS study, we found agreement for ", sum(isectPgc3prioritized_oursSczPrioritized_sameOrNot2), " (", round(sum(isectPgc3prioritized_oursSczPrioritized_sameOrNot2==T) / length(isectPgc3prioritized_oursSczPrioritized_sameOrNot2) * 100, 1), "%) of the ", length(isectPgc3prioritized_oursSczPrioritized_sameOrNot2), " loci that are covered by both studies. Notably, for ", length(unique(existingInAbc)), " out of ", sum(isectPgc3prioritized_oursSczPrioritized_sameOrNot2==F), " loci with different predicted causal gene, we confirmed the existence of E-P(ABC) (but not E-P(ABC-MAX)) link to the same gene as predicted by SMR/FINEMAP, indicating the existence of regulatory link in our data, even if it's not the top-ranked one."))
  
  # For SynGO: Run online with inputs generated here (all SCZ genes) vs background of all genes as exported here ;; then, save ZIP with results, extract it and convert XLS to CSV 
  write.csv(unique(unique(unlist(expressedGenes))), file=file.path(ROOT, "inputs", "SynGO_background.txt"), row.names=F, quote=F)
  
  SYNGO_SCZ_RESULTS = file.path(ROOT, "inputs", "syngo_scz_precalculated.tsv")
  syngoDf = read.csv(SYNGO_SCZ_RESULTS, sep="\t")
  syngoDf = syngoDf[(syngoDf$GO.domain == "BP"),]
  syngoGO_terms = unique(syngoDf$GO.term.ID)[startsWith(unique(syngoDf$GO.term.ID), prefix="GO")]
  TX$syngo = paste0("These prioritized genes are associated with ", length(syngoGO_terms), " unique biological terms within the synaptic hierarchy (Table SynGO), including several specific pathways and processes.")
  print(TX$syngo)
}

###################################################################################################################################################
##### SUMMARY OF GENE IMPLICATIONS BY TRAIT :: TABLE S9, FIG. 3F ##################################################################################

{
  outDir = "~/Desktop/3dg/abc_max/" # abcResultsFolder$ABC_MAX
  
  traitsByCtype = do.call("rbind.data.frame", lapply(TRAITS, function(trait) {
    table(finalAbc[finalAbc$trait == trait,"CellType_general"])[CELL_TYPES]
  }))
  colnames(traitsByCtype) = CELL_TYPES
  rownames(traitsByCtype) = TRAITS
  traitsByCtype$sumUnique = sapply(rownames(traitsByCtype), function(trait) { 
    length(unique(unlist(sapply(colnames(traitsByCtype), function(ctype) { finalAbc[(finalAbc$trait == trait) & (finalAbc$CellType_general == ctype),"gene"] }))))
  })
  traitsByCtypeDf = rbind.data.frame(traitsByCtype, sapply(colnames(traitsByCtype), function(ctype) { length(finalAbc[(finalAbc$CellType_general == ctype),"gene"]) }))
  rownames(traitsByCtypeDf)[nrow(traitsByCtypeDf)] = "SumNonUnique"
  write.csv(traitsByCtypeDf, file=file.path(ROOT, "outputs", "Table_S9_b.csv"))
  
  ####
  # Summary of covered loci (loci can be represented by multiple genes as we allow 2 links per loci in each cell type) in each cell type & trait
  traitsByCtypeAlt = do.call("rbind.data.frame", lapply(TRAITS, function(trait) {
    sapply(CELL_TYPES, function(ctype) length(unique(finalAbc[(finalAbc$CellType_general == ctype) & (finalAbc$trait == trait),"rs"])) )
  }))
  colnames(traitsByCtypeAlt) = CELL_TYPES
  rownames(traitsByCtypeAlt) = TRAITS
  traitsByCtypeAlt$sumUnique = sapply(rownames(traitsByCtypeAlt), function(trait) { 
    length(unique(unlist(sapply(colnames(traitsByCtypeAlt), function(ctype) { finalAbc[(finalAbc$trait == trait) & (finalAbc$CellType_general == ctype),"rs"] }))))
  })
  traitsByCtypeAlt$lociCoverage = sapply(TRAITS, function(trait) (traitsByCtypeAlt[trait,"sumUnique"] / nrow(lociTable_df[(lociTable_df$trait == trait),])) )
  traitsByCtypeDfAlt = rbind.data.frame(traitsByCtypeAlt, sapply(colnames(traitsByCtypeAlt), function(ctype) { length(unique(finalAbc[(finalAbc$CellType_general == ctype),"rs"])) }))
  rownames(traitsByCtypeDfAlt)[nrow(traitsByCtypeDfAlt)] = "SumNonUnique"
  write.csv(traitsByCtypeDfAlt, file=file.path(ROOT, "outputs", "Table_S9_a.csv"))
  
  ####
  # Figure 3f: Number of explained loci per cell type across all investigated GWAS traits
  df = data.frame(t(traitsByCtypeDfAlt[nrow(traitsByCtypeDfAlt),1:(length(traitsByCtypeDfAlt)-2)]))
  df$`Cell type` = rownames(df)
  colnames(df) = c("Explained loci", "Cell type")
  df$`Cell type` = ordered(df$`Cell type`, levels=df$`Cell type`[order(df$`Explained loci`, decreasing=T)] )
  explainedLociByCellPlot = ggplot(df, aes(x = `Cell type`, y = `Explained loci`)) + geom_bar(stat = "identity", color = "white", width=1) +  theme_classic() + xlab("Cell type") + 
    ylab("A total number of explained loci") + theme(panel.spacing=unit(1, "lines")) + theme(strip.background=element_rect(fill="#eeeeee",color="#eeeeee")) + 
    theme(strip.text=element_text(face="bold"), axis.text.x = element_text(angle = 45)) + scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0.5)) #+ facet_wrap(~ CellType, ncol=3)
  mpdf("Fig_3f", outDir=file.path(ROOT, "outputs"), height=4, width=4); print(explainedLociByCellPlot); dev.off()
}

###################################################################################################################################################
##### STATISTICS FOR E-P LINKS FROM ABC :: FIG. 3D-E ##############################################################################################

{
  # Mean / median distance for EP links
  medianEpDistance = median(na.omit(finalAbc$distance_rsToGene))
  meanEpDistance = mean(na.omit(finalAbc$distance_rsToGene))
  percUnder100kb = round(sum(na.omit(finalAbc$distance_rsToGene)<1e5) / length(na.omit(finalAbc$distance_rsToGene)), 2) * 100
  #sapply(unique(finalAbc$CellType_general), function(ctype) { round(mean(na.omit(finalAbc[finalAbc$CellType_general == ctype, "distance_rsToGene"]) / 10^3, 0)) })
  TX[["abcmax_epDistance"]] = paste0("The E-PABC-MAX distances varied significantly, with a median of ", round(medianEpDistance/1e3, 0), "kb, and the majority (", percUnder100kb, "%) were less than 100kb (Fig. X).")
  print(TX[["abcmax_epDistance"]])
  
  ####
  # Figure 3d: Histogram of the E-PABC-MAX distances between the GWAS risk variant and TSS of the putatively regulated gene
  distRsPeakLogPlot = ggplot(finalAbc, aes(x = distance_rsToGene)) + geom_density() + labs(x="E-P link distance (log10 scale)") + scale_x_continuous(trans='log10') + theme_bw() + 
    geom_vline(xintercept=meanEpDistance, linetype="dashed", color = "blue") + geom_vline(xintercept=medianEpDistance, linetype="dashed", color = "red") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)
  mpdf("Fig_3d", outDir=file.path(ROOT, "outputs"), width=3, height=3); print(distRsPeakLogPlot); dev.off()
  
  ####
  # Figure 3e: Histogram of the number of genes ‘skipped’ by E-PABC-MAX to reach the putatively regulated gene
  medianGenesInBetween = median(na.omit(finalAbc$genesInBetween)) # Distribution of the number of "skipped" genes by E-P(ABC-MAX) links 
  meanGenesInBetween = mean(na.omit(finalAbc$genesInBetween))
  
  df = finalAbc
  df[df$genesInBetween > 49,"genesInBetween"] = 50
  genesInBetweenHistogram = ggplot(df, aes(x=genesInBetween)) + geom_histogram() + labs(x="Genes skipped by E-P links") + theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1) + geom_vline(xintercept=medianGenesInBetween, linetype="dashed", color = "red") + 
    geom_vline(xintercept=meanGenesInBetween, linetype="dashed", color = "blue")
  mpdf("Fig_3e", outDir=file.path(ROOT, "outputs"), width=3, height=3); print(genesInBetweenHistogram); dev.off()
  
  skipped10orMoreGenes = round((1 - cumsum(table(finalAbc$genesInBetween) / nrow(finalAbc))[9]) * 100, 2)
  TX[["abcmax_genesInBetween"]] = paste0("Still, the vast majority of E-P(ABC-MAX) links (", round((1 - cumsum(table(finalAbc$genesInBetween) / nrow(finalAbc))[1]) * 100, 1), ") did not connect to the nearest gene, with more than half of them (", skipped10orMoreGenes, "%) having 10 or more genes in between.")
  print(TX[["abcmax_genesInBetween"]])
}

###################################################################################################################################################
##### MUTATIONAL CONSTRAINTS FROM GNOMAD :: FIG. 3B ###############################################################################################

{
  # Load mutational score from gnomAD
  mutScore = read.csv(GNOMAD_MUT_SCORE, sep="\t")
  mutScoreGr = makeGRangesFromDataFrame(mutScore, keep.extra.columns=T)
  
  # Calculate overlap of mutational score bins (1kb regions) and peaks
  qcPeakAnnoGr = makeGRangesFromDataFrame(qcPeakAnno, keep.extra.columns=T)
  hits = findOverlaps(qcPeakAnnoGr, mutScoreGr)
  hitsLong = findOverlapPairs(qcPeakAnnoGr, mutScoreGr)
  hitsCount = GenomicRanges::countOverlaps(qcPeakAnnoGr, mutScoreGr)
  table(hitsCount)
  
  calcMutability = function(qcPeakAnnoGr, mutScoreGr, variable) {
    sapply(1:nrow(qcPeakAnno), function(i) {
      ii = hits[(hits@from == i),]@to
      if(length(ii) == 0) {
        NA
      } else if(length(ii) == 1) {
        mutScoreGr[ii,]@elementMetadata[[variable]]
      } else {
        overlapWidth = sapply(ii, function(iii) {
          GenomicRanges::intersect(qcPeakAnnoGr[i,], mutScoreGr[iii,])@ranges@width
        })
        if(length(overlapWidth) == 0) {
          NA
        } else {
          sum(overlapWidth / sum(overlapWidth) * mutScoreGr[ii,]$z)
        }
      }
    })
  }
  
  qcPeakAnno$myZ = sapply(1:nrow(qcPeakAnno), function(i) {
    xx = unlist(qcPeakAnno$z[i], recursive=F)
    if(is.list(xx)) {
      data.frame(xx)$z
    } else {
      xx
    }
  })
  
  qcPeakAnno$z = calcMutability(qcPeakAnnoGr, mutScoreGr, "z")
  qcPeakAnno$observed = calcMutability(qcPeakAnnoGr, mutScoreGr, "observed")
  qcPeakAnno$expected = calcMutability(qcPeakAnnoGr, mutScoreGr, "expected")
  qcPeakAnno$oe = qcPeakAnno$observed / qcPeakAnno$expected
  
  # Textual sumary of differences in Z-score for all peak / peaks participating in EP links / peaks participating in EP links valid for GWAS loci
  a0 = qcPeakAnno[,"z"]
  print(paste0("> All peaks (n=", length(a0), "); mean=", round(mean(na.omit(a0)), 2), "; median=", round(median(na.omit(a0)), 2)))
  a1 = qcPeakAnno[unique(unique(unlist(lapply(epLinks, function(x) unique(x$PeakID))))),"z"]
  print(paste0("> All peaks connected by at least one ABC E-P (n=", length(a1), "); mean=", round(mean(na.omit(a1)), 2), "; median=", round(median(na.omit(a1)), 2)))
  a3 = qcPeakAnno[unique(finalAbc$peak),"z"]
  print(paste0("> All peaks connected by ABC E-P to GWAS loci (n=", length(a3), "); mean=", round(mean(na.omit(a3)), 2), "; median=", round(median(na.omit(a3)), 2)))
  
  # Create df for plotting
  linksFreqSq = table((unlist(lapply(epLinks, function(x) unique(x$PeakID)))))
  frqList = list()
  for(num in sort(unique(linksFreqSq))) {
    peakNames = names(linksFreqSq)[which(linksFreqSq == num)]
    frqList[[num]] = list(
      "peakNames" = peakNames,
      "n" = length(peakNames),
      "n_naOmit" = length(na.omit(qcPeakAnno[peakNames,"z"])),
      "distToTSS_mean" = mean(abs(qcPeakAnno[peakNames,"distanceToTSS"])),
      "distToTSS_median" = median(abs(qcPeakAnno[peakNames,"distanceToTSS"])),
      "z_mean" = round(mean(na.omit(qcPeakAnno[peakNames,"z"])), 2),
      "z_median" = round(median(na.omit(qcPeakAnno[peakNames,"z"])), 2),
      "oe_mean" = round(mean(na.omit(qcPeakAnno[peakNames,"oe"])), 2),
      "oe_median" = round(median(na.omit(qcPeakAnno[peakNames,"oe"])), 2)
    )
    print(paste0("> All ABC peaks with ", num," E-P links in cell types (n=", length(peakNames), "); mean(z)=", round(mean(na.omit(qcPeakAnno[peakNames,"z"])), 2), "; median(z)=", round(median(na.omit(qcPeakAnno[peakNames,"z"])), 2), ";  mean(oe)=", round(mean(na.omit(qcPeakAnno[peakNames,"oe"])), 2), "; median(oe)=", round(median(na.omit(qcPeakAnno[peakNames,"oe"])), 2), "; mean_distToTSS=", round(mean(abs(qcPeakAnno[peakNames,"distanceToTSS"]))), "; median_distToTSS=", round(median(abs(qcPeakAnno[peakNames,"distanceToTSS"]))) ))
  }
  
  frqTable = do.call("rbind.data.frame", lapply(frqList, function(x) { x[2:9]} ))
  frqTable$freq = rownames(frqTable)
  frqRows = frqTable[frqTable$freq %in% c(2, 3, 4, 5),]
  frqTable[frqTable$freq==2,] = c(sum(frqRows$n), sum(frqRows$n_naOmit), round(mean(abs(qcPeakAnno[unlist(sapply(frqList[2:5], function(x) x$peakNames)),"distanceToTSS"]))), 
                                  round(median(abs(qcPeakAnno[unlist(sapply(frqList[2:5], function(x) x$peakNames)),"distanceToTSS"]))), 
                                  sum(frqRows$n / sum(frqRows$n) * frqRows$z_mean), sum(frqRows$n / sum(frqRows$n) * frqRows$z_median),
                                  sum(frqRows$n / sum(frqRows$n) * frqRows$oe_mean), sum(frqRows$n / sum(frqRows$n) * frqRows$oe_median), 2)
  frqRows = frqTable[frqTable$freq %in% c(6, 7, 8, 9),]
  frqTable[frqTable$freq==3,] = c(sum(frqRows$n), sum(frqRows$n_naOmit), round(mean(abs(qcPeakAnno[unlist(sapply(frqList[6:9], function(x) x$peakNames)),"distanceToTSS"]))), 
                                  round(median(abs(qcPeakAnno[unlist(sapply(frqList[6:9], function(x) x$peakNames)),"distanceToTSS"]))), 
                                  sum(frqRows$n / sum(frqRows$n) * frqRows$z_mean), sum(frqRows$n / sum(frqRows$n) * frqRows$z_median),
                                  sum(frqRows$n / sum(frqRows$n) * frqRows$oe_mean), sum(frqRows$n / sum(frqRows$n) * frqRows$oe_median), 3)
  frqTable = frqTable[frqTable$freq <= 3,]
  peaksNotInAbc = qcPeakAnno[rownames(qcPeakAnno) %in% unique(unlist(lapply(epLinks, function(x) unique(x$PeakID)))),]
  peaksNotInAbc = peaksNotInAbc[(peaksNotInAbc$annotation %in% c("Distal Intergenic")),]
  
  frqTable = rbind.data.frame(frqTable, 
                              c(nrow(peaksNotInAbc), length(na.omit(peaksNotInAbc[,"z"])), round(mean(abs(peaksNotInAbc[,"distanceToTSS"]))),
                                round(median(abs(peaksNotInAbc[,"distanceToTSS"]))), 
                                mean(na.omit(peaksNotInAbc[,"z"])),
                                median(na.omit(peaksNotInAbc[,"z"])), 
                                mean(na.omit(peaksNotInAbc[,"oe"])), median(na.omit(peaksNotInAbc[,"oe"])), 0))
  head(frqTable)
  frqTable$desc = c("1", "2-5", "6-9", "0")
  frqTable$desc = ordered(frqTable$desc, levels=c("0", "1", "2-5", "6-9"))
  
  df = data.frame(t(traitsByCtypeDfAlt[nrow(traitsByCtypeDfAlt),1:(length(traitsByCtypeDfAlt)-2)]))
  df$`Cell type` = rownames(df)
  colnames(df) = c("Explained loci", "Cell type")
  df$`Cell type` = ordered(df$`Cell type`, levels=df$`Cell type`[order(df$`Explained loci`, decreasing=T)] )
  
  # Fig. 3b: Mean of mutational constraint Z-score for peaks grouped by the number of cell types for which the peaks have at least one E-PABC link
  mutability_oePlot = ggplot(frqTable, aes(x = frqTable$desc, y = frqTable$oe_mean)) +
    geom_bar(stat="identity") +  theme_classic() + xlab("Participation of peak(ABC-MAX) in cell types") + 
    ylab("Mean OE of mutability of peak(ABC-MAX)") + theme(panel.spacing=unit(1, "lines")) +
    theme(strip.background=element_rect(fill="#eeeeee",color="#eeeeee")) + 
    theme(strip.text=element_text(face="bold"), axis.text.x = element_text(angle = 45)) +
    scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0.5)) #+ facet_wrap(~ CellType, ncol=3)
  mpdf("Fig. 3b", outDir=file.path(ROOT, "outputs"), height=4, width=4); print(mutability_oePlot); dev.off()
}

# --------------------------------------------
# End of Script
# --------------------------------------------
