library(Seurat)
library(Signac)
library(ggplot2)
library(harmony)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)




###################### MULTI-MODAL SEURAT OBJECT & MERGE   ###################### 


##-----------  Seurat object on a single sample

# Set the working directory to the folder including its raw data
sampleID = "example_1"
setwd(paste("../data", sampleID, sep = "/"))
getwd()

# Read the 10x hdf5 file 
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
metadata <- read.csv('per_barcode_metrics.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# Generate Seurat object - RNA-seq data
pbmc <- CreateSeuratObject(counts = inputdata.10x$`Gene Expression`, 
                           assay = 'RNA',
                           project = sampleID,
                           meta.data = metadata)
# Calculate perc.mt
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Add in the ATAC-seq data (only use peaks in standard chromosomes)
atac_counts <- inputdata.10x$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Prepare gene annotations for hg38
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "atac_fragments.tsv.gz"
pbmc[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

# Calculate nucleosome signal and TSS enrichment
DefaultAssay(pbmc) <- 'ATAC'
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)
pbmc$pct_reads_in_peaks <- pbmc$atac_peak_region_fragments / pbmc$atac_fragments * 100
pbmc$blacklist_fraction <- FractionCountsInRegion(pbmc, assay = 'ATAC', regions = blacklist_hg38_unified)


# Call peaks separately using MACS2
peaks <- CallPeaks(pbmc, 
                   assay = 'ATAC',
                   macs2.path = '/home/terez/anaconda3/bin/macs2')
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = blacklist_hg38_unified, 
                          invert = TRUE)

# Quantify counts for each peak
macs_count <- FeatureMatrix(fragments = Fragments(pbmc),
                            features = peaks,
                            cells = colnames(pbmc))

pbmc[['ATAC']] <- CreateChromatinAssay(
  counts = macs_count,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

# Save the sample-specific Seurat object 
obj_example_1 = pbmc
save(obj_example_1, file = paste('obj', sampleID, 'processed.rda', sep = '_'))


#------- Merge samples

obj_all <- merge(obj_example_1, y = c(obj_example_2, obj_example_3, obj_example_4), 
                 project = "postnatal_brain")

save(obj_all, file = './data_preprocessed/ALL_snMultiome.rda')
saveRDS(obj_all, './data_processed/ALL_snMultiome.rds')






###################### FILTERING & CLUSTERING 	###################### 


options(future.globals.maxSize = 50000 * 1024^2)

OBJ_ALL <- readRDS('./data_processed/ALL_snMultiome.rds')


##---------------------------------- Filter cells
# RNA-seq filtering: nCount_RNA, percent.mt
# ATAC-seq filtering: nCount_ATAC, nucleosome_signal, TSS.enrichment
FILT_ALL <- subset(x = OBJ_ALL,
		 subset = nCount_RNA > 2e2 & 
                 nCount_RNA < 5e4 &
                 nCount_ATAC > 2e2 &
                 nCount_ATAC < 1e5 &
                 percent.mt < 5 &
                 nucleosome_signal < 3 &
                 TSS.enrichment > 1)



##--------------------------------- Precessing and Clustering 
N.PC = 30; N.LSI = 10
ClusterData <- function(object, n.pc=N.PC, n.lsi = N.LSI){
  
  
  ##-----Filter peaks/genes which are detected in < 10 cells
  
  tmp <- Matrix::rowSums(object[['RNA']]@counts > 0)
  object[['RNA']] <- subset(object[['RNA']], features = names(which(tmp >= 10)))
  tmp <- Matrix::rowSums(object[['ATAC']]@counts > 0)
  object[['ATAC']] <- subset(object[['ATAC']], features = names(which(tmp >= 10)))
  
  
  ##-----Normalization, dimensional reduction, and clustering on RNA-seq and ATAC-seq separately
  
  # RNA-seq
  DefaultAssay(object) <- 'RNA'
  object <- SCTransform(object)
  object <- RunPCA(object)
  object <- RunUMAP(object, reduction = 'pca', dims = 1:n.pc, assay = 'SCT',
                    reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  object <- FindNeighbors(object, reduction = 'pca', dims = 1:n.pc, assay = 'SCT')
  object <- FindClusters(object, graph.name = 'SCT_snn', algorithm = 3, resolution = 0.2)
  
  # ATAC-seq
  DefaultAssay(object) <- 'ATAC'
  object <- RunTFIDF(object, method = 3)
  object <- FindTopFeatures(object, min.cutoff = 'q75')
  object <- RunSVD(object)
  object <- RunUMAP(object, reduction = 'lsi', dims = 2:n.lsi, assay = 'ATAC',
                    reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  object <- FindNeighbors(object, reduction = 'lsi', dims = 2:n.lsi, assay = 'ATAC')
  object <- FindClusters(object, graph.name = 'ATAC_snn', algorithm = 3, resolution = 0.2)
  
  
  # Weighted nearest neighbor (WNN) analysis using both modalities
  object <- FindMultiModalNeighbors(object,
                                    reduction.list = list("pca", "lsi"),
                                    dims.list = list(1:n.pc, 2:n.lsi),
                                    modality.weight.name = 'RNA.weight')
  object <- RunUMAP(object, nn.name = "weighted.nn", 
                    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  object <- FindClusters(object, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.2)
  object <- FindMultiModalNeighbors(object,
                                    reduction.list = list("harmony.pca", "harmony.lsi"),
                                    dims.list = list(1:n.pc, 2:n.lsi),
                                    modality.weight.name = 'RNA.weight.harmony',
                                    knn.graph.name = 'wknn.harmony',
                                    snn.graph.name = 'wsnn.harmony',
                                    weighted.nn.name = 'weighted.nn.harmony')
  

  DefaultAssay(object) <- 'SCT'
  
  return(object)
  
}

OBJ_FILT = ClusterData(OBJ_FILT, N.PC, N.LSI)


