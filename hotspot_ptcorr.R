library(Seurat)
library(Signac)

#-------- hotspot trendline examination
library(Seurat)
library(Signac)
library(ggplot2)
opc_olig3 <- readRDS('OPCOligObject.rds')

annot.use <- read.csv('3DG_OPCOlig_annotation_Sep2023.csv',row.names = 1)
rownames(annot.use) <- annot.use$Colnames
opc_olig3 <- AddMetaData(opc_olig3,metadata = annot.use[Cells(opc_olig3),])


# gene modules from hotspot output
hs.modules <- read.csv('hs_results_modules.csv',row.names = 1)
hs.list <- lapply(setdiff(sort(unique(hs.modules$Module)),-1),function(x){
  rownames(hs.modules)[which(hs.modules == x)]
})

# single cell module score
opc_olig3 <- AddModuleScore(opc_olig3,features = hs.list)
colnames(opc_olig3@meta.data)[grep('^Cluster',colnames(opc_olig3@meta.data))] <- paste0('HS_',1:length(hs.list))

# merge gene modules with visually similar patterns
# hs1+3, hs2+6, hs4,5,7,8,11,12, hs9+10
hs.new <- list(c(1,3,8),c(2,6),c(4,5,7,11,12),c(9,10))
hs.list.2 <- lapply(hs.new,function(x){
  rownames(hs.modules)[hs.modules$Module %in% x]
})

# single cell module score with updated module assignment
opc_olig3 <- AddModuleScore(opc_olig3,features = hs.list.2)
colnames(opc_olig3@meta.data)[grep('^Cluster',colnames(opc_olig3@meta.data))] <- paste0('HS2_',1:length(hs.list.2))

opc_olig3$fgroup <- opc_olig3$anno_clus_dreamorigBB_bioarx_v2
opc_olig3$fgroup <- factor(as.character(opc_olig3$fgroup),levels = c('TypeI_OPC','TypeII_OPC','COP','TypeI_Oligo_OPALIN','TypeII_Oligo_RBFOX1'))
Idents(opc_olig3) <- 'fgroup'
plot.list <- lapply(1:length(hs.list.2),function(x){
  FeatureScatter(opc_olig3,feature1 = 'pseudotime_sc_OPC1',feature2 = paste0('HS2_',x),pt.size = 0.6)
})
cowplot::plot_grid(plotlist = plot.list)


# predicted trend lines of hotspot modules
pred.use <- c()
for (i in 1:length(hs.list)){
  data.use <- FetchData(opc_olig3,vars = c('pseudotime_sc_OPC1',paste0('HS_',i)))
  colnames(data.use) <- c('pt','feature')
  model.use <- loess(feature ~ pt,data = data.use)
  pred.use <- cbind(pred.use,model.use$fitted)
}
pred.use <- as.data.frame(pred.use)
colnames(pred.use) <- paste0('HS_',1:ncol(pred.use))
pred.use$pt <- data.use$pt
rownames(pred.use) <- Cells(opc_olig3)
# plot trend lines of the same category within the same plot
hs.new
plot.list <- lapply(1:length(hs.new),function(x){
  plot.data <- reshape2::melt(pred.use[,paste0('HS_',hs.new[[x]])])
  colnames(plot.data) <- c('Module','Value')
  plot.data$pt <- rep(pred.use$pt,times = length(hs.new[[x]]))
  g <- ggplot(plot.data,aes(x = pt,y = Value,color=Module))+
    geom_line(size=1.2)+theme_bw()+ggtitle(paste0('Trend_',x))
  return (g)
})
cowplot::plot_grid(plotlist = plot.list)


#-------- pseudotime vs eRegulon activity correlation

# pt correlation with eRegulons
opc_olig3 <- readRDS('OPCOligObject.rds')

ereg.use <- readRDS('scenic_sc.rds')
cells.use <- rownames(ereg.use)

opc_olig3 <- subset(opc_olig3,cells = cells.use)
opc_olig3[['scenic']] <- CreateAssayObject(data = t(ereg.use))
DefaultAssay(opc_olig3) <- 'scenic'
opc_olig3 <- ScaleData(opc_olig3,features = rownames(opc_olig3))
set.seed(42)
plot.cells <- lapply(levels(opc_olig3$pt_bin),function(x){
  sample(Cells(opc_olig3)[opc_olig3$pt_bin == x],size = min(1000,sum(opc_olig3$pt_bin == x)))
})
plot.cells <- unlist(plot.cells)
plot.cells <- plot.cells[order(opc_olig3@meta.data[plot.cells,'pseudotime_sc_OPC1'])]

# correlation: pseudotime vs eRegulon AUC scores
cor.ereg <- sapply(1:ncol(ereg.use),function(x) cor(opc_olig3$pseudotime_sc_OPC1,ereg.use[cells.use,x]))
names(cor.ereg) <- colnames(ereg.use)
cor.ereg <- sort(cor.ereg)
# top 10 correlated eRegulons
plot.features <- c(head(names(cor.ereg),10),tail(names(cor.ereg),10))
plot.features <- gsub('_','-',plot.features)
DoHeatmap(opc_olig3,features = plot.features,group.by = 'is_cell',cells = plot.cells,label = F)

# TF expression of the top pseudotime-correlated eRegulons
DefaultAssay(opc_olig3) <- 'SCT'
plot.features <- unname(sapply(plot.features,function(x) unlist(strsplit(x,'-'))[1]))
opc_olig3 <- GetResidual(opc_olig3,plot.features)
DoHeatmap(opc_olig3,features = plot.features,group.by = 'is_cell',cells = plot.cells,label = F)
