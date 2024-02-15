.libPaths()

# dreamlet universe
library(dreamlet)
library(crumblr)
library(zenith)

# data IO
library(SingleCellExperiment) 
library(zellkonverter)
library(tidyr)

# plotting
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(ggtree)
library(aplot)

# meta
library(muscat)
library(metafor)
library(broom)
library(tidyverse)

# gsea
library(topGO)
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

library(clusterProfiler)
library(enrichplot)




meta_analysis = function( tabList ){
    
    # set entry names of none
    if( is.null(names(tabList)) ){
        names(tabList) = as.character(seq(length(tabList)))
    }

    # define dataset
    for( key in names(tabList) ){
        tabList[[key]]$Dataset = key
    }

    # stack datasets
    df = do.call(rbind, tabList) 

    # meta-analysis for each matching gene and assay
    # compute se from logFC and t
    # Use the fact that t = logFC / se
    df %>%
    as_tibble %>%
    group_by(assay) %>%
    do(tidy(rma( yi = logFC, sei = logFC / t, data=., method = "FE"))) %>%
    select(-term, -type) %>%
    ungroup() %>%
    mutate(FDR = p.adjust(p.value, "fdr")) %>% 
    mutate('log10FDR' = -log10(FDR))
}




plotTree = function(tree, low="grey90", mid = "red", high="darkred", xmax.scale=1.5){

    fig = ggtree(tree, branch.length = "none") + 
        geom_tiplab(color = "black", size=4, hjust=0, offset=.2) +
        theme(legend.position="top left", plot.title = element_text(hjust = 0.5))

    # get default max value of x-axis
    xmax = layer_scales(fig)$x$range$range[2]

    # increase x-axis width
    fig + xlim(0, xmax*xmax.scale) 
}

plotCoef2 = function(tab, coef, fig.tree, low="grey90", mid = "red", high="darkred", ylab, tick_size){
    tab$logFC = tab$estimate
    tab$celltype = factor(tab$assay, rev(ggtree::get_taxa_name(fig.tree)))
    tab$se = tab$std.error
    fig.es = ggplot(tab, aes(celltype, logFC)) + 
        geom_hline(yintercept=0, linetype="dashed", color="grey", linewidth=1) +
        geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0) +
        # geom_point(color="dodgerblue") +
        geom_point2(aes(color=pmin(4,-log10(FDR)), size=pmin(4,-log10(FDR)))) + 
        scale_color_gradient2(name = bquote(-log[10]~FDR), limits=c(0,4), low=low, mid=mid, high=high, midpoint=-log10(0.01)) +
        scale_size_area(name = bquote(-log[10]~FDR), limits=c(0,4)) +
        geom_text2(aes(label = '+', subset=FDR < 0.05), color = "white", size=6, vjust=.3, hjust=.5) +
        theme_classic() +
        coord_flip() +
        xlab('') + 
        ylab(ylab) +
        theme(axis.text.y=element_blank(), axis.text=element_text(size = 12),
              axis.ticks.y=element_blank(), text = element_text(size = tick_size)) +
        scale_y_continuous(breaks = scales::breaks_pretty(3))
    return(fig.es)    
}


#--- load results from variance partition

load('./varpart_ALLsnMultiome_celltype.RDATA')
options(repr.plot.width=10, repr.plot.height=10)
plotVarPart(sortCols(vp.lst), label.angle=60, ncol=4) 


#--- crumblr on pseudobulk

cobj = crumblr(cellCounts(pb_celltype))
form <- ~ (Age) + (1|Brain.region) + (1|Sex) + (1|Brain.bank) + (1|Individual)

vp.c = fitExtractVarPartModel(cobj, form, colData(pb_celltype))
options(repr.plot.width=6, repr.plot.height=5)
plotVarPart(sortCols(vp.c), label.angle=60, ncol=4) 

fig.vp = plotPercentBars(sortCols(vp.c) )
fig.vp

#--- analysis with dream()
form <- ~ (Age) + (1|Brain.region) + (1|Sex) + (1|Brain.bank) + (1|Individual)
fit = dream( cobj, form, colData(pb_celltype))
fit = eBayes(fit)





## Meta-analysis
prefix = 'Age_analysis_'

topTable(fit,coef='Age',number=Inf)%>% rownames_to_column('assay')
res.age <- topTable(fit,coef='Age',number=Inf)%>% rownames_to_column('assay')
#head(res.age)

res.age = meta_analysis(list(topTable(fit,coef='Age',number=Inf)%>% rownames_to_column('assay')))
head(res.age)
unique(fit$coefficients)


hc = buildClusterTreeFromPB(pb_celltype)
fig.tree = plotTree(ape::as.phylo(hc), xmax.scale=2.2) + theme(legend.position="bottom")


### effect size
fig.es1 = plotCoef2(res.age, coef='Age', fig.tree, ylab='Age')


### combine plots
#options(repr.plot.width=13, repr.plot.height=5)
options(repr.plot.width=8, repr.plot.height=5)
fig.es1 %>% insert_left(fig.tree, width=1.3) 
ggsave(paste0(prefix,'_controls_age_sex.pdf'), width = 13, height = 5)





## GSEA analysis

# define functions
SelectGenes <- function(VP.LST,cova, topN_genes, analysis){
    DF <- data.frame()
    for(i in cova){
        print(i)
        VP.LST %>%
        as_tibble %>%
        arrange(desc(!!sym(i))) %>%
        slice_head(n=topN_genes)-> top100

        #print(top100)
        Gene <- top100$gene
       # print(Gene)
        df <- as.data.frame(Gene)
        df$covariate <- i
        df$analysis <- analysis
        df$topN_genes <- topN_genes
        DF <- rbind(DF, df)
    }

    DF$topN_genes <- DF$gene
    return(DF)
}


# get expression values for correlation
ListExpression <- function(dataframe, Gene_name){
    test <- subset(dataframe, subset=Gene==Gene_name)
    test.df <- test[head(seq_len(ncol(test)), -1)]
    test.list <- as.list(as.data.frame(t(test.df)))
    #head(test.list)
    return(test.list)
}




PrepforDotplot <- function(Df, Database, numGenes, numGenesCutoff){
    
    Df_bio_proc <- Df[which(Df$database==Database),]
    CNT <- c()
    for(i in 1:nrow(Df_bio_proc)){
        gens <- Df_bio_proc[i,]$Genes
        CNT <- append(CNT,count.fields(textConnection(gens), sep = ";"))
    }

    Df_bio_proc$Gene_count <- CNT
    D <- Df_bio_proc
    D0 <- data.frame()
    for(i in unique(D$Top_DE)){
        d<-D[which(D$Top_DE==i),]
        d %>% top_n(n = -numGenes, wt = Adjusted.P.value) -> d0
        print(nrow(d0))
        if(nrow(d0) > numGenes){
            d0 <- head(d0,numGenesCutoff)
        }
        D0 <- rbind(D0, d0)
    }
    return(D0)
}



runGOEnrichment <- function(df_corr){
    Df <- data.frame()

    for(i in unique(df_corr$covariate)){
        temp<- df_corr[which(df_corr$covariate==i),]
       #temp$X
        Cluster_DE = i #paste('HS_module',i,sep='')
        print(Cluster_DE)
        GENES <- temp$Gene



        #EnrichR and GO
        dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021",
                    "WikiPathway_2021_Human", "Elsevier_Pathway_Collection", "Reactome_2022")
        if (websiteLive) {
                enriched <- enrichr(GENES, dbs)
            }

        enriched[[1]]$Top_DE <- Cluster_DE
        enriched[[1]]$database <- 'GO_Molecular_Function_2021'
        TEXT <- c()
        GOT <- c()
        for(i in enriched[[1]]$Term){
            #print(i)
            A <- strsplit(i, split = "\\(")
            TEXT<- append(TEXT,A[[1]][1])
            B <- strsplit(A[[1]][2], split = "\\)")
            GOT <- append(GOT,B[[1]][1])
        }

        enriched[[1]]$GO_text <- TEXT
        enriched[[1]]$GO_number <- GOT
        #write.csv(enriched[[1]],'./GO_Molecular_Function_2021_3DG_ALL_GluDE.csv')
        #write.csv(enriched[[1]]$GO_number,'./GO_Molecular_Function_2021_3DG_ALL_GluDE_GO.csv')


        #----------
        enriched[[2]]$Top_DE <- Cluster_DE
        enriched[[2]]$database <- 'GO_Biological_Process_2021'
        TEXT <- c()
        GOT <- c()
        for(i in enriched[[2]]$Term){
            #print(i)
            A <- strsplit(i, split = "\\(")
            TEXT<- append(TEXT,A[[1]][1])
            B <- strsplit(A[[1]][2], split = "\\)")
            GOT <- append(GOT,B[[1]][1])
        }

        enriched[[2]]$GO_text <- TEXT
        enriched[[2]]$GO_number <- GOT
        #write.csv(enriched[[2]],'./GO_Biological_Process_2021_3DG_ALL_GluDE.csv')
        #write.csv(enriched[[2]]$GO_number,'./GO_Biological_Process_2021_3DG_ALL_GluDE_GO.csv')

        #----------
        enriched[[3]]$Top_DE <- Cluster_DE
        enriched[[3]]$database <- 'WikiPathway_2021_Human'
        TEXT <- c()
        GOT <- c()
        for(i in enriched[[3]]$Term){
            #print(i)
            A <- strsplit(i, split = "\\(")
            TEXT<- append(TEXT,A[[1]][1])
            B <- strsplit(A[[1]][2], split = "\\)")
            GOT <- append(GOT,B[[1]][1])
        }
        enriched[[3]]$GO_text <- TEXT
        enriched[[3]]$GO_number <- GOT
        #write.csv(enriched[[3]],'./WikiPathway_2021_Human_3DG_ALL_GluDE.csv')
        #write.csv(enriched[[3]]$GO_number,'./WikiPathway_2021_Human_3DG_ALL_GluDE_GO.csv')

        #----------
        enriched[[4]]$Top_DE <- Cluster_DE
        enriched[[4]]$database <- 'Elsevier_Pathway_Collection'
        TEXT <- c()
        GOT <- c()
        for(i in enriched[[4]]$Term){
            #print(i)
            A <- strsplit(i, split = "\\(")
            TEXT<- append(TEXT,A[[1]][1])
            B <- strsplit(A[[1]][2], split = "\\)")
            GOT <- append(GOT,B[[1]][1])
        }

        enriched[[4]]$GO_text <- TEXT
        enriched[[4]]$GO_number <- GOT
        #write.csv(enriched[[4]],'./Elsevier_Pathway_Collection_3DG_ALL_GluDE.csv')
        #write.csv(enriched[[4]]$GO_number,'./Elsevier_Pathway_Collection_3DG_ALL_GluDE_GO.csv')

        #----------
        enriched[[5]]$Top_DE <- Cluster_DE
        enriched[[5]]$database <- 'Reactome_2022'
        TEXT <- c()
        GOT <- c()
        for(i in enriched[[5]]$Term){
            #print(i)
            A <- strsplit(i, split = "\\(")
            TEXT<- append(TEXT,A[[1]][1])
            B <- strsplit(A[[1]][2], split = "\\)")
            GOT <- append(GOT,B[[1]][1])
        }
        enriched[[5]]$GO_text <- TEXT
        enriched[[5]]$GO_number <- GOT
        #write.csv(enriched[[5]],'./Reactome_2022_3DG_ALL_GluDE.csv')
        #write.csv(enriched[[5]]$GO_number,'./Reactome_2022_3DG_ALL_GluDE_GO.csv')


        Df <- rbind(Df, enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]],enriched[[5]])



        p0<- if (websiteLive) plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, 
                                              y = "Count", orderBy = "Adjusted.P.value")
        p1<- if (websiteLive) plotEnrich(enriched[[2]], showTerms = 20, numChar = 50, 
                                              y = "Count", orderBy = "Adjusted.P.value")
        p3<-if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 50, 
                                              y = "Count", orderBy = "Adjusted.P.value")
        p4<-if (websiteLive) plotEnrich(enriched[[4]], showTerms = 20, numChar = 50, 
                                              y = "Count", orderBy = "Adjusted.P.value")
        p5<-if (websiteLive) plotEnrich(enriched[[5]], showTerms = 20, numChar = 50, 
                                              y = "Count", orderBy = "Adjusted.P.value")
        options(repr.plot.width=20, repr.plot.height=8)

        print(p0+p1)
        print(p3+p4)
        print(p1+p5)



    }
    
    return(Df)
}










