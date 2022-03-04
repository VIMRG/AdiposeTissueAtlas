###############################################################################
######-----------------CAPILLARY PSEUDOTIME-------------------------#####
###############################################################################

# 1. SET UP___________________________________________________________________
##########################
# Set seed
##########################
set.seed(7612) # Reproducibility

##########################
# Set working directory
##########################
setwd('/data/p_koethe_lab/Atlas_AT/Analysis')

##########################
# Source function for conversion to h5ad
##########################
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/AnnData.R')

##########################
# Load Libraries
##########################
library(Seurat)
library(slingshot)
library(ggplot2)
library(dplyr)
library(gplots)
library(stats)
library(scater)
library(tradeSeq)
library(BiocParallel)
library(harmony)
library(future)
library(scater)
library(tidyverse)
library(seriation)
library(ComplexHeatmap)
library(phateR)
library(destiny)
library(TSCAN)
library(ElPiGraph.R)
library(monocle3)
library(SCORPIUS)
library(igraph)

plan("multiprocess", workers = 5)
options(future.globals.maxSize = 3000 * 1024^2)
options(future.seed = TRUE)

##########################
# Load Vascular
##########################
Endo <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Endo_RNA_1.3.rds')

##########################
# Function: Filter Genes
##########################
Genes.to.use.sce <- function(object, min.value = 1, min.cells = 200) {
    num.cells <- Matrix::rowSums(assays(object)$counts > min.value) # number of cells expressing > 1
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    return(genes.use)
}

# 2. CAPILLARY ENDO-MT___________________________________________________________________
##########################
# First Iteration
##########################
Idents(Endo) <- "seurat_clusters"
Cap <- subset(Endo, idents = c(8,16,26,2,0,11,20,17,5,1,15))
cap1 <- Cap # save UMAP

Cap <- Cap %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
Cap<- RunHarmony(Cap, "Project", assay.use = "RNA", max.iter.harmony = 30)

Cap <- Cap %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:25) %>%
        FindClusters(resolution = 1.5)

p1 <- DimPlot(Cap, label = T) + NoLegend()
p2 <- DimPlot(Cap, label = T, group.by = "Manual_Annotation") + NoLegend()
p1 + p2

##########################
# Second Iteration
##########################
Cap <- subset(Cap, idents = c(7,9,10,11,12,14), invert = T)

Cap <- Cap %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
Cap<- RunHarmony(Cap, "Project", assay.use = "RNA", max.iter.harmony = 30)

Cap <- Cap %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 1.5)

p1 <- DimPlot(Cap, label = T) + NoLegend()
p2 <- DimPlot(Cap, label = T, group.by = "Manual_Annotation") + NoLegend()
p1 + p2

##########################
# Third Iteration
##########################
Cap <- subset(Cap, idents = c(9,2), invert = T)

Cap <- Cap %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
Cap<- RunHarmony(Cap, "Project", assay.use = "RNA", max.iter.harmony = 30)

Cap <- Cap %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 1.5)

p1 <- DimPlot(Cap, label = T) + NoLegend()
p2 <- DimPlot(Cap, label = T, group.by = "Manual_Annotation") + NoLegend()
p1 + p2

##########################
# Final Iteration
##########################
Cap <- subset(Cap, idents = c(13), invert = T)

Cap <- Cap %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
Cap<- RunHarmony(Cap, "Project", assay.use = "RNA", max.iter.harmony = 30)

Cap <- Cap %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 0.5)

p1 <- DimPlot(Cap, label = T) + NoLegend()
p2 <- DimPlot(Cap, label = T, group.by = "Manual_Annotation") + NoLegend()
p1 + p2

Idents(Cap, cells = WhichCells(Cap, idents = c(3))) <- "Pericyte 1"
Idents(Cap, cells = WhichCells(Cap, idents = c(1))) <- "Pericyte 2"
Idents(Cap, cells = WhichCells(Cap, idents = c(5))) <- "EndoMT 1"
Idents(Cap, cells = WhichCells(Cap, idents = c(2))) <- "EndoMT 2"
Idents(Cap, cells = WhichCells(Cap, idents = c(4))) <- "Capillary 1"
Idents(Cap, cells = WhichCells(Cap, idents = c(0))) <- "Capillary 2"

Cap[['Manual_Annotation']] <- Idents(Cap)

# Going to recluster and remove likely contaminating doublet (19) as it expressed high level of CA4, CLDN5 and Pericyte Markers (ACTA2, COX4I2)
Cap <- FindClusters(Cap, resolution = 2.0)
Cap <- subset(Cap, idents = 19, invert = T)

Cap <- Cap %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 0.5)
        
saveRDS(Cap, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Capillary_EndoMT.rds')

# 3. Pseudotime___________________________________________________________________
##########################
# Perform Slingshot
##########################
# down sample
Cap.downsample = subset(Cap, cells = sample(Cells(Cap), 12000))

cap.sce <- as.SingleCellExperiment(Cap.downsample, assay = "RNA")

cap.slingshot <- slingshot(cap.sce, clusterLabels = "Manual_Annotation", reducedDim = "UMAP", start.clus = "Capillary 2", approx_points = 150)
pseudo.paths <- slingPseudotime(cap.slingshot)
head(pseudo.paths)
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

summary(cap.slingshot$slingPseudotime_1)
print(SlingshotDataSet(cap.slingshot))

##########################
# Plots 
##########################

color <- as.character(cap.slingshot$Manual_Annotation)
temp <- color
temp[temp == 'Pericyte 1'] <- "#03AF21"
temp[temp == 'Pericyte 2'] <- "#FE627D"
temp[temp == 'EndoMT 1'] <- "#EE7342"
temp[temp == 'EndoMT 2'] <- "#77D46E"
temp[temp == 'Capillary 1'] <- "#EBCE71"
temp[temp == 'Capillary 2'] <- "#DD63EB"
color <- temp

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Slingshot_curve_Capillary2.png", res = 600, width = 10, height = 5, units = 'in')
plot(reducedDims(cap.slingshot)$UMAP, col = color,
     pch = 16, asp = 1, xlab = "", ylab = "", cex=0.6, axes=FALSE, frame.plot=FALSE)   
legend(par(xpd = T), x= "topleft", pch = c(20), cex = 1.0, 
       legend = c("Pericyte 1", "Pericyte 2", "EndoMT 1", "EndoMT 2", "Capillary 1", "Capillary 2"), 
       col = c("#03AF21", "#FE627D", "#EE7342", "#77D46E", "#EBCE71", "#DD63EB"), bty = 'n')
lines(slingCurves(cap.slingshot)$curve1, lwd=3)
dev.off()

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Pseudotime_Capillary2.png", res = 600, width = 10, height = 5, units = 'in')
plotUMAP(cap.slingshot, colour_by = I(shared.pseudo)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14))
dev.off()

# 3. GAM___________________________________________________________________
##########################
# GAM Model
##########################
counts <- as.matrix(assay(cap.sce))
crv <- SlingshotDataSet(cap.slingshot)
genes.use <- Genes.to.use.sce(cap.sce, min.value = 2)

# Run GAM Model
set.seed(6)
BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 5

sce <- fitGAM(counts = counts, sds = crv, nknots = 6, genes = genes.use, verbose = T, parallel=TRUE, BPPARAM = BPPARAM)

##########################
# Association Test
##########################
assocRes <- as.data.frame(associationTest(sce), l2fc = log2(2))

assocRes$FDR <- p.adjust(assocRes$pvalue, method = 'fdr')
assocRes <- assocRes[assocRes$FDR <= 0.05,]

##########################
# Rolling Wave Plot
##########################
# Smooth Plot
Smooth <- predictSmooth(sce, gene = rownames(assocRes), tidy = F,  n = 100)
Smooth <- t(scale(t(Smooth))) # scale

# Seriate
Smooth <- Smooth[get_order(seriate(Smooth, method = "PCA_angle")),]

# Get Annotation
Pseudotime <- as.data.frame(pseudo.paths)
Pseudotime$Annotation <- cap.slingshot$Manual_Annotation
Pseudotime$Interval <- cut(Pseudotime$curve1, breaks = 100)
Pseudotime <- Pseudotime %>% mutate(lineage = factor(Interval, labels = c(1:100)))
Pseudotime <- Pseudotime %>% mutate(Lineage_Annotation = case_when(lineage %in% c(1:35) ~ "Capillary",
lineage %in% c(36:70) ~ "EndoMT", TRUE ~ "Pericyte"))

# Capillary: (6.62,6.82]
# EndoMT: (13.4,13.6]

# Create column names
col = list(CellType = c("Capillary" = "red", "EndoMT" = "blue", "Pericyte" = "brown"))
column_ha = HeatmapAnnotation(CellType = c(rep("Capillary", 35), rep("EndoMT", 35), rep("Pericyte", 30)), col = col, 
annotation_legend_param = list(CellType = list(title = "Cell Type", at = c("Capillary", "EndoMT", "Pericyte"), labels = c("Capillary", "EndoMT", "Pericyte"))))

# Create Heatmaps
png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Heatmap_Cap2_Recluster.png", width = 10, height = 10, units = 'in', res = 600)
Heatmap(Smooth, cluster_columns = F, cluster_rows = F, show_column_names = F, show_row_names = F, top_annotation = column_ha) # bitemporal
dev.off()

# Save TradeSeq
saveRDS(sce, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/GAM_Cap2_Recluster.rds')



