##############################################################################################
##---------------------------- MACROPHAGE PSEUDOTIME ANALYSIS-------------------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will evaluate Pseudotime and changes in gene expression
## along the trajectory, specifically for LAM and PVM.
##############################################################################################

# 1. SETTING UP ENVIRONMENT------------------------------------------------------------------

###############################
# set seed
###############################
set.seed(7612) # Reproducibility

###############################
# set directory
###############################
date <- "6.27"
Subset_dir <- paste0("../", date, "/SubsetAnalysis")
fig_dir <- "../Figures"

###############################
# Load Libraries
###############################
library(Seurat) # Seurat 4.0 version(note that several packages will also be attached)
library(tidyverse) # collection of packages for manipulating data (ggplot2, dplyr, tidyr, purr, stringr, tibble)
library(Matrix) # dense and sparse matrix manipulations
library(future) # Parallelize certain functions
library(slingshot) # Pseudotime Trajectory Analysis
library(scater)
library(tradeSeq) # GLM for Genes along lineage
library(BiocParallel) # Used by tradeseq for parallelization
library(harmony) # Integration algorithm
library(seriation) # Seriate pseudotime genes
library(phateR) # Trajectory analysis
library(destiny) # Trajectory analysis
library(TSCAN) # Trajectory Analysis
library(ElPiGraph.R) # Tree graph
library(monocle3) # Trajectory analysis
library(SCORPIUS) # Trajectory analysis
library(igraph)
library(ComplexHeatmap) # plotting
library(MAST) 
library(clusterProfiler) # GSEA analysis
library(org.Hs.eg.db) # Gene conversions
library(DOSE) # 


plan("multiprocess", workers = 10)
options(future.globals.maxSize = 3000 * 1024^2)
options(future.seed = TRUE)

Subset_dir <- paste0("../", date, "/SubsetAnalysis")
fig_dir <- "../Figures"
Pseudo_dir <- paste0("../Pseudotime")
dir.create(Pseudo_dir)

###############################
# Load Data
###############################
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage.rds"))

# 1. PVM PSEUDOTIME------------------------------------------------------------------

Macrophage <- FindClusters(Macrophage, resolution = 0.75)

# Subset on clusters 5 (Mo-Mac 2), 4 (Mo-Mac 1), 0 (PVM) and 6 (PVM)
PVM <- subset(Macrophage, idents = c(5,4,0,6))
Idents(PVM) <- "Annotation"
PVM <- subset(PVM, idents = c("IM", "LAM"), invert = T) # 78 cells

# Convert to SCE
PVM_sce <- as.SingleCellExperiment(PVM, assay = "RNA")

# Remove ribosomal and mitochondrial genes
MT_genes <- grep("^MT-", rownames(PVM_sce), ignore.case = F, value = T)
RP_genes <- grep("^RP", rownames(PVM_sce), ignore.case = F, value = T)
PVM_sce <- PVM_sce[!rownames(PVM_sce) %in% MT_genes,]
PVM_sce <- PVM_sce[!rownames(PVM_sce) %in% RP_genes,]

# Run Slingshot on UMAP with Mo-Mac 1 as starting cluster and PVM as ending cluster
PVM_slingshot <- slingshot(PVM_sce, clusterLabels = "Annotation", reducedDim = "UMAP", start.clus = "Mo-Mac 2", end.clus = "PVM", approx_points = 150)

pseudo.paths <- slingPseudotime(PVM_slingshot)
head(pseudo.paths)
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

summary(PVM_slingshot$slingPseudotime_1)
print(SlingshotDataSet(PVM_slingshot))

# Plot Pseudotime
plotUMAP(PVM_slingshot, colour_by = I(shared.pseudo))

# Remove Contaminating cells (LAM and IM)
DF <- data.frame(barcode = colnames(PVM_slingshot), pseudotime = PVM_slingshot$slingPseudotime_1, annotation = PVM_slingshot$Annotation)
DF <- DF %>% mutate(time_bin = cut(pseudotime, breaks=100))
prop.table(table(DF$time_bin, DF$annotation), margin = 1) * 100

IMcontaminant <- which(PVM_slingshot$slingPseudotime_1 > 15.0 | PVM_slingshot$slingPseudotime_1 < 4.3) # Exclude cells off the trajectory
PVM_sce <- PVM_sce[, -IMcontaminant,]

PVM_slingshot <- slingshot(PVM_sce, clusterLabels = "Annotation", reducedDim = "UMAP", start.clus = "Mo-Mac 2", end.clus = "PVM", approx_points = 150)

pseudo.paths <- slingPseudotime(PVM_slingshot)
head(pseudo.paths)
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Plot Pseudotime
png(paste0(fig_dir, "/", "SF7A.png"), width = 5, height = 5, units = 'in', res = 600)
plotUMAP(PVM_slingshot, colour_by = I(shared.pseudo))
dev.off()

##########################
# TradeSeq GAM Model
##########################
# Select PVM Trajectory
tmp_dir <- paste0("../", "GLM")
dir.create(tmp_dir)

# Set up Model
counts <- as.matrix(assay(PVM_sce))
crv <- SlingshotDataSet(PVM_slingshot)

Genes.to.use.sce <- function(object, min.value = 1, min.cells = 200) {
    num.cells <- Matrix::rowSums(assays(object)$counts > min.value) # number of cells expressing > 1
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    return(genes.use)
}

genes.use <- Genes.to.use.sce(PVM_sce, min.cells = 300) # 5-10% of cells expressing

# Run GAM Model
set.seed(6)
BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 10

# Find optimal nKots
icMat <- evaluateK(counts = counts, sds = crv, k = 3:7, nGenes = 100,
                   verbose = FALSE, plot = TRUE)
print(icMat[1:2, ])

# Build Model
sce <- fitGAM(counts = counts, sds = crv, nknots = 6, genes = genes.use, verbose = T, parallel=TRUE, BPPARAM = BPPARAM)
saveRDS(sce, file = paste0(tmp_dir, "/", "Macrophage_PVM_TradeSeq_4.1.rds"))

##########################
# Association Test
##########################
# Association Test
AssocRes <- as.data.frame(associationTest(sce), l2fc = log2(2))
AssocRes$FDR <- p.adjust(AssocRes$pvalue, method = "fdr")
AssocRes_Sig <- AssocRes[AssocRes$FDR <= 0.05,]

##########################
# Seriate for plotting
##########################
# Smooth
Smooth <- predictSmooth(sce, gene = rownames(AssocRes_Sig), tidy = F,  n = 100)
Smooth <- t(scale(t(Smooth))) # scale

# Seriate
Smooth <- Smooth[get_order(seriate(Smooth, method = "PCA_angle")),]

# Get Annotations
DF <- data.frame(barcode = colnames(PVM_slingshot), pseudotime = PVM_slingshot$slingPseudotime_1, annotation = PVM_slingshot$Annotation)
DF <- DF %>% mutate(time_bin = cut(pseudotime, breaks=100))
prop.table(table(DF$time_bin, DF$annotation), margin = 1) * 100
Ann <- c(rep("Mo-Mac 2", 29), rep("Mo-Mac 1", 28), rep("PVM", 43))
Arranged <- Smooth
genes <- c("LYZ", "HLA-CRB1", "CPVL", "CST3", "CEBPD", "CCL2", "LYVE1", "C1QB", "CD36", "HMOX1", "CD14", "VCAN", "ISG15", "S100A8", "S100A12")
index <- which(rownames(Arranged) %in% genes)
label <- rownames(Arranged)[index]
row_ha <- rowAnnotation(foo = anno_mark(at = index, labels = label, which = "row", side = "right"))
ha = HeatmapAnnotation(Cell = Ann, col = list(Cell = c("Mo-Mac 2" = 'magenta', "Mo-Mac 1" = "salmon",  "PVM" = "steelblue")))

png(paste0(fig_dir, "/", "SF7B.png"), width = 10, height = 10, units = 'in', res = 600)
Heatmap(Arranged, cluster_columns = F, cluster_rows = F, top_annotation = ha, right_annotation = row_ha, show_column_names = F, show_row_names = F)
dev.off() 

##########################
# Transcription Factor Regulation
##########################
TF <- read.delim("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/MetaData/TF.txt", header = F)
colnames(TF) <- "Gene"

AssocTF <- AssocRes_Sig[rownames(AssocRes_Sig) %in% TF$Gene,]
SmoothDiff <- predictSmooth(sce, gene = rownames(AssocTF), tidy = F, n = 100)
SmoothDiff <- t(scale(t(SmoothDiff))) # Scale
Smooth <- SmoothDiff[get_order(seriate(SmoothDiff, method = "PCA_angle")),]

DF <- data.frame(barcode = colnames(PVM_slingshot), pseudotime = PVM_slingshot$slingPseudotime_1, annotation = PVM_slingshot$Annotation)
DF <- DF %>% mutate(time_bin = cut(pseudotime, breaks=100))
prop.table(table(DF$time_bin, DF$annotation), margin = 1) * 100
#Ann <- c(rep("ISG+ Mo", 20), rep("Mo-Mac 1", 44), rep("PVM", 36))
Arranged <- Smooth

ha = HeatmapAnnotation(Cell = Ann, col = list(Cell = c("Mo-Mac 2" = 'magenta', "Mo-Mac 1" = "salmon", "PVM" = "steelblue")))
genes = c("EGR1", "EGR2", "ATF3", "JUN", "FOS", "MAF", "KLF4", "MAFB", "CEBPD", "JUND", "SON", "STAT1", "CEBPB")
index <- which(rownames(Arranged) %in% genes)
label <- rownames(Arranged)[index]
row_ha <- rowAnnotation(foo = anno_mark(at = index, labels = label, which = "row", side = "right"))

png(paste0(fig_dir, "/", "Figure6C.png"), width = 10, height = 10, units = 'in', res = 600)
Heatmap(Arranged, cluster_columns = F, cluster_rows = F, top_annotation = ha, right_annotation = row_ha, show_column_names = F, show_row_names = F)
dev.off()    

# 3. LAM PSEUDOTIME------------------------------------------------------------------
# Convert to SCE
LAM <- subset(Macrophage, idents = c(7,2))
Idents(LAM) <- "Annotation"
LAM <- subset(LAM, idents = c("PVM", "IM"), invert = T)

# Convert to SCE
LAM_sce <- as.SingleCellExperiment(LAM, assay = "RNA")

# Remove ribosomal and mitochondrial genes
MT_genes <- grep("^MT-", rownames(LAM_sce), ignore.case = F, value = T)
RP_genes <- grep("^RP", rownames(LAM_sce), ignore.case = F, value = T)
LAM_sce <- LAM_sce[!rownames(LAM_sce) %in% MT_genes,]
LAM_sce <- LAM_sce[!rownames(LAM_sce) %in% RP_genes,]

# Run Slingshot on UMAP with Mo-Mac 1 as starting cluster and PVM as ending cluster
LAM_slingshot <- slingshot(LAM_sce, clusterLabels = "Annotation", reducedDim = "UMAP", start.clus = "Mo-Mac 2", end.clus = "LAM", approx_points = 150)

pseudo.paths <- slingPseudotime(LAM_slingshot)
head(pseudo.paths)
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Plot Pseudotime
plotUMAP(LAM_slingshot, colour_by = I(shared.pseudo))

t <- which(LAM_slingshot$slingPseudotime_1 > 10)
LAM_sce <- LAM_sce[, -t,]

LAM_slingshot <- slingshot(LAM_sce, clusterLabels = "Annotation", reducedDim = "UMAP", start.clus = "Mo-Mac 2", end.clus = "LAM", approx_points = 150)

pseudo.paths <- slingPseudotime(LAM_slingshot)
head(pseudo.paths)
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Pseudotime
png(paste0(fig_dir, "/", "SF7C.png"), width = 5, height = 5, units = 'in', res = 600)
plotUMAP(LAM_slingshot, colour_by = I(shared.pseudo))
dev.off()

##########################
# TradeSeq GAM Model
##########################
# Select LAM Trajectory
tmp_dir <- paste0("../", "GLM")
dir.create(tmp_dir)

# Set up Model
counts <- as.matrix(assay(LAM_sce))
crv <- SlingshotDataSet(LAM_slingshot)

Genes.to.use.sce <- function(object, min.value = 1, min.cells = 200) {
    num.cells <- Matrix::rowSums(assays(object)$counts > min.value) # number of cells expressing > 1
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    return(genes.use)
}

genes.use <- Genes.to.use.sce(LAM_sce, min.cells = 400) # 5-10% of cells expressing

# Run GAM Model
set.seed(6)
BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 10

# Find optimal nKots
icMat <- evaluateK(counts = counts, sds = crv, k = 3:7, nGenes = 100,
                   verbose = FALSE, plot = TRUE)
print(icMat[1:2, ])

# Build Model
sce <- fitGAM(counts = counts, sds = crv, nknots = 6, genes = genes.use, verbose = T, parallel=TRUE, BPPARAM = BPPARAM)
saveRDS(sce, file = paste0(tmp_dir, "/", "Macrophage_LAM_TradeSeq_4.1.rds"))

##########################
# Association Test
##########################
# Association Test
AssocRes <- as.data.frame(associationTest(sce), l2fc = log2(2))
AssocRes$FDR <- p.adjust(AssocRes$pvalue, method = "fdr")
AssocRes_Sig <- AssocRes[AssocRes$FDR <= 0.05,]


##########################
# Seriate for plotting
##########################
# Smooth
Smooth <- predictSmooth(sce, gene = rownames(AssocRes_Sig), tidy = F,  n = 100)
Smooth <- t(scale(t(Smooth))) # scale

# Seriate
Smooth <- Smooth[get_order(seriate(Smooth, method = "PCA_angle")),]

# Get Annotations
DF <- data.frame(barcode = colnames(LAM_slingshot), pseudotime = LAM_slingshot$slingPseudotime_1, annotation = LAM_slingshot$Annotation)
DF <- DF %>% mutate(time_bin = cut(pseudotime, breaks=100))
prop.table(table(DF$time_bin, DF$annotation), margin = 1) * 100
Ann <- c(rep("Mo-Mac 2", 22), rep("LAM", 78))
Arranged <- Smooth
genes <- c("LYZ", "VCAN", "S100A8", "S100A12", "CXCL2", "CXCL3", "AREG", "CD36", "LIPA", "TREM2", "SPP1")
index <- which(rownames(Arranged) %in% genes)
label <- rownames(Arranged)[index]
row_ha <- rowAnnotation(foo = anno_mark(at = index, labels = label, which = "row", side = "right"))
ha = HeatmapAnnotation(Cell = Ann, col = list(Cell = c("Mo-Mac 2" = "peru",  "LAM" = "gold")))

png(paste0(fig_dir, "/", "SF7D.png"), width = 10, height = 10, units = 'in', res = 600)
Heatmap(Arranged, cluster_columns = F, cluster_rows = F, top_annotation = ha, right_annotation = row_ha, show_column_names = F, show_row_names = F)
dev.off()    

Heatmap(Smooth, cluster_columns = F, cluster_rows = F, show_column_names = F, show_row_names = F) # bitemporal


##########################
# Transcription Factor Regulation
##########################
TF <- read.delim("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/MetaData/TF.txt", header = F)
colnames(TF) <- "Gene"

AssocTF <- AssocRes_Sig[rownames(AssocRes_Sig) %in% TF$Gene,]
SmoothDiff <- predictSmooth(sce, gene = rownames(AssocTF), tidy = F, n = 100)
SmoothDiff <- t(scale(t(SmoothDiff))) # Scale
Smooth <- SmoothDiff[get_order(seriate(SmoothDiff, method = "PCA_angle")),]

DF <- data.frame(barcode = colnames(LAM_slingshot), pseudotime = LAM_slingshot$slingPseudotime_1, annotation = LAM_slingshot$Annotation)
DF <- DF %>% mutate(time_bin = cut(pseudotime, breaks=100))
prop.table(table(DF$time_bin, DF$annotation), margin = 1) * 100
Ann <- c(rep("Mo-Mac 2", 22), rep("LAM", 78))
Arranged <- Smooth

ha = HeatmapAnnotation(Cell = Ann, col = list(Cell = c("ISG+ Mo" = "coral", "Mo-Mac 2" = "peru",  "LAM" = "gold")))
genes = c("STAT1", "CEBPB", "HIF1A", "PPARG", "NR1H3", "NFKB1")
index <- which(rownames(Arranged) %in% genes)
label <- rownames(Arranged)[index]
row_ha <- rowAnnotation(foo = anno_mark(at = index, labels = label, which = "row", side = "right"))

png(paste0(fig_dir, "/", "Figure6D.png"), width = 10, height = 10, units = 'in', res = 600)
Heatmap(Arranged, cluster_columns = F, cluster_rows = F, top_annotation = ha, right_annotation = row_ha, show_column_names = F, show_row_names = F)
dev.off()    


