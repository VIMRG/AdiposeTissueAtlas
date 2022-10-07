##############################################################################################
##---------------------------- MACROPHAGE PSEUDOTIME ANALYSIS-------------------------------##
##------------------------- DATE: 6/27/21 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will evaluate Pseudotime and changes in gene expression
## along the trajectory, specifically for Mo-Mac and PVM.
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
Subset_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Figures/"

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

###############################
# Load Data
###############################
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage.rds"))

Pseudo_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Pseudotime")
dir.create(tmp_dir)

# 2. Macrophage Pseudotime Mo-Mac to PVM or to IM------------------------------------------------------------------
##########################
# IM/PVM Recluster
##########################
Idents(Macrophage) <- "Annotation"
Mac <- subset(Macrophage, idents = c("Mo-Mac", "PVM", "IM"))

DefaultAssay(Mac) <- "RNA"
Mac <- Mac %>%
        FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction (Note:Stochastic process)
Mac <- Mac %>% RunHarmony("Lane", assay.use = "RNA", max.iter.harmony = 20)

# Downstream Analysis
Mac <- Mac %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 0.5)

p1 <- DimPlot(Mac, label = T) + NoLegend()
p2 <- DimPlot(Mac, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

Idents(Mac) <- "seurat_clusters"
Idents(Mac, cells = WhichCells(Mac, idents = 4)) <- "Mo-Mac 1"
Idents(Mac, cells = WhichCells(Mac, idents = 3)) <- "Mo-Mac 2"
Idents(Mac, cells = WhichCells(Mac, idents = c(0,5))) <- "PVM"
Idents(Mac, cells = WhichCells(Mac, idents = c(1,2,6))) <- "IM"
Idents(Mac, cells = WhichCells(Mac, idents = c(7))) <- "Other"

Mac[['Annotation']] <- Idents(Mac)

# 7 is a small contaminanting cluster (LAM)
Mac <- subset(Mac, idents = c("Other"), invert = T)

##########################
# IM/PVM Round 2
##########################
DefaultAssay(Mac) <- "RNA"
Mac <- Mac %>%
        FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction (Note:Stochastic process)
Mac <- Mac %>% RunHarmony("Lane", assay.use = "RNA", max.iter.harmony = 20)

# Downstream Analysis
Mac <- Mac %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 0.7)

p1 <- DimPlot(Mac, label = T) + NoLegend()
p2 <- DimPlot(Mac, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

Idents(Mac) <- "seurat_clusters"
Idents(Mac, cells = WhichCells(Mac, idents = c(3))) <- "Mo-Mac 1"
Idents(Mac, cells = WhichCells(Mac, idents = c(2))) <- "Mo-Mac 2"
Idents(Mac, cells = WhichCells(Mac, idents = c(1,6))) <- "PVM"
Idents(Mac, cells = WhichCells(Mac, idents = c(4,5,0,8))) <- "IM"
Idents(Mac, cells = WhichCells(Mac, idents = c(7))) <- "Debris"

Mac[['Annotation']] <- Idents(Mac)

Mac <- subset(Mac, idents = "Debris", invert = T)

##########################
# IM/PVM Round 3
##########################
DefaultAssay(Mac) <- "RNA"
Mac <- Mac %>%
        FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction (Note:Stochastic process)
Mac <- Mac %>% RunHarmony("Lane", assay.use = "RNA", max.iter.harmony = 20)

# Downstream Analysis
Mac <- Mac %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 0.5)

p1 <- DimPlot(Mac, label = T) + NoLegend()
p2 <- DimPlot(Mac, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

Idents(Mac) <- "seurat_clusters"
Idents(Mac, cells = WhichCells(Mac, idents = c(4))) <- "Mo-Mac 1"
Idents(Mac, cells = WhichCells(Mac, idents = c(3))) <- "Mo-Mac 2"
Idents(Mac, cells = WhichCells(Mac, idents = c(0,6))) <- "PVM"
Idents(Mac, cells = WhichCells(Mac, idents = c(5,2,1,7))) <- "IM"
Mac[['Annotation']] <- Idents(Mac)

saveRDS(Mac, file = paste(Pseudo_dir, "Mac_seurat.rds", sep = "/"))

# 3. PVM Slingshot------------------------------------------------------------------
##########################
# Slingshot PVM 
##########################
Idents(Mac) <- "seurat_clusters"
PVM <- subset(Mac, idents = c(4,3,0,6))
Idents(PVM) <- "Annotation"

PVM <- PVM %>%
        FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction (Note:Stochastic process)
PVM <- PVM %>% RunHarmony("Lane", assay.use = "RNA", max.iter.harmony = 20)

# Downstream Analysis
PVM <- PVM %>% RunUMAP(reduction = 'harmony', dims = 1:15) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:15) %>%
        FindClusters(resolution = 0.35)

p1 <- DimPlot(PVM, label = T) + NoLegend()
p2 <- DimPlot(PVM, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

Idents(PVM) <- "seurat_clusters"
Idents(PVM, cells = WhichCells(PVM, idents = c(0,3))) <- "PVM"
Idents(PVM, cells = WhichCells(PVM, idents = c(1))) <- "Mo-Mac 2"
Idents(PVM, cells = WhichCells(PVM, idents = c(2))) <- "Mo-Mac 1"
Idents(PVM, cells = WhichCells(PVM, idents = c(4))) <- "Contaminant"
PVM[['Annotation']] <- Idents(PVM)

PVM <- subset(PVM, idents = "Contaminant", invert = T)

# Convert to SCE
PVM_sce <- as.SingleCellExperiment(PVM, assay = "RNA")

# Remove ribosomal and mitochondrial genes
MT_genes <- grep("^MT-", rownames(PVM_sce), ignore.case = F, value = T)
RP_genes <- grep("^RP", rownames(PVM_sce), ignore.case = F, value = T)
PVM_sce <- PVM_sce[!rownames(PVM_sce) %in% MT_genes,]
PVM_sce <- PVM_sce[!rownames(PVM_sce) %in% RP_genes,]

# Run Slingshot on UMAP with Mo-Mac 1 as starting cluster and PVM as ending cluster
PVM_slingshot <- slingshot(PVM_sce, clusterLabels = "Annotation", reducedDim = "UMAP", start.clus = "Mo-Mac 1", end.clus = "PVM", approx_points = 150)

pseudo.paths <- slingPseudotime(PVM_slingshot)
head(pseudo.paths)
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

summary(PVM_slingshot$slingPseudotime_1)
print(SlingshotDataSet(PVM_slingshot))

##########################
# Plots
##########################
# Pseudotime PVM
png(file = paste0(fig_dir, "/", "SF5A.png"), res = 300, units = "in", height = 10, width = 10)
plotUMAP(PVM_slingshot, colour_by = I(shared.pseudo))
dev.off()

col_shuffle <- c("burlywood", "salmon", "darkolivegreen", "red", "blue", "darkorange", "dodgerblue", "coral", "gold", "tan",
"magenta", "cyan", "steelblue", "peru", "green", "orchid")

color <- data.frame(color = c("salmon", "darkolivegreen", "steelblue"), Annotation = c("Mo-Mac 1", "Mo-Mac 2", "PVM"))
metadata <- data.frame(colData(PVM_slingshot))
metadata <- inner_join(metadata, color, by = "Annotation")
PVM_slingshot$color <- metadata$color

png(file = paste0(fig_dir, "/", "SF5B.png"), res = 300, units = "in", height = 10, width = 10)
plot(reducedDims(PVM_slingshot)$UMAP, col = PVM_slingshot$color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topright", pch = c(20), 
       legend = c("Mo-Mac 1", "Mo-Mac 2", "PVM"), 
       col = c("salmon", "darkolivegreen", "steelblue"), bty = 'n')
lines(slingCurves(PVM_slingshot)$curve1, lwd=3)
dev.off()

# 4. PVM GAM------------------------------------------------------------------
##########################
# TradeSeq GAM Model
##########################
# Select PVM Trajectory
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Pseudotime/", "GLM")
dir.create(tmp_dir)

# Set up Model
counts <- as.matrix(assay(PVM_sce))
crv <- SlingshotDataSet(PVM_slingshot)

Genes.to.use.sce <- function(object, min.value = 1, min.cells = 200) {
    num.cells <- Matrix::rowSums(assays(object)$counts > min.value) # number of cells expressing > 1
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    return(genes.use)
}

genes.use <- Genes.to.use.sce(PVM_sce, min.cells = 500) # roughly 10% of cells expressing

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
saveRDS(sce, file = paste0(tmp_dir, "/", "PVM_TradeSeq_9.8.rds"))

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
Ann <- c(rep("Mo-Mac 1", 20), rep("Mo-Mac 2", 26), rep("PVM", 47), rep("Mixed", 7))
Arranged <- Smooth[c(900:1090,1:899),]
genes <- c("LYZ", "HLA-CRB1", "CPVL", "CST3", "CEBPD", "CCL2", "LYVE1", "C1QB", "CD36", "HMOX1", "CD14", "LIPA", "CD9", "MIF", "FABP4", "LSP1")
index <- which(rownames(Arranged) %in% genes)
label <- rownames(Arranged)[index]
row_ha <- rowAnnotation(foo = anno_mark(at = index, labels = label, which = "row", side = "right"))
ha = HeatmapAnnotation(Cell = Ann, col = list(Cell = c("Mo-Mac 1" = "salmon", "Mo-Mac 2" = "darkolivegreen", "PVM" = "steelblue", "Mixed" = "gold")))

png(paste0(fig_dir, "/", "PVM_GLM.png"), width = 10, height = 10, units = 'in', res = 600)
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

DF <- data.frame(barcode = colnames(PVM_slingshot), pseudotime = PVM_slingshot$slingPseudotime_1, annotation = PVM_slingshot$Annotation)
DF <- DF %>% mutate(time_bin = cut(pseudotime, breaks=100))
prop.table(table(DF$time_bin, DF$annotation), margin = 1) * 100
Ann <- c(rep("Mo-Mac 1", 20), rep("Mo-Mac 2", 26), rep("PVM", 47), rep("Mixed", 7))
Arranged <- Smooth[c(19:32,33:54,1:18),]

ha = HeatmapAnnotation(Cell = Ann, col = list(Cell = c("Mo-Mac 1" = "salmon", "Mo-Mac 2" = "darkolivegreen", "PVM" = "steelblue", "Mixed" = "gold")))
genes = c("EGR1", "EGR2", "ATF3", "JUN", "FOS", "MAF", "KLF4", "MAFB", "HIF1A")
index <- which(rownames(Arranged) %in% genes)
label <- rownames(Arranged)[index]
row_ha <- rowAnnotation(foo = anno_mark(at = index, labels = label, which = "row", side = "right"))

png(paste0(fig_dir, "/", "Figure6C.png"), width = 10, height = 10, units = 'in', res = 600)
Heatmap(Arranged, col = viridis(100), cluster_columns = F, cluster_rows = F, top_annotation = ha, right_annotation = row_ha, show_column_names = F, show_row_names = F)
dev.off()    
