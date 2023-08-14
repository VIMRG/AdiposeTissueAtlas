###############################################################################
######----------------------------SUPPLEMENTARY FIGURE 11--------------------------------#####
###############################################################################

###-------------------####
# Import Libraries
###-------------------####
library(ggplot2)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(PResiduals)
library(rms)
library(gridExtra)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(stringi)
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(ggplot2)
library(scuttle)
library(circlize)
library(CellChat)

###-------------------####
# Source Code
###-------------------####
source('../Utils.R')

###-------------------####
# Directories
###-------------------####
fig_dir <- "../Figures"
cellchat_dir <- "../CellChat"

###-------------------####
# Load Data
###-------------------####
cellchat_HIVneg <- readRDS(paste0(cellchat_dir, "/", "cellchat_HIVneg.rds")) # non-diabetic
cellchat_HIVpos <- readRDS(paste0(cellchat_dir, "/", "cellchat_HIVpos.rds")) # Glucose Intolerant

object.list <- list(HIVneg = cellchat_HIVneg, HIVpos = cellchat_HIVpos)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

###-------------------####
# Supplementary Figure 11A: Overview of Interactions
###-------------------####
png(paste0(fig_dir, "/", "SF11A.png"), units = 'in',height = 15, width = 25, res = 300)
gg1 <- netVisual_heatmap(cellchat, font.size = 12, cluster.rows = TRUE, cluster.cols = TRUE, width = 20, height = 20)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", font.size = 12, cluster.rows = TRUE, cluster.cols = TRUE, width = 20, height = 20)
gg1 + gg2
dev.off()

###-------------------####
# Supplementary Figure 11B: Macrophage Target
###-------------------####
png(paste0(fig_dir, "/", "SF11B.png"), units = 'in',height = 11, width = 14, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, targets.use = c(21,24,28,29,37), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, targets.use = c(21,24,28,29,37), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###-------------------####
# Supplementary Figure 11C: Stromal Target
###-------------------####
png(paste0(fig_dir, "/", "SF11C.png"), units = 'in',height = 10, width = 12, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, targets.use = c(1,18,26,30,31,34,36), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, targets.use = c(1,18,26,30,31,34,36), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###-------------------####
# Supplementary Figure 11D: Endothelial to Monocyte
###-------------------####
png(paste0(fig_dir, "/", "SF11D.png"), units = 'in',height = 7, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("Capillary EC", "Venous EC", "Arterial EC"), targets.use = c(17,32,33), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,sources.use = c("Capillary EC", "Venous EC", "Arterial EC"),  targets.use = c(17,32,33), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###-------------------####
# Supplementary Figure 11E: T Cells Signaling
###-------------------####
png(paste0(fig_dir, "/", "SF11E.png"), units = 'in',height = 7, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c(5,6,7,8,9,11,12,13,14), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,sources.use = c(5,6,7,8,9,11,12,13,14), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()
