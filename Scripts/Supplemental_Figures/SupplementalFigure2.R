###############################################################################
######------------------SUPPLEMENTAL FIGURE 2-----------------------------#####
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

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Code/Utils.R')

###-------------------####
# Load Data
###-------------------####
tmp_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/SubsetAnalysis"
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Figures/"
markers_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/ClusterMarkers/"

Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Bcells <- readRDS(paste0(tmp_dir, "/", "Bcells.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage.rds"))
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8.rds"))

###-------------------####
# Supplmental Figure 2A: Myeloid ADT
###-------------------####
ADT <- c("CD11C", "CD71", "CD16", "CD1C")
Avg_myeloid <- CITESEQ_Avg(Myeloid)
myeloid_prep <- CITESEQ_PlotPrep(seurat_object = Myeloid, CITESeq_Averages = Avg_myeloid, ADT = ADT)

DF_long <- myeloid_prep %>% pivot_longer(cols = -c(UMAP_1, UMAP_2),
            names_to = "measure",
            values_to = "value")
            
DF_list <- DF_long %>% group_split(measure)
Plot_list <- lapply(DF_list, CITESEQ_Plot)

png(file = paste(fig_dir, "SF2A.png", sep = "/"), res = 300, units = "in", height = 8, width = 10)
marrangeGrob(Plot_list, ncol = 2, nrow = 2)
dev.off()

###-------------------####
# Supplmental Figure 2B: Macrophage UMAP
###-------------------####
Idents(Macrophage) <- "Annotation"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("Mo-Mac"))) <- 1
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("PVM"))) <- 2
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("IM"))) <- 3
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("LAM"))) <- 4
Macrophage[['Celltype_plot']] <- Idents(Macrophage)

col_shuffle <- c("coral", "dodgerblue", "darkolivegreen", "gold", "blue")

UMAP_FUN(Macrophage, path = paste(fig_dir, "SF2B.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:4), height = 5, width = 7.5,
label = c("Mo-Mac (1)", "PVM (2)", "IM (3)", "LAM (4)"))

###-------------------####
# Supplmental Figure 2C: Macrophage Pathway Analysis
###-------------------####
markers_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/ClusterMarkers/"
Mac_markers <- read.csv(paste(markers_dir, "Macrophage.Markers.csv", sep = "/"))
Mac.KEGG.Plot <- KEGG_ORA(Mac_markers, celltype = c("IM", "LAM", "Mo-Mac", "PVM"), filename = paste(fig_dir, "SF2C.png", sep = "/"))

###-------------------####
# Supplmental Figure 2D: Lymphoid ADT
###-------------------####
Avg_Lymphoid <- CITESEQ_Avg(Lymphoid)
Lymphoid_prep <- CITESEQ_PlotPrep(seurat_object = Lymphoid, CITESeq_Averages = Avg_Lymphoid, ADT = c("CD4", "CD8", "CD27", "CD45RA", "CD57", "CD16"))

DF_long <- Lymphoid_prep %>% pivot_longer(cols = -c(UMAP_1, UMAP_2),
            names_to = "measure",
            values_to = "value")
            
DF_list <- DF_long %>% group_split(measure)
Plot_list <- lapply(DF_list, CITESEQ_Plot)
png(file = paste(fig_dir, "SF2D.png", sep = "/"), res = 300, units = "in", height = 5, width = 8.5)
marrangeGrob(Plot_list, ncol = 3, nrow = 2)
dev.off()

###-------------------####
# Supplemental Figure 2E: CD4 UMAP
###-------------------####
Idents(CD4) <- "Annotation"
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Naive"))) <- 1
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 TCM"))) <- 2
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 TEM"))) <- 3
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Regulatory"))) <- 4
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Cytotoxic"))) <- 5
CD4[['Celltype_plot']] <- Idents(CD4)

col_shuffle <- c("coral", "dodgerblue", "darkolivegreen", "gold", "magenta")

UMAP_FUN(CD4, path = paste(fig_dir, "SF2E.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:5), height = 5, width = 8,
label = c("CD4 Naive (1)", "CD4 TCM (2)", "CD4 TEM (3)", "CD4 Regulatory (4)", "CD4 Cytotoxic (5)"))

###-------------------####
# Supplmental Figure 2F: CD4 ADT
###-------------------####
Avg_CD4 <- CITESEQ_Avg(CD4)
CD4_prep <- CITESEQ_PlotPrep(seurat_object = CD4, CITESeq_Averages = Avg_CD4, ADT = c("CD27", "CD45RA", "CD57", "CD69"))

DF_long <- CD4_prep %>% pivot_longer(cols = -c(UMAP_1, UMAP_2),
            names_to = "measure",
            values_to = "value")
            
DF_list <- DF_long %>% group_split(measure)
Plot_list <- lapply(DF_list, CITESEQ_Plot)
png(file = paste(fig_dir, "SF2F.png", sep = "/"), res = 300, units = "in", height = 8, width = 10)
marrangeGrob(Plot_list, ncol = 2, nrow = 2)
dev.off()

###-------------------####
# Supplemental Figure 2G: CD8 UMAP
###-------------------####
Idents(CD8) <- "Annotation"
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 Naive"))) <- 1
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 TCM"))) <- 2
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 TEM"))) <- 3
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 Cytotoxic"))) <- 4
CD8[['Celltype_plot']] <- Idents(CD8)

col_shuffle <- c("coral", "dodgerblue", "darkolivegreen", "gold")

UMAP_FUN(CD8, path = paste(fig_dir, "SF2G.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:4), height = 5, width = 8,
label = c("CD8 Naive (1)", "CD8 TCM (2)", "CD8 TEM (3)", "CD8 Cytotoxic (4)"))

###-------------------####
# Supplmental Figure 2H: CD8 ADT
###-------------------####
Avg_CD8 <- CITESEQ_Avg(CD8)
CD8_prep <- CITESEQ_PlotPrep(seurat_object = CD8, CITESeq_Averages = Avg_CD8, ADT = c("CD27", "CD45RA", "CD57", "CD69"))

DF_long <- CD8_prep %>% pivot_longer(cols = -c(UMAP_1, UMAP_2),
            names_to = "measure",
            values_to = "value")
            
DF_list <- DF_long %>% group_split(measure)
Plot_list <- lapply(DF_list, CITESEQ_Plot)
png(file = paste(fig_dir, "SF2H.png", sep = "/"), res = 300, units = "in", height = 8, width = 10)
marrangeGrob(Plot_list, ncol = 2, nrow = 2)
dev.off()
