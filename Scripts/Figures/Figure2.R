###############################################################################
######----------------------------FIGURE 2--------------------------------#####
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
# Figure 2A: Myeloid UMAP
###-------------------####
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cMo"))) <- 1
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Other Mo"))) <- 2
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("ISG+ Mo"))) <- 3
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("nMo"))) <- 4
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Other Mac"))) <- 5
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Mo-Mac"))) <- 6
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("PVM"))) <- 7
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("IM"))) <- 8
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("LAM"))) <- 9
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Cycling Myeloid"))) <- 10
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Migratory DC"))) <- 11
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cDC2B"))) <- 12
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cDC1"))) <- 13
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("pDC"))) <- 14
Myeloid[['Celltype_plot']] <- Idents(Myeloid)

col_shuffle <- c("#f78b8b","#cc3b3b","#9242c2", "#14e380", "#eda151", "#88b4f2", "#7273ba", "#58b555", "#e63c3c", "#c940c9", "#ab3f7c", "#79a38d",  "#e63c3c","#3041b0")
UMAP_FUN(Myeloid, path = paste(fig_dir, "Figure2A_Collapsed.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:14), height = 5, width = 10,
label = c("cMo (1)", "Other Mo (2)", "ISG+ Mo (3)", "nMo (4)", "Other Mac (5)", "Mo-Mac (6)", "PVM (7)", "IM (8)", "LAM (9)", "Cycling Myeloid (10)", "Migratory DC (11)", "cDC2B (12)", "cDC1 (13)",
"pDC (14)"))

###-------------------####
# Figure 2B: Myeloid Cannonical Markers
###-------------------####
Myeloid[['Annotation']] <- factor(Myeloid$Annotation, levels = c('cMo', 'nMo', 'Other Mo', 'ISG+ Mo', 'Mo-Mac', 'Other Mac', 'PVM', 'IM', 'LAM', 'cDC2B', 'Migratory DC', 'cDC1', 'pDC', 'Cycling Myeloid'))

DOTPLOT_FUN(Myeloid, path = paste(fig_dir, 'Figure2B.png', sep = "/"), features = c("S100A8", "VCAN", "CD14", "FCGR3A", "SERPINA1", "WARS", "APOBEC3A", "AREG", "ISG15",
'C1QB', 'CD68', "CD163", 'HLA-DRA', 'FCN1', "RNASE1", "LYVE1", 'CXCL2', 'CXCL3', 'TREM2', 'CD9', 'CD1C', 'CLEC10A', 'CCR7', 'CLEC9A', 'LILRA4', 'MKI67'))

###-------------------####
# Figure 2C: Lymphoid UMAP
###-------------------####
Idents(Lymphoid) <- "Annotation"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD4 Naive"))) <- 1
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD4 TCM"))) <- 2
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD4 TEM"))) <- 3
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD4 Regulatory"))) <- 4
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD4 Cytotoxic"))) <- 5
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD8 Naive"))) <- 6
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD8 TCM"))) <- 7
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD8 TEM"))) <- 8
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD8 Cytotoxic"))) <- 9
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("MAIT"))) <- 10
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("Gamma Delta"))) <- 11
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD57 mNK"))) <- 12
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("CD16 mNK"))) <- 13
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("Immature NK"))) <- 14
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("ILC"))) <- 15
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("Cycling T/NK"))) <- 16
Lymphoid[['Celltype_plot']] <- Idents(Lymphoid)

col_shuffle <- c("burlywood", "salmon", "darkolivegreen", "red", "blue", "darkorange", "dodgerblue", "coral", "gold", "tan",
"magenta", "cyan", "steelblue", "peru", "green", "orchid")

UMAP_FUN(Lymphoid, path = paste(fig_dir, "Figure2D.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:16), height = 5, width = 10,
label = c("CD4 Naive (1)", "CD4 TCM (2)", "CD4 TEM (3)", "CD4 Regulatory (4)", "CD4 Cytotoxic (5)", "CD8 Naive (6)", "CD8 TCM (7)", "CD8 TEM (8)", "CD8 Cytotoxic (9)", "MAIT (10)", 
"Gamma Delta (11)", "CD57 mNK (12)", "CD16 mNK (13)", "Immature NK (14)", "ILC (15)", "Cycling Lymphoid (16)"))

###-------------------####
# Figure 2D: Lymphoid DotPlot
###-------------------####
Lymphoid[['Annotation']] <- factor(Lymphoid$Annotation, levels = c('CD4 Naive', 'CD4 TCM', 'CD4 TEM', 'CD4 Cytotoxic', 'CD4 Regulatory', 'CD8 Naive', 'CD8 TCM', 
'CD8 TEM', 'CD8 Cytotoxic', 'MAIT', 'Gamma Delta', 'CD16 mNK', 'CD57 mNK', 'Immature NK', 'ILC', "Cycling T & NK"))

DOTPLOT_FUN(Lymphoid, path = paste(fig_dir, 'Figure2D.png', sep = "/"), height = 10, width = 15, features = c("IL7R", "SELL", "LEF1",
"CCR7", "LTB", "LDHB", "GPR183", "COTL1", "ALOX5AP", "CCL5", "NKG7", "GNLY", "PRF1", "FOXP3", "CD8A", "CCR6", "TRDV1", "KLRC2", "FCER1G", "FCGR3A", "XCL1", "KIT", "GATA2", "SPINK2", "MKI67"))

###-------------------####
# Figure 2E: Stromal UMAP
###-------------------####
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1"))) <- 1
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 2"))) <- 2
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Early Preadipocyte"))) <- 3
Idents(Stromal, cells = WhichCells(Stromal, idents = c("ECM-Producing Early Preadipocyte"))) <- 4
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 1"))) <- 5
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 2"))) <- 6
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Metallothionein+ Preadipocyte"))) <- 7
Idents(Stromal, cells = WhichCells(Stromal, idents = c("ISG+ Preadipocyte"))) <- 8
Idents(Stromal, cells = WhichCells(Stromal, idents = c("PCOLCE+ Fibroblast"))) <- 9
Idents(Stromal, cells = WhichCells(Stromal, idents = c("MYOC+ Fibroblast"))) <- 10
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Myofibroblast"))) <- 11
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Cycling Myofibroblast"))) <- 12
Stromal[['Celltype_plot']] <- Idents(Stromal)

col_shuffle <- c("burlywood", "salmon", "darkolivegreen", "red", "blue", "darkorange", "dodgerblue", "coral", "gold", "tan",
"magenta", "cyan", "steelblue", "peru", "green", "orchid")

UMAP_FUN(Stromal, path = paste(fig_dir, "Figure2E.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:12), height = 5, width = 10,
label = c("Adipose Progenitor Cell 1 (1)", "Adipose Progenitor Cell 2 (2)", "Early Preadipocyte (3)", "ECM-Producing Early Preadipocyte (4)",
"Mature Preadipocyte 1 (5)", "Mature Preadipocyte 2 (6)", "Metallothionein+ Preadipocyte (7)", "ISG+ Preadipocyte (8)", "PCOLCE+ Fibroblast (9)", "MYOC+ Fibroblast (10)",
"Myofibroblast (11)","Cycling Myofibroblast (12)"))

###-------------------####
# Figure 2F: Stromal DotPlot
###-------------------####
Stromal[['Annotation']] <- factor(Stromal$Annotation, levels = c("PCOLCE+ Fibroblast", "MYOC+ Fibroblast", "Myofibroblast", "Cycling Myofibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", 
"Early Preadipocyte", "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "ISG+ Preadipocyte", "Metallothionein+ Preadipocyte"))

DOTPLOT_FUN(Stromal, path = paste(fig_dir, "Figure2F.png", sep = "/"), height = 10, width = 15, features = c("CD55", "DPP4", "PCOLCE2", "MYOC", "ITM2A",
"IGFBP7","POSTN", "TIMP1", "APOE", "MKI67", "DCN", "CLU", "LUM", "GSN", "CCN5", "PI16", "ZFP36", "MYC", "FOS", "CEBPD", "CXCL12", "CXCL14", "EGR1", "ATF3", "CD9", "THY1", "COL3A1", "COL6A1", "FABP4", "C7", "CIDEC", "LPL", "ADIPOQ",
"ISG15", "IFI6", "MT1X", "MT2A"))

###-------------------####
# Figure 2G: Vascular UMAP
###-------------------####
Idents(Vascular) <- "Annotation"
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Arterial EC"))) <- 1
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Intermediate Capillary EC"))) <- 2
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Capillary EC"))) <- 3
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Ven-Cap EC"))) <- 4
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Venous EC"))) <- 5
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Arterial EndoMT-like"))) <- 6
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Capillary EndoMT-like"))) <- 7
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Venous EndoMT-like"))) <- 8
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Pericyte"))) <- 9
Idents(Vascular, cells = WhichCells(Vascular, idents = c("VSMC"))) <- 10
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Cycling Vascular Cells"))) <- 11
Vascular[['Celltype_plot']] <- Idents(Vascular)

col_shuffle <- c("burlywood", "salmon", "darkolivegreen", "red", "blue", "darkorange", "dodgerblue", "coral", "gold", "tan",
"magenta", "cyan", "steelblue", "peru", "green", "orchid")

UMAP_FUN(Vascular, path = paste(fig_dir, "Figure2G.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:11), height = 5, width = 10,
label = c("Arterial EC (1)", "Intermediate Capillary EC (2)", "Capillary EC (3)", "Ven-Cap (4)", "Venous EC (5)", "Arterial EndoMT-like (6)", "Capillary EndoMT-like (7)", 
"Venous EndoMT-like (8)", "Pericyte (9)", "VSMC (10)", "Cycling Vascular (11)"))

###-------------------####
# Figure 2H: Vascular DotPlot
###-------------------####
Vascular[['Annotation']] <- factor(Vascular$Annotation, levels = c("Capillary EC", "Intermediate Capillary EC", "Arterial EC", "Venous EC", "Ven-Cap EC", "Capillary EndoMT-like", "Arterial EndoMT-like", "Venous EndoMT-like",
"Pericyte", "VSMC", "Cycling Vascular Cells"))

DOTPLOT_FUN(Vascular, path = paste(fig_dir, "Figure2H.png", sep = "/"), height = 10, width = 15, features = c("CLDN5", "GNG11", "CA4", "VWF", "GJA5", "GJA4", "HEY1", "SEMA3G", "ACKR1", "RGS5", "COX4I2",
"CSPG4", "SNAI1", "ACTA2", "TAGLN", "NOTCH3", "MYH11", "MKI67"))