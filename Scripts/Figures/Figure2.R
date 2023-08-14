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
source('../Utils.R')

###-------------------####
# Load Data
###-------------------####
tmp_dir <- "../SubsetAnalysis"
fig_dir <- "../Figures"

Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Bcells <- readRDS(paste0(tmp_dir, "/", "Bcells.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))

###-------------------####
# Figure 2A: Myeloid UMAP
###-------------------####
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cMo"))) <- 1
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Other Mo"))) <- 2
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("ISG+ Mo"))) <- 3
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("nMo"))) <- 4
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Mo-Mac 1"))) <- 5
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Mo-Mac 2"))) <- 6
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
UMAP_FUN(Myeloid, path = paste(fig_dir, "Figure2A.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:14), height = 5, width = 10,
label = c("cMo (1)", "Other Mo (2)", "ISG+ Mo (3)", "nMo (4)", "Mo-Mac 1 (5)", "Mo-Mac 2 (6)", "PVM (7)", "IM (8)", "LAM (9)", "Cycling Myeloid (10)", "Migratory DC (11)", "cDC2B (12)", "cDC1 (13)",
"pDC (14)"))

###-------------------####
# Figure 2B: Myeloid Cannonical Markers
###-------------------####
Myeloid[['Annotation']] <- factor(Myeloid$Annotation, levels = c('cMo', 'nMo', 'Other Mo', 'ISG+ Mo', 'Mo-Mac 1', 'Mo-Mac 2', 'PVM', 'IM', 'LAM', 'cDC2B', 'Migratory DC', 'cDC1', 'pDC', 'Cycling Myeloid'))

DOTPLOT_FUN(Myeloid, path = paste(fig_dir, 'Figure2B.png', sep = "/"), features = c("S100A8", "VCAN", "CD14", "FCGR3A", "SERPINA1", "WARS", "APOBEC3A", "ISG15", "AREG",
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
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c("Cycling T & NK"))) <- 16
Lymphoid[['Celltype_plot']] <- Idents(Lymphoid)

col_shuffle <- c("burlywood", "salmon", "darkolivegreen", "red", "blue", "darkorange", "dodgerblue", "coral", "gold", "tan",
"magenta", "cyan", "steelblue", "peru", "green", "orchid")

UMAP_FUN(Lymphoid, path = paste(fig_dir, "Figure2C.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:16), height = 5, width = 10,
label = c("CD4 Naive (1)", "CD4 TCM (2)", "CD4 TEM (3)", "CD4 Regulatory (4)", "CD4 Cytotoxic (5)", "CD8 Naive (6)", "CD8 TCM (7)", "CD8 TEM (8)", "CD8 Cytotoxic (9)", "MAIT (10)", 
"Gamma Delta (11)", "CD57 mNK (12)", "CD16 mNK (13)", "Immature NK (14)", "ILC (15)", "Cycling Lymphoid (16)"))

###-------------------####
# Figure 2D: Lymphoid DotPlot
###-------------------####
Lymphoid[['Annotation']] <- factor(Lymphoid$Annotation, levels = c('CD4 Naive', 'CD4 TCM', 'CD4 TEM', 'CD4 Cytotoxic', 'CD4 Regulatory', 'CD8 Naive', 'CD8 TCM', 
'CD8 TEM', 'CD8 Cytotoxic', 'MAIT', 'Gamma Delta', 'CD16 mNK', 'CD57 mNK', 'Immature NK', 'ILC', "Cycling T & NK"))

DOTPLOT_FUN(Lymphoid, path = paste(fig_dir, 'Figure2D.png', sep = "/"), features = c("IL7R", "SELL", "LEF1",
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
Idents(Stromal, cells = WhichCells(Stromal, idents = c("PCOLCE+ Fibroblast"))) <- 8
Idents(Stromal, cells = WhichCells(Stromal, idents = c("MYOC+ Fibroblast"))) <- 9
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Myofibroblast"))) <- 10
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Cycling Myofibroblast"))) <- 11
Stromal[['Celltype_plot']] <- Idents(Stromal)

col_shuffle <- c("burlywood", "salmon", "darkolivegreen", "red", "blue", "darkorange", "dodgerblue", "coral", "gold", "tan",
"magenta", "cyan", "steelblue", "peru", "green", "orchid")

UMAP_FUN(Stromal, path = paste(fig_dir, "Figure2E.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:11), height = 5, width = 10,
label = c("Progenitor 1 (1)", "Progenitor 2 (2)", "Early PreAd (3)", "ECM Early PreAd (4)",
"Mature PreAd 1 (5)", "Mature PreAd 2 (6)", "MT+ PreAd (7)", "PCOLCE+ FIB (8)", "MYOC+ FIB (9)",
"MyoFIB (10)","Cycling MyoFIB (11)"))

###-------------------####
# Figure 2F: Stromal DotPlot
###-------------------####
Stromal[['Annotation']] <- factor(Stromal$Annotation, levels = c("PCOLCE+ Fibroblast", "MYOC+ Fibroblast", "Myofibroblast", "Cycling Myofibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", 
"Early Preadipocyte", "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "Metallothionein+ Preadipocyte"), labels = c("PCOLCE+ FIB", "MYOC+ FIB", "MyoFIB", 
"Cycling MyoFIB", "Progenitor 1", "Progenitor 2", "Early PreAd", "ECM Early PreAd", "Mature PreAd 1", "Mature PreAd 2", "MT+ PreAd"))

DOTPLOT_FUN(Stromal, path = paste(fig_dir, "Figure2F.png", sep = "/"), features = c("CD55", "DPP4", "PCOLCE2", "MYOC", "IGFBP7","POSTN", "TIMP1", "APOE", "MKI67", "CLU", "LUM", "PI16", 
"ZFP36", "EGR1", "CEBPD", "CXCL12", "CD9", "THY1", "COL6A1", "PTGDS", "C7", "CIDEC", "ADIPOQ","MT1X"))

###-------------------####
# Figure 2G: Vascular UMAP
###-------------------####
Idents(Vascular) <- "Annotation"
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Arterial EC"))) <- 1
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Capillary EC"))) <- 2
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Venous EC"))) <- 3
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Pericyte"))) <- 4
Idents(Vascular, cells = WhichCells(Vascular, idents = c("VSMC 1"))) <- 5
Idents(Vascular, cells = WhichCells(Vascular, idents = c("VSMC 2"))) <- 6
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Cycling Vascular"))) <- 7
Vascular[['Celltype_plot']] <- Idents(Vascular)

col_shuffle <- c("burlywood", "salmon", "darkolivegreen", "red", "blue", "darkorange", "dodgerblue", "coral", "gold", "tan",
"magenta", "cyan", "steelblue", "peru", "green", "orchid")

UMAP_FUN(Vascular, path = paste(fig_dir, "Figure2G.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:7), height = 5, width = 10,
label = c("Arterial EC (1)", "Capillary EC (2)", "Venous EC (3)", "Pericyte (4)", "VSMC 1 (5)", "VSMC 2 (6)", "Cycling Vascular (7)"))

###-------------------####
# Figure 2H: Vascular DotPlot
###-------------------####
Vascular[['Annotation']] <- factor(Vascular$Annotation, levels = c("Capillary EC", "Arterial EC", "Venous EC", "Pericyte", "VSMC 1", "VSMC 2", "Cycling Vascular"))

DOTPLOT_FUN(Vascular, path = paste(fig_dir, "Figure2H.png", sep = "/"), features = c("CLDN5", "GNG11", "CA4", "VWF", "GJA5", "GJA4", "HEY1", "SEMA3G", "ACKR1", "RGS5", "COX4I2",
"CSPG4", "ACTA2", "TAGLN", "NOTCH3", "MYH11", "MUSTN1", "PHLDA2", "MKI67"))
