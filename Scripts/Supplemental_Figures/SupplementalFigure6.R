###############################################################################
######------------------SUPPLEMENTAL FIGURE 6-----------------------------#####
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
library(circlize)

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Code/Utils.R')

###-------------------####
# Directories
###-------------------####
Subset_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/SubsetAnalysis"
AllCells_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Dialogue/All/AllCells" 
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Figures/"

###-------------------####
# Load Seurat Object
###-------------------####
Stromal <- readRDS(paste0(Subset_dir, "/", "Stromal.rds"))
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage.rds"))
CD4 <- readRDS(paste0(Subset_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(Subset_dir, "/", "CD8.rds"))

# Process Seurat Objects to Match DIALOGUE
# CD4 T cells All (>30 cells)
# Get samples to use
metadata <- CD4@meta.data
HATIMID_CD4 <- metadata %>% dplyr::count(HATIMID) %>% filter(n > 30) %>% dplyr::select(HATIMID)
metadata <- CD8@meta.data
HATIMID_CD8 <- metadata %>% dplyr::count(HATIMID) %>% filter(n > 30) %>% dplyr::select(HATIMID)

CD4[["Annotation"]] <- "CD4"
CD4 <- subset(CD4, cells = colnames(CD4)[CD4$HATIMID %in% HATIMID_CD4$HATIMID])
CD4[['HATIMID']] <- droplevels(CD4$HATIMID)
CD4[['HATIMID']] <- as.character(CD4$HATIMID)

metadata <- CD4@meta.data
metadata <- left_join(metadata, data_hatim[,c("HATIMID", "hba1c")], by = "HATIMID")
hba1c <- metadata %>% dplyr::select(hba1c)
rownames(hba1c) <- colnames(CD4)
CD4 <- AddMetaData(CD4, hba1c)

# CD8 T cells
CD8[["Annotation"]] <- "CD8"
CD8 <- subset(CD8, cells = colnames(CD8)[CD8$HATIMID %in% HATIMID_CD8$HATIMID])
CD8[['HATIMID']] <- droplevels(CD8$HATIMID)
CD8[['HATIMID']] <- as.character(CD8$HATIMID)

metadata <- CD8@meta.data
metadata <- left_join(metadata, data_hatim[,c("HATIMID", "hba1c")], by = "HATIMID")
hba1c <- metadata %>% dplyr::select(hba1c)
rownames(hba1c) <- colnames(CD8)
CD8 <- AddMetaData(CD8, hba1c)

# Macrophage Collapse
metadata <- Macrophage@meta.data
HATIMID_Macrophage <- metadata %>% dplyr::count(HATIMID) %>% filter(n > 30) %>% dplyr::select(HATIMID)
Macrophage <- subset(Macrophage, cells = colnames(Macrophage)[Macrophage$HATIMID %in% HATIMID_Macrophage$HATIMID])
Macrophage[['HATIMID']] <- droplevels(Macrophage$HATIMID)
Macrophage[['HATIMID']] <- as.character(Macrophage$HATIMID)

metadata <- Macrophage@meta.data
metadata <- left_join(metadata, data_hatim[,c("HATIMID", "hba1c")], by = "HATIMID")
hba1c <- metadata %>% dplyr::select(hba1c)
rownames(hba1c) <- colnames(Macrophage)
Macrophage <- AddMetaData(Macrophage, hba1c)

# Stromal
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 1", "Mature Preadipocyte 2", "Early Preadipocyte", "ECM-Producing Early Preadipocyte"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Progenitor"
Stromal[['Annotation']] <- Idents(Stromal)
Stromal[['HATIMID']] <- droplevels(Stromal$HATIMID)
Stromal[['HATIMID']] <- as.character(Stromal$HATIMID)

metadata <- Stromal@meta.data
metadata <- left_join(metadata, data_hatim[,c("HATIMID", "hba1c")], by = "HATIMID")
hba1c <- metadata %>% dplyr::select(hba1c)
rownames(hba1c) <- colnames(Stromal)
Stromal <- AddMetaData(Stromal, hba1c)

###-------------------####
# Load Data
###-------------------####
R <- readRDS(paste(AllCells_dir, "Results.rds", sep = "/"))

###-------------------####
# Supplemental Figure 6A: PVM
###-------------------####
plot <- Heatmap_fun(seurat_object = Macrophage, subset = T, metadata = data_hatim, cluster = "PVM", genes = c("ATF4", "CCL3", "CCL4", "CCL2", "CCL8", "EGR1", "EGR2", "FOSB", "KLF6", "KLF4", "KLF10", 
"TNFAIP2", "TNFRSF12A", "TNF", "LIPA", "IFI27", "TREM2", "IL1RN", "ZNF331", "TIMP3", "STAT1"), 
filename = "SF6A.png",tmp_dir = fig_dir, R = R, MCP_dim = 1, height = 10, width = 15)

###-------------------####
# Supplemental Figure 6B: LAM
###-------------------####
plot <- Heatmap_fun(seurat_object = Macrophage, subset = T, metadata = data_hatim, cluster = "LAM", genes = c("KLF4", "MAFF", "TRIB1", "ADM", "HMOX1", "APOE", "JUNB", "AREG", 'ID3'), 
filename = "SF6B.png", tmp_dir = fig_dir, R = R, height = 10, width = 15)

###-------------------####
# Supplemental Figure 6C:IM
###-------------------####
plot <- Heatmap_fun(seurat_object = Macrophage, subset = T, metadata = data_hatim, cluster = "IM", genes = c("EGR1", "EGR2", "CEBPB", "HMOX1", "KLF2", "KLF4", "KLF6", "MAFF", "PLAU", "NR4A1", "TRIB1",
"APOE", "LY96", "TREM2", "TMEM140", "ZNF331"), filename = "SF6C.png", tmp_dir = fig_dir, R = R, height = 10, width = 15)

###-------------------####
#  Supplemental Figure 6D: Mo-Mac
###-------------------####
plot <- Heatmap_fun(seurat_object = Macrophage, subset = T, metadata = data_hatim, cluster = "Mo-Mac", genes = c("ATF4", "CCL2", "EGR1", "EGR2", "IER3", "JUN", "KLF10", "KLF4", "KLF6",
"MAFB", "MAPK6", "TNF", "TNFRSF12A", "S100A6", "NFKBIA", "APOE", "FABP4", "FTL", "IFI27", "TREM2", "JUNB", "MMP19"),tmp_dir = fig_dir, filename = "SF6D.png", R = R, MCP_dim = 1,height = 10, width = 15)

###-------------------####
# Supplemental Figure 6E: Preadipocyte
###-------------------####
plot <- Heatmap_fun(seurat_object = Stromal, subset = T, metadata = data_hatim, cluster = "Preadipocyte", genes = c("CEBPB", "GPX3", "GSN", "KLF2", "KLF6", "MYC", "HMOX1", "EGR1", "EGR2", 
"ZFP36", "COL14A1", "MIF","PCOLCE", "PDGFD", "PLAAT4", "POSTN", "FABP4", "TIMP1"),
filename = "SF6E.png", tmp_dir = fig_dir, R = R,height = 10, width = 15)
