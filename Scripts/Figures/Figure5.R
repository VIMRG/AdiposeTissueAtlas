###############################################################################
######----------------------------FIGURE 5--------------------------------#####
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
# Directories
###-------------------####
Subset_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/SubsetAnalysis"
Prop_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Proportion" 
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Figures/"

###-------------------####
# Load Seurat Object
###-------------------####
Vascular <- readRDS(paste0(Subset_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(Subset_dir, "/", "Stromal.rds"))
Myeloid <- readRDS(paste0(Subset_dir, "/", "Myeloid.rds"))
Lymphoid <- readRDS(paste0(Subset_dir, "/", "Lymphoid.rds"))
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage.rds"))
CD4 <- readRDS(paste0(Subset_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(Subset_dir, "/", "CD8.rds"))

###-------------------####
# Load Cell Proportions
###-------------------####
Lymphoid.Prop <- read.table(paste(Prop_dir, "Lymphoid.Prop.txt", sep = "/"), check.names = FALSE)
Myeloid.Prop <- read.table(paste(Prop_dir, "Myeloid.Prop.txt", sep = "/"), check.names = F)
Stromal.Prop <- read.table(paste(Prop_dir, "Stromal.Prop.txt", sep = "/"), check.names = F)
Vascular.Prop <- read.table(paste(Prop_dir, "Vascular.Prop.txt", sep = "/"), check.names = F)
CD4.Prop <- read.table(paste(Prop_dir, "CD4.Prop.txt", sep = "/"), check.names = F)
CD8.Prop <- read.table(paste(Prop_dir, "CD8.Prop.txt", sep = "/"), check.names = F)
Macrophage.Prop <- read.table(paste(Prop_dir, "Macrophage.Prop.txt", sep = "/"), check.names = F)
CD4.Num <- read.table(paste(Prop_dir, "CD4.Num.txt", sep = "/"), check.names = F)
CD8.Num <- read.table(paste(Prop_dir, "CD8.Num.txt", sep = "/"), check.names = F)
Macrophage.Num <- read.table(paste(Prop_dir, "Macrophage.Num.txt", sep = "/"), check.names = F)

# Remove NAs
Lymphoid.Prop <- Lymphoid.Prop %>% na.omit()
Myeloid.Prop <- Myeloid.Prop %>% na.omit()
Stromal.Prop <- Stromal.Prop %>% na.omit()
Vascular.Prop <- Vascular.Prop %>% na.omit()
CD4.Prop <- CD4.Prop %>% na.omit()
CD8.Prop <- CD8.Prop %>% na.omit()
Macrophage.Prop <- Macrophage.Prop %>% na.omit()

###-------------------####
# Merge with Metadata
###-------------------####
# Main Subsets
Myeloid.prop <- Prop_Merge(Myeloid.Prop, data_hatim)
Lymphoid.prop <- Prop_Merge(Lymphoid.Prop, data_hatim)
Stromal.prop <- Prop_Merge(Stromal.Prop, data_hatim)
Vascular.prop <- Prop_Merge(Vascular.Prop, data_hatim)

# Macrophage Subset
Macrophage.prop <- Prop_Merge(Macrophage.Prop, data_hatim)
colnames(Macrophage.Num) <- c("HATIMID", "Macnum")
Macrophage.Num$HATIMID <- as.factor(as.character(Macrophage.Num$HATIMID))

Macrophage.prop <- left_join(Macrophage.prop, Macrophage.Num, by = "HATIMID")
Macrophage.prop <- Macrophage.prop %>% dplyr::filter(Macnum > 30) # Filter > 30
Macrophage.prop <- Macrophage.prop %>% dplyr::select(-c(Macnum))

# CD4 Subset
CD4.prop <- Prop_Merge(CD4.Prop, data_hatim)
colnames(CD4.Num) <- c("HATIMID", "CD4num")
CD4.Num$HATIMID <- as.factor(as.character(CD4.Num$HATIMID))

CD4.prop <- left_join(CD4.prop, CD4.Num, by = "HATIMID")
CD4.prop <- CD4.prop %>% dplyr::filter(CD4num > 30) # Filter > 30
CD4.prop <- CD4.prop %>% dplyr::select(-c(CD4num))

# CD8 Subset
CD8.prop <- Prop_Merge(CD8.Prop, data_hatim)
colnames(CD8.Num) <- c("HATIMID", "CD8num")
CD8.Num$HATIMID <- as.factor(as.character(CD8.Num$HATIMID))

CD8.prop <- left_join(CD8.prop, CD8.Num, by = "HATIMID")
CD8.prop <- CD8.prop %>% dplyr::filter(CD8num > 30) # Filter > 30
CD8.prop <- CD8.prop %>% dplyr::select(-c(CD8num))

###-------------------####
# Figure 5A: Inter-Correlation: CD4 TEM & Myeloid
###-------------------####
# CD4 TEM: cDC1, Mo-Mac, Migratory DC, Other Mac, ISG+ Mo, Cycling Myeloid, LAM, PVM, IM, Other Mo
Cor_res <- Cor_rho_subset(df1 = Myeloid.prop, df2 = CD4.prop, subset = "CD4 TEM", HIV = T)
CD4TEM <- Cor_plot(Cor_res, Myeloid, title = "")

CD4TEM
ggsave(paste0(fig_dir, "/", "Figure5A.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 5B: Inter-Correlation: CD8 TEM & Myeloid
###-------------------####
# CD8 TEM (cDC1, cDC2B, Mo-Mac, mreg, Other Mac, cMo, ISG+ Mo, Cycling, LAM, IM, Other Mo
Cor_res <- Cor_rho_subset(df1 <- Myeloid.prop, df2 <- CD8.prop, subset = "CD8 TEM", HIV = T)
CD8TEM <- Cor_plot(Cor_res, Myeloid, title = "")

CD8TEM
ggsave(paste0(fig_dir, "/", "Figure5B.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 5C: Inter-Correlation: CD4 TEM & Adipogenic
###-------------------####
# CD4 TEM:ECM-Producing Early Preadipocyte, Early Preadipocyte, Myofibroblast
Cor_res <- Cor_rho_subset(df1 <- Stromal.prop, df2 <- CD4.prop, subset = "CD4 TEM", HIV = T)
CD4TEM_Stromal <- Cor_plot(Cor_res, Stromal, title = "")

CD4TEM_Stromal
ggsave(paste0(fig_dir, "/", "Figure5C.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 5D: Inter-Correlation: CD8 TEM & Adipogenic
###-------------------####
# CD8 TEM (ECM-Producing Early Preadipocyte, Early Preadipocyte, Adipose Progenitor Cell 1)
Cor_res <- Cor_rho_subset(df1 <- Stromal.prop, df2 <- CD8.prop, subset = "CD8 TEM", HIV = T)
CD8TEM_Stromal <- Cor_plot(Cor_res, Stromal, title = "")

CD8TEM_Stromal
ggsave(paste0(fig_dir, "/", "Figure5D.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 5E: Inter-Correlation: IM & Adipogenic
###-------------------####
# IM (Adipose Progenitor Cell 2, Mature Preadipocyte 1, Myofibroblast, PCOLCE+ Fibroblast)
Cor_res <- Cor_rho_subset(df1 <- Stromal.prop, df2 <- Macrophage.prop, subset = "IM", HIV = T)
IM<- Cor_plot(Cor_res, Stromal, title = "")

IM
ggsave(paste0(fig_dir, "/", "Figure5E.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 5F: Inter-Correlation: LAM & Adipogenic
###-------------------####
# LAM (MMature Preadipocyte 1, Adipose Progenitor Cell 2, Myofibroblast, PCOLCE+ Fibroblast)
Cor_res <- Cor_rho_subset(df1 <- Stromal.prop, df2 <- Macrophage.prop, subset = "LAM", HIV = T)
LAM<- Cor_plot(Cor_res, Stromal, title = "")

LAM
ggsave(paste0(fig_dir, "/", "Figure5F.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 5G: Inter-Correlation: PVM & Adipogenic
###-------------------####
# PVM (PCOLCE+, Myofibroblast, Adipose Progenitor Cell 2)
Cor_res <- Cor_rho_subset(df1 <- Stromal.prop, df2 <- Macrophage.prop, subset = "PVM", HIV = T)
PVM <- Cor_plot(Cor_res, Stromal, title = "")

PVM
ggsave(paste0(fig_dir, "/", "Figure5G.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)
