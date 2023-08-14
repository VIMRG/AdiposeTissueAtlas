###############################################################################
######----------------------------FIGURE 4--------------------------------#####
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
# Directories
###-------------------####
Subset_dir <- "../SubsetAnalysis"
Prop_dir <- "../Proportion"
fig_dir <- "../Figures"

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
# Figure 4A: BMI and Age Spearman's: Myeloid
###-------------------####
Cor_res <- Cor_rho(Myeloid.prop[,c(5:9,15:24)], Measures = c("BMI", "Age"), nonDM = F)
plot <- Heatmap_circle(Cor_res, limits = c(-.4, 0.4), legend_breaks = c(-0.3, 0, 0.3), order = c("LAM", "IM", "PVM", "Mo-Mac 1", "Mo-Mac 2"))
ggsave(file = paste(fig_dir, "Figure4A.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 4B: BMI and Age Spearman's: CD4
###-------------------####
Cor_res <- Cor_rho(CD4.prop, Measures = c("BMI", "Age"), nonDM = F)
plot <- Heatmap_circle(Cor_res, limits = c(-.45, .45), legend_breaks = c(-0.4, 0, 0.4), order = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Cytotoxic", "CD4 Regulatory"),
                       p_breaks = c(0.05, 0.1, 0.4, 0.8), range = c(11,5))
ggsave(paste0(fig_dir, "/", "Figure4B.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 4C: BMI and Age Spearman's: CD8
###-------------------####
Cor_res <- Cor_rho(CD8.prop, Measures = c("BMI", "Age"), nonDM = F)
plot <- Heatmap_circle(Cor_res, limits = c(-.45, .45), legend_breaks = c(-0.4, 0, 0.4), order = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Cytotoxic"),
                       p_breaks = c(0.05, 0.1, 0.4, 0.8), range = c(11,5))
ggsave(paste0(fig_dir, "/", "Figure4C.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 4D: Sex OLR: Macrophage
###-------------------####
OR <- ORM_Model(Macrophage.prop) # PVM, IM, LAM
plot_Macrophage <- ORM_Plot(OR, title = "Macrophage")
ggsave(paste0(fig_dir, "/", "Figure4D.png"), device = "png", height = 5, width = 7, units = "in", dpi = 300)

###-------------------####
# Figure 4E: Sex OLR: CD4
###-------------------####
OR <- ORM_Model(CD4.prop) # CD4 Cytotoxic, CD4 TCM
plot_CD4 <- ORM_Plot(OR, title = "CD4")
ggsave(paste0(fig_dir, "/", "Figure4E.png"), device = "png", height = 5, width = 7, units = "in", dpi = 300)

###-------------------####
# Figure 4F: Sex OLR: CD8
###-------------------####
OR <- ORM_Model(CD8.prop) # CD8 Cytotoxic
plot_CD8 <- ORM_Plot(OR, title = "CD8")
ggsave(paste0(fig_dir, "/", "Figure4F.png"), device = "png", height = 5, width = 7, units = "in", dpi = 300)

###-------------------####
# Figure 4G: Sex OLR: Stromal
###-------------------####
OR <- ORM_Model(Stromal.prop) # Adipose Progenitor Cell 2, Early Preadipocyte, Mature Preadipocyte 1, Mature Preadipocyte 2, PCLOCE+ Fibroblast
plot_Stromal <- ORM_Plot(OR, title = "Stromal")
ggsave(paste0(fig_dir, "/", "Figure4G.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

