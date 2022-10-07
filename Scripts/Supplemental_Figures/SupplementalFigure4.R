###############################################################################
######------------------SUPPLEMENTAL FIGURE 4-----------------------------#####
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
# Supplemental Figure 4A: BMI and Age Correlations: Lymphoid
###-------------------####
# BMI Adjusted (MAIT)
Cor_res <- Cor_rho(Lymphoid.prop, adjusted = T, variable = "BMI")
BMI_Lymph <- Cor_plot(Cor_res, Lymphoid, title = "BMI")

# Age Adjusted (MAIT)
Cor_res <- Cor_rho(Lymphoid.prop, adjusted = T, variable = "Age")
Age_Lymph <- Cor_plot(Cor_res, Lymphoid, title = "Age")

plot_grid(BMI_Lymph, Age_Lymph)
ggsave(paste0(fig_dir, "/", "SF4A.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

###-------------------####
# Supplemental Figure 4B: BMI and Age Correlations: Vascular
###-------------------####
# BMI Adjusted (IM Cap, Venous EndoMT, Venous EC)
Cor_res <- Cor_rho(Vascular.prop, adjusted = T, variable = "BMI")
BMI_Vascular <- Cor_plot(Cor_res, Vascular, title = "BMI")

# Age Adjusted
Cor_res <- Cor_rho(Vascular.prop, adjusted = T, variable = "Age")
Age_Vascular <- Cor_plot(Cor_res, Vascular, title = "Age")

plot_grid(BMI_Vascular, Age_Vascular)
ggsave(paste0(fig_dir, "/", "SF4B.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

###-------------------####
# Supplemental Figure 4C: OLR Sex: Myeloid
###-------------------####
OR <- ORM_Model(Myeloid.prop) # mreg, PVM
plot_Myeloid <- ORM_Plot(OR, Myeloid, title = "")

plot_Myeloid
ggsave(paste0(fig_dir, "/", "SF4C.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Supplemental Figure 4D: OLR Sex: Lymphoid
###-------------------####
OR <- ORM_Model(Lymphoid.prop) # CD4 Cytotoxic, CD4 Regulatory, Immature NK, CD57 mNK
plot_Lymphoid <- ORM_Plot(OR, Lymphoid, title = "")

plot_Lymphoid
ggsave(paste0(fig_dir, "/", "SF4D.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Supplemental Figure 4E: OLR Sex: Vascular
###-------------------####
OR <- ORM_Model(Vascular.prop) # none
plot_Vascular <- ORM_Plot(OR, Vascular, title = "")

plot_Vascular
ggsave(paste0(fig_dir, "/", "SF4E.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)


