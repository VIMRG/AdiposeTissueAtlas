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
# Figure 4A: BMI and Age Spearman's: Myeloid
###-------------------####
# BMI Adjusted (Mo-Mac, pDC, LAM, IM)
Cor_res <- Cor_rho(Myeloid.prop, adjusted = T, variable = "BMI")
BMI_Myeloid <- Cor_plot(Cor_res, Myeloid, title = "BMI")

# Age Adjusted (LAM)
Cor_res <- Cor_rho(Myeloid.prop, adjusted = T, variable = "Age")
Age_Myeloid <- Cor_plot(Cor_res, Myeloid, title = "Age")

plot_grid(BMI_Myeloid, Age_Myeloid)
ggsave(paste0(fig_dir, "/", "Figure4A.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 4B: BMI and Age Spearman's: CD4
###-------------------####
# BMI Adjusted
Cor_res <- Cor_rho(CD4.prop, adjusted = T, variable = "BMI")
BMI_CD4 <- Cor_plot(Cor_res, CD4, title = "BMI")

# Age Adjusted (CD4 Regulatory, CD4 Naive)
Cor_res <- Cor_rho(CD4.prop, adjusted = T, variable = "Age")
Age_CD4 <- Cor_plot(Cor_res, CD4, title = "Age")

plot_grid(BMI_CD4, Age_CD4)
ggsave(paste0(fig_dir, "/", "Figure4B.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 4C: BMI and Age Spearman's: CD8
###-------------------####
# BMI Adjusted (CD8 TCM)
Cor_res <- Cor_rho(CD8.prop, adjusted = T, variable = "BMI")
BMI_CD8 <- Cor_plot(Cor_res, CD8, title = "BMI")

# Age Adjusted ()
Cor_res <- Cor_rho(CD8.prop, adjusted = T, variable = "Age")
Age_CD8 <- Cor_plot(Cor_res, CD8, title = "Age")

plot_grid(BMI_CD8, Age_CD8)
ggsave(paste0(fig_dir, "/", "Figure4C.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 4D: BMI and Age Spearman's: Stromal
###-------------------####
# BMI Adjusted
Cor_res <- Cor_rho(Stromal.prop, adjusted = T, variable = "BMI")
BMI_Stromal <- Cor_plot(Cor_res, Stromal, title = "BMI")

# Age Adjusted (APC 1)
Cor_res <- Cor_rho(Stromal.prop, adjusted = T, variable = "Age")
Age_Stromal <- Cor_plot(Cor_res, Stromal, title = "Age")

plot_grid(BMI_Stromal, Age_Stromal)
ggsave(paste0(fig_dir, "/", "Figure4D.png"), device = "png", height = 5, width = 10, units = "in", dpi = 300)

#plot_grid(BMI_Myeloid, Age_Myeloid, BMI_CD4, Age_CD4, BMI_CD8, Age_CD8, BMI_Stromal, Age_Stromal, ncol = 4)
#ggsave(paste0(fig_dir, "/", "Figure4A_D.png"), device = "png", height = 10, width = 20, units = "in", dpi = 300)

###-------------------####
# Figure 4E: Sex OLR: Macrophage
###-------------------####
OR <- ORM_Model(Macrophage.prop) # PVM, IM, LAM
plot_Macrophage <- ORM_Plot(OR, Macrophage, title = "")

plot_Macrophage
ggsave(paste0(fig_dir, "/", "Figure4E.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 4F: Sex OLR: CD4
###-------------------####
OR <- ORM_Model(CD4.prop) # CD4 Cytotoxic, CD4 TCM
plot_CD4 <- ORM_Plot(OR, CD4, title = "")

plot_CD4
ggsave(paste0(fig_dir, "/", "Figure4F.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 4G: Sex OLR: CD8
###-------------------####
OR <- ORM_Model(CD8.prop) # CD8 Cytotoxic
plot_CD8 <- ORM_Plot(OR, CD8, title = "")

plot_CD8
ggsave(paste0(fig_dir, "/", "Figure4G.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

###-------------------####
# Figure 4H: Sex OLR: Stromal
###-------------------####
OR <- ORM_Model(Stromal.prop) # Adipose Progenitor Cell 2, Early Preadipocyte, Mature Preadipocyte 1, Mature Preadipocyte 2, PCLOCE+ Fibroblast
plot_Stromal <- ORM_Plot(OR, Stromal, title = "")

plot_Stromal
ggsave(paste0(fig_dir, "/", "Figure4H.png"), device = "png", height = 5, width = 6.5, units = "in", dpi = 300)

