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
# Figure 5A: Inter-Correlation: T cells with Macrophages
###-------------------####
Merged.prop <- Stromal.prop %>% left_join(Myeloid.prop[,c(1:15)], by = "HATIMID") %>% left_join(CD4.prop[1:6], by = "HATIMID") %>% left_join(CD8.prop[1:5], by = "HATIMID")
Merged.prop <- Merged.prop[, c(1:11,22:44,12:21)]
Cor_res <- Cor_rho_subset(df = Merged.prop, subset1 = c("CD4 TEM", "CD8 TEM"), subset2 = c("IM", "LAM", "PVM", "Mo-Mac 1", 'Mo-Mac 2', "cMo", "Other Mo", "nMo", "ISG+ Mo"))

plot <- Heatmap_subset(Cor_res, legend_breaks = c(-0.5, 0, 0.5), order = c("IM", "LAM", "PVM", "Mo-Mac 1", "Mo-Mac 2", "cMo", "Other Mo", "nMo", "ISG+ Mo"),
                       p_breaks = c(0.01, 0.05, 0.1, 0.5)) 
ggsave(paste0(fig_dir, "/", "Figure5A.png"), device = "png", height = 6, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 5B: Inter-Correlation: T cells with all Myeloid
###-------------------####
Cor_res <- Cor_rho_subset(df = Merged.prop, subset1 = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Regulatory", "CD4 Cytotoxic", "CD8 TCM", "CD8 TEM", "CD8 Cytotoxic", "CD8 Naive"), 
subset2 = c("IM", "LAM", "PVM", "Mo-Mac 1", 'Mo-Mac 2', "cMo", "Other Mo", "nMo", "ISG+ Mo", "cDC1", "cDC2B", "Migratory DC", "pDC"))

plot <- Heatmap_subset(Cor_res, legend_breaks = c(-0.5, 0, 0.5), order = c("IM", "LAM", "PVM", "Mo-Mac 1", 'Mo-Mac 2', "cMo", "Other Mo", "nMo", "ISG+ Mo", "cDC1", "cDC2B", "Migratory DC", "pDC"),
                       p_breaks = c(0.01, 0.05, 0.1, 0.5)) 
ggsave(paste0(fig_dir, "/", "Figure5B.png"), device = "png", height = 10, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 5C: Inter-Correlation: TEM & Adipogenic
###-------------------####
Cor_res <- Cor_rho_subset(df = Merged.prop, subset1 = c("CD4 TEM", "CD8 TEM"), subset2 = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", "Early Preadipocyte", 
"ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast"))

plot <- Heatmap_subset(Cor_res, legend_breaks = c(-0.4, 0, 0.4), plot.margin = margin(5.5, 25, 5.5, 5.5, "pt"), order = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", "Early Preadipocyte", 
"Mature Preadipocyte 1", "Mature Preadipocyte 2", "ECM-Producing Early Preadipocyte", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast")) 
ggsave(paste0(fig_dir, "/", "Figure5C.png"), device = "png", height = 7, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 5D: Inter-Correlation: Macrophage & Adipogenic
###-------------------####
Merged.mac <- Stromal.prop %>% left_join(Macrophage.prop[,c(1:6)], by = "HATIMID")
Merged.mac <- Merged.mac[,c(1:11, 22:26, 12:21)]

Cor_res <- Cor_rho_subset(df = Merged.mac, subset1 = c("LAM", "IM", "PVM", "Mo-Mac 1", "Mo-Mac 2"), subset2 = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", "Early Preadipocyte", 
                                                                                                                 "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast"))

plot <- Heatmap_subset(Cor_res, p_breaks = c(0.01, 0.05, 0.1, 0.5), legend_breaks = c(-0.5, 0, 0.5), plot.margin = margin(5.5, 25, 5.5, 5.5, "pt"), order = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", "Early Preadipocyte", 
        "Mature Preadipocyte 1", "Mature Preadipocyte 2", "ECM-Producing Early Preadipocyte", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast")) 
ggsave(paste0(fig_dir, "/", "Figure5D.png"), device = "png", height = 7, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 5E: Scatter Plot IM vs Myofibroblast
###-------------------####
plot <- scatter_plot(Merged.mac, x = "IM", y = "Myofibroblast", xlabel = "IM (% Macrophage)", ylabel = "Myofibroblast (% Stromal)", dotsize = 4)
ggsave(paste0(fig_dir, "/", "Figure5E.png"), device = "png", height = 8, width = 8, units = "in", dpi = 300)

