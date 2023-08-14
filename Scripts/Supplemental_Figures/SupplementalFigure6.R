###############################################################################
######-------------------SUPPLEMENTAL FIGURE 6----------------------------#####
###############################################################################

###-------------------####
# Import Libraries
###-------------------####
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
library(Hmisc)


###-------------------####
# Source Code
###-------------------####
source('../Utils.R')

###-------------------####
# Directories
###-------------------####
tmp_dir <- "../SubsetAnalysis"
Prop_dir <- "../Proportion"
fig_dir <- "../Figures"

###-------------------####
# Load Seurat Object
###-------------------####
Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage.rds"))
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8.rds"))

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
# Supplemental Figure 6A: Myeloid Correlation: Age and BMI
###-------------------####
Cor_res <- Cor_rho(Myeloid.prop[,-c(5:9)], Measures = c("BMI", "Age"), nonDM = F)
plot <- Heatmap_circle(Cor_res, legend_breaks = c(-0.2, 0, 0.2), order = c("cMo", "Other Mo", "nMo", "ISG+ Mo", "cDC1", "cDC2B", "Migratory DC", "pDC", "Cycling Myeloid"),
                       limit = c(-.25, .25), range(9,5))
ggsave(file = paste(fig_dir, "SF6A.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)


###-------------------####
# Supplemental Figure 6B: Stromal Correlation: Age and BMI
###-------------------####
Cor_res <- Cor_rho(Stromal.prop, Measures = c("BMI", "Age"), nonDM = F)
plot <- Heatmap_circle(Cor_res, legend_breaks = c(-0.3, 0, 0.3), plot.margin = margin(5.5, 25, 5.5, 5.5, "pt"), order = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2",
"Early Preadipocyte", "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast", "Cycling Myofibroblast"),
        limit = c(-.35, .35), range(9,5))
ggsave(file = paste(fig_dir, "SF6B.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 6, width = 10)

###-------------------####
# Supplemental Figure 6C: Vascular Correlation: Age and BMI
###-------------------####
Cor_res <- Cor_rho(Vascular.prop, Measures = c("BMI", "Age"), nonDM = F)
plot <- Heatmap_circle(Cor_res, legend_breaks = c(-0.4, 0, 0.4), order = c("VSMC 1", "VSMC 2", "Pericyte", "Arterial EC", "Capillary EC", "Venous EC"),
                       limit = c(-.5, .5))
ggsave(file = paste(fig_dir, "SF6C.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Supplemental Figure 6D: Inter-Correlation: T Cells & Adipogenic
###-------------------####
Cor_res <- Cor_rho_subset(df = Merged.prop, subset1 = c("CD4 Naive", "CD4 TCM", "CD4 Cytotoxic", "CD4 Regulatory", "CD8 Naive", "CD8 TCM", "CD8 Cytotoxic"), subset2 = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", "Early Preadipocyte", 
                                                                                                                                                                         "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast"))
plot <- Heatmap_subset(Cor_res, range = c(8,5), legend_breaks = c(-0.4, 0, 0.4),limits = c(-.45, .45), plot.margin = margin(5.5, 25, 5.5, 5.5, "pt"), order = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", "Early Preadipocyte", 
       "Mature Preadipocyte 1", "Mature Preadipocyte 2", "ECM-Producing Early Preadipocyte", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast")) 
ggsave(paste0(fig_dir, "/", "SF6D.png"), device = "png", height = 10, width = 10, units = "in", dpi = 300)

###-------------------####
# Supplemental Figure 6E: Inter-Correlation: Myeloid & Adipogenic
###-------------------####
Cor_res <- Cor_rho_subset(df = Merged.prop, subset1 = c("cMo", "nMo", "Other Mo", "ISG+ Mo", "cDC1", "cDC2B", "pDC", "Migratory DC"), subset2 = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", "Early Preadipocyte", 
          "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast"))

plot <- Heatmap_subset(Cor_res, legend_breaks = c(-0.4, 0, 0.4), range = c(13,5), limits = c(-.5, .5), plot.margin = margin(5.5, 25, 5.5, 5.5, "pt"), order = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2", "Early Preadipocyte", 
        "Mature Preadipocyte 1", "Mature Preadipocyte 2", "ECM-Producing Early Preadipocyte", "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast")) 
ggsave(paste0(fig_dir, "/", "SF6E.png"), device = "png", height = 10, width = 10, units = "in", dpi = 300)



