###############################################################################
######------------------SUPPLEMENTAL FIGURE 3-----------------------------#####
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
library(Hmisc)

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
# Supplemental Figure 3A: Continuous Measures Lymphoid
###-------------------####
Cor_res <- Cor_rho(Lymphoid.prop, adjusted = T, variable = "meta_fbg", nonDM = T)
FBG_Lymph <- Cor_plot(Cor_res, Lymphoid, title = "FBG")

# HbA1c Adjusted ()
Cor_res <- Cor_rho(Lymphoid.prop, adjusted = T, variable = "meta_hba1c", nonDM = T)
HbA1c_Lymph <- Cor_plot(Cor_res, Lymphoid, title = "HbA1c")

plot_grid(FBG_Lymph, HbA1c_Lymph)

ggsave(file = paste(fig_dir, "SF3A.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Supplemental Figure 3B: Continuous Measures of Glucose Adipogenic
###-------------------####
Cor_res <- Cor_rho(Stromal.prop, adjusted = T, variable = "meta_fbg", nonDM = T)
FBG_Stromal <- Cor_plot(Cor_res, Stromal, title = "FBG")

# HbA1c Adjusted ()
Cor_res <- Cor_rho(Stromal.prop, adjusted = T, variable = "meta_hba1c", nonDM = T)
HbA1c_Stromal <- Cor_plot(Cor_res, Stromal, title = "HbA1c")

plot_grid(FBG_Stromal, HbA1c_Stromal)
ggsave(file = paste(fig_dir, "SF3B.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Supplemental Figure 3C: Box Plot Fibroblast
###-------------------####
ECM.Prop <- Stromal.Prop %>% dplyr::select(`ECM-Producing Early Preadipocyte`, `PCOLCE+ Fibroblast`, `MYOC+ Fibroblast`, `Myofibroblast`)
ECM.prop <- Prop_Merge(ECM.Prop, data_hatim)
colnames(ECM.prop)[c(1:4)] <- c("ECM Early PreAd", "PCOLCE+ FIB", "MYOC+ FIB", "MyoFIB")
Prop_long <- ECM.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

Prop_list <- Prop_long %>% group_split(measure)
plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "Stromal", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 6)
                    
ggsave(file = paste(fig_dir, "SF3C.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 15)

###-------------------####
# Supplemental Figure 3D: Continuous Measures of Glucose Vascular
###-------------------####
Cor_res <- Cor_rho(Vascular.prop, adjusted = T, variable = "meta_fbg", nonDM = T)
FBG_Vascular <- Cor_plot(Cor_res, Vascular, title = "FBG")

# HbA1c Adjusted ()
Cor_res <- Cor_rho(Vascular.prop, adjusted = T, variable = "meta_hba1c", nonDM = T)
HbA1c_Vascular <- Cor_plot(Cor_res, Vascular, title = "HbA1c")

plot_grid(FBG_Vascular, HbA1c_Vascular)
ggsave(file = paste(fig_dir, "SF3D.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)
