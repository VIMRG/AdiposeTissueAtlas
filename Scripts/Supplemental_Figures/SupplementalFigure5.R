###############################################################################
######-------------------SUPPLEMENTAL FIGURE 5----------------------------#####
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
# Supplemental Figure 5A: Fibroblast Box Plots
###-------------------####
ECM.Prop <- Stromal.Prop %>% dplyr::select(`ECM-Producing Early Preadipocyte`, `PCOLCE+ Fibroblast`, `MYOC+ Fibroblast`, `Myofibroblast`)
ECM.prop <- Prop_Merge(ECM.Prop, data_hatim)
colnames(ECM.prop)[c(1:4)] <- c("ECM PreAd", "PCOLCE+ FIB", "MYOC+ FIB", "MyoFIB")
Prop_long <- ECM.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "Stromal", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 6)
                    
ggsave(file = paste(fig_dir, "SF5A.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Supplemental Figure 5B: Heatmap Stromal and glucose tolerance measures
###-------------------####
Cor_res <- Cor_rho(Stromal.prop, Measures = c("meta_hba1c", "meta_fbg"))
plot <- Heatmap_circle(df = Cor_res, order = c("PCOLCE+ Fibroblast", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2",
                                               "Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "ECM-Producing Early Preadipocyte",
                                               "Metallothionein+ Preadipocyte", "MYOC+ Fibroblast", "Myofibroblast", "Cycling Myofibroblast"),
                       legend_breaks = c(-0.2, 0, 0.2), p_breaks = c(0.05, 0.1, 0.5, 0.9), range = c(6,5))

ggsave(file = paste(fig_dir, "SF5B.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Supplemental Figure 5C: Vascular Box Plots
###-------------------####
Vascular.prop <- Vascular.prop[,-c(1)]
Prop_long <- Vascular.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "Vascular", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 7, legend.pos = "none")
                    
ggsave(file = paste(fig_dir, "SF5C.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 7, width = 20)

###-------------------####
# Supplemental Figure 5D: Vascular Continuous Measures Glucose
###-------------------####
Cor_res <- Cor_rho(Vascular.prop, Measures = c("meta_hba1c", "meta_fbg"))
plot <- Heatmap_circle(df = Cor_res, order = c("VSMC 1", "VSMC 2", "Pericyte", "Capillary EC", "Venous EC", "Arterial EC", "Cycling Vascular"),
                       legend_breaks = c(-0.2, 0, 0.2), p_breaks = c(0.05, 0.1, 0.4, 0.8), range = c(9,5))

ggsave(file = paste(fig_dir, "SF5D.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

