###############################################################################
######----------------------------SUPPLEMENTAL FIGURE 9--------------------------------#####
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
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(ggplot2)
library(scuttle)

###-------------------####
# Source Code
###-------------------####
source('../Utils.R')

###-------------------####
# Directories
###-------------------####
Subset_dir <- "../SubsetAnalysis_HIVneg"
Final_dir <- "../Final_Integration_HIVneg/"
fig_dir <- "../FiguresRevised"
Prop_dir <- "../Proportion" 

###-------------------####
# Load Seurat Object
###-------------------####
Integrated <- readRDS(paste0(Final_dir, "/", "Integrated_HIVneg.rds"))
Stromal <- readRDS(paste0(Subset_dir, "/", "Stromal_HIVneg.rds"))
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage_HIVneg.rds"))
Myeloid <- readRDS(paste0(Subset_dir, "/", "Myeloid_HIVneg.rds"))
Lymphoid <- readRDS(paste0(Subset_dir, "/", "Lymphoid_HIVneg.rds"))
Vascular <- readRDS(paste0(Subset_dir, "/", "Vascular_HIVneg.rds"))
Bcells <- readRDS(paste0(Subset_dir, "/", "Bcells_HIVneg.rds"))
CD4 <- readRDS(paste0(Subset_dir, "/", "CD4_HIVneg.rds"))
CD8 <- readRDS(paste0(Subset_dir, "/", "CD8_HIVneg.rds"))

###-------------------####
# Load Cell Proportions
###-------------------####
Overall.prop <- read.table(paste(Prop_dir, "Major.CellType.prop.HIVneg.txt", sep = "/"), check.names = F)
Stromal.Prop <- read.table(paste(Prop_dir, "Stromal.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD4.Prop <- read.table(paste(Prop_dir, "CD4.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD8.Prop <- read.table(paste(Prop_dir, "CD8.Prop.HIVneg.txt", sep = "/"), check.names = F)
Vascular.Prop <- read.table(paste(Prop_dir, "Vascular.Prop.HIVneg.txt", sep = "/"), check.names = F)
Myeloid.Prop <- read.table(paste(Prop_dir, "Myeloid.Prop.HIVneg.txt", sep = "/"), check.names = F)
Lymphoid.Prop <- read.table(paste(Prop_dir, "Lymphoid.Prop.HIVneg.txt", sep = "/"), check.names = F)
Macrophage.Prop <- read.table(paste(Prop_dir, "Macrophage.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD4.Num <- read.table(paste(Prop_dir, "CD4.Num.HIVneg.txt", sep = "/"), check.names = F)
CD8.Num <- read.table(paste(Prop_dir, "CD8.Num.HIVneg.txt", sep = "/"), check.names = F)
Macrophage.Num <- read.table(paste(Prop_dir, "Macrophage.Num.HIVneg.txt", sep = "/"), check.names = F)

# Remove NAs
Overall.prop <- Overall.prop %>% na.omit()
Stromal.Prop <- Stromal.Prop %>% na.omit()
Myeloid.Prop <- Myeloid.Prop %>% na.omit()
Lymphoid.Prop <- Lymphoid.Prop %>% na.omit()
Vascular.Prop <- Vascular.Prop %>% na.omit()
CD4.Prop <- CD4.Prop %>% na.omit()
CD8.Prop <- CD8.Prop %>% na.omit()
Macrophage.Prop <- Macrophage.Prop %>% na.omit()

###-------------------####
# Merge with Metadata
###-------------------####
# Main Subsets
Overall.prop <- Prop_Merge(Overall.prop, data_hatim)
Stromal.prop <- Prop_Merge(Stromal.Prop, data_hatim)
Myeloid.prop <- Prop_Merge(Myeloid.Prop, data_hatim)
Vascular.prop <- Prop_Merge(Vascular.Prop, data_hatim)
Lymphoid.prop <- Prop_Merge(Lymphoid.Prop, data_hatim)

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
# Supplementary Figure 10A: Overall Integrated UMAP
###-------------------####
DefaultAssay(Integrated) <- "RNA"
Integrated[['StudyGroup']] <- factor(Integrated$StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic", "HIV- diabetic"))
png(paste(fig_dir, "SF10A.png", sep = "/"), width = 5, height = 5, units = 'in', res = 600)
DimPlot(Integrated, label = F, cols = c("#91C46C", "#7ABAF2", "#FFBE00", "coral"), group.by = 'StudyGroup', split.by = "StudyGroup", raster = F, ncol = 2) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          plot.title = element_text(face = "bold", size = 16)) + NoLegend()
dev.off()

###-------------------####
# Supplementary Figure 10B: CD4 T cells and Age/BMI
###-------------------####
Cor_rho_neg <- function(df, Measures) {
    Results <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(Results) <- c("rho", "pvalue", "Cluster", "Measure")
    numb <- dim(df)[2] - 10
    
    for (measure in Measures) {
        if (measure == "BMI") {
        for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,measure] ~ Sex  + Age, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i], Measure = measure)
                Results <- rbind(Results, res)
            }
        } else if (measure == "Age") {
             for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,measure] ~ Sex + BMI, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i], Measure = measure)
                Results <- rbind(Results, res)
            }
        }
    } 
    Results$padj <- p.adjust(Results$pvalue, method = "BH")
    return(Results)
}

df <- CD4.prop %>% filter(StudyGroup == "HIV- diabetic")
Cor_res <- Cor_rho_neg(df, Measures = c("BMI", "Age"))

plot <- Heatmap_circle(Cor_res, limits = c(-.65, .65), legend_breaks = c(-0.5, 0, 0.5), order = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Cytotoxic", "CD4 Regulatory"),
                       p_breaks = c(0.05, 0.1, 0.4), range = c(15,5))
ggsave(paste0(fig_dir, "/", "SF10B.png"), device = "png", height = 7, width = 10, units = "in", dpi = 300)

###-------------------####
# Supplementary Figure 10C: CD8 T cells and Age/BMI
###-------------------####
df <- CD8.prop %>% filter(StudyGroup == "HIV- diabetic")
Cor_res <- Cor_rho_neg(df, Measures = c("BMI", "Age"))

plot <- Heatmap_circle(Cor_res, limits = c(-.30, .30), legend_breaks = c(-0.25, 0, 0.25), order = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Cytotoxic"),
                       p_breaks = c(0.4, 0.9), range = c(9,5))
ggsave(paste0(fig_dir, "/", "SF10C.png"), device = "png", height = 7, width = 10, units = "in", dpi = 300)

###-------------------####
# Supplementary Figure 10D: Myeloid and Age/BMI
###-------------------####
df <- Myeloid.prop %>% filter(StudyGroup == "HIV- diabetic")
Cor_res <- Cor_rho_neg(df, Measures = c("BMI", "Age"))

plot <- Heatmap_circle(Cor_res, limits = c(-.50, .50), legend_breaks = c(-0.30, 0, 0.30), order = c("cMo", "nMo", "Other Mo", "ISG+ Mo", "Mo-Mac 1", "Mo-Mac 2", "PVM", "IM", "LAM", "cDC2B", "cDC1", "pDC", "Migratory DC", "Cycling Myeloid"),
                       p_breaks = c(0.4, 0.8), range = c(10,5))
ggsave(paste0(fig_dir, "/", "SF10D.png"), device = "png", height = 7, width = 10, units = "in", dpi = 300)

###-------------------####
# Supplementary Figure 10E: Lymphoid and Age/BMI
###-------------------####
df <- Lymphoid.prop %>% filter(StudyGroup == "HIV- diabetic")
Cor_res <- Cor_rho_neg(df, Measures = c("BMI", "Age"))

plot <- Heatmap_circle(Cor_res, limits = c(-.55, .55), legend_breaks = c(-0.30, 0, 0.30), order = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Cytotoxic", "CD4 Regulatory", "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Cytotoxic", "MT+ T cells",
"Gamma Delta", "MAIT", "CD57 mNK", "CD16 mNK", "Immature NK", "ILC", "Cycling T & NK"),
                       p_breaks = c(0.1, 0.4, 0.8), range = c(10,5))
ggsave(paste0(fig_dir, "/", "SF10E.png"), device = "png", height = 7, width = 10, units = "in", dpi = 300)
