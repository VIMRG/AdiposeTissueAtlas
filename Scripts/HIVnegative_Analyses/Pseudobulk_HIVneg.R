###############################################################################
######----------------------------FIGURE 8--------------------------------#####
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
source('/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/CodeRevised/Utils.R')

###-------------------####
# Directories
###-------------------####
Subset_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/SubsetAnalysis_HIVneg"
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/FiguresRevised"
Prop_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/RevisedProportion" 
Pseudo_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/Pseudobulk"
HIV_dir <- paste(Pseudo_dir, "HIVneg", sep = "/")
dir.create(HIV_dir)

###-------------------####
# Load Seurat Object
###-------------------####
Stromal <- readRDS(paste0(Subset_dir, "/", "Stromal_HIVneg_4.1.rds"))
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage_HIVneg_4.1.rds"))
Myeloid <- readRDS(paste0(Subset_dir, "/", "Myeloid_HIVneg_4.1.rds"))
Lymphoid <- readRDS(paste0(Subset_dir, "/", "Lymphoid_HIVneg_4.1.rds"))
Vascular <- readRDS(paste0(Subset_dir, "/", "Vascular_HIVneg_4.1.rds"))
Bcells <- readRDS(paste0(Subset_dir, "/", "Bcells_HIVneg_4.1.rds"))
CD4 <- readRDS(paste0(Subset_dir, "/", "CD4_HIVneg_4.1.rds"))
CD8 <- readRDS(paste0(Subset_dir, "/", "CD8_HIVneg_4.1.rds"))


###-------------------####
# Macrophage DGE
###-------------------####
# HIV+ DM vs HIV- diabetic
Macrophage[['Pseudobulk']] <- "Macrophage"
agg.macrophage <- Aggregate(Macrophage)

# HIV+ DM vs HIV+ nonDM
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Macrophage", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Macrophage", StudyGroup2 = "HIV- diabetic")

Mac.HIV.HIVneg <- Pseudobulk(agg.macrophage, cluster = "Macrophage", exclude = c("3013"), StudyGroup2 = "HIV- diabetic")
write.csv(Mac.HIV.HIVneg, file = paste(HIV_dir, "Macrophage.HIVDMvsHIVnegDM.csv", sep = "/"))

###-------------------####
# Not Figure (Table S12)
###-------------------####
# HIV+ DM vs HIV- diabetic
CD4[['Pseudobulk']] <- "CD4"
agg.CD4 <- Aggregate(CD4)

# HIV+ DM vs HIV+ nonDM
Outlier <- Outlier_Eval(agg.CD4, cluster = "CD4", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.CD4, vst = Outlier, cluster = "CD4", StudyGroup2 = "HIV- diabetic")

CD4.HIV.HIVneg <- Pseudobulk(agg.CD4, cluster = "CD4", exclude = c(""), StudyGroup2 = "HIV- diabetic")
write.csv(CD4.HIV.HIVneg, file = paste(HIV_dir, "CD4.HIVDMvsHIVnegDM.csv", sep = "/"))

###-------------------####
# CD8
###-------------------####
# HIV+ DM vs HIV- diabetic
CD8[['Pseudobulk']] <- "CD8"
agg.CD8 <- Aggregate(CD8)

# HIV+ DM vs HIV+ nonDM
Outlier <- Outlier_Eval(agg.CD8, cluster = "CD8", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.CD8, vst = Outlier, cluster = "CD8", StudyGroup2 = "HIV- diabetic")

CD8.HIV.HIVneg <- Pseudobulk(agg.CD8, cluster = "CD8", exclude = c("5017", "5036"), StudyGroup2 = "HIV- diabetic")
write.csv(CD8.HIV.HIVneg, file = paste(HIV_dir, "CD8.HIVDMvsHIVnegDM.csv", sep = "/"))

###-------------------####
# Stromal
###-------------------####
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte", "Early Preadipocyte", "ECM-Producing Early Preadipocyte"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Progenitor"
Stromal[['Pseudobulk']] <- Idents(Stromal)
agg.Stromal <- Aggregate(Stromal)

# HIV+ DM vs HIV- diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Preadipocyte", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Preadipocyte", StudyGroup2 = "HIV- diabetic")

Stromal.HIV.HIVneg <- Pseudobulk(agg.Stromal, cluster = "Preadipocyte", exclude = c("5018", "5020", "5021", "5026"), StudyGroup2 = "HIV- diabetic") # Batch effect notable for this lane for unclear reasons
write.csv(Stromal.HIV.HIVneg, file = paste(HIV_dir, "Preadipocyte.HIVDM.HIVnegDM.csv", sep = "/"))

# Progenitor
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Progenitor", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Progenitor", StudyGroup2 = "HIV- diabetic")

Progenitor.HIV.HIVneg <- Pseudobulk(agg.Stromal, cluster = "Progenitor", exclude = c("3037", "3013", "5026"), StudyGroup2 = "HIV- diabetic")
write.csv(Progenitor.HIV.HIVneg, file = paste(HIV_dir, "Progenitor.HIVDMvsHIVnegDM.csv", sep = "/"))

# PCOLCE
Outlier <- Outlier_Eval(agg.Stromal, cluster = "PCOLCE+ Fibroblast", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "PCOLCE+ Fibroblast", StudyGroup2 = "HIV- diabetic")

PCOLCE.HIV.HIVneg <- Pseudobulk(agg.Stromal, cluster = "PCOLCE+ Fibroblast", exclude = c(""), StudyGroup2 = "HIV- diabetic")
write.csv(PCOLCE.HIV.HIVneg, file = paste(HIV_dir, "PCOLCE.HIVDMvsHIVnegDM.csv", sep = "/"))

# Myofibroblast
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Myofibroblast", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Myofibroblast", StudyGroup2 = "HIV- diabetic")

Myofibroblast.HIV.HIVneg <- Pseudobulk(agg.Stromal, cluster = "Myofibroblast", exclude = c("5034"), StudyGroup2 = "HIV- diabetic")
write.csv(Myofibroblast.HIV.HIVneg, file = paste(HIV_dir, "Myofibroblast.HIVDMvsHIVnegDM.csv", sep = "/"))

###-------------------####
# Not Figure (Table S12)
###-------------------####
# HIV+ DM vs HIV- diabetic
Idents(Vascular) <- "Annotation"
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Arterial EC", "Venous EC", "Capillary EC"))) <- "EC"
Vascular[['Pseudobulk']] <- Idents(Vascular)
agg.Vascular <- Aggregate(Vascular)

# EC
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EC", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EC", StudyGroup2 = "HIV- diabetic")

EC.HIV.HIVneg <- Pseudobulk(agg.Vascular, cluster = "EC", exclude = c("5026"), StudyGroup2 = "HIV- diabetic")
write.csv(EC.HIV.HIVneg, file = paste(HIV_dir, "EC.HIVDMvsHIVnegDM.csv", sep = "/"))

###-------------------####
# Not Figure (Table S12)
###-------------------####
# Monocytes
Idents(Myeloid) <- "Annotation"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("nMo", "cMo", "Other Mo"))) <- "Monocyte"
Myeloid[['Pseudobulk']] <- Idents(Myeloid)
agg.Myeloid <- Aggregate(Myeloid)

Outlier <- Outlier_Eval(agg.Myeloid, cluster = "Monocyte", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Myeloid, vst = Outlier, cluster = "Monocyte", StudyGroup2 = "HIV- diabetic")

Monocyte.HIV.HIVneg <- Pseudobulk(agg.Myeloid, cluster = "Monocyte", exclude = c("5019"), StudyGroup2 = "HIV- diabetic")
write.csv(Monocyte.HIV.HIVneg, file = paste(HIV_dir, "Monocyte.HIVDMvsHIVnegDM.csv", sep = "/"))


