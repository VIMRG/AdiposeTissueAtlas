#############################################################################################
##------------------------- PSEUDOBULK ANALYSIS ---------------------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will evaluate DGE using Pseudobulk.
##############################################################################################

# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(tidyverse)
library(scuttle)
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(ggplot2)

set.seed(7412)

date = "1.7.23"

source("../Utils.R")

###############################
# Load seurat objects
###############################
tmp_dir <- paste0("../SubsetAnalysis")

Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage.rds"))
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8.rds"))

tmp_dir <- paste0("../", date, "/Pseudobulk")
dir.create(tmp_dir)
DM_dir <- paste(tmp_dir, "DiabeticVsNonDM", sep = "/")
dir.create(DM_dir)
PreDM_dir <- paste(tmp_dir, "PreDMvsNonDM", sep = "/")
dir.create(PreDM_dir)

# 2. MACROPHAGE PSEUDOBULK_____________________________________________________
###############################
# OVERALL MACROPHAGE
###############################
# HIV+ DM vs HIV+ non-diabetic
Macrophage[['Pseudobulk']] <- "Macrophage"
agg.macrophage <- Aggregate(Macrophage)

# HIV+ DM vs HIV+ nonDM
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Macrophage")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Macrophage")

Mac.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Macrophage", exclude = c("1109"))
write.csv(Mac.DM.NonDM, file = paste(DM_dir, "Macrophage.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Macrophage", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Macrophage", StudyGroup1 = "HIV+ prediabetic")

Mac.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Macrophage", exclude = c("1109"), StudyGroup1 = "HIV+ prediabetic")
write.csv(Mac.PreDM.NonDM, file = paste(PreDM_dir, "Macrophage.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# PVM
###############################
# Pseudobulk
Macrophage[['Pseudobulk']] <- Idents(Macrophage)
agg.macrophage <- Aggregate(Macrophage)

# HIV+ DM vs HIV+ non-DM
Outlier <- Outlier_Eval(agg.macrophage, cluster = "PVM")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "PVM")

PVM.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "PVM", exclude = c("1132"))
write.csv(PVM.DM.NonDM, file = paste(DM_dir, "PVM.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "PVM", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "PVM", StudyGroup1 = "HIV+ prediabetic")

PVM.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "PVM", exclude = c("1132", "1165"), StudyGroup1 = "HIV+ prediabetic")
write.csv(PVM.PreDM.NonDM, file = paste(PreDM_dir, "PVM.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# IM
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "IM")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "IM")

IM.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "IM", exclude = c("1132", "3013")) # Remove substantial outliers
write.csv(IM.DM.NonDM, file = paste(DM_dir, "IM.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "IM", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "IM", StudyGroup1 = "HIV+ prediabetic")

IM.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "IM", exclude = c("1132"), StudyGroup1 = "HIV+ prediabetic")
write.csv(IM.PreDM.NonDM, file = paste(PreDM_dir, "IM.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# LAM
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "LAM")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "LAM")

LAM.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "LAM", exclude = c("3014"))
write.csv(LAM.DM.NonDM, file = paste(DM_dir, "LAM.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "LAM", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "LAM", StudyGroup1 = "HIV+ prediabetic")

LAM.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "LAM", exclude = c("1137"), StudyGroup1 = "HIV+ prediabetic")
write.csv(LAM.PreDM.NonDM, file = paste(PreDM_dir, "LAM.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# Mo-Mac 1
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Mo-Mac 1")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Mo-Mac 1")

MoMac.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Mo-Mac 1", exclude = c(""))
write.csv(MoMac.DM.NonDM, file = paste(DM_dir, "MoMac1.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Mo-Mac 1", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Mo-Mac 1", StudyGroup1 = "HIV+ prediabetic")

MoMac.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Mo-Mac 1", exclude = c(""), StudyGroup1 = "HIV+ prediabetic")
write.csv(MoMac.PreDM.NonDM, file = paste(PreDM_dir, "MoMac1.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# Mo-Mac 2
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Mo-Mac 2")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Mo-Mac 2")

MoMac.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Mo-Mac 2", exclude = c("1109"))
write.csv(MoMac.DM.NonDM, file = paste(DM_dir, "MoMac2.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Mo-Mac 2", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Mo-Mac 2", StudyGroup1 = "HIV+ prediabetic")

MoMac.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Mo-Mac 2", exclude = c("1109"), StudyGroup1 = "HIV+ prediabetic")
write.csv(MoMac.PreDM.NonDM, file = paste(PreDM_dir, "MoMac2.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 3. CD4 PSEUDOBULK_____________________________________________________
###############################
# OVERALL CD4
###############################
CD4[['Pseudobulk']] <- "CD4"
agg.CD4 <- Aggregate(CD4)

# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.CD4, cluster = "CD4")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.CD4, vst = Outlier, cluster = "CD4")

CD4.DM.NonDM <- Pseudobulk(agg.CD4, cluster = "CD4", exclude = c("1167"))
write.csv(CD4.DM.NonDM, file = paste(DM_dir, "CD4.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.CD4, cluster = "CD4", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.CD4, vst = Outlier, cluster = "CD4", StudyGroup1 = "HIV+ prediabetic")

CD4.PreDM.NonDM <- Pseudobulk(agg.CD4, cluster = "CD4", exclude = "", StudyGroup1 = "HIV+ prediabetic")
write.csv(CD4.PreDM.NonDM, file = paste(PreDM_dir, "CD4.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 4. MONOCYTE PSEUDOBULK_____________________________________________________
Idents(Myeloid) <- "Annotation"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cMo", "nMo", "Other Mo", "ISG+ Mo"))) <- "Monocyte"
Myeloid[['Pseudobulk']] <- Idents(Myeloid)
agg.Myeloid <- Aggregate(Myeloid)

# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Myeloid, cluster = "Monocyte")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Myeloid, vst = Outlier, cluster = "Monocyte")

Mo.DM.NonDM <- Pseudobulk(agg.Myeloid, cluster = "Monocyte", exclude = c("1167", '3015', '3031'))
write.csv(Mo.DM.NonDM, file = paste(DM_dir, "Monocyte.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Myeloid, cluster = "Monocyte", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Myeloid, vst = Outlier, cluster = "Monocyte", StudyGroup1 = "HIV+ prediabetic")

Mo.PreDM.NonDM <- Pseudobulk(agg.Myeloid, cluster = "Monocyte", exclude = "1165", StudyGroup1 = "HIV+ prediabetic")
write.csv(Mo.PreDM.NonDM, file = paste(PreDM_dir, "Monocyte.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 5. CD8 PSEUDOBULK_____________________________________________________
###############################
# OVERALL CD8
###############################
# HIV+ DM vs HIV+ non-diabetic
CD8[['Pseudobulk']] <- "CD8"
agg.CD8 <- Aggregate(CD8)

Outlier <- Outlier_Eval(agg.CD8, cluster = "CD8")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.CD8, vst = Outlier, cluster = "CD8")

CD8.DM.NonDM <- Pseudobulk(agg.CD8, cluster = "CD8", exclude = c("3015", "3037"))
write.csv(CD8.DM.NonDM, file = paste(DM_dir, "CD8.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.CD8, cluster = "CD8", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.CD8, vst = Outlier, cluster = "CD8", StudyGroup1 = "HIV+ prediabetic")

CD8.PreDM.NonDM <- Pseudobulk(agg.CD8, cluster = "CD8", exclude = c("2110"), StudyGroup1 = "HIV+ prediabetic")
write.csv(CD8.PreDM.NonDM, file = paste(PreDM_dir, "CD8.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 6. VASCULAR PSEUDOBULK_____________________________________________________
###############################
# OVERALL ENDOTHELIAL
###############################
Idents(Vascular) <- "Annotation"
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Venous EC", "Arterial EC", "Capillary EC"))) <- "EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c("VSMC 1", "VSMC 2"))) <- "VSMC"
Vascular[['Pseudobulk']] <- Idents(Vascular)
agg.Vascular <- Aggregate(Vascular)

# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EC")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EC")

EC.DM.NonDM <- Pseudobulk(agg.Vascular, cluster = "EC", exclude = c("1113", "1182", "1141"))
write.csv(EC.DM.NonDM, file = paste(DM_dir, "EC.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EC", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EC", StudyGroup1 = "HIV+ prediabetic")

EC.PreDM.NonDM <- Pseudobulk(agg.Vascular, cluster = "EC", exclude = c("1113", "1141", "1137", "1182"), StudyGroup1 = "HIV+ prediabetic")
write.csv(EC.PreDM.NonDM, file = paste(PreDM_dir, "EC.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 7. STROMAL PSEUDOBULK_____________________________________________________
###############################
# Preadipocyte
###############################
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 1", "Mature Preadipocyte 2", "Early Preadipocyte", "ECM-Producing Early Preadipocyte"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Progenitor"
Stromal[['Pseudobulk']] <- Idents(Stromal)
agg.Stromal <- Aggregate(Stromal)

# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Preadipocyte")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Preadipocyte")

PreAd.DM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Preadipocyte", exclude = "")
write.csv(PreAd.DM.NonDM, file = paste(DM_dir, "Preadipocyte.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Preadipocyte", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Preadipocyte", StudyGroup1 = "HIV+ prediabetic")

Preadipocyte.PreDM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Preadipocyte", exclude = "", StudyGroup1 = "HIV+ prediabetic")
write.csv(Preadipocyte.PreDM.NonDM, file = paste(PreDM_dir, "Preadipocyte.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# Progenitor
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Progenitor")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Progenitor")

Progenitor.DM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Progenitor", exclude = "3013")
write.csv(Progenitor.DM.NonDM, file = paste(DM_dir, "Progenitor.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Progenitor", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Progenitor", StudyGroup1 = "HIV+ prediabetic")

Progenitor.PreDM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Progenitor", exclude = "1165", StudyGroup1 = "HIV+ prediabetic")
write.csv(Progenitor.PreDM.NonDM, file = paste(PreDM_dir, "Progenitor.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# PCOLCE+ Fibroblast
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "PCOLCE+ Fibroblast")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "PCOLCE+ Fibroblast")

PCOLCE.DM.NonDM <- Pseudobulk(agg.Stromal, cluster = "PCOLCE+ Fibroblast", exclude = "3031")
write.csv(PCOLCE.DM.NonDM, file = paste(DM_dir, "PCOLCE.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "PCOLCE+ Fibroblast", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "PCOLCE+ Fibroblast", StudyGroup1 = "HIV+ prediabetic")

PCOLCE.PreDM.NonDM <- Pseudobulk(agg.Stromal, cluster = "PCOLCE+ Fibroblast", exclude = c("1109", "1151"), StudyGroup1 = "HIV+ prediabetic")
write.csv(PCOLCE.PreDM.NonDM, file = paste(PreDM_dir, "PCOLCE.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# Myofibroblast
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Myofibroblast")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Myofibroblast")

Myofibroblast.DM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Myofibroblast", exclude = "1167")
write.csv(Myofibroblast.DM.NonDM, file = paste(DM_dir, "Myofibroblast.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Myofibroblast", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Myofibroblast", StudyGroup1 = "HIV+ prediabetic")

Myofibroblast.PreDM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Myofibroblast", exclude = c("2110", "1167"), StudyGroup1 = "HIV+ prediabetic")
write.csv(Myofibroblast.PreDM.NonDM, file = paste(PreDM_dir, "Myofibroblast.HIVpreDM.HIVnonDM.csv", sep = "/"))