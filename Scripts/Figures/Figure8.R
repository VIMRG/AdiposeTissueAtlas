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
source('/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Code/Utils.R')

###-------------------####
# Directories
###-------------------####
Subset_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/SubsetAnalysis_HIVneg"
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Figures/"
Prop_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Proportion" 
Pseudo_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Pseudobulk/"
HIV_dir <- paste(Pseudo_dir, "HIVneg", sep = "/")
dir.create(HIV_dir)

###-------------------####
# Load Seurat Object
###-------------------####
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
Stromal.Prop <- read.table(paste(Prop_dir, "Stromal.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD4.Prop <- read.table(paste(Prop_dir, "CD4.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD8.Prop <- read.table(paste(Prop_dir, "CD8.Prop.HIVneg.txt", sep = "/"), check.names = F)
Macrophage.Prop <- read.table(paste(Prop_dir, "Macrophage.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD4.Num <- read.table(paste(Prop_dir, "CD4.Num.HIVneg.txt", sep = "/"), check.names = F)
CD8.Num <- read.table(paste(Prop_dir, "CD8.Num.HIVneg.txt", sep = "/"), check.names = F)
Macrophage.Num <- read.table(paste(Prop_dir, "Macrophage.Num.HIVneg.txt", sep = "/"), check.names = F)

# Remove NAs
Stromal.Prop <- Stromal.Prop %>% na.omit()
CD4.Prop <- CD4.Prop %>% na.omit()
CD8.Prop <- CD8.Prop %>% na.omit()
Macrophage.Prop <- Macrophage.Prop %>% na.omit()

###-------------------####
# Merge with Metadata
###-------------------####
# Main Subsets
Stromal.prop <- Prop_Merge(Stromal.Prop, data_hatim)

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
# Figure 9A: Major Cell Type Proportions
###-------------------####
Vascular[['Major_CellType']] <- "Vascular"
Stromal[['Major_CellType']] <- "Fibro-adipogenic"
Bcells[['Major_CellType']] <- "Lymphoid"
Myeloid[['Major_CellType']] <- "Myeloid"
Lymphoid[['Major_CellType']] <- "Lymphoid"

Vasc_meta <- Vascular@meta.data %>% dplyr::select(HATIMID, Major_CellType)
Strom_meta <- Stromal@meta.data  %>% dplyr::select(HATIMID, Major_CellType)
Bcell_meta <- Bcells@meta.data  %>% dplyr::select(HATIMID, Major_CellType)
Myel_meta <- Myeloid@meta.data  %>% dplyr::select(HATIMID, Major_CellType)
Lymph_meta <- Lymphoid@meta.data  %>% dplyr::select(HATIMID, Major_CellType)

metadata <- Strom_meta %>% rbind(Vasc_meta) %>% rbind(Bcell_meta) %>% rbind(Myel_meta) %>% rbind(Lymph_meta)
write.table(prop.table(table(metadata$HATIMID, metadata$Major_CellType), margin = 1) * 100, file = paste(Prop_dir, "Major.CellType.prop.txt", sep = "/"))
Major.prop <- read.table(paste(Prop_dir, "Major.CellType.prop.txt", sep = "/"), check.names = F)
Major.prop <- Major.prop %>% na.omit()

Major.prop <- Prop_Merge(Major.prop, data_hatim)

Major.prop <- Major.prop %>% filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))
Prop_long <- Major.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

Prop_list <- Prop_long %>% group_split(measure)

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "All Cells", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ diabetic", "HIV- diabetic")), specify_x_order = NA, ncol = 5, plot_cols = c("#FFBE00", "red"))
                    
ggsave(file = paste(fig_dir, "Figure8A.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 7.5, width = 15)

###-------------------####
# Figure 8B: CD4 Proportion BoxPlot
###-------------------####
CD4.prop <- CD4.prop %>% filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))
Prop_long <- CD4.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

Prop_list <- Prop_long %>% group_split(measure)

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "CD4", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ diabetic", "HIV- diabetic")), specify_x_order = NA, ncol = 5, plot_cols = c("#FFBE00", "red"))
                    
ggsave(file = paste(fig_dir, "Figure8B.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 7.5, width = 15)

###-------------------####
# Figure 8C: Preadipocyte Transcriptional Pattern
###-------------------####
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 1", "Mature Preadipocyte 2", "Early Preadipocyte", "ECM-Producing Early Preadipocyte"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Progenitor"
Stromal[['Pseudobulk']] <- Idents(Stromal)
agg.Stromal <- Aggregate(Stromal)

# HIV+ DM vs HIV- diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Preadipocyte", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Preadipocyte", StudyGroup2 = "HIV- diabetic")

Stromal.HIV.HIVneg <- Pseudobulk(agg.Stromal, cluster = "Preadipocyte", exclude = c("5018", "5020", "5021", "5026"), StudyGroup2 = "HIV- diabetic") # Batch effect notable for this lane
write.csv(Stromal.HIV.HIVneg, file = paste(HIV_dir, "Stromal.HIVDM.HIVnegDM.csv", sep = "/"))

Stromal_volcano <- VOLCANO(resOrderedDF = Stromal.HIV.HIVneg, filename = paste(fig_dir, "Figure8C.png", sep = "/"), genes = c("BAG3", "ELN", "HSP90AA1", "HSPH1", "NID2", "HSPA8", "FSTL1", "SMOC2", "H19",
"CDH11", "MOCS2", "FBN1", "COL5A1", "IGF2", "VAMP2", "CSF1", "C3"), ylim = c(0,4))

###-------------------####
# Figure 8D: Preadipocyte Pathway Analysis
###-------------------####
pathways_kegg <- KEGG_GSEA(Stromal.HIV.HIVneg)
pathways_kegg <- GSEA_order(pathways_kegg)
Stromal_gsea_plot <- GSEA_plot(df = data.frame(pathways_kegg), filename = paste(fig_dir, "Figure8D.png", sep = "/"), num_path = 8)

###-------------------####
# Figure 9E: Dialogue
###-------------------####
# See Dialogue_HIVneg.R

###-------------------####
# Not Figure (Table S12)
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
# Not Figure (Table S12)
###-------------------####
# HIV+ DM vs HIV- diabetic
CD8[['Pseudobulk']] <- "CD8"
agg.CD8 <- Aggregate(CD8)

# HIV+ DM vs HIV+ nonDM
Outlier <- Outlier_Eval(agg.CD8, cluster = "CD8", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Sex')
Plot <- Outlier_Plot(agg.CD8, vst = Outlier, cluster = "CD8", StudyGroup2 = "HIV- diabetic")

CD8.HIV.HIVneg <- Pseudobulk(agg.CD8, cluster = "CD8", exclude = c("5017", "5036"), StudyGroup2 = "HIV- diabetic")
write.csv(CD8.HIV.HIVneg, file = paste(HIV_dir, "CD8.HIVDMvsHIVnegDM.csv", sep = "/"))

###-------------------####
# Not Figure (Table S12)
###-------------------####
# Progenitor
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Progenitor", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Progenitor", StudyGroup2 = "HIV- diabetic")

Progenitor.HIV.HIVneg <- Pseudobulk(agg.Stromal, cluster = "Progenitor", exclude = c(""), StudyGroup2 = "HIV- diabetic")
write.csv(Progenitor.HIV.HIVneg, file = paste(HIV_dir, "Progenitor.HIVDMvsHIVnegDM.csv", sep = "/"))

# PCOLCE
Outlier <- Outlier_Eval(agg.Stromal, cluster = "PCOLCE+ Fibroblast", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "PCOLCE+ Fibroblast", StudyGroup2 = "HIV- diabetic")

PCOLCE.HIV.HIVneg <- Pseudobulk(agg.Stromal, cluster = "PCOLCE+ Fibroblast", exclude = c("3031", "5017"), StudyGroup2 = "HIV- diabetic")
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
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Arterial EC", "Venous EC", "Ven-Cap EC", "Capillary EC", "Intermediate Capillary EC"))) <- "EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Arterial EndoMT-like", "Capillary EndoMT-like", "Venous EndoMT-like"))) <- "EndoMT"
Vascular[['Pseudobulk']] <- Idents(Vascular)
agg.Vascular <- Aggregate(Vascular)

# EC
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EC", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EC", StudyGroup2 = "HIV- diabetic")

EC.HIV.HIVneg <- Pseudobulk(agg.Vascular, cluster = "EC", exclude = c("3014"), StudyGroup2 = "HIV- diabetic")
write.csv(EC.HIV.HIVneg, file = paste(HIV_dir, "EC.HIVDMvsHIVnegDM.csv", sep = "/"))

# Endo
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EndoMT", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EndoMT", StudyGroup2 = "HIV- diabetic")

EndoMT.HIV.HIVneg <- Pseudobulk(agg.Vascular, cluster = "EndoMT", exclude = c("5015"), StudyGroup2 = "HIV- diabetic")
write.csv(EndoMT.HIV.HIVneg, file = paste(HIV_dir, "EndoMT.HIVDMvsHIVnegDM.csv", sep = "/"))

###-------------------####
# Not Figure (Table S12)
###-------------------####
# Monocytes
Idents(Myeloid) <- "Annotation"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("nMo", "cMo"))) <- "Monocyte"
Myeloid[['Pseudobulk']] <- Idents(Myeloid)
agg.Myeloid <- Aggregate(Myeloid)

Outlier <- Outlier_Eval(agg.Myeloid, cluster = "Monocyte", StudyGroup2 = "HIV- diabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.Myeloid, vst = Outlier, cluster = "Monocyte", StudyGroup2 = "HIV- diabetic")

Monocyte.HIV.HIVneg <- Pseudobulk(agg.Myeloid, cluster = "Monocyte", exclude = c("5019"), StudyGroup2 = "HIV- diabetic")
write.csv(Monocyte.HIV.HIVneg, file = paste(HIV_dir, "Monocyte.HIVDMvsHIVnegDM.csv", sep = "/"))

