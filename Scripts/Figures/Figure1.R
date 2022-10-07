###############################################################################
######----------------------------FIGURE 1--------------------------------#####
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
library(DiagrammeR)
library(rsvg)
library(DiagrammeRsvg)
library(anndata)

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Code/Utils.R')

###-------------------####
# Load Data
###-------------------####
dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Final_Integration"
Integrated <- readRDS(paste(dir, "Integrated.rds", sep = "/"))

###-------------------####
# ReAnnotate
###-------------------####
Idents(Integrated) <- "Annotation"
# Stromal
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Early Preadipocyte"))) <- 1
Idents(Integrated, cells = WhichCells(Integrated, idents = c("ECM-Producing Early Preadipocyte"))) <- 2
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Mature Preadipocyte 1", "Mature Preadipocyte 2"))) <- 3
Idents(Integrated, cells = WhichCells(Integrated, idents = c("ISG+ Preadipocyte"))) <- 4
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Metallothionein+ Preadipocyte"))) <- 5
Idents(Integrated, cells = WhichCells(Integrated, idents = c("MYOC+ Fibroblast"))) <- 6
Idents(Integrated, cells = WhichCells(Integrated, idents = c("PCOLCE+ Fibroblast"))) <- 7
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- 8
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Myofibroblast"))) <- 9
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling Myofibroblast"))) <- 10

# Vascular
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Venous EC"))) <- 11
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Ven-Cap EC"))) <- 12
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Capillary EC"))) <- 13
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Intermediate Capillary EC"))) <- 14
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Arterial EC"))) <- 15
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Arterial EndoMT-like"))) <- 16
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Capillary EndoMT-like"))) <- 17
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Venous EndoMT-like"))) <- 18
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Pericyte"))) <- 19
Idents(Integrated, cells = WhichCells(Integrated, idents = c("VSMC"))) <- 20
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling Vascular Cells"))) <- 21

# Monocyte
Idents(Integrated, cells = WhichCells(Integrated, idents = c("cMo"))) <- 22
Idents(Integrated, cells = WhichCells(Integrated, idents = c("nMo"))) <- 23
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Other Mo", "ISG+ Mo"))) <- 24
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Other Mac", "Mo-Mac"))) <- 25
Idents(Integrated, cells = WhichCells(Integrated, idents = c("PVM"))) <- 26
Idents(Integrated, cells = WhichCells(Integrated, idents = c("LAM"))) <- 27
Idents(Integrated, cells = WhichCells(Integrated, idents = c("IM"))) <- 28
Idents(Integrated, cells = WhichCells(Integrated, idents = c("cDC2B"))) <- 29
Idents(Integrated, cells = WhichCells(Integrated, idents = c("cDC1"))) <- 30
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Migratory DC"))) <- 31
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling Myeloid"))) <- 32
Idents(Integrated, cells = WhichCells(Integrated, idents = c("pDC"))) <- 33

# B cells
Idents(Integrated, cells = WhichCells(Integrated, idents = c("B Cell"))) <- 34
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Plasmablast"))) <- 35

# CD4 T Cells
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 Naive"))) <- 36
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 TCM"))) <- 37
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 TEM"))) <- 38
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 Regulatory"))) <- 39
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 Cytotoxic"))) <- 40
# CD8 T Cells
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD8 Naive"))) <- 41
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD8 TCM"))) <- 42
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD8 TEM"))) <- 43
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD8 Cytotoxic"))) <- 44
Idents(Integrated, cells = WhichCells(Integrated, idents = c("MAIT"))) <- 45
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Gamma Delta"))) <- 46

# NK
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD57 mNK"))) <- 47
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD16 mNK"))) <- 48
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Immature NK"))) <- 49
Idents(Integrated, cells = WhichCells(Integrated, idents = c("ILC"))) <- 50
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling T/NK"))) <- 51
Integrated[['Celltype_plot']] <- Idents(Integrated)

###-------------------####
# Figure 1B: Annotated UMAP Overall
###-------------------####
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Figures"
UMAP_FUN(Integrated, path = paste(fig_dir, "Figure1B_Collapsed.png", sep = "/"), color = col_shuffle, breaks = c(1:51), height = 8, width = 18,
label = c("Preadipocytes (1)", "ECM-Producing Preadipocyte (2)", "Mature Preadipocyte (3)", "ISG+ Preadipocyte (4)", "Metallothionein+ Preadipocyte (5)","MYOC+ Fibroblasts (6)", "PCOLCE+ Fibroblast (7)", "Adipose Progenitor Cell (8)",
"Myofibroblast (9)", "Cycling fibroblast (10)", "Venous EC (11)", "Ven-Cap EC (12)", "Capillary EC (13)", "Intermediate Capillary EC (14)", "Arterial EC (15)", "Arterial EndoMT-like (16)", "Capillary EndoMT-like (17)",
"Venous EndoMT-like (18)", "Pericyte (19)", "VSMC (20)", "Cycling Vascular Cells (21)", "cMo (22)", "nMo (23)", "Other Mo (24)", "Mo-Mac (25)", "PVM (26)", "LAM (27)", "IM (28)", "cDC2B (29)", "cDC1 (30)",
"Migratory DC (31)", "Cycling Myeloid (32)", "pDC (33)",
"B Cell (34)", "Plasmablast (35)", "CD4 Naive (36)", "CD4 TCM (37)", "CD4 TEM (38)", "CD4 Regulatory (39)", "CD4 Cytotoxic (40)", "CD8 Naive (41)", "CD8 TCM (42)", "CD8 TEM (43)",
"CD8 Cytotoxic (44)", "MAIT (45)", "Gamma Delta (46)", "CD57+ mNK (47)", "CD16+ mNK (48)", "Immature NK (49)", "ILC (50)", "Cycling T & NK (51)"))

###-------------------####
# Figure 1C: Harmony Integration
###-------------------####
DefaultAssay(Integrated) <- "RNA"
Integrated[['StudyGroup']] <- factor(Integrated$StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic"))
png(paste(fig_dir, "Figure1C.png", sep = "/"), width = 7, height = 5, units = 'in', res = 600)
DimPlot(Integrated, label = F, cols = c("#91C46C", "#7ABAF2", "#FFBE00"), group.by = 'StudyGroup', raster = F) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          plot.title = element_text(face = "bold", size = 16)) + NoLegend() + ggtitle("Harmony Integration")
dev.off()

###-------------------####
# Figure 1D: Gene Markers for Major Cell Types
###-------------------####
Genes <- c("CLDN5", "ACTA2", "CD68", "CD1C", "COL1A2", "CCDC80", "CD3E", "NKG7")

Feature_Plot <- function(seurat_object, genes) {
    Plot <- FeaturePlot(seurat_object, features = genes, raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) +
          theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14),
          plot.title = element_text(face = "bold", size = 20)) + ggtitle(genes)
    return(Plot)
}

Plot_list <- lapply(Genes, Feature_Plot, seurat_object = Integrated)
png(file = paste(fig_dir, "Figure1D.png", sep = "/"), units = 'in', height = 10, width = 20, res = 300)
marrangeGrob(Plot_list, nrow=2, ncol=4)
dev.off()

###-------------------####
# Figure 1E: Proportion Differences Major Groups
###-------------------####
Prop_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Proportion/"
Major.prop <- read.table(paste(Prop_dir, "All.Prop.txt", sep = "/"), check.names = FALSE)
Major.prop <- Major.prop %>% na.omit()
Major.prop <- Prop_Merge(Major.prop, data_hatim)
Prop_long <- Major.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

use_colors <- c("#91C46C", "#7ABAF2", "#FFBE00")

plot <- box_plot_fun_paper(Prop_long, plot_cols = use_colors, parent_cell_type = "All Cell Types", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 4)
ggsave(file = paste(fig_dir, "Figure1E.png", sep = "/"), dpi = 300, units = "in", height = 8, width = 16)

# Miscellaneous
FDR_Fun <- function(df1, df2, clusters) {
    Results_DM <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(Results_DM) <- c("Cluster", "pvalue", "StudyGroup")
    for (i in clusters) {
        new <- data.frame(Cluster = i, pvalue = df1$results$.ALL.$testresults[[i]]$P, StudyGroup = "HIV+ diabetic", NonDM.prop = round(median(df1$results$.ALL.$data[[i]]$`HIV+ non-diabetic`), digits = 1), prop2 = round(median(df1$results$.ALL.$data[[i]]$`HIV+ diabetic`), digits = 1))
        Results_DM <- rbind(Results_DM, new)
    }    
    Results_preDM <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(Results_preDM) <- c("Celltype", "pvalue", "StudyGroup")
    for (i in clusters) {
        new <- data.frame(Cluster = i, pvalue = df2$results$.ALL.$testresults[[i]]$P, StudyGroup = "HIV+ prediabetic", NonDM.prop = round(median(df2$results$.ALL.$data[[i]]$`HIV+ non-diabetic`), digits = 1), prop2 = round(median(df2$results$.ALL.$data[[i]]$`HIV+ prediabetic`), digits = 1))
        Results_preDM <- rbind(Results_preDM, new)
    }
    Results <- rbind(Results_DM, Results_preDM)
    Results$p.adj <- p.adjust(Results$pvalue, method = "BH")
    return(Results)
}

## Proportion Multiple Comparison Adjustment
# Macrophage
DM <- Hmisc::summaryM(Stromal + Vascular + Myeloid + Lymphoid ~ StudyGroup, data = Major.prop[Major.prop$StudyGroup %in% c("HIV+ diabetic", "HIV+ non-diabetic"),], overall = T, test = T)
preDM <- Hmisc::summaryM(Stromal + Vascular + Myeloid + Lymphoid ~ StudyGroup, data = Major.prop[Major.prop$StudyGroup %in% c("HIV+ prediabetic", "HIV+ non-diabetic"),], overall = T, test = T)
Results <- FDR_Fun(df1 = DM, df2 = preDM, clusters = c("Stromal", "Vascular", "Myeloid", "Lymphoid"))
write.csv(Results_Mac, file = paste(Prop_dir, "Macrophage.Proportion.Adjusted.csv", sep = "/"))


