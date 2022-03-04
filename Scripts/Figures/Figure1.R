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


###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Figure_Functions.R')

###-------------------####
# Load Data
###-------------------####
Integrated <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Final_Integrated.1.3.rds')

###-------------------####
# Proportion
###-------------------####
#write.table(prop.table(table(Integrated$HATIMID, Integrated$Major_CellType), margin = 1) * 100, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Integrated.HATIM.by.Major.CellType.1.3.txt')

Prop_Long <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Proportions_Long.csv')
Prop <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Merged_Proportions.csv')

###-------------------####
# Set ggplot theme
###-------------------####
theme_set(theme_bw() +
            theme(axis.text = (element_text(size = 16)),
                  plot.caption = element_text(face = "italic", size = 9),
                  legend.position = "bottom",
                  legend.background = element_rect(fill = "transparent"), 
                  legend.text = element_text(size = 18),
                  panel.background = element_rect(fill = "transparent"),
                  panel.border = element_blank(),
                  axis.line = element_line(color="black"),
                  strip.text = element_text(size = 24), 
                  axis.title = element_text(size = 18),
                  #plot.margin=unit(c(0,0,0,2.5),"cm"),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)
          ))

###-------------------####
# Environmental Parameters
###-------------------####
set.seed(7612)

#########################
# Figure 1B: Integrated UMAP (Consolidated)
#########################
Idents(Integrated) <- "Manual_Annotation"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Early Preadipocyte 1", "Early Preadipocyte 2", "Early Preadipocyte 3", "Early Preadipocyte 4", "Preadipocyte", "PTGDS+ Preadipocyte"))) <- 1
Idents(Integrated, cells = WhichCells(Integrated, idents = c("PTGDS+ Extracellular Preadipocyte", "Extracellular Preadipocyte", "CD9hi Preadipocytes"))) <- 2
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Mature Preadipocyte 1", "Mature Preadipocyte 2"))) <- 3
Idents(Integrated, cells = WhichCells(Integrated, idents = c("ISG+ Preadipocyte"))) <- 4
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Metallothionein+ Preadipocyte/ASC"))) <- 5
Idents(Integrated, cells = WhichCells(Integrated, idents = c("MYOC+ Fibroblast"))) <- 6
Idents(Integrated, cells = WhichCells(Integrated, idents = c("PCOLCE+ Fibroblast"))) <- 7
Idents(Integrated, cells = WhichCells(Integrated, idents = c("ASC 1", "ASC 2", "ASC 3"))) <- 8
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Myofibroblast", "Lipomyofibroblast"))) <- 9
Idents(Integrated, cells = WhichCells(Integrated, idents = c("DKK3+ Fibroblast"))) <- 10
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Proliferating Myofibroblast"))) <- 11

Idents(Integrated, cells = WhichCells(Integrated, idents = c("Venous"))) <- 12
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Capillary"))) <- 13
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Intermediate Capillary"))) <- 14
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Arterial"))) <- 15
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Arterial-EndoMT-like"))) <- 16
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Capillary-EndoMT-like"))) <- 17
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Venous-EndoMT-like"))) <- 18
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Pericyte"))) <- 19
Idents(Integrated, cells = WhichCells(Integrated, idents = c("VSMC"))) <- 20

Idents(Integrated, cells = WhichCells(Integrated, idents = c("cMo"))) <- 21
Idents(Integrated, cells = WhichCells(Integrated, idents = c("nMo"))) <- 22
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Other Mo"))) <- 23
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Other Mac"))) <- 24
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Mo-Mac"))) <- 25
Idents(Integrated, cells = WhichCells(Integrated, idents = c("PVM"))) <- 26
Idents(Integrated, cells = WhichCells(Integrated, idents = c("LAM"))) <- 27
Idents(Integrated, cells = WhichCells(Integrated, idents = c("cDC2B"))) <- 28
Idents(Integrated, cells = WhichCells(Integrated, idents = c("cDC1"))) <- 29
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CCR7+ DC"))) <- 30
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling Myeloid"))) <- 31

Idents(Integrated, cells = WhichCells(Integrated, idents = c("Naive B cell", "Memory B cell"))) <- 32
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Plasmablast"))) <- 33
Idents(Integrated, cells = WhichCells(Integrated, idents = c("pDC"))) <- 34

Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 Naive"))) <- 35
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 TCM", "CD4 TEM"))) <- 36
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 Regulatory"))) <- 37
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 TEMRA & Senescent"))) <- 38

Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD8 Naive"))) <- 39
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD8 TCM", "CD8 TEM"))) <- 40
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD8 TEMRA & Senescent"))) <- 41
Idents(Integrated, cells = WhichCells(Integrated, idents = c("MAIT"))) <- 42
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Gamma Delta"))) <- 43

Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD57+ mNK"))) <- 44
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD16+ mNK"))) <- 45
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Immature NK"))) <- 46
Idents(Integrated, cells = WhichCells(Integrated, idents = c("ILC"))) <- 47
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling T & NK"))) <- 48

Integrated[['Celltype_plot']] <- Idents(Integrated)

UMAP_FUN(Integrated, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1B_Collapsed.png', color = col_shuffle, breaks = c(1:48), height = 8, width = 18,
label = c("Preadipocytes (1)", "ECM-Producing Preadipocyte (2)", "Mature Preadipocyte (3)", "ISG+ Preadipocyte (4)", "Metallothionein+ Preadipocyte/ASC (5)","MYOC+ Fibroblasts (6)", "PCOLCE+ Fibroblast (7)", "Adipose Stem Cell (8)",
"Myofibroblast (9)", "DKK3+ Fibroblast (10)", "Proliferating fibroblast (11)", "Venous EC (12)", "Capillary EC (13)", "Intermediate Capillary EC (14)", "Arterial EC (15)", "Arterial-EndoMT-like (16)", "Capillary-EndoMT-like (17)",
"Venous-EndoMT-like (18)", "Pericyte (19)", "VSMC (20)", "cMo (21)", "nMo (22)", "Other Mo (23)", "Other Mac (24)", "Mo-Mac (25)", "PVM (26)", "LAM (27)", "cDC2B (28)", "cDC1 (29)", "CCR7+ DC (30)", "Cycling Myeloid (31)",
"B Cell (32)", "Plasmablast (33)", "pDC (34)", "CD4 Naive (35)", "CD4 Memory (36)", "CD4 Regulatory (37)", "CD4 TEMRA & Senescent (38)", "CD8 Naive (39)", "CD8 Memory (40)", "CD8 TEMRA & Senescent (41)",
"MAIT (42)", "Gamma Delta (43)", "CD57+ mNK (44)", "CD16+ mNK (45)", "Immature NK (46)", "ILC (47)", "Cycling T & NK (48)"))

# With labels
UMAP_FUN(Integrated, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1B_Collapsed_labeled.png', color = col_shuffle, plot_label = T, breaks = c(1:48), height = 16, width = 20,
label = c("Preadipocytes (1)", "ECM-Producing Preadipocyte (2)", "Mature Preadipocyte (3)", "ISG+ Preadipocyte (4)", "Metallothionein+ Preadipocyte/ASC (5)","MYOC+ Fibroblasts (6)", "PCOLCE+ Fibroblast (7)", "Adipose Stem Cell (8)",
"Myofibroblast (9)", "DKK3+ Fibroblast (10)", "Proliferating fibroblast (11)", "Venous EC (12)", "Capillary EC (13)", "Intermediate Capillary EC (14)", "Arterial EC (15)", "Arterial-EndoMT-like (16)", "Capillary-EndoMT-like (17)",
"Venous-EndoMT-like (18)", "Pericyte (19)", "VSMC (20)", "cMo (21)", "nMo (22)", "Other Mo (23)", "Other Mac (24)", "Mo-Mac (25)", "PVM (26)", "LAM (27)", "cDC2B (28)", "cDC1 (29)", "CCR7+ DC (30)", "Cycling Myeloid (31)",
"B Cell (32)", "Plasmablast (33)", "pDC (34)", "CD4 Naive (35)", "CD4 Memory (36)", "CD4 Regulatory (37)", "CD4 TEMRA & Senescent (38)", "CD8 Naive (39)", "CD8 Memory (40)", "CD8 TEMRA & Senescent (41)",
"MAIT (42)", "Gamma Delta (43)", "CD57+ mNK (44)", "CD16+ mNK (45)", "Immature NK (46)", "ILC (47)", "Cycling T & NK (48)"))

#########################
# Figure 1C: Integrated UMAP By StudyGroup
#########################
DefaultAssay(Integrated) <- "RNA"
Integrated[['StudyGroup']] <- factor(Integrated$StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic", "HIV- diabetic"))
png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1C.png', width = 7, height = 5, units = 'in', res = 600)
DimPlot(Integrated, label = F, cols = c("#F04832", "#F57636", "#EBAA3B", "#4287F5"), group.by = 'StudyGroup', raster = F) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) + NoLegend()
dev.off()


#########################
# Figure 1D: Integrated UMAP with Gene Expression
#########################
png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1D_CLDN5.png", width = 5, height = 5, units = "in", res = 600)
FeaturePlot(Integrated, features = "CLDN5", raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1D_ACTA2.png", width = 5, height = 5, units = "in", res = 600)
FeaturePlot(Integrated, features = "ACTA2", raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1D_CD68.png", width = 5, height = 5, units = "in", res = 600)
FeaturePlot(Integrated, features = "CD68", raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1D_CD1C.png", width = 5, height = 5, units = "in", res = 600)
FeaturePlot(Integrated, features = "CD1C", raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1D_COL1A1.png", width = 5, height = 5, units = "in", res = 600)
FeaturePlot(Integrated, features = "COL1A1", raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1D_CCDC80.png", width = 5, height = 5, units = "in", res = 600)
FeaturePlot(Integrated, features = "CCDC80", raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1D_CD3E.png", width = 5, height = 5, units = "in", res = 600)
FeaturePlot(Integrated, features = "CD3E", raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1D_NKG7.png", width = 5, height = 5, units = "in", res = 600)
FeaturePlot(Integrated, features = "NKG7", raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

#########################
# Figure 1E: Integrated Major Cell Proportions by Study Group
#########################
Plot_prop(Prop_Long, filter = "Integrated", sort = "FALSE", y.axis = "all", filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure1/Figure1E.png', width = 18, height = 6, Margins = c(0,0,0,4), plot_unit = 'cm')

#########################
# Other Statistics (Results Section 1)
#########################
# Total Number of Cells per Major Cell Type
table(Integrated$Major_CellType)
#Endothelium    Lymphoid     Myeloid     Stromal 
#      93453       33733       68631       87674 

# Median Proportion of Major Cell Types by Disease State
Prop %>% group_by(StudyGroup) %>% 
    summarise(Lymphoid = median(Lymphoid),
    Vascular = median(Endothelium),
    Stromal = median(Stromal),
    Myeloid = median(Myeloid))

#  StudyGroup        Lymphoid Vascular Stromal Myeloid
#  <chr>                <dbl>    <dbl>   <dbl>   <dbl>
#1 HIV- diabetic         4.74    30.2     38.7    23.9
#2 HIV+ diabetic        15.1     24.4     29.3    23.7
#3 HIV+ non-diabetic    24.7      7.71    21.5    34.1
#4 HIV+ prediabetic      8.82    34.1     22.6    19.3

# Median Proportion of Major Cell Types by Glucose tolerance
Prop %>% filter(HIV == "HIVpos") %>%
group_by(Glucose) %>%
summarise(Lymphoid = median(Lymphoid),
Myeloid = median(Myeloid))

#  Glucose            Lymphoid Myeloid
#1 Glucose Intolerant     9.93    23.0
#2 Glucose Tolerant      24.7     34.1

# Wilcoxon Sum Testing
wilcox.test(Lymphoid~StudyGroup, data = Prop, subset = StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic")) # p = 0.002
wilcox.test(Myeloid~StudyGroup, data = Prop, subset = StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic")) # p = 0.52
wilcox.test(Myeloid~Glucose, data = Prop, subset = HIV %in% c("HIVpos")) # p = 0.03
wilcox.test(Lymphoid~Glucose, data = Prop, subset = HIV %in% c("HIVpos")) # p = 0.02

# Kruskal-Walis
kruskal.test(Lymphoid~StudyGroup, data = Prop, subset = StudyGroup %in% c("HIV+ diabetic", "HIV+ prediabetic", "HIV+ diabetic"))
kruskal.test(Myeloid~StudyGroup, data = Prop, subset = StudyGroup %in% c("HIV+ diabetic", "HIV+ prediabetic", "HIV+ diabetic"))


