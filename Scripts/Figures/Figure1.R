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

###-------------------####
# Source Code
###-------------------####
source('../Utils.R')

###-------------------####
# Load Data
###-------------------####
dir <- "../Final_Integration"
Integrated <- readRDS(paste(dir, "Integrated_3.20.rds", sep = "/"))

fig_dir <- "../Figures"
Prop_dir <- "../Proportion"

###-------------------####
# Figure 1B: UMAP Overall
###-------------------####
Idents(Integrated) <- "Annotation"

# Stromal
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Early Preadipocyte", "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2",
"Metallothionein+ Preadipocyte", "Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Preadipocyte"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("PCOLCE+ Fibroblast"))) <- "Adipose Progenitor"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Myofibroblast"))) <- "Myofibroblast"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("MYOC+ Fibroblast"))) <- "MYOC+ Fibroblast"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling Myofibroblast"))) <- "Cycling"

# Vascular
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Venous EC", "Capillary EC", "Arterial EC"))) <- "Endothelial" 
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Pericyte"))) <- "Pericyte"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("VSMC 1", "VSMC 2"))) <- "Smooth Muscle"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling Vascular"))) <- "Cycling"

# Monocyte
Idents(Integrated, cells = WhichCells(Integrated, idents = c("cMo", "nMo", "Other Mo", "ISG+ Mo"))) <- "Monocyte"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Mo-Mac 2", "Mo-Mac 1", "PVM", "LAM", "IM"))) <- "Macrophage"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("cDC2B", "cDC1", "Migratory DC", "pDC"))) <- "Dendritic"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cycling Myeloid"))) <- "Cycling"

# B cells
Idents(Integrated, cells = WhichCells(Integrated, idents = c("B Cell", "Plasmablast"))) <- "B cell"

# T Cells
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD4 Naive", "CD4 TEM", "CD4 TCM", "CD8 TCM", "CD4 Regulatory", "CD4 Cytotoxic", "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Cytotoxic", 
"MAIT", "Gamma Delta"))) <- "T cell"

# NK
Idents(Integrated, cells = WhichCells(Integrated, idents = c("CD57 mNK", "CD16 mNK", "Immature NK"))) <- "NK cell"
Idents(Integrated, cells = WhichCells(Integrated, idents = c("ILC"))) <- "Innate Lymphoid"
Integrated[['Celltype_plot']] <- Idents(Integrated)

plot <- DimPlot(Integrated, label = F, repel = F, raster = F) +
    scale_colour_manual(breaks = c("Preadipocyte", "Adipose Progenitor", "Myofibroblast", "Antiadipogenic", "Endothelial", "Pericyte", "Smooth Muscle", "Monocyte", "Macrophage", "Dendritic", "B cell",
    "T cell", "NK cell", "Innate Lymphoid", "Cycling"), values = c("gold1", "chocolate4", "slategrey", "peru", "firebrick", "olivedrab", "lightpink2", "darkolivegreen", "steelblue", "tomato",
    "dodgerblue", "darkseagreen4", "darkorchid", "burlywood2", "coral4", "darkorange3")) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.position = "none")
ggsave(filename = paste(fig_dir, "Figure1B.png", sep = "/"), dpi = 600, height = 5, width = 5, units = 'in', device = 'png')  

###-------------------####
# Figure 1C: Harmony Integration
###-------------------####
DefaultAssay(Integrated) <- "RNA"
Integrated[['StudyGroup']] <- factor(Integrated$StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic"))
png(paste(fig_dir, "Figure1C.png", sep = "/"), width = 5, height = 5, units = 'in', res = 600)
DimPlot(Integrated, label = F, cols = c("#91C46C", "#7ABAF2", "#FFBE00"), group.by = 'StudyGroup', raster = F) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          plot.title = element_text(face = "bold", size = 16)) + NoLegend()
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
Major.prop <- read.table(paste(Prop_dir, "All.Prop.txt", sep = "/"), check.names = FALSE)
Major.prop <- Major.prop %>% na.omit()
Major.prop <- Prop_Merge(Major.prop, data_hatim)
Prop_long <- Major.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

use_colors <- c("#91C46C", "#7ABAF2", "#FFBE00")

plot <- box_plot_fun_paper(Prop_long, plot_cols = use_colors, parent_cell_type = "All Cell Types", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 4, legend.pos = "bottom")
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
DM <- Hmisc::summaryM(Stromal + Vascular + Myeloid + Lymphoid ~ StudyGroup, data = Major.prop[Major.prop$StudyGroup %in% c("HIV+ diabetic", "HIV+ non-diabetic"),], overall = T, test = T)
preDM <- Hmisc::summaryM(Stromal + Vascular + Myeloid + Lymphoid ~ StudyGroup, data = Major.prop[Major.prop$StudyGroup %in% c("HIV+ prediabetic", "HIV+ non-diabetic"),], overall = T, test = T)
Results <- FDR_Fun(df1 = DM, df2 = preDM, clusters = c("Stromal", "Vascular", "Myeloid", "Lymphoid"))
write.csv(Results, file = paste(Prop_dir, "MajorCellType.Proportion.Adjusted.csv", sep = "/"))

