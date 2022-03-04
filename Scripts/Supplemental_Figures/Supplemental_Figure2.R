###############################################################################
######--------------------SUPPLEMENTAL FIGURE 2---------------------------#####
###############################################################################

# 1. SET UP___________________________________________________________________
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
library(clusterProfiler)
library(rWikiPathways)
library(org.Hs.eg.db)
library(DOSE)

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Figure_Functions.R')

###-------------------####
# Import Data
###-------------------####
Myeloid <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_RNA_11.1.rds') # Seurat Object
Prop_Long <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Proportions_Long.csv') # Proportions Cell Types Long
Prop <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Merged_Proportions.csv')  # Proportions Cell Types Wide
Myeloid.markers <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_cluster_RNAmarkers.csv') # Myeloid Markers
Hildreth <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Merged_Hildreth.rds')

###-------------------####
# Environment
###-------------------####
set.seed(7612)

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
          
# 2. SUPPLEMENTAL FIGURE 2___________________________________________________________________
#########################
# Figure S2A: Myeloid UMAP by Groups
#########################
DefaultAssay(Myeloid) <- "RNA"
Myeloid[['StudyGroup']] <- factor(Myeloid$StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic", "HIV- diabetic"))
png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2A.png', width = 7, height = 5, units = 'in', res = 600)
DimPlot(Myeloid, label = F, cols = c("#F04832", "#F57636", "#EBAA3B", "#4287F5"), group.by = 'StudyGroup') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

#########################
# Figure S2B: Myeloid Marker Heatmap
#########################
HEATMAP_FUN(Myeloid, markers = Myeloid.markers, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2B.png')

#########################
# Figure S2C: Myeloid CITE-Seq Antibody
#########################
DefaultAssay(Myeloid) <- "ADT"
png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2C_CD14.png', res = 600, height = 5, width = 5, units = 'in')
FeaturePlot(Myeloid, features = c("CD14"), min.cutoff = 'q20', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2C_CD16.png', res = 600, height = 5, width = 5, units = 'in')
FeaturePlot(Myeloid, features = c("CD16"), min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2C_CD11C.png', res = 600, height = 5, width = 5, units = 'in')
FeaturePlot(Myeloid, features = c("CD11C"), min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2C_CD1C.png', res = 600, height = 5, width = 5, units = 'in')
FeaturePlot(Myeloid, features = c("CD1C"), min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

#########################
# Figure S2D: Correlation Between LAM & cMo
#########################
Prop %>% ggplot(aes(x = LAM, y = cMo)) + 
        geom_point(size = 2) +
        xlab("Lipid-Associated Macrophages (% Total Myeloid)") +
        ylab("Classical Monocytes (% Total Myeloid)") + 
        theme(axis.title.x = element_text(size = 16, face = "bold", margin(t = 20)),
        axis.title.y = element_text(size = 16, face = "bold", margin(r = 20)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) + 
        ylim(0,70) +
        geom_smooth(color = "black", method=lm, se=FALSE)

ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2D.png', dpi = 600, units = 'in', width = 8, height = 8, device = 'png')

#########################
# Figure S2E: Intermediate Macrophage Scoring
#########################
DefaultAssay(Myeloid) <- "RNA"
IM_list <- list(c('CXCL3', 'CXCL2', 'IER3', 'IL1B', 'CCL20', 'CXCL8', 'EREG', 'SOD2', 'CCL3L1', 'BCL2A1', 'THBS1', 'OLR1', 'VCAN', 'C5AR1', 'PLAUR', 'IL10', 'IL1RN', 'AQP3',
'PTGS2', 'LYZ', 'CD163', 'NLRP3', 'NAMPT', 'ABL2', 'EGR1', 'ICAM1', 'FCN1', 'MARCKS', 'MMP19', 'PLIN2', 'CCL3', 'NFKBIA', 'TNFAIP2', 'CLEC4E', 'CTSB', 'FCGR2A', 'MAFB', 'IL1A',
'GK', 'MIR3945HG', 'ACSL1', 'PLEK', 'CYBB', 'SELENOP', 'MS4A6A', 'KYNU', 'TLR2', 'C1QA', 'TNF'))

Myeloid <- AddModuleScore(Myeloid, features = IM_list, name = 'IM')

png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2E_1.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Myeloid, features = "IM1", min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

Hildreth <- AddModuleScore(Hildreth, features = IM_list, name = 'IM')

png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2E_2.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Hildreth, features = "IM1", min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "#f59d0f", "#f52d05")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

#########################
# Figure S2F: Anchor-based Reference Mapping
#########################
# See IM_Anchor_Based_Reference.R script for processing code
Idents(Hildreth) <- "predicted.id"
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("cMo"))) <- 1
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("Other Mo"))) <- 2
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("nMo"))) <- 3
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("Mo-Mac"))) <- 4
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("PVM"))) <- 5
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("LAM"))) <- 6
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("Other Mac"))) <- 7
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("cDC2B"))) <- 8
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("CCR7+ DC"))) <- 9
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("cDC1"))) <- 10
Idents(Hildreth, cells = WhichCells(Hildreth, idents = c("Cycling Myeloid"))) <- 11

Hildreth[['Celltype_plot']] <- Idents(Hildreth)

UMAP_FUN(Hildreth, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2F.png', color = col_shuffle, breaks = c(1:11),  
label = c("cMo (1)","Other Mo (2)", "nMo (3)", "Mo-Mac (4)", "PVM (5)", "LAM (6)", "Other Macrophage (7)", "cDC2B (8)", "CCR7+ DC (9)", "cDC1 (10)", "Cycling Myeloid (11)"))
