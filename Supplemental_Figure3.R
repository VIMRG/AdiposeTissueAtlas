###############################################################################
######--------------------SUPPLEMENTAL FIGURE 3---------------------------#####
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

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Figure_Functions.R')

###-------------------####
# Import Data
###-------------------####
CD4 <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/CD4_RNA_11.1.rds')
CD8 <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/CD8_RNA_11.1.rds')
Tcells <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_10.30.rds')

Prop_Long <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Proportions_Long.csv') # Proportions Cell Types Long
Prop <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Merged_Proportions.csv')  # Proportions Cell Types Wide

Tcell.markers <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcell_cluster_RNAmarkers.csv')

###-------------------####
# Environment
###-------------------####
set.seed(7612)

###-------------------####
# Set ggplot theme
###-------------------####
theme_set(theme_bw() +
            theme(axis.text = (element_text(size = 6)),
                  plot.caption = element_text(face = "italic", size = 9),
                  legend.position = "bottom",
                  legend.background = element_rect(fill = "transparent"), 
                  legend.text = element_text(size = 18),
                  panel.background = element_rect(fill = "transparent"),
                  panel.border = element_blank(),
                  axis.line = element_line(color="black"),
                  strip.text = element_text(size = 24), 
                  axis.title = element_text(size = 18),
                  plot.margin=unit(c(1,1,0,0),"cm"),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)
          ))
          
# 2. SUPPLEMENTAL FIGURE 3 (Tcells) ___________________________________________________________________
#########################
# Figure S3A: T Cell Heatmap
#########################
HEATMAP_FUN(Tcells, markers = Tcell.markers, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3A.png', height = 12, Margins = c(0, 1, 0, 0), plot_unit = "cm")

# 3. SUPPLEMENTAL FIGURE 3 (CD4 T Cells) ___________________________________________________________________
#########################
# Figure S3B: CD4 T cell UMAP
#########################
Idents(CD4) <- "Fine_Annotation"
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Naive"))) <- 1
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 TCM"))) <- 2
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Transitional"))) <- 3
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 TEM"))) <- 4
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 TEMRA & Senescent"))) <- 5
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Regulatory"))) <- 6

CD4[['Celltype_plot']] <- Idents(CD4)

UMAP_FUN(CD4, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3B.png', color = col_shuffle[24:29], breaks = c(1:6),
label = c("CD4 Naive (1)", "CD4 TCM (2)", "CD4 Transitional (3)", "CD4 TEM (4)", "CD4 TEMRA & Senescent (5)", "CD4 Regulatory (6)"))

#########################
# Figure S3C: CD4 T cell DotPlot
#########################
CD4[['Fine_Annotation']] <- factor(CD4$Fine_Annotation, levels = c("CD4 Naive", "CD4 TCM", "CD4 Transitional", "CD4 TEM", "CD4 TEMRA & Senescent", "CD4 Regulatory"))

DOTPLOT_FUN(CD4, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3C.png', features = c("SELL", "LEF1", "CCR7", "LTB", "LDHB", "GPR183",
"COTL1", "DUSP2", "ALOX5AP", "IL32", "CCL5", "GNLY", "PRF1", "FOXP3"), ident = "Fine_Annotation", width = 10, height = 5)

#########################
# Figure S3D: CD4 T cell Antibody
#########################
DefaultAssay(CD4) <- "ADT"
Idents(CD4) <- "Fine_Annotation"
png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3D.png', width = 10, height = 5, units = 'in', res = 600)
VlnPlot(CD4, features = c("CD45RA", "CD27", "CD57"), pt.size = 0, col = col_shuffle[24:29], stack = TRUE, flip = TRUE, fill.by = 'ident') + 
NoLegend() + 
theme(axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold')) +
xlab("")
dev.off()

#########################
# Figure S3E: CD4 T cell Proportion
#########################
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
          
Plot_prop(Prop_Long, filter = "CD4", sort = "FALSE", y.axis = "CD4", filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3E.png')

# 3. SUPPLEMENTAL FIGURE 3 (CD8 T Cells) ___________________________________________________________________
#########################
# Figure S3F: CD8 T cell UMAP
#########################
Idents(CD8) <- "Fine_Annotation"
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 Naive"))) <- 1
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 TCM"))) <- 2
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 TEM"))) <- 3
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 TEMRA & Senescent"))) <- 4

CD8[['Celltype_plot']] <- Idents(CD8)

UMAP_FUN(CD8, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3F.png', color = col_shuffle[30:33], breaks = c(1:4),
label = c("CD8 Naive (1)", "CD8 TCM (2)", "CD8 TEM (3)", "CD8 TEMRA & Senescent (4)"), height = 5, width = 10)

#########################
# Figure S3G: CD8 T cell DotPlot
#########################
CD8[['Fine_Annotation']] <- factor(CD8$Fine_Annotation, levels = c("CD8 Naive", "CD8 TCM", "CD8 Transitional", "CD8 TEM", "CD8 TEMRA & Senescent"))

DOTPLOT_FUN(CD8, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3G.png', features = c("SELL", "NELL2", "LEF1", "CCR7", "LTB", "LDHB", "GPR183",
"COTL1", "DUSP2", "ALOX5AP", "IL32", "CCL5", "GNLY", "PRF1"), ident = "Fine_Annotation", width = 10, height = 5)

#########################
# Figure S3H: CD8 T cell Antibody (Violin Plot)
#########################
DefaultAssay(CD8) <- "ADT"
Idents(CD8) <- "Fine_Annotation"
png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3H.png', width = 10, height = 5, units = 'in', res = 600)
VlnPlot(CD8, features = c("CD45RA", "CD27", "CD57"), pt.size = 0, col = col_shuffle[30:33], stack = TRUE, flip = TRUE, fill.by = 'ident') + 
NoLegend() + 
theme(axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold')) +
xlab("")
dev.off()

#########################
# Figure S3I: CD8 T cell Proportion
#########################
Plot_prop(Prop_Long, filter = "CD8", sort = "FALSE", y.axis = "CD8", filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure3/SF3I.png', Margins = c(0,0,0,0.6), plot_unit = "in")




