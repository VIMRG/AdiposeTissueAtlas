###############################################################################
######---------------------SUPPLEMENTAL FIGURE 1----------------------------#####
###############################################################################

###-------------------####
# Import Libraries
###-------------------####
library(ggplot2)
library(Seurat)
library(dplyr)
library(patchwork)
library(DiagrammeR)
library(rsvg)
library(DiagrammeRsvg)
library(tidyverse)

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Figure_Functions.R')

###-------------------####
# Load Datasets
###-------------------####
Integrated <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Final_Integrated.1.3.rds')
Tcells <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_10.30.rds')

###-------------------####
# SF1A: Algorithm
###-------------------####
# Homotypic Doublet Removal
graph <- grViz("
digraph boxes_and_circles {
graph [overlap = true, fontsize = 10, fontname = Arial]
node [shape = box, fontsize = 10, fontname = Arial]
'10x lane'; 'Fastq files'; 'Cellranger output'; 'Genotypic clusters'; 'Seurat object'
node [shape = circle, fixedsize = false, width = 0.9]
'Demultiplex'; 'Cellranger count'; 'Souporcell'; 'SoupX'; 'Read 10x data'
'10x lane'->'Demultiplex' 'Demultiplex'-> 'Fastq files' 'Fastq files'->'Cellranger count'
'Cellranger count'->'Cellranger output' 'Cellranger output'->'SoupX' 'SoupX'->'Read 10x data'
'Cellranger output'->'Souporcell' 'Souporcell'->'Genotypic clusters' 'Read 10x data'->'Seurat object'
'Genotypic clusters'->'Seurat object' 'Seurat object'->'Genotypic doublet removal' 
}")

tmp <- DiagrammeRsvg::export_svg(graph)
tmp <- charToRaw(tmp)
rsvg::rsvg_png(tmp, "/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/SF1A_HomotypicDoublet.png")

# Heterotypic Doublets
graph1 <- grViz("
digraph boxes_and_circles {
graph [overlap = true, fontsize = 10, fontnames = Arial]
node [shape = box, fontsize = 10, fontnames = Arial]
'Seurat object'; 'Filtered Seurat object'; 'Seurat clusters'; 'Heterotypic doublet identification'; 'Processed object'
node [shape = circle, fixedsize = false, width = 0.9]
'Quality control'; 'PCA/UMAP/Neighbor/Cluster'; 'Visual inspection'; 'DoubletFinder'
'Seurat object'->'Quality control' 'Quality control'->'Filtered Seurat object' 'Filtered Seurat object'->'PCA/UMAP/Neighbor/Cluster' 
'PCA/UMAP/Neighbor/Cluster'->'Seurat clusters' 'Seurat clusters'->'Visual inspection'
'Seurat clusters'->'DoubletFinder' 'Visual inspection'->'Heterotypic doublet identification'
'DoubletFinder'->'Heterotypic doublet identification'
'Heterotypic doublet identification'->'Processed object'
}")

tmp <- DiagrammeRsvg::export_svg(graph1)
tmp <- charToRaw(tmp)
rsvg::rsvg_png(tmp, "/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/SF1A_HeterotypicDoublet.png")

# Integration
graph2 <- grViz("
digraph boxes_and_circles {
graph [overlap = true, fontsize = 10, fontnames = Arial]
node [shape = box, fontsize = 10, fontnames = Arial]
'Singlet lane 1'; 'Singlet lane 2'; 'Singlet lane 3';
'Integrated Seurat object'; 'Subcluster on major cell types'; 'Suspicious for doublet'; 'Clean'
node [shape = circle, fixedsize = false, width = 0.9]
'Harmony'; 'PCA/UMAP/Neighbors/Cluster'; 'Recluster'; 'Remove cluster'; 'Integrate cleaned subclusters'
'Singlet lane 1'->'Harmony' 'Singlet lane 2'->'Harmony' 'Singlet lane 3'->'Harmony' 'Harmony'->'Integrated Seurat object' 
'Integrated Seurat object'->'PCA/UMAP/Neighbors/Cluster' 'PCA/UMAP/Neighbors/Cluster'->'Subcluster on major cell types' 
'Subcluster on major cell types'->'Recluster' 'Recluster'->'Suspicious for doublet' 'Recluster'->'Clean'
'Suspicious for doublet'->'Remove cluster' 'Remove cluster'->'Recluster' 'Clean'->'Integrate cleaned subclusters'}")

tmp <- DiagrammeRsvg::export_svg(graph2)
tmp <- charToRaw(tmp)
rsvg::rsvg_png(tmp, "/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/SF1A_Integration.png")


###-------------------####
# SF1B: Processing Counts
###-------------------####
DoubletGraph<- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/DoubletGraph.csv")
colnames(DoubletGraph) <- c("Project", "Cellranger Total Barcodes", "Homotypic Doublet Removal", "QC and Heterotypic Doublet Removal")

# Reshape
Long_graph <- reshape2::melt(DoubletGraph, id.vars = 'Project')

# Plot Graph
ggplot(Long_graph, aes(x = Project, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme(axis.text = (element_text(size = 16)),
                  plot.caption = element_text(face = "italic", size = 9),
                  legend.position = "right",
                  legend.background = element_rect(fill = "transparent"), 
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 24, face = "bold"),
                  panel.background = element_rect(fill = "transparent"),
                  panel.border = element_blank(),
                  axis.line = element_line(color="black"),
                  strip.text = element_text(size = 24), 
                  axis.title = element_text(size = 20, face = "bold"),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)
          ) +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) + 
  scale_fill_hue(l=45, name = 'Filter') +
  ylab('Number of Cells') + xlab("")

# Save Graph
ggsave(filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/SF1B.png', dpi = 600, device = 'png', width = 15, height = 8, units = 'in')


###-------------------####
# SF1C: Counts per Cell
###-------------------####
Idents(Integrated) <- "Project"
VlnPlot(Integrated, features = 'nCount_RNA', pt.size = 0, col = col_shuffle) + NoLegend() +
theme(axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold')) +
xlab("") + ylab("Reads Per Cell") + ggtitle("")
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/SF1C.png', dpi = 600, width = 10, height = 5, units = 'in', device = 'png')

###-------------------####
# SF1D: Features per Cell
###-------------------####
Idents(Integrated) <- "Project"
VlnPlot(Integrated, features = 'nFeature_RNA', pt.size = 0, col = col_shuffle) + NoLegend() +
theme(axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold')) +
xlab("") + ylab("Features Per Cell") + ggtitle("")
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/SF1D.png', dpi = 600, width = 10, height = 5, units = 'in', device = 'png')

###-------------------####
# SF1E: T cells as example of batch effect
###-------------------####
# Batch Uncorrected
Tcells_uncorrected <- Tcells %>% RunUMAP(dims = 1:30) %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 2.0)

Tcells[['StudyGroup']] <- factor(Tcells$StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic", "HIV- diabetic"))
Tcells_uncorrected[['StudyGroup']] <- factor(Tcells_uncorrected$StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic", "HIV- diabetic"))

png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/SF1F.png', width = 10, height = 5, units = 'in', res = 600)
p1 <- DimPlot(Tcells, label = F, cols = c("#F04832", "#F57636", "#EBAA3B", "#4287F5"), group.by = 'StudyGroup', raster = F) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) + ggtitle("")
p2 <- DimPlot(Tcells_uncorrected, label = F, cols = c("#F04832", "#F57636", "#EBAA3B", "#4287F5"), group.by = 'StudyGroup', raster = F) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) + ggtitle("")
p1 + p2
dev.off()

###-------------------####
# SF1F: DotPlot
###-------------------####
Integrated[['Manual_Annotation']] <- factor(Integrated$Manual_Annotation, levels = c("PCOLCE+ Fibroblast", "MYOC+ Fibroblast", "ASC 3", "ASC 2", "ASC 1",
"Early Preadipocyte 1", "Early Preadipocyte 2", "Early Preadipocyte 3", "Early Preadipocyte 4", "PTGDS+ Preadipocyte", "Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "CD9hi Preadipocyte", "PTGDS+ ECM-Producing Preadipocyte",
"ECM-Producing Preadipocyte", "ISG+ Preadipocyte", "Metallothionein+ Preadipocyte/ASC", "Myofibroblast", "Lipomyofibroblast", "Proliferating Myofibroblast", "DKK3+ Fibroblast",
"VSMC", "Pericyte", "Arterial-EndoMT-like", "Venous-EndoMT-like", "Capillary-EndoMT-like", "Arterial", "Intermediate Capillary", "Capillary", "Venous", "LAM", "PVM", "Mo-Mac", "Other Mac", "cMo",
"nMo", "Other Mo", "cDC2B", "cDC1", "CCR7+ DC", "pDC", "Cycling Myeloid", "Naive B cell", "Memory B cell", "Plasmablast", "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Regulatory", "CD4 TEMRA & Senescent", "CD8 Naive",
"CD8 TCM", "CD8 TEM", "CD8 TEMRA & Senescent", "Gamma Delta", 'MAIT', "CD57+ mNK", "CD16+ mNK", "Immature NK", "ILC", "Cycling T & NK"))

DOTPLOT_FUN(Integrated, ident = 'Manual_Annotation', features = c("PCOLCE2", "CD55", "MYOC", "PI16", "CLU", "DCN", "LUM", "ZFP36", "MYC", "CEBPD", "CXCL14", "EGR1", "PTGDS", "LPL", "CIDEC", "CD9", "THY1", "COL6A1", "ISG15", "IFI6",
"MT1X", "MT2A", "POSTN", "TIMP1", "APOE", "DKK3", "MKI67", "ACTA2", "TAGLN","COX4I2", "RGS5", "SNAI1", "ETS2", "ZEB2", "GJA4", "ACKR1", "CA4", "HEY1", "VWF", "CD68", "TREM2", "RNASE1", "LYVE1", "S100A12", "VCAN", 
"FCGR3A", "CD1C", "CLEC9A", "CCR7", "LILRA4", "MS4A1", "MZB1", "JCHAIN", "IL7", "SELL", "COTL1", "DUSP2", "GPR183", "ANXA5", "ALOX5AP", "CCL5", "NKG7", "PRF1", "CD8A", "TRDV1", "KLRG1", "XCL1", "KIT"), 
path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/SF1E.png', width = 25, height = 15)


