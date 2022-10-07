###############################################################################
######------------------SUPPLEMENTAL FIGURE 1-----------------------------#####
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
# Supplemental Figure 1A: Flow Chart
###-------------------####
# Preprocessing
graph <- grViz("
digraph boxes_and_circles {
graph [overlap = true, fontsize = 10, fontname = Arial]
node [shape = box, fontsize = 10, fontname = Arial]
'Cellranger output'; 'Genotypic clusters'; 'Seurat object'; 'Doublet-tagged object'
node [shape = circle, fixedsize = false, width = 0.9]
'Souporcell'; 'SoupX'; 'Doublet & QC';
'Cellranger output'->'SoupX' 'SoupX'->'Seurat object' 'Cellranger output'->'Souporcell' 'Souporcell'->'Genotypic clusters'
'Genotypic clusters'->'Seurat object' 'Seurat object'->'Doublet & QC'
'Doublet & QC'->'Doublet-tagged object'
}")

tmp <- DiagrammeRsvg::export_svg(graph)
tmp <- charToRaw(tmp)
rsvg::rsvg_png(tmp, paste(fig_dir, "SF1A_preprocessing.png", sep = "/"))

# Downstream
graph1 <- grViz("
digraph boxes_and_circles {
graph [overlap = true, fontsize = 10, fontnames = Arial]
node [shape = box, fontsize = 10, fontnames = Arial]
'Merged object'; 'Singlet object'; 'Subset analysis'; 'Cleaned object';
node [shape = circle, fixedsize = false, width = 0.9]
'Harmony integration'; 'Doublet removal'; 'QC & annotation'
'Merged object'->'Harmony integration' 'Harmony integration'->'Doublet removal' 'Doublet removal'->'Singlet object'
'Singlet object'->'Subset analysis' 'Subset analysis'->'QC & annotation' 'QC & annotation'->'Cleaned object'
}")

tmp <- DiagrammeRsvg::export_svg(graph1)
tmp <- charToRaw(tmp)
rsvg::rsvg_png(tmp, paste(fig_dir, "SF1A_downstream.png", sep = "/"))

###-------------------####
# Supplemental Figure 1B: Doublet Plots
###-------------------####
DoubletGraph<- read.csv("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/MetaData/DoubletGraph_HIV.csv")
colnames(DoubletGraph) <- c("Project", "Cellranger Total Barcodes", "QC and Doublet Removal")

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
        panel.grid.major.y = element_line(color = "grey90", linetype = 2),
        axis.text.x = element_text(size = 14, face = "bold")
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) + 
  scale_fill_hue(l=45, name = 'Filter') +
  ylab('Number of Cells') + xlab("")


# Save Graph
ggsave(filename = paste(fig_dir, 'SF1B.png', sep = "/"), dpi = 600, device = 'png', width = 15, height = 8, units = 'in')

###-------------------####
# Supplemental Figure 1C: nCount/Cell Plot
###-------------------####
Idents(Integrated) <- "Lane"
VlnPlot(Integrated, features = 'nCount_RNA', pt.size = 0, col = col_shuffle) + NoLegend() + theme_bw() +
theme(axis.text.x = element_text(size = 16, face = 'bold', angle = 45, vjust = 0.95, hjust = 1), axis.title.y = element_text(size = 18, face = 'bold')) +
xlab("") + ylab("Reads Per Cell") + ylim(0, 25000) + ggtitle("") + NoLegend()
ggsave(paste(fig_dir, 'SF1C.png', sep = "/"), dpi = 600, width = 10, height = 5, units = 'in', device = 'png')

###-------------------####
# Supplemental Figure 1D: nFeature/Cell Plot
###-------------------####
Idents(Integrated) <- "Lane"
VlnPlot(Integrated, features = 'nFeature_RNA', pt.size = 0, col = col_shuffle) + NoLegend() +
theme(axis.text.x = element_text(size = 16, face = 'bold', angle = 45, vjust = 0.95, hjust = 1), axis.title.y = element_text(size = 18, face = 'bold')) +
xlab("") + ylab("Reads Per Cell") + ggtitle("") + NoLegend()
ggsave(paste(fig_dir, 'SF1D.png', sep = "/"), dpi = 600, width = 10, height = 5, units = 'in', device = 'png')

###-------------------####
# Supplemental Figure 1E: Integration metrics
###-------------------####
# Will subsample Atlas to 50K cells (from 178K) to conserve memory. Then feed this as input to the SCIB pipeline
Idents(Integrated) <- "orig.ident"
SubInt <- Integrated[, sample(colnames(Integrated), size = 50000, replace = F)]

# Will recluster to ensure this captures variability in dataset before feeding into the SCIB pipeline
SubInt <- SubInt %>%
        FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Downstream Analysis
SubInt <- SubInt %>% RunUMAP(dims = 1:30) %>%
        FindNeighbors(dims = 1:30) %>%
        FindClusters(resolution = 2.0)

p1 <- DimPlot(PVM, label = T) + NoLegend()
p2 <- DimPlot(PVM, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

saveRDS(SubInt, file = paste(dir, "DownSampled50K.rds", sep = "/"))

# Convert seurat rds to h5ad format for SCIB. This will create an anndata object with the normalized and log transformed gene counts (X), raw counts (layers) along with metadata (batch) needed for SCIB
writeAnnData(SubInt, filename = "/data/p_koethe_lab/Integration_Pipeline/scib-pipeline/data/DownSampled.h5ad")

###-------------------####
# Supplemental Figure 1F: Dot Plot
###-------------------####
Integrated[['Annotation']] <- factor(Integrated$Annotation, levels = c("PCOLCE+ Fibroblast", "MYOC+ Fibroblast", "Adipose Progenitor Cell 2", "Adipose Progenitor Cell 1",
"Early Preadipocyte", "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", 
"ISG+ Preadipocyte", "Metallothionein+ Preadipocyte", "Myofibroblast", "Cycling Myofibroblast",
"VSMC", "Pericyte", "Arterial EndoMT-like", "Venous EndoMT-like", "Capillary EndoMT-like", "Arterial EC", "Intermediate Capillary EC", "Capillary EC", "Venous EC", "Ven-Cap EC", "Cycling Vascular Cells", "LAM", "IM", "PVM", "Mo-Mac", "Other Mac", "cMo",
"nMo", "Other Mo", "ISG+ Mo", "cDC2B", "cDC1", "Migratory DC", "pDC", "Cycling Myeloid", "B Cell", "Plasmablast", "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Regulatory", "CD4 Cytotoxic", "CD8 Naive",
"CD8 TCM", "CD8 TEM", "CD8 Cytotoxic", "Gamma Delta", 'MAIT', "CD57 mNK", "CD16 mNK", "Immature NK", "ILC", "Cycling T & NK"))

DOTPLOT_FUN(Integrated, ident = 'Annotation', features = c("PCOLCE2", "CD55", "MYOC", "PI16", "CLU", "DCN", "LUM", "ZFP36", "MYC", "CEBPD", "CXCL14", "EGR1", "CD9", "THY1", "COL6A1", "LPL", "CIDEC", "ISG15", "IFI6",
"MT1X", "MT2A", "POSTN", "TIMP1", "APOE", "MKI67", "ACTA2", "TAGLN","COX4I2", "RGS5", "SNAI1", "ETS2", "ZEB2", "GJA4", "ACKR1", "CA4", "HEY1", "VWF", "CD68", "TREM2", "CXCL2", "CXCL3", "RNASE1", "LYVE1", "S100A12", "VCAN", 
"FCGR3A", "CD1C", "CLEC9A", "LAMP3", "LILRA4", "MS4A1", "MZB1", "JCHAIN", "IL7", "SELL", "COTL1", "DUSP2", "GPR183", "ANXA5", "ALOX5AP", "CCL5", "NKG7", "PRF1", "CD8A", "TRDV1", "KLRG1", "XCL1", "KIT"), 
path = paste(fig_dir, 'SF1F.png', sep = "/"), width = 25, height = 15)

