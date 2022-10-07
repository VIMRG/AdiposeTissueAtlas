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
Idents(Integrated, cells = WhichCells(Integrated, idents = c("mreg DC"))) <- 31
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
Idents(Integrated, cells = WhichCells(Integrated, idents = c("Cyclcing T/NK"))) <- 51
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


