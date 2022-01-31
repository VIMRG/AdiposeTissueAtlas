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
# Figure S2C: Expression Plots Macrophages
#########################
Idents(Myeloid) <- "Manual_Annotation"
Macrophage <- subset(Myeloid, idents = c("LAM", "PVM", "Mo-Mac", "Other Mac"))
VlnPlot(Macrophage, features = c("TREM2", "CD9", "RNASE1", "LYVE1", "HLA-DPB1", "CXCL3", "CXCL2", "TIMP1", "ISG15"), pt.size = 0, col = c("#F53722", "#40E66F", "#F9C99B", "#607C3B"), stack = TRUE, flip = TRUE, fill.by = 'ident') + 
NoLegend() + 
theme(axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold')) +
xlab("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2C.png", device = 'png', width = 10, height = 10, units = 'in', dpi = 600)

#########################
# Figure S2D: GO ORA
#########################
Macrophage_clusters <- Myeloid.markers %>% filter(cluster %in% c("LAM", "PVM", "Other Mac", "Mo-Mac")) # Pull out Macrophage clusters
Entrez <- clusterProfiler::bitr(Macrophage_clusters$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Convert to ENTREZID
Macrophage_Markers <- left_join(Entrez, Macrophage_clusters, by = c("SYMBOL" = "gene")) # Add ENTREZID
Macrophage_Markers <- Macrophage_Markers %>% dplyr::filter(p_val_adj < 0.05) # Filter to genes < 0.05

Macrophage_list <- Macrophage_Markers %>% group_by(cluster) %>% group_split() # Make List
names(Macrophage_list) <- c("LAM", "Mo-Mac", "Other Mac", "PVM")

# Loop through list of DGE for each Macrophage cluster
DF <- data.frame()
DF <- lapply(Macrophage_list, function(x) {
    bp <- enrichGO(gene = x$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
    #kk2 <- as.data.frame(DOSE::setReadable(kk2, org.Hs.eg.db, keyType = "ENTREZID"))
    results <- data.frame(bp2)
    #results$cluster <- names(x)
    DF <- rbind(DF, results)
} )

# Plot Results
Pathways <- bind_rows(DF, .id = "Cluster") # Combine in data frame

Pathways %>% group_by(Cluster) %>%
        do(head(.,5)) %>% 
        mutate(Description = factor(Description, levels = Description)) %>%
        ggplot() + geom_point(aes(x = Cluster, y = Description, size = -log10(pvalue), fill = Count), shape = 21) +
        xlab("") + ylab("") + ggtitle("") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1, size = 20, face = "bold"), 
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 20, face = "bold")) + scale_y_discrete(labels = function(x) str_wrap(x, width = 50))

ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2D.png", device = 'png', width = 10, height = 5, units = 'in', dpi = 600)

#########################
# Figure S2E: Intermediate Macrophage Scoring
#########################
IM_list <- list(c('CXCL3', 'CXCL2', 'IER3', 'IL1B', 'CCL20', 'CXCL8', 'EREG', 'SOD2', 'CCL3L1', 'BCL2A1', 'THBS1', 'OLR1', 'VCAN', 'C5AR1', 'PLAUR', 'IL10', 'IL1RN', 'AQP3',
'PTGS2', 'LYZ', 'CD163', 'NLRP3', 'NAMPT', 'ABL2', 'EGR1', 'ICAM1', 'FCN1', 'MARCKS', 'MMP19', 'PLIN2', 'CCL3', 'NFKBIA', 'TNFAIP2', 'CLEC4E', 'CTSB', 'FCGR2A', 'MAFB', 'IL1A',
'GK', 'MIR3945HG', 'ACSL1', 'PLEK', 'CYBB', 'SELENOP', 'MS4A6A', 'KYNU', 'TLR2', 'C1QA', 'TNF'))

Myeloid <- AddModuleScore(Myeloid, features = IM_list, name = 'IM')

png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2E.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Myeloid, features = "IM1", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
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
Merged <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Merged_Hildreth.rds') # Generated using IM_Anchor_Based_Reference.R script
Idents(Merged) <- "predicted.id"
Idents(Merged, cells = WhichCells(Merged, idents = c("cMo"))) <- 1
Idents(Merged, cells = WhichCells(Merged, idents = c("Other Mo"))) <- 2
Idents(Merged, cells = WhichCells(Merged, idents = c("nMo"))) <- 3
Idents(Merged, cells = WhichCells(Merged, idents = c("Mo-Mac"))) <- 4
Idents(Merged, cells = WhichCells(Merged, idents = c("PVM"))) <- 5
Idents(Merged, cells = WhichCells(Merged, idents = c("LAM"))) <- 6
Idents(Merged, cells = WhichCells(Merged, idents = c("Other Mac"))) <- 7
Idents(Merged, cells = WhichCells(Merged, idents = c("cDC2B"))) <- 8
Idents(Merged, cells = WhichCells(Merged, idents = c("CCR7+ DC"))) <- 9
Idents(Merged, cells = WhichCells(Merged, idents = c("cDC1"))) <- 10
Idents(Merged, cells = WhichCells(Merged, idents = c("Cycling Myeloid"))) <- 11

Merged[['Celltype_plot']] <- Idents(Merged)

UMAP_FUN(Merged, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2F.png', color = col_shuffle, breaks = c(1:11),  
label = c("cMo (1)","Other Mo (2)", "nMo (3)", "Mo-Mac (4)", "PVM (5)", "LAM (6)", "Other Macrophage (7)", "cDC2B (8)", "CCR7+ DC (9)", "cDC1 (10)", "Cycling Myeloid (11)"))


#########################
# Figure S2G: Intermediate Macrophage Scoring (Hildreth dataset)
#########################
# See IM_Anchor_Based_Reference.R script for preprocessing code

png('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2G.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Merged, features = "IM1", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

#########################
# Figure S2G: LAM score
#########################
LAMgenes <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/LAM_Genes_1.5FC.csv")
LAMgenesList <- list(c(as.character(LAMgenes$gene)))
Myeloid <- AddModuleScore(Myeloid, features = LAMgenesList, name = "LAM")
Myeloid$LAM1 <- scale(Myeloid$LAM1)

FeaturePlot(Stromal, features = "LAM1", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/LAM_Module.png", dpi = 600, height = 5, width = 5, units = 'in', device = 'png')
