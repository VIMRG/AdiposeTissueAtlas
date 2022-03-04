###############################################################################
######----------------------------FIGURE 2--------------------------------#####
###############################################################################

###-------------------####
# Environment
###-------------------####
set.seed(7612)

###-------------------####
# Import Libraries
###-------------------####
library(ggplot2)
library(Seurat)
library(tidyverse)
library(rms)
library(MAST)
library(future)
library(ggrepel)
plan("multiprocess", workers = 5)

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Figure_Functions.R')
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/DGE/DGE_Function.R')

###-------------------####
# Load Data
###-------------------####
Myeloid <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_RNA_11.1.rds')
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
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)
          ))

#########################
# Figure 2A: Myeloid UMAP
#########################
Idents(Myeloid) <- "Manual_Annotation"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cMo"))) <- 1
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Other Mo"))) <- 2
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("nMo"))) <- 3
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Mo-Mac"))) <- 4
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("PVM"))) <- 5
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("LAM"))) <- 6
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Other Mac"))) <- 7
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cDC2B"))) <- 8
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("CCR7+ DC"))) <- 9
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("pDC"))) <- 10
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cDC1"))) <- 11
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("Cycling Myeloid"))) <- 12

Myeloid[['Celltype_plot']] <- Idents(Myeloid)

UMAP_FUN(Myeloid, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure2/Figure2A.png', color = col_shuffle, breaks = c(1:12),  
label = c("cMo (1)","Other Mo (2)", "nMo (3)", "Mo-Mac (4)", "PVM (5)", "LAM (6)", "Other Macrophage (7)", "cDC2B (8)", "CCR7+ DC (9)", "pDC (10)", "cDC1 (11)", "Cycling Myeloid (12)"))

#########################
# Figure 2B: Myeloid Dot Plot
#########################
Myeloid[['Manual_Annotation']] <- factor(Myeloid$Manual_Annotation, levels = c('cMo', 'nMo', 'Other Mo', 'Mo-Mac', 'Other Mac', 'PVM', 'LAM', 'cDC2B', 'CCR7+ DC', 'cDC1', 'pDC', 'Cycling Myeloid'))

DOTPLOT_FUN(Myeloid, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure2/Figure2B.png', features = c("S100A8", "VCAN", "CD14", "FCGR3A", "SERPINA1", "WARS", "APOBEC3A", "AREG",
'C1QB', 'CD68', "CD163", 'HLA-DRA', 'FCN1', "RNASE1", "LYVE1", 'TREM2', 'CD9', 'CD1C', 'CLEC10A', 'CCR7', 'CLEC9A', 'LILRA4', 'MKI67'), ident = 'Manual_Annotation')

#########################
# Figure 2C: Expression Plots Macrophages
#########################
Idents(Myeloid) <- "Manual_Annotation"
Macrophage <- subset(Myeloid, idents = c("LAM", "PVM", "Mo-Mac", "Other Mac"))
VlnPlot(Macrophage, features = c("TREM2", "CD9", "RNASE1", "LYVE1", "HLA-DPB1", "CXCL3", "CXCL2", "TIMP1", "ISG15"), pt.size = 0, col = c("#F53722", "#40E66F", "#F9C99B", "#607C3B"), stack = TRUE, flip = TRUE, fill.by = 'ident') + 
NoLegend() + 
theme(axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold')) +
xlab("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure2/SF2C.png", device = 'png', width = 10, height = 10, units = 'in', dpi = 600)

#########################
# Figure 2D: GO ORA
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
# Figure 2E: Myeloid Proportions
#########################
Plot_prop(Prop_Long, filter = "Myeloid", sort = "TRUE", y.axis = "Myeloid", filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure2/Figure2D.png')

#########################
# Figure 2E: Correlation Between LAM & cMo
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

ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure2/Figure2E.png', dpi = 600, units = 'in', width = 8, height = 8, device = 'png')

#########################
# Figure 2F: Macrophage DGE
#########################
Idents(Myeloid) <- "Manual_Annotation"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("LAM", "PVM", "Other Mac", "Mo-Mac"))) <- "Macrophage" # Aggregate Macrophage
Myeloid[['Manual_Annotation']] <- Idents(Myeloid)
Myeloid <- CellType_Fun(Myeloid) # From DGE functions

Myeloid.DGE.HIVDMvHIVnonDM <- HIV_comparison.5(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Macrophage.HIVDM.HIVnonDM.csv")
volcano.HIVDMvHIVnonDM.plot <- VOLCANO(Myeloid.DGE.HIVDMvHIVnonDM, filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/HIV.DM.NonDM.LAM.volcano.png', 
genes = c("HLA-DRB5", "CCL2", "MMP19", "FABP5", "APOE", "FABP4", "STAB1", "KLF2", "KLF4", "TREM2", "CEBPD","PLIN2", "CD9", "LPL", "CD36", "LIPA", "TGFBI", "LYVE1", "SPP1"), filter = "Macrophage")

Myeloid.DGE.HIVDMvHIVnegDM <- DM_comparison.5(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Macrophage.HIVDM.HIVnegDM.csv")
volcano.HIVDMvHIVnegDM.plot <- VOLCANO(Myeloid.DGE.HIVDMvHIVnegDM, filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/HIV.DM.NonDM.LAM.volcano.png', 
genes = c("CCL2", "MMP19", "FABP5", "APOE", "FABP4", "STAB1", "KLF2", "KLF4", "TREM2", "CEBPD","PLIN2", "CD9", "LPL", "CD36", "LIPA", "TGFBI", "LYVE1", "SPP1"), filter = "Macrophage")

#########################
# Other Statistics (Results Section 2)
#########################
# Median Proportion of Myeloid Cells by Study Group
Prop %>% group_by(StudyGroup) %>% 
    summarise(LAM = median(LAM),
    PVM = median(PVM),
    cMo = median(cMo))

# Median Proportion of Myeloid Cells by Glucose tolerance
Prop %>% group_by(Glucose) %>%
summarise(LAM = median(LAM))

# Median BMI by Study Group
Prop %>% group_by(StudyGroup) %>%
summarise(BMI = median(meta_bmi))

# Median Macrophage Proportion by Study Group
Prop %>% group_by(StudyGroup) %>%
summarise(Macrophages = median(LAM) + median(PVM))

# Median cDC1 Proportion by Study Group
Prop %>% group_by(Glucose) %>%
summarise(cDC1 = median(cDC1))

# Wilcoxon Sum Testing
wilcox.test(LAM~Glucose, data = Prop) # p = 0.0001
wilcox.test(PVM~Glucose, data = Prop) # p = 0.27
wilcox.test(cDC1~Glucose, data = Prop) # p < 0.001

# Spearman's Correlation LAM vs cMo
cor.test(Prop$LAM, Prop$cMo, method = "spearman")

# Ordinal Linear Regression LAM with glucose intolerance adjusting for age, bmi, HIV serostatus
rms::orm(LAM ~ Glucose + HIV + meta_bmi + meta_age, data = Prop, x = TRUE, y = TRUE)







