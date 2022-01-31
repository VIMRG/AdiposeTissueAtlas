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
# Figure 2C: Myeloid CITE-Seq Antibody
#########################
DefaultAssay(Myeloid) <- "ADT"
png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure2/Figure2C_CD14.png', res = 600, height = 5, width = 5, units = 'in')
FeaturePlot(Myeloid, features = c("CD14"), min.cutoff = 'q20', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure2/Figure2C_CD16.png', res = 600, height = 5, width = 5, units = 'in')
FeaturePlot(Myeloid, features = c("CD16"), min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure2/Figure2C_CD11C.png', res = 600, height = 5, width = 5, units = 'in')
FeaturePlot(Myeloid, features = c("CD11C"), min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure2/Figure2C_CD1C.png', res = 600, height = 5, width = 5, units = 'in')
FeaturePlot(Myeloid, features = c("CD1C"), min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

#########################
# Figure 2D: Myeloid Proportions
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







