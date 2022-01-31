###############################################################################
######----------------------------FIGURE 3--------------------------------#####
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

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Figure_Functions.R')

###-------------------####
# Load Data
###-------------------####
Tcells <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_10.30.rds')
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
                  plot.margin=unit(c(0,0,0,2.5),"cm"),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)
          ))

#########################
# Figure 3A: T Cell UMAP
#########################
Idents(Tcells) <- "Manual_Annotation"
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD4 Naive"))) <- 1
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD4 TCM"))) <- 2
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD4 TEM"))) <- 3
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD4 TEMRA & Senescent"))) <- 4
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD4 Regulatory"))) <- 5
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD8 Naive"))) <- 6
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD8 TCM"))) <- 7
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD8 TEM"))) <- 8
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD8 TEMRA & Senescent"))) <- 9
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD57+ mNK"))) <- 10
Idents(Tcells, cells = WhichCells(Tcells, idents = c("CD16+ mNK"))) <- 11
Idents(Tcells, cells = WhichCells(Tcells, idents = c("Immature NK"))) <- 12
Idents(Tcells, cells = WhichCells(Tcells, idents = c("MAIT"))) <- 13
Idents(Tcells, cells = WhichCells(Tcells, idents = c("Gamma Delta"))) <- 14
Idents(Tcells, cells = WhichCells(Tcells, idents = c("ILC"))) <- 15
Idents(Tcells, cells = WhichCells(Tcells, idents = c("Cycling T & NK"))) <- 16

Tcells[['Celltype_plot']] <- Idents(Tcells)

UMAP_FUN(Tcells, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3A.png', color = col_shuffle, breaks = c(1:16), 
label = c("CD4 Naive (1)","CD4 TCM (2)", "CD4 TEM (3)", "CD4 TEMRA & Senescent (4)", "CD4 Regulatory (5)", "CD8 Naive (6)", "CD8 TCM (7)", "CD8 TEM (8)", "CD8 TEMRA & Senescent (9)",
"CD57+ mNK (10)", "CD16+ mNK (11)", "Immature NK (12)", "MAIT (13)", "Gamma Delta (14)", "ILC (15)", "Cycling T & NK (16)"))

#########################
# Figure 3B: T Cell DotPlot
#########################
Tcells[['Manual_Annotation']] <- factor(Tcells$Manual_Annotation, levels = c('CD4 Naive', 'CD4 TCM', 'CD4 TEM', 'CD4 TEMRA & Senescent', 'CD4 Regulatory', 'CD8 Naive', 'CD8 TCM', 
'CD8 TEM', 'CD8 TEMRA & Senescent', 'MAIT', 'Gamma Delta', 'CD16+ mNK', 'CD57+ mNK', 'Immature NK', 'ILC', "Cycling T & NK"))

DOTPLOT_FUN(Tcells, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3B.png', height = 8, width = 14, features = c("IL7R", "SELL", "LEF1",
"CCR7", "LTB", "LDHB", "GPR183", "COTL1", "ALOX5AP", "CCL5", "NKG7", "GNLY", "PRF1", "FOXP3", "CD8A", "CCR6", "TRDV1", "KLRC2", "FCER1G", "FCGR3A", "XCL1", "KIT", "GATA2", "SPINK2", "MKI67"))

#########################
# Figure 3C: T Cell CITE-Seq Ab
#########################
DefaultAssay(Tcells) <- "ADT"
png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3C_CD4.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Tcells, features = "CD4", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3C_CD8.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Tcells, features = "CD8", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3C_CD27.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Tcells, features = "CD27", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3C_CD45RA.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Tcells, features = "CD45RA", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3C_CD57.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Tcells, features = "CD57", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

png('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3C_CD16.png', width = 5, height = 5, units = 'in', res = 600)
FeaturePlot(Tcells, features = "CD16", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()

#########################
# Figure 3D: T Cell Proportions
#########################
Plot_prop(Prop_Long, filter = "Tcell", sort = "TRUE", y.axis = "T & NK", filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3D.png', width = 20, height = 8, Margins = c(0,0,0,4), plot_unit = 'cm')

#########################
# Figure 3E: CD4 TEM & Macrophage Correlation by HIV status
#########################
# Filter out individuals contributing < 30 CD4 T cells
CD4.LAM <- Prop %>% 
        dplyr::filter(CD4_number >= 30)

# HIVpos Correlation
HIV.CD4.LAM <- CD4.LAM %>% dplyr::filter(HIV == 'HIVpos')
cor.test(HIV.CD4.LAM$CD4.TEM_CD4, HIV.CD4.LAM$LAM, method = 'spearman') #rho = 0.67, p-value = 2.69e-6

# HIV+ Plot
HIV.CD4.LAM %>% ggplot(aes(x = CD4.TEM_CD4, y = LAM, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#F04832", "#F57636", "#EBAA3B")) +
                xlab("CD4 T Effector Memory (% Total CD4)") +
                ylab("Lipid-Associated Macrophages (% Total Myeloid)") + 
                geom_smooth(method = lm, se = FALSE) +
                theme(legend.position = "bottom", axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)))
                
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3E_HIVpos.png', dpi = 600, units = 'in', width = 10, height = 10, device = 'png')


# Diabetic Plot
Diabetic.CD4.LAM <- CD4.LAM %>% dplyr::filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))

Diabetic.CD4.LAM %>% ggplot(aes(x = CD4.TEM_CD4, y = LAM, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#4287F5", "#EBAA3B")) +
                xlab("CD4 T Effector Memory (% Total CD4)") +
                ylab("Lipid-Associated Macrophages (% Total Myeloid)") + 
                geom_smooth(method = lm, se = FALSE) +
                stat_cor(aes(color = StudyGroup), method = 'spearman') +
                theme(legend.position = "right", axis.title.x = element_text(margin = margin(t = 20), face = "bold"), axis.title.y = element_text(margin = margin(r = 20), face = "bold"))

ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3E_HIVneg.png', dpi = 600, units = 'in', width = 10, height = 10, device = 'png')

#########################
# Figure 3F: CD8 TEM & Macrophage Correlation by HIV status
#########################
CD8.LAM <- Prop %>% 
        dplyr::filter(CD8_number >= 30)

# HIVpos Correlation
HIV.CD8.LAM <- CD8.LAM %>% dplyr::filter(HIV == 'HIVpos')
cor.test(HIV.CD8.LAM$CD8.TEM_CD8, HIV.CD8.LAM$LAM, method = 'spearman') #rho = 0.63, p-value = 3e-6

# HIV+ Plot
HIV.CD8.LAM %>% ggplot(aes(x = CD8.TEM_CD8, y = LAM, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#F04832", "#F57636", "#EBAA3B")) +
                xlab("CD8 T Effector Memory (% Total CD8)") +
                ylab("Lipid-Associated Macrophages (% Total Myeloid)") + 
                geom_smooth(method = lm, se = FALSE) +
                theme(legend.position = "bottom", axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)))
                
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3F_HIVpos.png', dpi = 600, units = 'in', width = 10, height = 10, device = 'png')

# Diabetic Plot
Diabetic.CD8.LAM <- CD8.LAM %>% dplyr::filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))

Diabetic.CD8.LAM %>% ggplot(aes(x = CD8.TEM_CD8, y = LAM, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#4287F5", "#EBAA3B")) +
                xlab("CD8 T Effector Memory (% Total CD8)") +
                ylab("Lipid-Associated Macrophages (% Total Myeloid)") + 
                geom_smooth(method = lm, se = FALSE) +
                stat_cor(aes(color = StudyGroup), method = 'spearman') +
                theme(legend.position = "bottom", axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)))
                
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure3/Figure3F_HIVneg.png', dpi = 600, units = 'in', width = 10, height = 10, device = 'png')

#########################
# Other Statistics (Results Section 3)
#########################
# Proportion CD4 and CD8 T cells by Study Group
Merged <- merge(CD4, CD8)
prop.table(table(Merged$StudyGroup, Merged$Tcell_subcluster), margin = 1)*100

# Grab only CD4 and CD8 with >= 30 cells
CD4 <- Prop %>% filter(CD4_number >= 30)
CD8 <- Prop %>% filter(CD8_number >= 30)

# Proportion Gamma Delta
Prop %>% group_by(HIV) %>%
summarise(GD = median(Gamma.Delta))

# Proportion iNK
Prop %>% group_by(StudyGroup) %>%
summarise(iNK = median(Immature.NK))

# Proportion CD57+ mNK
Prop %>% group_by(HIV) %>%
summarise(mNK = median(CD57..mNK))

# Proportion MAIT
Prop %>% group_by(StudyGroup) %>%
summarise(Mait = median(MAIT))

# Proportion CD4 T cells
CD4 %>% group_by(StudyGroup) %>%
summarise(CD4TEM = median(CD4.TEM_CD4),
CD4TEMRA = median(CD4.TEMRA...Senescent_CD4),
CD4Reg = median(CD4.Regulatory_CD4))

# Proportion CD4 T cells by Glucose Tolerance
CD4 %>% filter(HIV == 'HIVpos') %>% group_by(Glucose) %>%
summarise(CD4TEM = median(CD4.TEM_CD4),
CD4TEMRA = median(CD4.TEMRA...Senescent_CD4),
CD4Reg = median(CD4.Regulatory_CD4))

# Proportion CD8 T cells
CD8 %>% group_by(StudyGroup) %>%
summarise(CD8TEM = median(CD8.TEM_CD8),
CD8TEMRA = median(CD8.TEMRA...Senescent_CD8))

# Proportion CD8 T cells by Glucose intolerance
CD8 %>% filter(HIV == 'HIVpos') %>% group_by(Glucose) %>%
summarise(CD8TEM = median(CD8.TEM_CD8),
CD8TEMRA = median(CD8.TEMRA...Senescent_CD8))

# Wilcoxon sum test
wilcox.test(Gamma.Delta~StudyGroup, data = Prop, subset = StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))
wilcox.test(Immature.NK~StudyGroup, data = Prop, subset = StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))
wilcox.test(CD57..mNK~HIV, data = Prop)
wilcox.test(Gamma.Delta~HIV, data = Prop)

wilcox.test(CD4.TEM_CD4~StudyGroup, data = CD4, subset = StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))
wilcox.test(CD4.TEMRA...Senescent_CD4~StudyGroup, data = CD4, subset = StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))
wilcox.test(CD4.Regulatory_CD4~StudyGroup, data = CD4, subset = StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))
wilcox.test(CD4.TEM_CD4~Glucose, data = CD4, subset = HIV %in% c("HIVpos"))
wilcox.test(CD4.TEMRA...Senescent_CD4~Glucose, data = CD4, subset = HIV %in% c("HIVpos"))

wilcox.test(CD8.TEMRA...Senescent_CD8~StudyGroup, data = CD8, subset = StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))
wilcox.test(CD8.TEM_CD8~Glucose, data = CD8, subset = HIV %in% c("HIVpos"))
wilcox.test(CD8.TEM_CD8~StudyGroup, data = Prop, subset = StudyGroup %in% c("HIV+ diabetic", "HIV+ non-diabetic"))


# HIV+ CD4 TEM & LAM adjusted for Age, BMI, Glucose intolerance
PResiduals::partial_Spearman(LAM | CD4.TEM_CD4 ~ meta_age + meta_bmi + Glucose, data = HIV.CD4.LAM, link.x = "logit", link.y = "logit")
# rho = 0.53 (0.24, 0.73), p = 0.0008

# HIV+ CD8 TEM & LAM adjusted for Age, BMI, Glucose intolerance
PResiduals::partial_Spearman(LAM | CD8.TEM_CD8 ~ meta_age + meta_bmi + Glucose, data = HIV.CD8.LAM, link.x = "logit", link.y = "logit")
# rho = 0.42 (0.10), p = 0.67

