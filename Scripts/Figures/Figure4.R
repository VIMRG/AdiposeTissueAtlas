###############################################################################
######----------------------------FIGURE 4--------------------------------#####
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

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Figure_Functions.R')

###-------------------####
# Load Data
###-------------------####
Stromal <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_1.3.rds')
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
                  #plot.margin=unit(c(0,0,0,2.5),"cm"),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)
          ))

###-------------------####
# Environmental Parameters
###-------------------####
set.seed(7612)

#########################
# Figure 4A: Stromal UMAP
#########################
Idents(Stromal) <- "Manual_Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("ASC 1"))) <- 1
Idents(Stromal, cells = WhichCells(Stromal, idents = c("ASC 2"))) <- 2
Idents(Stromal, cells = WhichCells(Stromal, idents = c("ASC 3"))) <- 3
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Early Preadipocyte 1"))) <- 4
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Early Preadipocyte 2"))) <- 5
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Early Preadipocyte 3"))) <- 6
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Early Preadipocyte 4"))) <- 7
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Preadipocyte"))) <- 8
Idents(Stromal, cells = WhichCells(Stromal, idents = c("PTGDS+ Preadipocyte"))) <- 9
Idents(Stromal, cells = WhichCells(Stromal, idents = c("CD9hi Preadipocyte"))) <- 10
Idents(Stromal, cells = WhichCells(Stromal, idents = c("PTGDS+ ECM-Producing Preadipocyte"))) <- 11
Idents(Stromal, cells = WhichCells(Stromal, idents = c("ECM-Producing Preadipocyte"))) <- 12
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 1"))) <- 13
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 2"))) <- 14
Idents(Stromal, cells = WhichCells(Stromal, idents = c("ISG+ Preadipocyte"))) <- 15
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Metallothionein+ Preadipocyte/ASC"))) <- 16
Idents(Stromal, cells = WhichCells(Stromal, idents = c("MYOC+ Fibroblast"))) <- 17
Idents(Stromal, cells = WhichCells(Stromal, idents = c("PCOLCE+ Fibroblast"))) <- 18
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Myofibroblast"))) <- 19
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Lipomyofibroblast"))) <- 20
Idents(Stromal, cells = WhichCells(Stromal, idents = c("DKK3+ Fibroblast"))) <- 21
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Proliferating Myofibroblast"))) <- 22

Stromal[['Celltype_plot']] <- Idents(Stromal)

UMAP_FUN(Stromal, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure4/Figure4A.png', color = col_shuffle, breaks = c(1:22), width = 15, height = 10, 
label = c("Adipose Stem Cell 1 (1)", "Adipose Stem Cell 2 (2)", "Adipose Stem Cell 3 (3)", "Early Preadipocyte 1 (4)", "Early Preadipocyte 2 (5)", "Early Preadipocyte 3 (6)", "Early Preadipocyte 4 (7)", "Preadipocyte (8)",
"PTGDS+ Preadipocyte (9)", "CD9hi Preadipocyte (10)", "PTGDS+ ECM-Producing Preadipocyte (11)", "ECM-Producing Preadipocyte (12)", "Mature Preadipocyte 1 (13)", "Mature Preadipocyte 2 (14)",
"ISG+ Preadipocyte (15)", "Metallothionein+ Preadipocyte/ASC (16)", "MYOC+ Fibroblast (17)", 
"PCOLCE+ Fibroblast (18)", "Myofibroblast (19)", "Lipomyofibroblast (20)", "DKK3+ Fibroblast (21)", "Proliferating Myofibroblast (22)"))


#########################
# Figure 4B: Stromal DotPlot
#########################
Stromal[['Manual_Annotation']] <- factor(Stromal$Manual_Annotation, levels = c("PCOLCE+ Fibroblast", "MYOC+ Fibroblast", "ASC 1", "ASC 2", "ASC 3", "Early Preadipocyte 1", "Early Preadipocyte 2", "Early Preadipocyte 3",
"Early Preadipocyte 4", "PTGDS+ Preadipocyte", "CD9hi Preadipocyte", "PTGDS+ ECM-Producing Preadipocyte", "ECM-Producing Preadipocyte", "Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "ISG+ Preadipocyte", 
"Metallothionein+ Preadipocyte/ASC", "Myofibroblast", "Lipomyofibroblast", "DKK3+ Fibroblast", "Proliferating Myofibroblast"))

DOTPLOT_FUN(Stromal, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure4/Figure4B.png', height = 8, width = 18, features = c("CD55", "DPP4", "PCOLCE2", "MYOC", "ITM2A",
"IGFBP7", "DCN", "CLU", "LUM", "GSN", "CCN5", "PI16", "ZFP36", "MYC", "FOS", "CEBPD", "CXCL12", "CXCL14", "EGR1", "ATF3", "PTGDS", "CD9", "THY1", "COL3A1", "COL6A1", "FABP4", "C7", "CIDEC", "LPL", "ADIPOQ",
"ISG15", "IFI6", "MT1X", "MT2A", "POSTN", "TIMP1", "APOE", "DKK3", "MKI67"))

#########################
# Figure 4C: Stromal Proportion
#########################
Plot_prop(Prop_Long, filter = "Stromal", sort = "TRUE", y.axis = "Stromal", filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure4/Figure4C.png', width = 18, height = 9, Margins = c(0,0,0,0), plot_unit = "cm")

#########################
# Figure 4D: Correlation between LAM and Myofibroblasts
#########################
# Correlation between LAMs and Myofibroblasts
cor.test(Prop$LAM, Prop$Myofibroblast, method = "spearman") # rho = 0.74, p-value = 2.1e-15

# LAM vs Myofibroblast Plot
Prop <- Prop %>% mutate(StudyGroup = factor(StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic", "HIV- diabetic")))

Prop %>% ggplot(aes(x = LAM, y = Myofibroblast, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#F04832", "#F57636", "#EBAA3B", "#4287F5")) +
                xlab("Lipid-Associated Macrophages (% Total Myeloid)") +
                ylab("Myofibroblast (% Total Stromal)") + 
                geom_smooth(method = lm, se = FALSE) +
                theme(legend.position = "right", axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)))
                
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure4/Figure4D.png', dpi = 600, units = 'in', width = 10, height =10, device = 'png')

#########################
# Figure 4E: Correlation between CD4 TEM and Myofibroblasts
#########################
# CD4 vs LAM
CD4 <- Prop %>% dplyr::filter(CD4_number >= 30)

# Correlation Test
cor.test(CD4$CD4.TEM_CD4, CD4$Myofibroblast, method = "spearman") # rho = 0.22, p-value = 0.10

# Split by HIV status
CD4.HIV <- CD4 %>% filter(HIV == "HIVpos")

CD4.HIV %>% ggplot(aes(x = CD4.TEM_CD4, y = Myofibroblast, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#F04832", "#F57636", "#EBAA3B")) +
                xlab("CD4 T Effector Memory (% Total CD4)") +
                ylab("Myofibroblast (% Total Stromal)") + 
                geom_smooth(method = lm, se = FALSE) +
                theme(legend.position = "right", axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)))
                
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure4/Figure4E_1.png', dpi = 600, units = 'in', width = 10, height = 6, device = 'png')

# Diabetes Status Correlation
CD4.Diabetes <- CD4 %>% filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))

CD4.Diabetes %>% ggplot(aes(x = CD4.TEM_CD4, y = Myofibroblast, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#4287F5", "#EBAA3B")) +
                xlab("CD4 T Effector Memory (% Total CD4)") +
                ylab("Myofibroblast (% Total Stromal)") + 
                geom_smooth(method = lm, se = FALSE) +
                stat_cor(aes(color = StudyGroup), method = 'spearman') +
                theme(legend.position = "right", axis.title.x = element_text(margin = margin(t = 20), face = "bold"), axis.title.y = element_text(margin = margin(r = 20), face = "bold"))
                
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure4/Figure4E_2.png', dpi = 600, units = 'in', width = 10, height = 6, device = 'png')

#########################
# Figure 4F: Correlation between CD8 TEM and Myofibroblasts
#########################
CD8 <- Prop %>% 
        dplyr::filter(CD8_number >= 30)

# HIVpos Correlation
HIV.CD8 <- CD8 %>% dplyr::filter(HIV == 'HIVpos')
cor.test(HIV.CD8$CD8.TEM_CD8, HIV.CD8$Myofibroblast, method = 'spearman') #rho = 0.42, p-value = 0.004

HIV.CD8 %>% ggplot(aes(x = CD8.TEM_CD8, y = Myofibroblast, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#F04832", "#F57636", "#EBAA3B")) +
                xlab("CD8 T Effector Memory (% Total CD8)") +
                ylab("Myofibroblast (% Total Stromal)") + 
                geom_smooth(method = lm, se = FALSE) +
                theme(legend.position = "right", axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(r = 20)))

ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure4/Figure4F_1.png', dpi = 600, units = 'in', width = 10, height = 6, device = 'png')
 
# Diabetic Plot
Diabetic.CD8 <- CD8 %>% dplyr::filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic"))

# Diabetic Plot
Diabetic.CD8 %>% ggplot(aes(x = CD8.TEM_CD8, y = Myofibroblast, color = StudyGroup)) +
                geom_point(size = 2, l = 45) +
                scale_colour_manual(values = c("#4287F5", "#EBAA3B")) +
                xlab("CD8 T Effector Memory (% Total CD8)") +
                ylab("Myofibroblast (% Total Stromal)") + 
                geom_smooth(method = lm, se = FALSE) +
                stat_cor(aes(color = StudyGroup), method = 'spearman') +
                theme(legend.position = "right", axis.title.x = element_text(margin = margin(t = 20), face = "bold"), axis.title.y = element_text(margin = margin(r = 20), face = "bold"))

ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Figures/Figure4/Figure4F_2.png', dpi = 600, units = 'in', width = 10, height = 6, device = 'png')

#########################
# Other Statistics
#########################
# Correlation Between LAMS and Myofibroblasts (HIV+ only)
HIVpos <- Prop %>% filter(HIV == 'HIVpos')
cor.test(HIVpos$LAM, HIVpos$Myofibroblast, method = "spearman") # rho = 0.71, p-value = 2.4e-10

# Correlation Between LAMS and Myofibroblasts (HIV+ only)
HIVneg <- Prop %>% filter(HIV == 'HIVneg')
cor.test(HIVneg$LAM, HIVneg$Myofibroblast, method = "spearman") # rho = 0.79, p-value = 9.6e-06

# Partial Spearman LAM and Myofibroblasts
PResiduals::partial_Spearman(Myofibroblast | LAM ~ meta_age + meta_bmi + Glucose + HIV, data = Prop, link.x = "logit", link.y = "logit") # rho = 0.71 (0.55, 0.82), p = 3.6e-10

# Partial Spearman CD4 TEM and Myofibroblasts (HIV+ only)
PResiduals::partial_Spearman(Myofibroblast | CD4.TEM_CD4 ~ meta_age + meta_bmi + Glucose, data = CD4.HIV, link.x = "logit", link.y = "logit") # rho = 0.28, p = 0.08

# Partial Spearman CD4 TEM and Myofibroblasts (HIV+ only)
PResiduals::partial_Spearman(Myofibroblast | CD8.TEM_CD8 ~ meta_age + meta_bmi + Glucose, data = HIV.CD8, link.x = "logit", link.y = "logit") # rho = 0.27 (-0.5, 0.54), p = 0.09

# Partial Spearman LAM and CD9hi
PResiduals::partial_Spearman(CD9hi.Preadipocyte | LAM ~ meta_age + meta_bmi + Glucose + HIV, data = Prop, link.x = "logit", link.y = "logit") # rho = 0.27, p = 0.02
rms::orm(CD9hi.Preadipocyte~ LAM + meta_age + meta_bmi + Glucose + HIV, data = Prop, x = TRUE, y = TRUE)


