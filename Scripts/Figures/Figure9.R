###############################################################################
######----------------------------FIGURE 9--------------------------------#####
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
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(stringi)
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(ggplot2)
library(scuttle)

###-------------------####
# Source Code
###-------------------####
source('../Utils.R')

###-------------------####
# Directories
###-------------------####
Subset_dir <- "../SubsetAnalysis_HIVneg"
Final_dir <- "../Final_Integration_HIVneg/"
fig_dir <- "../FiguresRevised"
Prop_dir <- "../Proportion" 

###-------------------####
# Load Seurat Object
###-------------------####
Integrated <- readRDS(paste0(Final_dir, "/", "Integrated_HIVneg.rds"))
Stromal <- readRDS(paste0(Subset_dir, "/", "Stromal_HIVneg.rds"))
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage_HIVneg.rds"))
Myeloid <- readRDS(paste0(Subset_dir, "/", "Myeloid_HIVneg.rds"))
Lymphoid <- readRDS(paste0(Subset_dir, "/", "Lymphoid_HIVneg.rds"))
Vascular <- readRDS(paste0(Subset_dir, "/", "Vascular_HIVneg.rds"))
Bcells <- readRDS(paste0(Subset_dir, "/", "Bcells_HIVneg.rds"))
CD4 <- readRDS(paste0(Subset_dir, "/", "CD4_HIVneg.rds"))
CD8 <- readRDS(paste0(Subset_dir, "/", "CD8_HIVneg.rds"))

###-------------------####
# Load Cell Proportions
###-------------------####
Overall.prop <- read.table(paste(Prop_dir, "Major.CellType.prop.HIVneg.txt", sep = "/"), check.names = F)
Stromal.Prop <- read.table(paste(Prop_dir, "Stromal.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD4.Prop <- read.table(paste(Prop_dir, "CD4.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD8.Prop <- read.table(paste(Prop_dir, "CD8.Prop.HIVneg.txt", sep = "/"), check.names = F)
Vascular.Prop <- read.table(paste(Prop_dir, "Vascular.Prop.HIVneg.txt", sep = "/"), check.names = F)
Myeloid.Prop <- read.table(paste(Prop_dir, "Myeloid.Prop.HIVneg.txt", sep = "/"), check.names = F)
Lymphoid.Prop <- read.table(paste(Prop_dir, "Lymphoid.Prop.HIVneg.txt", sep = "/"), check.names = F)
Macrophage.Prop <- read.table(paste(Prop_dir, "Macrophage.Prop.HIVneg.txt", sep = "/"), check.names = F)
CD4.Num <- read.table(paste(Prop_dir, "CD4.Num.HIVneg.txt", sep = "/"), check.names = F)
CD8.Num <- read.table(paste(Prop_dir, "CD8.Num.HIVneg.txt", sep = "/"), check.names = F)
Macrophage.Num <- read.table(paste(Prop_dir, "Macrophage.Num.HIVneg.txt", sep = "/"), check.names = F)

# Remove NAs
Overall.prop <- Overall.prop %>% na.omit()
Stromal.Prop <- Stromal.Prop %>% na.omit()
Myeloid.Prop <- Myeloid.Prop %>% na.omit()
Lymphoid.Prop <- Lymphoid.Prop %>% na.omit()
Vascular.Prop <- Vascular.Prop %>% na.omit()
CD4.Prop <- CD4.Prop %>% na.omit()
CD8.Prop <- CD8.Prop %>% na.omit()
Macrophage.Prop <- Macrophage.Prop %>% na.omit()

###-------------------####
# Merge with Metadata
###-------------------####
# Main Subsets
Overall.prop <- Prop_Merge(Overall.prop, data_hatim)
Stromal.prop <- Prop_Merge(Stromal.Prop, data_hatim)
Myeloid.prop <- Prop_Merge(Myeloid.Prop, data_hatim)
Vascular.prop <- Prop_Merge(Vascular.Prop, data_hatim)
Lymphoid.prop <- Prop_Merge(Lymphoid.Prop, data_hatim)

# Macrophage Subset
Macrophage.prop <- Prop_Merge(Macrophage.Prop, data_hatim)
colnames(Macrophage.Num) <- c("HATIMID", "Macnum")
Macrophage.Num$HATIMID <- as.factor(as.character(Macrophage.Num$HATIMID))

Macrophage.prop <- left_join(Macrophage.prop, Macrophage.Num, by = "HATIMID")
Macrophage.prop <- Macrophage.prop %>% dplyr::filter(Macnum > 30) # Filter > 30
Macrophage.prop <- Macrophage.prop %>% dplyr::select(-c(Macnum))

# CD4 Subset
CD4.prop <- Prop_Merge(CD4.Prop, data_hatim)
colnames(CD4.Num) <- c("HATIMID", "CD4num")
CD4.Num$HATIMID <- as.factor(as.character(CD4.Num$HATIMID))

CD4.prop <- left_join(CD4.prop, CD4.Num, by = "HATIMID")
CD4.prop <- CD4.prop %>% dplyr::filter(CD4num > 30) # Filter > 30
CD4.prop <- CD4.prop %>% dplyr::select(-c(CD4num))

# CD8 Subset
CD8.prop <- Prop_Merge(CD8.Prop, data_hatim)
colnames(CD8.Num) <- c("HATIMID", "CD8num")
CD8.Num$HATIMID <- as.factor(as.character(CD8.Num$HATIMID))

CD8.prop <- left_join(CD8.prop, CD8.Num, by = "HATIMID")
CD8.prop <- CD8.prop %>% dplyr::filter(CD8num > 30) # Filter > 30
CD8.prop <- CD8.prop %>% dplyr::select(-c(CD8num))

###-------------------####
# Figure 9A: Major Cell Type Proportions
###-------------------####
Major.prop <- Overall.prop %>% filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic", "HIV+ non-diabetic"))
Prop_long <- Major.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "All Cells", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ diabetic", "HIV- diabetic"), c("HIV+ diabetic", "HIV+ non-diabetic"), c("HIV- diabetic", "HIV+ non-diabetic")), specify_x_order = NA, ncol = 5, plot_cols = c("#91C46C", "#FFBE00", "#e69f00"))
                    
ggsave(file = paste(fig_dir, "Figure9A.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 9B: Macrophage Proportion BoxPlot
###-------------------####
Macrophage.prop <- Macrophage.prop %>% filter(StudyGroup %in% c("HIV+ non-diabetic", "HIV+ diabetic", "HIV- diabetic"))
Macrophage.prop$StudyGroup <- droplevels(Macrophage.prop$StudyGroup)
Prop_long <- Macrophage.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "Macrophage", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ diabetic", "HIV- diabetic"), c("HIV+ non-diabetic", "HIV- diabetic"), c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 5, plot_cols = c("#91C46C", "#FFBE00", "#e69f00"), legend.pos = "none")
                    
ggsave(file = paste(fig_dir, "Figure9B.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 9C: CD4 Proportion BoxPlot
###-------------------####
CD4.prop <- CD4.prop %>% filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic", "HIV+ non-diabetic"))
Prop_long <- CD4.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "CD4", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ diabetic", "HIV- diabetic"), c("HIV+ non-diabetic", "HIV- diabetic"), c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 5, plot_cols = c("#91C46C", "#FFBE00", "#e69f00"), legend.pos = "none")
                    
ggsave(file = paste(fig_dir, "Figure9C.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 9D: CD8 Proportion BoxPlot
###-------------------####
CD8.prop <- CD8.prop %>% filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic", "HIV+ non-diabetic"))
Prop_long <- CD8.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "CD8", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ diabetic", "HIV- diabetic"), c("HIV+ non-diabetic", "HIV- diabetic"), c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 4, plot_cols = c("#91C46C", "#FFBE00", "#e69f00"), legend.pos = "none")
                    
ggsave(file = paste(fig_dir, "Figure9D.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 9E: Vascular Proportion BoxPlot
###-------------------####
Vascular.prop <- Vascular.prop %>% filter(StudyGroup %in% c("HIV+ diabetic", "HIV- diabetic", "HIV+ non-diabetic"))
Vascular.prop <- Vascular.prop[,-c(3)]
Prop_long <- Vascular.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "Vascular", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ diabetic", "HIV- diabetic"), c("HIV+ non-diabetic", "HIV- diabetic"), c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 7, plot_cols = c("#91C46C", "#FFBE00", "#e69f00"), legend.pos = "none")
                    
ggsave(file = paste(fig_dir, "Figure9E.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 9F: Stromal vs BMI/Age Heatmap
###-------------------####
Cor_rho_neg <- function(df, Measures) {
    Results <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(Results) <- c("rho", "pvalue", "Cluster", "Measure")
    numb <- dim(df)[2] - 10
    
    for (measure in Measures) {
        if (measure == "BMI") {
        for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,measure] ~ Sex  + Age, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i], Measure = measure)
                Results <- rbind(Results, res)
            }
        } else if (measure == "Age") {
             for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,measure] ~ Sex + BMI, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i], Measure = measure)
                Results <- rbind(Results, res)
            }
        }
    } 
    Results$padj <- p.adjust(Results$pvalue, method = "BH")
    return(Results)
}

df <- Stromal.prop %>% filter(StudyGroup == "HIV- diabetic")
Cor_res <- Cor_rho_neg(df, Measures = c("BMI", "Age"))

plot <- Heatmap_circle(Cor_res, limits = c(-.60, .60), legend_breaks = c(-0.5, 0, 0.5), order = c(),
                       p_breaks = c(0.01, 0.05, 0.1, 0.4, 0.8), range = c(15,5))
ggsave(paste0(fig_dir, "/", "Figure9F.png"), device = "png", height = 7, width = 10, units = "in", dpi = 300)

###-------------------####
# Figure 9G: Vascular vs BMI/Age Heatmap
###-------------------####
df <- Vascular.prop %>% filter(StudyGroup == "HIV- diabetic")
Cor_res <- Cor_rho_neg(df, Measures = c("BMI", "Age"))

plot <- Heatmap_circle(Cor_res, limits = c(-.70, .70), legend_breaks = c(-0.5, 0, 0.5), order = c("Arterial EC", "Capillary EC", "Venous EC", "Pericyte", "VSMC 1", "VSMC 2", "Cycling Vascular Cells"),
                       p_breaks = c(0.01, 0.05, 0.1, 0.4), range = c(15,5))
ggsave(paste0(fig_dir, "/", "Figure9G.png"), device = "png", height = 7, width = 10, units = "in", dpi = 300)


## Proportion Multiple Comparison Adjustment
###-------------------####
# Miscellaneous Data For Paper
###-------------------####
# Multiple Comparisons Adjustment
FDR_Fun <- function(df1, df2, df3, clusters) {
    Results_HIV <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(Results_HIV) <- c("Cluster", "pvalue", "StudyGroup")
    for (i in clusters) {
        new <- data.frame(Cluster = i, pvalue = df1$results$.ALL.$testresults[[i]]$P, StudyGroup = "HIVDMvHIVnonDM", prop1 = round(median(df1$results$.ALL.$data[[i]]$`HIV+ non-diabetic`), digits = 1), prop2 = round(median(df1$results$.ALL.$data[[i]]$`HIV+ diabetic`), digits = 1))
        Results_HIV <- rbind(Results_HIV, new)
    }    
    Results_nonDM <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(Results_nonDM) <- c("Celltype", "pvalue", "StudyGroup")
    for (i in clusters) {
        new <- data.frame(Cluster = i, pvalue = df2$results$.ALL.$testresults[[i]]$P, StudyGroup = "HIVnonDMvHIVneg", prop1 = round(median(df2$results$.ALL.$data[[i]]$`HIV+ non-diabetic`), digits = 1), prop2 = round(median(df2$results$.ALL.$data[[i]]$`HIV- diabetic`), digits = 1))
        Results_nonDM <- rbind(Results_nonDM, new)
    }
    Results_DM <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(Results_DM) <- c("Celltype", "pvalue", "StudyGroup")
    for (i in clusters) {
        new <- data.frame(Cluster = i, pvalue = df3$results$.ALL.$testresults[[i]]$P, StudyGroup = "HIVvHIVneg", prop1 = round(median(df3$results$.ALL.$data[[i]]$`HIV+ diabetic`), digits = 1), prop2 = round(median(df3$results$.ALL.$data[[i]]$`HIV- diabetic`), digits = 1))
        Results_DM <- rbind(Results_DM, new)
    }
    
    Results <- rbind(Results_DM, Results_nonDM, Results_HIV)
    Results$p.adj <- p.adjust(Results$pvalue, method = "BH")
    return(Results)
}


# Overall
DM <- Hmisc::summaryM(Myeloid + Lymphoid + Stromal + Vascular ~ StudyGroup, data = Major.prop[Major.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Neg <- Hmisc::summaryM(Myeloid + Lymphoid + Stromal + Vascular ~ StudyGroup, data = Major.prop[Major.prop$StudyGroup != "HIV- diabetic",], overall = TRUE, test = TRUE)
Non <- Hmisc::summaryM(Myeloid + Lymphoid + Stromal + Vascular ~ StudyGroup, data = Major.prop[Major.prop$StudyGroup != "HIV+ non-diabetic",], overall = TRUE, test = TRUE)
Results <- FDR_Fun(df1 = Neg, df2 = DM, df3 = Non, clusters = c("Myeloid", "Lymphoid", "Stromal", "Vascular"))
write.csv(Results_Mac, file = paste(Prop_dir, "Overall.Proportion.HIVneg.Adjusted.csv", sep = "/"))

# Macrophage
Mac_DM <- Hmisc::summaryM(LAM + PVM + IM + `Mo-Mac 1` + `Mo-Mac 2` ~ StudyGroup, data = Macrophage.prop[Macrophage.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Mac_Neg <- Hmisc::summaryM(LAM + PVM + IM + `Mo-Mac 1` + `Mo-Mac 2` ~ StudyGroup, data = Macrophage.prop[Macrophage.prop$StudyGroup != "HIV- diabetic",], overall = TRUE, test = TRUE)
Mac_Non <- Hmisc::summaryM(LAM + PVM + IM + `Mo-Mac 1` + `Mo-Mac 2` ~ StudyGroup, data = Macrophage.prop[Macrophage.prop$StudyGroup != "HIV+ non-diabetic",], overall = TRUE, test = TRUE)
Results_Mac <- FDR_Fun(df1 = Mac_Neg, df2 = Mac_DM, df3 = Mac_Non, clusters = c("LAM", "IM", "PVM", "Mo-Mac 1", "Mo-Mac 2"))
write.csv(Results_Mac, file = paste(Prop_dir, "Macrophage.Proportion.HIVneg.Adjusted.csv", sep = "/"))

# CD4
CD4_DM <- Hmisc::summaryM(`CD4 Naive` + `CD4 TCM` + `CD4 TEM` + `CD4 Cytotoxic` + `CD4 Regulatory` ~ StudyGroup, data = CD4.prop[CD4.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
CD4_Neg <- Hmisc::summaryM(`CD4 Naive` + `CD4 TCM` + `CD4 TEM` + `CD4 Cytotoxic` + `CD4 Regulatory` ~ StudyGroup, data = CD4.prop[CD4.prop$StudyGroup != "HIV- diabetic",], overall = TRUE, test = TRUE)
CD4_Non <- Hmisc::summaryM(`CD4 Naive` + `CD4 TCM` + `CD4 TEM` + `CD4 Cytotoxic` + `CD4 Regulatory` ~ StudyGroup, data = CD4.prop[CD4.prop$StudyGroup != "HIV+ non-diabetic",], overall = TRUE, test = TRUE)
Results_CD4 <- FDR_Fun(df1 = CD4_Neg, df2 = CD4_DM, df3 = CD4_Non, clusters = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Cytotoxic", "CD4 Regulatory"))
write.csv(Results_CD4, file = paste(Prop_dir, "CD4.Proportion.HIVneg.Adjusted.csv", sep = "/"))

CD8_DM <- Hmisc::summaryM(`CD8 Naive` + `CD8 TCM` + `CD8 TEM` + `CD8 Cytotoxic` ~ StudyGroup, data = CD8.prop[CD8.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
CD8_Neg <- Hmisc::summaryM(`CD8 Naive` + `CD8 TCM` + `CD8 TEM` + `CD8 Cytotoxic` ~ StudyGroup, data = CD8.prop[CD8.prop$StudyGroup != "HIV- diabetic",], overall = TRUE, test = TRUE)
CD8_Non <- Hmisc::summaryM(`CD8 Naive` + `CD8 TCM` + `CD8 TEM` + `CD8 Cytotoxic` ~ StudyGroup, data = CD8.prop[CD8.prop$StudyGroup != "HIV+ non-diabetic",], overall = TRUE, test = TRUE)
Results_CD8 <- FDR_Fun(df1 = CD8_Neg, df2 = CD8_DM, df3 = CD8_Non, clusters = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Cytotoxic"))
write.csv(Results_CD8, file = paste(Prop_dir, "CD8.Proportion.HIVneg.Adjusted.csv", sep = "/"))

Vascular_DM <- Hmisc::summaryM(`Capillary EC` + `Venous EC` + `VSMC 1` + `Cycling Vascular Cells` +  `Pericyte` + `VSMC 2` + `Arterial EC` ~ StudyGroup, data = Vascular.prop[Vascular.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Vascular_Neg <- Hmisc::summaryM(`Capillary EC` + `Venous EC` + `VSMC 1` + `Cycling Vascular Cells` +  `Pericyte` + `VSMC 2` + `Arterial EC` ~ StudyGroup, data = Vascular.prop[Vascular.prop$StudyGroup != "HIV- diabetic",], overall = TRUE, test = TRUE)
Vascular_Non <- Hmisc::summaryM(`Capillary EC` + `Venous EC` + `VSMC 1` + `Cycling Vascular Cells` +  `Pericyte` + `VSMC 2` + `Arterial EC` ~ StudyGroup, data = Vascular.prop[Vascular.prop$StudyGroup != "HIV+ non-diabetic",], overall = TRUE, test = TRUE)
Results_Vascular <- FDR_Fun(df1 = Vascular_Neg, df2 = Vascular_DM, df3 = Vascular_Non, clusters = c("Capillary EC", "Venous EC", "VSMC 1", "Cycling Vascular Cells", "Pericyte", "VSMC 2", "Arterial EC"))
write.csv(Results_Vascular, file = paste(Prop_dir, "Vascular.Proportion.HIVneg.Adjusted.csv", sep = "/"))

Stromal_DM <- Hmisc::summaryM(`PCOLCE+ Fibroblast` + `MYOC+ Fibroblast` + `Adipose Progenitor Cell 1` + `Adipose Progenitor Cell 2` + `Early Preadipocyte` + `ECM-Producing Early Preadipocyte` + `Mature Preadipocyte` +
`Myofibroblast` + `Metallothionein+ Preadipocyte` + `Cycling Myofibroblast` ~ StudyGroup, data = Stromal.prop[Stromal.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Stromal_Neg <- Hmisc::summaryM(`PCOLCE+ Fibroblast` + `MYOC+ Fibroblast` + `Adipose Progenitor Cell 1` + `Adipose Progenitor Cell 2` + `Early Preadipocyte` + `ECM-Producing Early Preadipocyte` + `Mature Preadipocyte` +
`Myofibroblast` + `Metallothionein+ Preadipocyte` + `Cycling Myofibroblast` ~ StudyGroup, data = Stromal.prop[Stromal.prop$StudyGroup != "HIV- diabetic",], overall = TRUE, test = TRUE)
Stromal_Non <- Hmisc::summaryM(`PCOLCE+ Fibroblast` + `MYOC+ Fibroblast` + `Adipose Progenitor Cell 1` + `Adipose Progenitor Cell 2` + `Early Preadipocyte` + `ECM-Producing Early Preadipocyte` + `Mature Preadipocyte` +
`Myofibroblast` + `Metallothionein+ Preadipocyte` + `Cycling Myofibroblast` ~ StudyGroup, data = Stromal.prop[Stromal.prop$StudyGroup != "HIV+ non-diabetic",], overall = TRUE, test = TRUE)
Results_Stromal <- FDR_Fun(df1 = Stromal_Neg, df2 = Stromal_DM, df3 = Stromal_Non, clusters = c('PCOLCE+ Fibroblast', 'MYOC+ Fibroblast', 'Adipose Progenitor Cell 1', 'Adipose Progenitor Cell 2', 'Early Preadipocyte', 'ECM-Producing Early Preadipocyte', 'Mature Preadipocyte',
'Myofibroblast', 'Metallothionein+ Preadipocyte', 'Cycling Myofibroblast'))
write.csv(Results_Stromal, file = paste(Prop_dir, "Stromal.Proportion.HIVneg.Adjusted.csv", sep = "/"))
