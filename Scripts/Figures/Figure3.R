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
library(rms)
library(gridExtra)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Hmisc)

###-------------------####
# Source Code
###-------------------####
source('../Utils.R')

###-------------------####
# Directories
###-------------------####
tmp_dir <- "../SubsetAnalysis"
Prop_dir <- "../Proportion"
fig_dir <- "../Figures"

###-------------------####
# Load Seurat Object
###-------------------####
Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage.rds"))
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8.rds"))

###-------------------####
# Load Cell Proportions
###-------------------####
Lymphoid.Prop <- read.table(paste(Prop_dir, "Lymphoid.Prop.txt", sep = "/"), check.names = FALSE)
Myeloid.Prop <- read.table(paste(Prop_dir, "Myeloid.Prop.txt", sep = "/"), check.names = F)
Stromal.Prop <- read.table(paste(Prop_dir, "Stromal.Prop.txt", sep = "/"), check.names = F)
Vascular.Prop <- read.table(paste(Prop_dir, "Vascular.Prop.txt", sep = "/"), check.names = F)
CD4.Prop <- read.table(paste(Prop_dir, "CD4.Prop.txt", sep = "/"), check.names = F)
CD8.Prop <- read.table(paste(Prop_dir, "CD8.Prop.txt", sep = "/"), check.names = F)
Macrophage.Prop <- read.table(paste(Prop_dir, "Macrophage.Prop.txt", sep = "/"), check.names = F)
CD4.Num <- read.table(paste(Prop_dir, "CD4.Num.txt", sep = "/"), check.names = F)
CD8.Num <- read.table(paste(Prop_dir, "CD8.Num.txt", sep = "/"), check.names = F)
Macrophage.Num <- read.table(paste(Prop_dir, "Macrophage.Num.txt", sep = "/"), check.names = F)

# Remove NAs
Lymphoid.Prop <- Lymphoid.Prop %>% na.omit()
Myeloid.Prop <- Myeloid.Prop %>% na.omit()
Stromal.Prop <- Stromal.Prop %>% na.omit()
Vascular.Prop <- Vascular.Prop %>% na.omit()
CD4.Prop <- CD4.Prop %>% na.omit()
CD8.Prop <- CD8.Prop %>% na.omit()
Macrophage.Prop <- Macrophage.Prop %>% na.omit()

###-------------------####
# Merge with Metadata
###-------------------####
# Main Subsets
Myeloid.prop <- Prop_Merge(Myeloid.Prop, data_hatim)
Lymphoid.prop <- Prop_Merge(Lymphoid.Prop, data_hatim)
Stromal.prop <- Prop_Merge(Stromal.Prop, data_hatim)
Vascular.prop <- Prop_Merge(Vascular.Prop, data_hatim)

# Macrophage Subset (Include individuals > 30 cells)
Macrophage.prop <- Prop_Merge(Macrophage.Prop, data_hatim)
colnames(Macrophage.Num) <- c("HATIMID", "Macnum")
Macrophage.Num$HATIMID <- as.factor(as.character(Macrophage.Num$HATIMID))

Macrophage.prop <- left_join(Macrophage.prop, Macrophage.Num, by = "HATIMID")
Macrophage.prop <- Macrophage.prop %>% dplyr::filter(Macnum > 30) # Filter > 30
Macrophage.prop <- Macrophage.prop %>% dplyr::select(-c(Macnum))

# CD4 Subset (Include individuals > 30 cells)
CD4.prop <- Prop_Merge(CD4.Prop, data_hatim)
colnames(CD4.Num) <- c("HATIMID", "CD4num")
CD4.Num$HATIMID <- as.factor(as.character(CD4.Num$HATIMID))

CD4.prop <- left_join(CD4.prop, CD4.Num, by = "HATIMID")
CD4.prop <- CD4.prop %>% dplyr::filter(CD4num > 30) # Filter > 30
CD4.prop <- CD4.prop %>% dplyr::select(-c(CD4num))

# CD8 Subset (Include individuals > 30 cells)
CD8.prop <- Prop_Merge(CD8.Prop, data_hatim)
colnames(CD8.Num) <- c("HATIMID", "CD8num")
CD8.Num$HATIMID <- as.factor(as.character(CD8.Num$HATIMID))

CD8.prop <- left_join(CD8.prop, CD8.Num, by = "HATIMID")
CD8.prop <- CD8.prop %>% dplyr::filter(CD8num > 30) # Filter > 30
CD8.prop <- CD8.prop %>% dplyr::select(-c(CD8num))

###-------------------####
# Figure 3A: Macrophage BoxPlot
###-------------------####
Prop_long <- Macrophage.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "Macrophage", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 5)
                    
ggsave(file = paste(fig_dir, "Figure3A.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 3B: Continuous Glucose Measures Macrophage
###-------------------####
Cor_res <- Cor_rho(Macrophage.prop, Measures = c("meta_hba1c", "meta_fbg"))

plot <- Heatmap_circle(df = Cor_res, order = c("LAM", "IM", "PVM", "Mo-Mac 1", "Mo-Mac 2"))
ggsave(file = paste(fig_dir, "Figure3B.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 3C: Continuous Glucose Measures Myeloid
###-------------------####
Cor_res <- Cor_rho(Myeloid.prop, Measures = c("meta_hba1c", "meta_fbg"))

plot <- Heatmap_circle(df = Cor_res, order = c("cMo", "nMo", "Other Mo", "ISG+ Mo", "Mo-Mac 1", "Mo-Mac 2", "PVM", "IM", "LAM", "cDC2B", "cDC1", "Migratory DC", "pDC", "Cycling Myeloid"))
ggsave(file = paste(fig_dir, "Figure3C.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 3D: CD8 Boxplots
###-------------------####
Prop_long <- CD8.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "CD8 T Cells", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 4)
                    
ggsave(file = paste(fig_dir, "Figure3D.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 3E: Continuous Glucose Measures CD8 T cells
###-------------------####
Cor_res <- Cor_rho(CD8.prop, Measures = c("meta_hba1c", "meta_fbg"))

plot <- Heatmap_circle(df = Cor_res, order = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Cytotoxic"), legend_breaks = c(-0.4, 0, 0.4))
ggsave(file = paste(fig_dir, "Figure3E.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 6, width = 10)

###-------------------####
# Figure 3F: CD4 Proportion Boxplot
###-------------------####
Prop_long <- CD4.prop %>% pivot_longer(cols = -c(HATIMID, Sex, Glucose, HIV, Age, BMI, StudyGroup, meta_fbg, meta_hba1c, homa2_ir),
            names_to = "measure",
            values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "CD4 T Cells", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 5, legend.pos = "none")
                    
ggsave(file = paste(fig_dir, "Figure3F.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 5, width = 10)

###-------------------####
# Figure 3G: Continuous Glucose Measures CD4 T cells
###-------------------####
Cor_res <- Cor_rho(CD4.prop, Measures = c("meta_hba1c", "meta_fbg"))

plot <- Heatmap_circle(df = Cor_res, order = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Cytotoxic", "CD4 Regulatory"), legend_breaks = c(-0.4, 0, 0.4), p_breaks = c(0.01, 0.05, 0.1, 0.5))
ggsave(file = paste(fig_dir, "Figure3G.png", sep = "/"), device = "png", dpi = 300, units = "in", height = 6, width = 10)

###-------------------####
# Miscellaneous Data For Paper
###-------------------####
# Multiple Comparisons Adjustment
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
Mac_DM <- Hmisc::summaryM(LAM + PVM + IM + `Mo-Mac 1` + `Mo-Mac 2` ~ StudyGroup, data = Macrophage.prop[Macrophage.prop$StudyGroup != "HIV+ prediabetic",], overall = TRUE, test = TRUE)
Mac_preDM <- Hmisc::summaryM(LAM + PVM + IM + `Mo-Mac 1` + `Mo-Mac 2` ~ StudyGroup, data = Macrophage.prop[Macrophage.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Results_Mac <- FDR_Fun(df1 = Mac_DM, df2 = Mac_preDM, clusters = c("LAM", "IM", "PVM", "Mo-Mac 1", "Mo-Mac 2"))
write.csv(Results_Mac, file = paste(Prop_dir, "Macrophage.Proportion.Adjusted.csv", sep = "/"))

# CD4
CD4_DM <- Hmisc::summaryM(`CD4 Naive` + `CD4 TCM` + `CD4 TEM` + `CD4 Cytotoxic` + `CD4 Regulatory` ~ StudyGroup, data = CD4.prop[CD4.prop$StudyGroup != "HIV+ prediabetic",], overall = TRUE, test = TRUE)
CD4_preDM <- Hmisc::summaryM(`CD4 Naive` + `CD4 TCM` + `CD4 TEM` + `CD4 Cytotoxic` + `CD4 Regulatory` ~ StudyGroup, data = CD4.prop[CD4.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Results_CD4 <- FDR_Fun(df1 = CD4_DM, df2 = CD4_preDM, clusters = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Cytotoxic", "CD4 Regulatory"))
write.csv(Results_CD4, file = paste(Prop_dir, "CD4.Proportion.Adjusted.csv", sep = "/"))

# CD8
CD8_DM <- Hmisc::summaryM(`CD8 Naive` + `CD8 TCM` + `CD8 TEM` + `CD8 Cytotoxic` ~ StudyGroup, data = CD8.prop[CD8.prop$StudyGroup != "HIV+ prediabetic",], overall = TRUE, test = TRUE)
CD8_preDM <- Hmisc::summaryM(`CD8 Naive` + `CD8 TCM` + `CD8 TEM` + `CD8 Cytotoxic` ~ StudyGroup, data = CD8.prop[CD8.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Results_CD8 <- FDR_Fun(df1 = CD8_DM, df2 = CD8_preDM, clusters = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Cytotoxic"))
write.csv(Results_CD8, file = paste(Prop_dir, "CD8.Proportion.Adjusted.csv", sep = "/"))

# Myeloid
Myeloid_DM <- Hmisc::summaryM(cMo + nMo + `ISG+ Mo` + `Other Mo` + `Mo-Mac 1` + LAM + IM + PVM + cDC1 + cDC2B + pDC + `Migratory DC` + `Cycling Myeloid` + `Mo-Mac 2`  ~ StudyGroup, data = Myeloid.prop[Myeloid.prop$StudyGroup != "HIV+ prediabetic",], overall = TRUE, test = TRUE)
Myeloid_preDM <- Hmisc::summaryM(cMo + nMo + `ISG+ Mo` + `Other Mo` + `Mo-Mac 1` + LAM + IM + PVM + cDC1 + cDC2B + pDC + `Migratory DC` + `Cycling Myeloid` + `Mo-Mac 2`  ~ StudyGroup, data = Myeloid.prop[Myeloid.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Results_Myeloid <- FDR_Fun(df1 = Myeloid_DM, df2 = Myeloid_preDM, clusters = c("cMo", "nMo", "ISG+ Mo", "Other Mo", "Mo-Mac 1", "LAM", "PVM", "IM", "cDC1", "cDC2B", "pDC", "Migratory DC", "Cycling Myeloid", "Mo-Mac 2"))
write.csv(Results_Myeloid, file = paste(Prop_dir, "Myeloid.Proportion.Adjusted.csv", sep = "/"))

# Stromal
Stromal_DM <- Hmisc::summaryM(`PCOLCE+ Fibroblast` + `MYOC+ Fibroblast` + `Adipose Progenitor Cell 1` + `Adipose Progenitor Cell 2` + `Early Preadipocyte` + `ECM-Producing Early Preadipocyte` + `Mature Preadipocyte 1` +
`Mature Preadipocyte 2` + `Myofibroblast` + `Metallothionein+ Preadipocyte` + `Cycling Myofibroblast` ~ StudyGroup, data = Stromal.prop[Stromal.prop$StudyGroup != "HIV+ prediabetic",], overall = TRUE, test = TRUE)
Stromal_preDM <- Hmisc::summaryM(`PCOLCE+ Fibroblast` + `MYOC+ Fibroblast` + `Adipose Progenitor Cell 1` + `Adipose Progenitor Cell 2` + `Early Preadipocyte` + `ECM-Producing Early Preadipocyte` + `Mature Preadipocyte 1` +
`Mature Preadipocyte 2` + `Myofibroblast` + `Metallothionein+ Preadipocyte` + `Cycling Myofibroblast` ~ StudyGroup, data = Stromal.prop[Stromal.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Results_Stromal <- FDR_Fun(df1 = Stromal_DM, df2 = Stromal_preDM, clusters = as.character(colnames(Stromal.prop)[c(1:11)]))
write.csv(Results_Stromal, file = paste(Prop_dir, "Stromal.Proportion.Adjusted.csv", sep = "/"))

# Lymphoid
Lymphoid_DM <- Hmisc::summaryM(`CD4 Naive` + `CD4 TCM` + `CD4 TEM` + `CD4 Regulatory` + `CD4 Cytotoxic` + `CD8 Naive` + `CD8 TCM` + `CD8 TEM` + `CD8 Cytotoxic` +
`CD57 mNK` + `CD16 mNK` + `Immature NK` + ILC + `Cycling T & NK` + MAIT + `Gamma Delta` ~ StudyGroup, data = Lymphoid.prop[Lymphoid.prop$StudyGroup != "HIV+ prediabetic",], overall = TRUE, test = TRUE)
Lymphoid_preDM <- Hmisc::summaryM(`CD4 Naive` + `CD4 TCM` + `CD4 TEM` + `CD4 Regulatory` + `CD4 Cytotoxic` + `CD8 Naive` + `CD8 TCM` + `CD8 TEM` + `CD8 Cytotoxic` +
`CD57 mNK` + `CD16 mNK` + `Immature NK` + ILC + `Cycling T & NK` + MAIT + `Gamma Delta` ~ StudyGroup, data = Lymphoid.prop[Lymphoid.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Results_Lymphoid <- FDR_Fun(df1 = Lymphoid_DM, df2 = Lymphoid_preDM, clusters = as.character(colnames(Lymphoid.prop)[c(1:16)]))
write.csv(Results_Lymphoid, file = paste(Prop_dir, "Lymphoid.Proportion.Adjusted.csv", sep = "/"))

# Vascular
Vascular_DM <- Hmisc::summaryM(`Capillary EC` + `Venous EC` + `VSMC 1` + `Cycling Vascular` +  `Pericyte` + `VSMC 2` + `Arterial EC` ~ StudyGroup, data = Vascular.prop[Vascular.prop$StudyGroup != "HIV+ prediabetic",], overall = TRUE, test = TRUE)
Vascular_preDM <- Hmisc::summaryM(`Capillary EC` + `Venous EC` + `VSMC 1` + `Cycling Vascular` +  `Pericyte` + `VSMC 2` + `Arterial EC` ~ StudyGroup, data = Vascular.prop[Vascular.prop$StudyGroup != "HIV+ diabetic",], overall = TRUE, test = TRUE)
Results_Vascular <- FDR_Fun(df1 = Vascular_DM, df2 = Vascular_preDM, clusters = as.character(colnames(Vascular.prop)[c(1:7)]))
write.csv(Results_Vascular, file = paste(Prop_dir, "Vascular.Proportion.Adjusted.csv", sep = "/"))


# Correlations
prop_list <- list(Myeloid.prop, Macrophage.prop, CD4.prop, CD8.prop, Lymphoid.prop)

cor_result_padj <- function(df, adjusted = T, variables = c("meta_fbg", "meta_hba1c"), nonDM = T) {
    Results <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(Results) <- c("rho", "pvalue", "Cluster", "Variable")
    numb <- dim(df)[2]-10
    if (nonDM == T) {
        df <- df %>% filter(StudyGroup != "HIV+ diabetic")
    }
    
    for(variable in variables){
    if (adjusted == T) {
        if (variable == "BMI") {
            for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,variable] ~ Sex  + StudyGroup + Age, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i], Variable = variable)
                Results <- rbind(Results, res)
            }
        } else if (variable == "Age") {
             for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,variable] ~ Sex  + StudyGroup + BMI, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i], Variable = variable)
                Results <- rbind(Results, res)
            }
        } else if (variable %in% c("homa2_ir", "meta_hba1c", "meta_fbg")) {
            for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,variable] ~ Sex + BMI + Age, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i], Variable = variable)
                Results <- rbind(Results, res)
            }
        } 
    } else {
        for (i in 1:numb) {
            cor_res <- cor.test(df[,i], df[,variable], method = "spearman")
            res <- data.frame(rho = cor_res$estimate, pvalue = cor_res$p.value, colnames(df)[i], Variable = variable)
            Results <- rbind(Results, res)
        }
    }
    }
    Results$padj = p.adjust(Results$pvalue, method = "BH")  
    return(Results)
}

Cor_list = lapply(prop_list, cor_result_padj)
names(Cor_list) <- c("Myeloid", "Macrophage", "CD4 T Cells", "CD8 T Cells", "Lymphoid")


