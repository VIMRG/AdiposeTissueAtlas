###############################################################################
######------------------SUPPLEMENTAL FIGURE 3-----------------------------#####
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

###-------------------####
# Source Code
###-------------------####
source('../Utils.R')

###-------------------####
# Load Data
###-------------------####
tmp_dir <- "../SubsetAnalysis"
fig_dir <- "../Figures"
markers_dir <- "../Markers"

Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8.rds"))

###-------------------####
# Supplmental Figure 3A: Lymphoid ADT
###-------------------####
Feature_Plot <- function(seurat_object, ADT) {
    DefaultAssay(seurat_object) <- "ADT"
    Plot <- FeaturePlot(seurat_object, features = ADT, raster = F, min.cutoff = 'q10', max.cutoff = 'q90', cols = c("lightgrey", "darkgreen")) +
          theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14),
          plot.title = element_text(face = "bold", size = 20)) + ggtitle(ADT)
    return(Plot)
}

ADT = c("CD4", "CD8", "CD27", "CD45RA", "CD57", "CD16")

Plot_list <- lapply(ADT, Feature_Plot, seurat_object = Lymphoid)

png(file = paste(fig_dir, "SF3A.png", sep = "/"), units = 'in', height = 5, width = 8.5, res = 300)
marrangeGrob(Plot_list, nrow=2, ncol=3)
dev.off()

###-------------------####
# Supplmental Figure 3B: CD4 UMAP
###-------------------####
Idents(CD4) <- "Annotation"
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Naive"))) <- 1
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 TCM"))) <- 2
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 TEM"))) <- 3
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Regulatory"))) <- 4
Idents(CD4, cells = WhichCells(CD4, idents = c("CD4 Cytotoxic"))) <- 5
CD4[['Celltype_plot']] <- Idents(CD4)

col_shuffle <- c("coral", "dodgerblue", "darkolivegreen", "gold", "magenta")

UMAP_FUN(CD4, path = paste(fig_dir, "SF3B.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:5), height = 5, width = 8,
label = c("CD4 Naive (1)", "CD4 TCM (2)", "CD4 TEM (3)", "CD4 Regulatory (4)", "CD4 Cytotoxic (5)"))

###-------------------####
# Supplmental Figure 3C: CD4 ADT
###-------------------####
ADT = c("CD27", "CD45RA", "CD57", "CD69")

Plot_list <- lapply(ADT, Feature_Plot, seurat_object = CD4)

png(file = paste(fig_dir, "SF3C.png", sep = "/"), units = "in", height = 8, width = 10, res = 300)
marrangeGrob(Plot_list, nrow=2, ncol=2)
dev.off()

###-------------------####
# Supplmental Figure 3D: CD8 UMAP
###-------------------####
Idents(CD8) <- "Annotation"
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 Naive"))) <- 1
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 TCM"))) <- 2
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 TEM"))) <- 3
Idents(CD8, cells = WhichCells(CD8, idents = c("CD8 Cytotoxic"))) <- 4
CD8[['Celltype_plot']] <- Idents(CD8)

col_shuffle <- c("coral", "dodgerblue", "darkolivegreen", "gold")

UMAP_FUN(CD8, path = paste(fig_dir, "SF3D.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:4), height = 5, width = 8,
label = c("CD8 Naive (1)", "CD8 TCM (2)", "CD8 TEM (3)", "CD8 Cytotoxic (4)"))

###-------------------####
# Supplemental Figure 3E: CD8 ADT
###-------------------####
ADT = c("CD27", "CD45RA", "CD57", "CD69")

Plot_list <- lapply(ADT, Feature_Plot, seurat_object = CD8)

png(file = paste(fig_dir, "SF3E.png", sep = "/"), units = "in", height = 8, width = 10, res = 300)
marrangeGrob(Plot_list, nrow=2, ncol=2)
dev.off()

###-------------------####
# Supplemental Figure 3F: Correlation scRNA-seq and Flow
###-------------------####
Prop_dir <- "../Proportion"
CD4.prop <- read.table(paste(Prop_dir, "CD4.Prop.txt", sep = "/"), check.names = FALSE)
CD4.Num <- read.table(paste(Prop_dir, "CD4.Num.txt", sep = "/"), check.names = FALSE)
colnames(CD4.Num) <- c("HATIMID", "Freq")
CD4.Num$HATIMID <- as.factor(CD4.Num$HATIMID)

CD4.prop <- CD4.prop %>% na.omit()
CD4.prop <- Prop_Merge(CD4.prop, data_hatim)
CD4.prop <- left_join(CD4.prop, data_hatim[, c("HATIMID", "at_cd4_cd69_mem1")], by = "HATIMID")
CD4.prop <- left_join(CD4.prop, CD4.Num, by = "HATIMID")
CD4.prop <- CD4.prop %>% dplyr::filter(!is.na(at_cd4_cd69_mem1) & Freq > 30)

plot <- scatter_plot(DF = CD4.prop, x = "`CD4 TEM`", y = "at_cd4_cd69_mem1", xlabel = "scRNA-seq CD4 TEM Proportion", ylabel = "Flow CD4 CD69+ (as % CD4)", dotsize = 4)
ggsave(filename = paste(fig_dir, "SF3F.png", sep = "/"), units = "in", height = 7, width = 7, dpi = 300, device = "png")

###-------------------####
# Supplemental Figure 3G: Comparisons
###-------------------####
comp_dir = "..Comparisons/"
Stromal_Comparisons <- read.csv(paste(comp_dir, "Stromal_Comparisons.csv", sep = "/"))

###############################
# Gather list and Remove genes that are blank or NA (don't have human orthologous gene)
###############################
GeneList <- list(as.character(Stromal_Comparisons[,1]), as.character(Stromal_Comparisons[,2]), as.character(Stromal_Comparisons[,3]), as.character(Stromal_Comparisons[,4]), as.character(Stromal_Comparisons[,5]), as.character(Stromal_Comparisons[,6]),
as.character(Stromal_Comparisons[,7]), as.character(Stromal_Comparisons[,8]), as.character(Stromal_Comparisons[,9]), as.character(Stromal_Comparisons[,10]), as.character(Stromal_Comparisons[,11]), as.character(Stromal_Comparisons[,12]), as.character(Stromal_Comparisons[,13]),
as.character(Stromal_Comparisons[,14]), as.character(Stromal_Comparisons[,15]), as.character(Stromal_Comparisons[,16]), as.character(Stromal_Comparisons[,17]), as.character(Stromal_Comparisons[,18]), as.character(Stromal_Comparisons[,19]), as.character(Stromal_Comparisons[,20]),
as.character(Stromal_Comparisons[,21]))
   
GeneList <- lapply(GeneList, function(x) {
    x[x != "N/A" & x != ""]
})

names(GeneList) <- gsub(".*_", "", colnames(Stromal_Comparisons))

Stromal <- AddModuleScore(Stromal, features = GeneList, name = "Comparison")
colnames(Stromal@meta.data)[49:69] <- names(GeneList)

# Scale the module scores
for (i in 49:69) {Stromal@meta.data[,i] <- scale(Stromal@meta.data[,i])}

###-------------------####
# Schwalie Comparisons
###-------------------####
Stromal$Schwalie <- "None"
Stromal@meta.data[Stromal@meta.data$P1 > Stromal@meta.data$P2 & Stromal@meta.data$P1 > Stromal@meta.data$P3, "Schwalie"] <- "P1"
Stromal@meta.data[Stromal@meta.data$P2 > Stromal@meta.data$P1 & Stromal@meta.data$P2 > Stromal@meta.data$P3, "Schwalie"] <- "P2"
Stromal@meta.data[Stromal@meta.data$P3 > Stromal@meta.data$P1 & Stromal@meta.data$P3 > Stromal@meta.data$P2, "Schwalie"] <- "P3"

color_list <- c("#F3A312","#70A494FF", "#12A5F3")
DimPlot(Stromal, group.by = "Schwalie", label = F, cols = color_list) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Schwalie et. al.")
ggsave(filename= paste(fig_dir, "SF3G_Schwalie.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')          

###-------------------####
# Merrick Comparisons
###-------------------####
Stromal$Merrick <- "None"
Stromal@meta.data[Stromal@meta.data$CD142 > Stromal@meta.data$DPP4 & Stromal@meta.data$CD142 > Stromal@meta.data$ICAM, "Merrick"] <- "CD142"
Stromal@meta.data[Stromal@meta.data$DPP4 > Stromal@meta.data$CD142 & Stromal@meta.data$DPP4 > Stromal@meta.data$ICAM, "Merrick"] <- "DPP4"
Stromal@meta.data[Stromal@meta.data$ICAM > Stromal@meta.data$CD142 & Stromal@meta.data$ICAM > Stromal@meta.data$DPP4, "Merrick"] <- "ICAM"

color_list <- c("#F3A312","#70A494FF", "#12A5F3")
DimPlot(Stromal, group.by = "Merrick", label = F, cols = color_list) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Merrick et. al.")
ggsave(filename= paste(fig_dir, "SF3G_Merrick.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')          

###-------------------####
# Hepler Comparisons
###-------------------####
Stromal$Hepler <- "None"
Stromal@meta.data[Stromal@meta.data$APC > Stromal@meta.data$CAP & Stromal@meta.data$APC > Stromal@meta.data$FIP, "Hepler"] <- "APC"
Stromal@meta.data[Stromal@meta.data$CAP > Stromal@meta.data$APC & Stromal@meta.data$CAP > Stromal@meta.data$FIP, "Hepler"] <- "CAP"
Stromal@meta.data[Stromal@meta.data$FIP > Stromal@meta.data$APC & Stromal@meta.data$FIP > Stromal@meta.data$CAP, "Hepler"] <- "FIP"

color_list <- c("#F3A312","#70A494FF", "#12A5F3")
DimPlot(Stromal, group.by = "Hepler", label = F, cols = color_list) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Hepler et. al.")
ggsave(filename= paste(fig_dir, "SF3G_Hepler.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')          


###-------------------####
# Sarvari et al.
###-------------------####
Stromal$Sarvari <- "None"
Stromal@meta.data[Stromal@meta.data$FAP1 > Stromal@meta.data$FAP2 & Stromal@meta.data$FAP1 > Stromal@meta.data$FAP3 & Stromal@meta.data$FAP1 > Stromal@meta.data$FAP4, "Sarvari"] <- "FAP1"
Stromal@meta.data[Stromal@meta.data$FAP2 > Stromal@meta.data$FAP1 & Stromal@meta.data$FAP2 > Stromal@meta.data$FAP3 & Stromal@meta.data$FAP2 > Stromal@meta.data$FAP4, "Sarvari"] <- "FAP2"
Stromal@meta.data[Stromal@meta.data$FAP3 > Stromal@meta.data$FAP1 & Stromal@meta.data$FAP3 > Stromal@meta.data$FAP2 & Stromal@meta.data$FAP3 > Stromal@meta.data$FAP4, "Sarvari"] <- "FAP3"
Stromal@meta.data[Stromal@meta.data$FAP4 > Stromal@meta.data$FAP1 & Stromal@meta.data$FAP4 > Stromal@meta.data$FAP2 & Stromal@meta.data$FAP4 > Stromal@meta.data$FAP3, "Sarvari"] <- "FAP4"

color_list <- c("#F3A312","#70A494FF", "#12A5F3", "#F54CE3")
DimPlot(Stromal, group.by = "Sarvari", label = F, cols = color_list) + 
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Sarvari et. al.")
ggsave(filename= paste(fig_dir, "SF3G_Sarvari.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')          

###-------------------####
# Emont et al.
###-------------------####
Stromal$Emont <- "None"
Stromal@meta.data[Stromal@meta.data$hASPC1 > Stromal@meta.data$hASPC2 & Stromal@meta.data$hASPC1 > Stromal@meta.data$hASPC3 & Stromal@meta.data$hASPC1 > Stromal@meta.data$hASPC4 & Stromal@meta.data$hASPC1 > Stromal@meta.data$hASPC5 & Stromal@meta.data$hASPC1 > Stromal@meta.data$hASPC6 , "Emont"] <- "hASPC1"
Stromal@meta.data[Stromal@meta.data$hASPC2 > Stromal@meta.data$hASPC1 & Stromal@meta.data$hASPC2 > Stromal@meta.data$hASPC3 & Stromal@meta.data$hASPC2 > Stromal@meta.data$hASPC4 & Stromal@meta.data$hASPC2 > Stromal@meta.data$hASPC5 & Stromal@meta.data$hASPC2 > Stromal@meta.data$hASPC6 , "Emont"] <- "hASPC2"
Stromal@meta.data[Stromal@meta.data$hASPC3 > Stromal@meta.data$hASPC1 & Stromal@meta.data$hASPC3 > Stromal@meta.data$hASPC2 & Stromal@meta.data$hASPC3 > Stromal@meta.data$hASPC4 & Stromal@meta.data$hASPC3 > Stromal@meta.data$hASPC5 & Stromal@meta.data$hASPC3 > Stromal@meta.data$hASPC6 , "Emont"] <- "hASPC3"
Stromal@meta.data[Stromal@meta.data$hASPC4 > Stromal@meta.data$hASPC1 & Stromal@meta.data$hASPC4 > Stromal@meta.data$hASPC2 & Stromal@meta.data$hASPC4 > Stromal@meta.data$hASPC3 & Stromal@meta.data$hASPC4 > Stromal@meta.data$hASPC5 & Stromal@meta.data$hASPC4 > Stromal@meta.data$hASPC6 , "Emont"] <- "hASPC4"
Stromal@meta.data[Stromal@meta.data$hASPC5 > Stromal@meta.data$hASPC1 & Stromal@meta.data$hASPC5 > Stromal@meta.data$hASPC2 & Stromal@meta.data$hASPC5 > Stromal@meta.data$hASPC3 & Stromal@meta.data$hASPC5 > Stromal@meta.data$hASPC4 & Stromal@meta.data$hASPC5 > Stromal@meta.data$hASPC6 , "Emont"] <- "hASPC5"
Stromal@meta.data[Stromal@meta.data$hASPC6 > Stromal@meta.data$hASPC1 & Stromal@meta.data$hASPC6 > Stromal@meta.data$hASPC2 & Stromal@meta.data$hASPC6 > Stromal@meta.data$hASPC3 & Stromal@meta.data$hASPC6 > Stromal@meta.data$hASPC4 & Stromal@meta.data$hASPC6 > Stromal@meta.data$hASPC5 , "Emont"] <- "hASPC6"

color_list <- c("#F3A312","#70A494FF", "#12A5F3", "#F54CE3", "#4C8920", "#DCC928")
DimPlot(Stromal, group.by = "Emont", label = F, cols = color_list) + 
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Emont et. al.")
ggsave(filename= paste(fig_dir, "SF3G_Emont.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')      
