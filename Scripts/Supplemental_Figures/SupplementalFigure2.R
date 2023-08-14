###############################################################################
######------------------SUPPLEMENTAL FIGURE 2-----------------------------#####
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
library(paletteer) 

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

Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage.rds"))

###-------------------####
# Supplmental Figure 2A: Myeloid ADT
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

ADT <- c("CD11C", "CD1C", "CD14", "CD16")

Plot_list <- lapply(ADT, Feature_Plot, seurat_object = Myeloid)

png(file = paste(fig_dir, "SF2A.png", sep = "/"), units = 'in', height = 10, width = 10, res = 300)
marrangeGrob(Plot_list, nrow=2, ncol=2)
dev.off()

###-------------------####
# Supplmental Figure 2B: Macrophage UMAP
###-------------------####
Idents(Macrophage) <- "Annotation"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("Mo-Mac 1"))) <- 1
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("Mo-Mac 2"))) <- 2
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("PVM"))) <- 3
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("IM"))) <- 4
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c("LAM"))) <- 5
Macrophage[['Celltype_plot']] <- Idents(Macrophage)

col_shuffle <- c("coral", "dodgerblue", "darkolivegreen", "gold", "blue")

UMAP_FUN(Macrophage, path = paste(fig_dir, "SF2B.png", sep = "/"), color = col_shuffle, plot_label = F, breaks = c(1:5), height = 5, width = 7.5,
label = c("Mo-Mac 1 (1)", "Mo-Mac 2 (2)", "PVM (3)", "IM (4)", "LAM (5)"))

###-------------------####
# Supplmental Figure 2C: Macrophage Pathway Analysis
###-------------------####
Mac_markers <- read.csv(paste(markers_dir, "Macrophage.csv", sep = "/"))
Mac.KEGG.Plot <- KEGG_ORA(Mac_markers, celltype = c("IM", "LAM", "Mo-Mac 1", "Mo-Mac 2", "PVM"), filename = paste(fig_dir, "SF2C.png", sep = "/"))

###-------------------####
# Supplmental Figure 2D: Macrophage ADT
###-------------------####
# Load CITE-seq Macrophage-specific data
Subset_dir <- "../Integrated_Filtered/"
Macrophage2 <- readRDS(paste(Subset_dir, "Macrophage.rds", sep = "/"))

# Save model
Macrophage <- RunUMAP(Macrophage, reduction = "harmony", dims = 1:20, return.model = T)

anchors <- FindTransferAnchors(
  reference = Macrophage,
  query = Macrophage2,
  reference.reduction = "pca",
  dims = 1:20
)

Macrophage2 <- MapQuery(
  anchorset = anchors,
  query = Macrophage2,
  reference = Macrophage,
  refdata = list(
    celltype.l1 = "Annotation"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

p1 <- DimPlot(Macrophage2, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 <- DimPlot(Macrophage, label = T) + NoLegend()
p3 <- DimPlot(Macrophage2, reduction = "ref.umap", group.by = "Annotation", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p1 + p2
p1 + p3

Idents(Macrophage2) <- "Visit"
Mac_ADT <- subset(Macrophage2, idents = "Second")

ADT <- c("CD9", "CD206", "CD163", "CD64")

Plot_list <- lapply(ADT, Feature_Plot, seurat_object = Mac_ADT, reduction = 'ref.umap')

png(file = paste(fig_dir, "SF2D.png", sep = "/"), units = 'in', height = 10, width = 15, res = 300)
marrangeGrob(Plot_list, nrow=2, ncol=2)
dev.off()

###-------------------####
# Supplemental Figure 2E: Vijay
###-------------------####
comp_dir = "../Comparisons"

Vijay <- read.csv(paste(comp_dir, "Vijay.csv", sep = "/"))
Emont <- read.csv(paste(comp_dir, "Emont.csv", sep = "/"))
Hildreth <- read.csv(paste(comp_dir, "Hildreth.csv", sep = "/"))
Dick <- read.csv(paste(comp_dir, "Dick.csv", sep = "/"))

###############################
# Modify Gene Names
###############################
Vijay <- Vijay %>% mutate_all(str_to_upper)
Emont <- Emont %>% mutate_all(str_to_upper)
Hildreth <- Hildreth %>% mutate_all(str_to_upper)
Dick <- Dick %>% mutate_all(str_to_upper)

###############################
# Generate Gene Lists for Module Scoring
###############################
GeneList <- list(as.character(Vijay[,1]), as.character(Vijay[,2]), as.character(Vijay[,3]), as.character(Vijay[,4]), as.character(Emont[,1]), as.character(Emont[,2]),
as.character(Emont[,3]), as.character(Hildreth[,1]), as.character(Hildreth[,2]), as.character(Hildreth[,3]), as.character(Dick[,1]), as.character(Dick[,2]), as.character(Dick[,3]))

names(GeneList) <- c("IS2", "IS3", "IS9", "IS12", "hMac3", "hMac1", "hMac2", "IM", "LAM", "PVM", "TLF", "CCR2",
"MHCII")

Macrophage <- AddModuleScore(Macrophage, features = GeneList, name = "Comparison")

colnames(Macrophage@meta.data)[49:61] <- names(GeneList)

# Scale the module scores
for (i in 49:61) {Macrophage@meta.data[,i] <- scale(Macrophage@meta.data[,i])}

Macrophage$Vijay <- "None"
Macrophage@meta.data[Macrophage@meta.data$IS2 > Macrophage@meta.data$IS3 & Macrophage@meta.data$IS2 > Macrophage@meta.data$IS9 & Macrophage@meta.data$IS2 > Macrophage@meta.data$IS12, "Vijay"] <- "IS2"
Macrophage@meta.data[Macrophage@meta.data$IS3 > Macrophage@meta.data$IS2 & Macrophage@meta.data$IS3 > Macrophage@meta.data$IS9 & Macrophage@meta.data$IS3 > Macrophage@meta.data$IS12, "Vijay"] <- "IS3"
Macrophage@meta.data[Macrophage@meta.data$IS9 > Macrophage@meta.data$IS3 & Macrophage@meta.data$IS9 > Macrophage@meta.data$IS2 & Macrophage@meta.data$IS9 > Macrophage@meta.data$IS12, "Vijay"] <- "IS9"
Macrophage@meta.data[Macrophage@meta.data$IS12 > Macrophage@meta.data$IS3 & Macrophage@meta.data$IS12 > Macrophage@meta.data$IS9 & Macrophage@meta.data$IS12 > Macrophage@meta.data$IS2, "Vijay"] <- "IS12"

color_list <- c("#F6EDBDFF","#B4C8A8FF", "#DE8A5AFF", "#70A494FF")
DimPlot(Macrophage, group.by = "Vijay", label = T, cols = color_list) + NoLegend() +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Vijay et. al.")
ggsave(filename= paste(fig_dir, "SF2E.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')          
          
###-------------------####
# Supplemental Figure 2F: Hildreth et al.
###-------------------####
Macrophage$Hildreth <- "None"
Macrophage@meta.data[Macrophage@meta.data$IM > Macrophage@meta.data$PVM & Macrophage@meta.data$IM > Macrophage@meta.data$LAM, "Hildreth"] <- "IM"
Macrophage@meta.data[Macrophage@meta.data$PVM > Macrophage@meta.data$IM & Macrophage@meta.data$PVM > Macrophage@meta.data$LAM, "Hildreth"] <- "PVM"
Macrophage@meta.data[Macrophage@meta.data$LAM > Macrophage@meta.data$IM & Macrophage@meta.data$LAM > Macrophage@meta.data$PVM, "Hildreth"] <- "LAM"

color_list <- c("#F6EDBDFF", "#DE8A5AFF", "#70A494FF")
DimPlot(Macrophage, group.by = "Hildreth", label = T, cols = color_list) + NoLegend() +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Hildreth et. al.")
ggsave(filename= paste(fig_dir, "SF2F.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')          

###-------------------####
# Supplemental Figure 2G: Emont et. al.
###-------------------####
Macrophage$Emont <- "None"
Macrophage@meta.data[Macrophage@meta.data$hMac1 > Macrophage@meta.data$hMac2 & Macrophage@meta.data$hMac1 > Macrophage@meta.data$hMac3, "Emont"] <- "hMac1"
Macrophage@meta.data[Macrophage@meta.data$hMac2 > Macrophage@meta.data$hMac1 & Macrophage@meta.data$hMac2 > Macrophage@meta.data$hMac3, "Emont"] <- "hMac2"
Macrophage@meta.data[Macrophage@meta.data$hMac3 > Macrophage@meta.data$hMac1 & Macrophage@meta.data$hMac3 > Macrophage@meta.data$hMac2, "Emont"] <- "hMac3"

DimPlot(Macrophage, group.by = "Emont", label = T, cols = color_list) + NoLegend() +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Emont et. al.")
ggsave(filename= paste(fig_dir, "SF2G.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')          

###-------------------####
# Supplemental Figure 2H: Dick et. al.
###-------------------####
Macrophage$Emont <- "None"
Macrophage@meta.data[Macrophage@meta.data$TLF > Macrophage@meta.data$CCR2 & Macrophage@meta.data$TLF > Macrophage@meta.data$MHCII, "Dick"] <- "TLF"
Macrophage@meta.data[Macrophage@meta.data$CCR2 > Macrophage@meta.data$TLF & Macrophage@meta.data$CCR2 > Macrophage@meta.data$MHCII, "Dick"] <- "CCR2"
Macrophage@meta.data[Macrophage@meta.data$MHCII > Macrophage@meta.data$TLF & Macrophage@meta.data$MHCII > Macrophage@meta.data$CCR2, "Dick"] <- "MHCII"

DimPlot(Macrophage, group.by = "Dick", label = T, cols = color_list) + NoLegend() +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("Dick et. al.")
ggsave(filename= paste(fig_dir, "SF2H.png", sep = "/"), dpi = 300, units = "in", height = 5, width = 5, device = 'png')          
          