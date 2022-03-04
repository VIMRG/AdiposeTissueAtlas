##############################################################################################
##------------------------- SINGLE LANE PROCESSING ---------------------------------##
##------------------------- DATE: 8/18/2021 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: The following code will assess each lane independently after applying QC filtering
# and removing homotypic doublets via genetic demultiplexing. DoubletFinder results will be loaded
# into a common metadata slot called "DoubletFinder." It is difficult to accurately annotate true
# doublets and remove at this step so instead, the following code will grossly annotate each dataset
# including possible transcriptional doublets and then remove after integration of the datasets so as
# to preserve any small cell populations that may otherwise be thrown out.
##############################################################################################

###-------------------------------####
# Environment
###------------------------------####
set.seed(7612) # random seed
setwd('/data/p_koethe_lab/Atlas_AT/Analysis') # working directory

###-------------------####
# Load Libraries
###-------------------####
library(Seurat) # Seurat 4.0 version(note that several packages will also be attached)
library(tidyverse) # collection of packages for manipulating data (ggplot2, dplyr, tidyr, purr, stringr, tibble)
library(patchwork) # tool to assist with plotting
library(Matrix) # dense and sparse matrix manipulations
library(cowplot) # tool to assist with plotting
library(SeuratDisk) # Stores objects for reference assignments
library(stringr) # grab strings
library(data.table) # convert table to dataframe
library(mgsub) # multiple substitutions
library(DoubletFinder) # find heterotypic doublets
library(stringr) # find strings
library(MAST) # DGE

###-------------------####
# Load Data
###-------------------####

# Read in files
path <- "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_objects/"

# Pull all files in path
file_list <- list.files(path = path, pattern = ".rds")

names(file_list) <- c('P5344_CW1', 'P5344_CW2', 'P5344_CW3', 'P5544_CW1', 'P5544_CW2', 'P5544_CW3', 'P5573_CW1',
                      'P5573_CW2', 'P5573_CW3', 'P5657_CW1', 'P5657_CW3', 'P5836_CW1', 'P5836_CW2',
                      'P5836_CW3', 'P5877_CW1', 'P5877_CW2', 'P5877_CW3', 'P5903_CW1', 'P5903_CW2', 'P5903_CW3',
                      'P5963_CW1', 'P5963_CW2', 'P5963_CW3')


for (i in 1:length(file_list)){
  print(file_list[[i]])
  print(names(file_list[i]))
  seurat_object <- readRDS(paste(path, file_list[[i]], sep = ""))
  assign(names(file_list[i]), seurat_object)
}

# 1. CREATE SINGLE METADATA SLOT FOR DOUBLETFINDER RESULTS--------------------------------------
# For ease of comparison after integration, will create metadata slot with common name for all lanes
# identifying singlet and doublet populations that resulted from DoubletFinder

###-------------------####
# Add DoubletFinder Metadata
###-------------------####

P5344_CW1[['DoubletFinder']] <- P5344_CW1[["DF.classifications_0.25_0.03_852"]]
P5344_CW2[['DoubletFinder']] <- P5344_CW2[["DF.classifications_0.25_0.08_917"]]
P5344_CW3[['DoubletFinder']] <- P5344_CW3[["DF.classifications_0.25_0.24_882"]]
P5544_CW1[['DoubletFinder']] <- P5544_CW1[["DF.classifications_0.25_0.005_605"]]
P5544_CW2[['DoubletFinder']] <- P5544_CW2[["DF.classifications_0.25_0.19_527"]]
P5544_CW3[['DoubletFinder']] <- P5544_CW3[["DF.classifications_0.25_0.25_770"]]
P5573_CW1[['DoubletFinder']] <- P5573_CW1[["DF.classifications_0.25_0.09_769"]]
P5573_CW2[['DoubletFinder']] <- P5573_CW2[["DF.classifications_0.25_0.23_705"]]
P5573_CW3[['DoubletFinder']] <- P5573_CW3[['DF.classifications_0.25_0.04_565']]
P5657_CW1[['DoubletFinder']] <- P5657_CW1[['DF.classifications_0.25_0.24_755']]
P5657_CW3[['DoubletFinder']] <- P5657_CW3[['DF.classifications_0.25_0.02_1059']]
P5836_CW1[['DoubletFinder']] <- P5836_CW1[['DF.classifications_0.25_0.04_453']]
P5836_CW2[['DoubletFinder']] <- P5836_CW2[['DF.classifications_0.25_0.03_664']]
P5836_CW3[['DoubletFinder']] <- P5836_CW3[['DF.classifications_0.25_0.3_804']]
P5877_CW1[['DoubletFinder']] <- P5877_CW1[['DF.classifications_0.25_0.06_615']]
P5877_CW2[['DoubletFinder']] <- P5877_CW2[['DF.classifications_0.25_0.05_630']]
P5877_CW3[['DoubletFinder']] <- P5877_CW3[['DF.classifications_0.25_0.21_481']]
P5903_CW1[['DoubletFinder']] <- P5903_CW1[['DF.classifications_0.25_0.1_614']]
P5903_CW2[['DoubletFinder']] <- P5903_CW2[['DF.classifications_0.25_0.03_1037']]
P5903_CW3[['DoubletFinder']] <- P5903_CW3[['DF.classifications_0.25_0.13_1232']]
P5963_CW1[['DoubletFinder']] <- P5963_CW1[['DF.classifications_0.25_0.005_665']]
P5963_CW2[['DoubletFinder']] <- P5963_CW2[['DF.classifications_0.25_0.27_1000']]
P5963_CW3[['DoubletFinder']] <- P5963_CW3[['DF.classifications_0.25_0.03_1168']]

# 2. COARSE MANUAL ANNOTATION BEFORE INTEGRATING LANES-------------------------------------

###-------------------------------####
# P5344_CW1 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5344_CW1, label = T) + NoLegend()

p2 <- FeaturePlot(P5344_CW1, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(14, 15, 18))) <- "NK/T Cells"

p2 <- FeaturePlot(P5344_CW1, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(16))) <- "Dendritic Cells"

p2 <- FeaturePlot(P5344_CW1, features = c("CD14", "VCAN", "LYZ", "S100A8"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(6, 21))) <- "Monocytes"

p2 <- FeaturePlot(P5344_CW1, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(24, 10, 4, 29))) <- "Macrophage"

p2 <- FeaturePlot(P5344_CW1, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(5, 13, 3))) <- "Capillary"

p2 <- FeaturePlot(P5344_CW1, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(8, 22, 27))) <- "Venous"

p2 <- FeaturePlot(P5344_CW1, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(20, 9, 2))) <- "Vascular Smooth Muscle"

p2 <- FeaturePlot(P5344_CW1, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(31))) <- "Mature Adipocytes"

p2 <- FeaturePlot(P5344_CW1, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(0, 12, 1, 11, 23, 7, 17, 19))) <- "Stromal cells"

p2 <- FeaturePlot(P5344_CW1, features = c("MS4A1", "CD74", "CD79A", "BANK1"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(25))) <- "B Cells"

p2 <- FeaturePlot(P5344_CW1, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(28, 30))) <- "Transcriptional Doublets"

p2 <- FeaturePlot(P5344_CW1, features = c("HBB"))
p1 + p2
Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(26))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5344_CW1, label = T) + NoLegend()
p2 <- DimPlot(P5344_CW1, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5344_CW1/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5344_CW1[['Gross_Annotation']] <- Idents(P5344_CW1)


saveRDS(P5344_CW1, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5344_CW1_singlet.rds")

###-------------------------------####
# P5344_CW2 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5344_CW2, label = T) + NoLegend()

p2 <- FeaturePlot(P5344_CW2, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(10, 21, 18, 26))) <- "NK/T Cells"

p2 <- FeaturePlot(P5344_CW2, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(27))) <- "Dendritic Cells"

p2 <- FeaturePlot(P5344_CW2, features = c("CD14", "VCAN", "LYZ", "S100A8"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(22, 24))) <- "Monocytes"

p2 <- FeaturePlot(P5344_CW2, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(23, 28, 25, 9, 30))) <- "Macrophage"

p2 <- FeaturePlot(P5344_CW2, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(5, 15))) <- "Venous"

p2 <- FeaturePlot(P5344_CW2, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(16, 11, 6, 20))) <- "Vascular Smooth Muscle"

p2 <- FeaturePlot(P5344_CW2, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(1, 13, 14, 17))) <- "Capillary"

p2 <- FeaturePlot(P5344_CW2, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2
#Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(31))) <- "Mature Adipocytes"

p2 <- FeaturePlot(P5344_CW2, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(0, 3, 19, 12, 8, 7, 4, 2))) <- "Stromal cells"

p2 <- FeaturePlot(P5344_CW2, features = c("MS4A1", "CD74", "CD79A", "BANK1"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(31))) <- "B Cells"

p2 <- FeaturePlot(P5344_CW2, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2
Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(29))) <- "Transcriptional Doublets"

p2 <- FeaturePlot(P5344_CW2, features = c("HBB"))
p1 + p2
#Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(26))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5344_CW2, label = T) + NoLegend()
p2 <- DimPlot(P5344_CW2, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5344_CW2/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5344_CW2[['Gross_Annotation']] <- Idents(P5344_CW2)


saveRDS(P5344_CW2, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5344_CW2_singlet.rds")

###-------------------------------####
# P5344_CW3 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5344_CW3, label = T) + NoLegend()

p2 <- FeaturePlot(P5344_CW3, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(16, 23, 35))) <- "NK/T Cells"

p2 <- FeaturePlot(P5344_CW3, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(25))) <- "Dendritic Cells"

p2 <- FeaturePlot(P5344_CW3, features = c("CD14", "VCAN", "LYZ", "S100A8"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(21, 28))) <- "Monocytes"

p2 <- FeaturePlot(P5344_CW3, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(6, 29, 10, 24, 20))) <- "Macrophage"

p2 <- FeaturePlot(P5344_CW3, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(22, 30))) <- "Venous"

p2 <- FeaturePlot(P5344_CW3, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(13, 27, 1))) <- "Vascular Smooth Muscle"

p2 <- FeaturePlot(P5344_CW3, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(18, 4, 15, 31))) <- "Capillary"

p2 <- FeaturePlot(P5344_CW3, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(32))) <- "Mature Adipocytes"

p2 <- FeaturePlot(P5344_CW3, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(26, 11, 12, 9, 19, 7, 8, 2, 14, 0, 5, 17, 3))) <- "Stromal cells"

p2 <- FeaturePlot(P5344_CW3, features = c("MS4A1", "CD74", "CD79A", "BANK1"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(33, 36))) <- "B Cells"

p2 <- FeaturePlot(P5344_CW3, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(35, 34))) <- "Transcriptional Doublets"

p2 <- FeaturePlot(P5344_CW3, features = c("HBB"))
p1 + p2
#Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(26))) <- "Erythrocytes"

Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(16, 23))) <- "NK/T Cells"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(25))) <- "Dendritic Cells"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(21, 28))) <- "Monocytes"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(6, 29, 10, 24, 20))) <- "Macrophage"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(22, 30))) <- "Venous"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(13, 27, 1))) <- "Vascular Smooth Muscle"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(18, 4, 15, 31))) <- "Capillary"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(32))) <- "Mature Adipocytes"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(26, 11, 12, 9, 19, 7, 8, 2, 14, 0, 5, 17, 3))) <- "Stromal cells"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(33, 36))) <- "B Cells"
Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(35, 34))) <- "Transcriptional Doublets"
#Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(26))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5344_CW3, label = T) + NoLegend()
p2 <- DimPlot(P5344_CW3, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5344_CW3/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5344_CW3[['Gross_Annotation']] <- Idents(P5344_CW3)


saveRDS(P5344_CW3, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5344_CW3_singlet.rds")

###-------------------------------####
# P5544_CW1 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5544_CW1, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5544_CW1, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5544_CW1, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5544_CW1, features = c("CD14", "VCAN", "LYZ", "S100A8"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5544_CW1, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5544_CW1, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5544_CW1, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5544_CW1, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5544_CW1, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5544_CW1, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5544_CW1, features = c("MS4A1", "CD74", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5544_CW1, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5544_CW1, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(6, 17, 9, 25, 21, 3, 22, 13, 8, 27))) <- "NK/T Cells"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(28))) <- "Dendritic Cells"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(14, 19, 7, 4, 18, 12))) <- "Monocytes"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(20, 10, 5))) <- "Macrophage"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(23))) <- "Venous"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(26))) <- "Vascular Smooth Muscle"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(16, 24))) <- "Capillary"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(29))) <- "Mature Adipocytes"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(15, 11, 1, 0))) <- "Stromal cells"
Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(2))) <- "B Cells"
#Idents(P5544_CW1, cells = WhichCells(P5544_CW1, idents = c(35, 34))) <- "Transcriptional Doublets"
#Idents(P5344_CW1, cells = WhichCells(P5344_CW1, idents = c(26))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5544_CW1, label = T) + NoLegend()
p2 <- DimPlot(P5544_CW1, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5544_CW1/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5544_CW1[['Gross_Annotation']] <- Idents(P5544_CW1)


saveRDS(P5544_CW1, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5544_CW1_singlet.rds")

###-------------------------------####
# P5544_CW2 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5544_CW2, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5544_CW2, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5544_CW2, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5544_CW2, features = c("CD14", "VCAN", "LYZ", "S100A8"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5544_CW2, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5544_CW2, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5544_CW2, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5544_CW2, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5544_CW2, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5544_CW2, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5544_CW2, features = c("MS4A1", "CD74", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5544_CW2, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5544_CW2, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(10, 2, 3, 17, 23, 22, 8, 29))) <- "NK/T Cells"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(30, 13))) <- "Dendritic Cells"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(15, 4, 14, 25))) <- "Monocytes"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(6, 5, 12))) <- "Macrophage"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(28))) <- "Venous"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(26, 21))) <- "Vascular Smooth Muscle"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(18, 19, 16))) <- "Capillary"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(31))) <- "Mature Adipocytes"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(9,1, 0, 7, 24, 11, 20))) <- "Stromal cells"
Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(27))) <- "B Cells"
#Idents(P5544_CW2, cells = WhichCells(P5544_CW2, idents = c(35, 34))) <- "Transcriptional Doublets"
#Idents(P5344_CW2, cells = WhichCells(P5344_CW2, idents = c(26))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5544_CW2, label = T) + NoLegend()
p2 <- DimPlot(P5544_CW2, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5544_CW2/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5544_CW2[['Gross_Annotation']] <- Idents(P5544_CW2)


saveRDS(P5544_CW2, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5544_CW2_singlet.rds")

###-------------------------------####
# P5544_CW3 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5544_CW3, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5544_CW3, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5544_CW3, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5544_CW3, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5544_CW3, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5544_CW3, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5544_CW3, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5544_CW3, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5544_CW3, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5544_CW3, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5544_CW3, features = c("MS4A1", "CD74", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5544_CW3, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5544_CW3, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(18, 11, 22, 29))) <- "NK/T Cells"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(28))) <- "Dendritic Cells"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(26, 20))) <- "Monocytes"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(24, 12, 14))) <- "Macrophage"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(15, 30))) <- "Venous"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(23, 10, 6, 5, 21))) <- "Vascular Smooth Muscle"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(16, 7, 19, 8, 27, 17))) <- "Capillary"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(32))) <- "Mature Adipocytes"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(0, 1, 9, 3, 4, 2, 13, 25))) <- "Stromal cells"
Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(31))) <- "B Cells"
#Idents(P5544_CW3, cells = WhichCells(P5544_CW3, idents = c(35, 34))) <- "Transcriptional Doublets"
#Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(26))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5544_CW3, label = T) + NoLegend()
p2 <- DimPlot(P5544_CW3, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5544_CW3/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5544_CW3[['Gross_Annotation']] <- Idents(P5544_CW3)


saveRDS(P5544_CW3, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5544_CW3_singlet.rds")

###-------------------------------####
# P5573_CW1 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5573_CW1, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5573_CW1, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5573_CW1, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5573_CW1, features = c("CD14", "VCAN", "LYZ", "S100A8"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5573_CW1, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5573_CW1, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5573_CW1, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5573_CW1, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5573_CW1, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5573_CW1, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5573_CW1, features = c("MS4A1", "CD74", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5573_CW1, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5573_CW1, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(0, 14, 23, 27, 26, 11))) <- "NK/T Cells"
Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(33, 28))) <- "Dendritic Cells"
Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(10, 7, 15, 29, 19))) <- "Monocytes"
Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(3, 22, 17, 2, 24, 21, 32, 12))) <- "Macrophage"
#Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c())) <- "Venous"
Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(30, 9))) <- "Vascular Smooth Muscle"
Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(25, 20, 6))) <- "Capillary"
Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(34))) <- "Mature Adipocytes"
Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(8, 31, 5, 13, 16, 4, 1))) <- "Stromal cells"
Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(18))) <- "B Cells"
#Idents(P5573_CW1, cells = WhichCells(P5573_CW1, idents = c(35, 34))) <- "Transcriptional Doublets"
#Idents(P5344_CW3, cells = WhichCells(P5344_CW3, idents = c(26))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5573_CW1, label = T) + NoLegend()
p2 <- DimPlot(P5573_CW1, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5573_CW1/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5573_CW1[['Gross_Annotation']] <- Idents(P5573_CW1)


saveRDS(P5573_CW1, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5573_CW1_singlet.rds")

###-------------------------------####
# P5573_CW2 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5573_CW2, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5573_CW2, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5573_CW2, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5573_CW2, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5573_CW2, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5573_CW2, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5573_CW2, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5573_CW2, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5573_CW2, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5573_CW2, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5573_CW2, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5573_CW2, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5573_CW2, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(31, 18, 8, 12, 14, 27, 28, 5, 24, 21))) <- "NK/T Cells"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(38))) <- "Dendritic Cells"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(25, 9, 6, 29))) <- "Monocytes"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(2, 19, 23, 32, 34))) <- "Macrophage"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(11))) <- "Venous"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(16, 3, 30))) <- "Vascular Smooth Muscle"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(0, 26, 36, 15, 13))) <- "Capillary"
#Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c())) <- "Mature Adipocytes"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(10, 1, 7, 4, 22))) <- "Stromal cells"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(17, 20))) <- "B Cells"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(37))) <- "Transcriptional Doublets"
Idents(P5573_CW2, cells = WhichCells(P5573_CW2, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5573_CW2, label = T) + NoLegend()
p2 <- DimPlot(P5573_CW2, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5573_CW2/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5573_CW2[['Gross_Annotation']] <- Idents(P5573_CW2)


saveRDS(P5573_CW2, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5573_CW2_singlet.rds")

###-------------------------------####
# P5573_CW3 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5573_CW3, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5573_CW3, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5573_CW3, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5573_CW3, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5573_CW3, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5573_CW3, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5573_CW3, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5573_CW3, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5573_CW3, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5573_CW3, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5573_CW3, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5573_CW3, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5573_CW3, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(13, 3, 21))) <- "NK/T Cells"
Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(16))) <- "Dendritic Cells"
Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(10, 15))) <- "Monocytes"
Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(9, 25, 5))) <- "Macrophage"
Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(19))) <- "Venous"
Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(24, 20, 14, 23, 7))) <- "Vascular Smooth Muscle"
Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(18, 27, 11, 8))) <- "Capillary"
#Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c())) <- "Mature Adipocytes"
Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(2, 1, 0, 4, 22, 6, 12, 17))) <- "Stromal cells"
Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(26))) <- "B Cells"
#Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(37))) <- "Transcriptional Doublets"
#Idents(P5573_CW3, cells = WhichCells(P5573_CW3, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5573_CW3, label = T) + NoLegend()
p2 <- DimPlot(P5573_CW3, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5573_CW3/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5573_CW3[['Gross_Annotation']] <- Idents(P5573_CW3)


saveRDS(P5573_CW3, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5573_CW3_singlet.rds")

###-------------------------------####
# P5657_CW1 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5657_CW1, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5657_CW1, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5657_CW1, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5657_CW1, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5657_CW1, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5657_CW1, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5657_CW1, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5657_CW1, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5657_CW1, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5657_CW1, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5657_CW1, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5657_CW1, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5657_CW1, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(0, 26, 5, 28, 7, 14, 24))) <- "NK/T Cells"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(29, 15))) <- "Dendritic Cells"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(11, 16, 10, 25))) <- "Monocytes"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(9, 23, 12))) <- "Macrophage"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(22, 21))) <- "Venous"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(20, 8, 18))) <- "Vascular Smooth Muscle"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(6, 13))) <- "Capillary"
#Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c())) <- "Mature Adipocytes"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(2, 4, 19, 1, 3))) <- "Stromal cells"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(17))) <- "B Cells"
Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(27))) <- "Transcriptional Doublets"
#Idents(P5657_CW1, cells = WhichCells(P5657_CW1, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5657_CW1, label = T) + NoLegend()
p2 <- DimPlot(P5657_CW1, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5657_CW1/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5657_CW1[['Gross_Annotation']] <- Idents(P5657_CW1)


saveRDS(P5657_CW1, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5657_CW1_singlet.rds")

###-------------------------------####
# P5657_CW3 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5657_CW3, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5657_CW3, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5657_CW3, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5657_CW3, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5657_CW3, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5657_CW3, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5657_CW3, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5657_CW3, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5657_CW3, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5657_CW3, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5657_CW3, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5657_CW3, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5657_CW3, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(27, 24))) <- "NK/T Cells"
Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(17))) <- "Dendritic Cells"
Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(7, 32))) <- "Monocytes"
Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(8, 28, 14, 13))) <- "Macrophage"
Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(21, 19))) <- "Venous"
Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(15, 6, 29, 0, 5))) <- "Vascular Smooth Muscle"
Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(30, 10, 3, 4, 16))) <- "Capillary"
#Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c())) <- "Mature Adipocytes"
Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(9, 11, 12, 1, 25, 22, 18, 2, 23, 26, 20))) <- "Stromal cells"
#Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(17))) <- "B Cells"
Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(31))) <- "Transcriptional Doublets"
#Idents(P5657_CW3, cells = WhichCells(P5657_CW3, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5657_CW3, label = T) + NoLegend()
p2 <- DimPlot(P5657_CW3, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5657_CW3/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5657_CW3[['Gross_Annotation']] <- Idents(P5657_CW3)


saveRDS(P5657_CW3, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5657_CW3_singlet.rds")

###-------------------------------####
# P5836_CW1 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5836_CW1, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5836_CW1, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5836_CW1, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5836_CW1, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5836_CW1, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5836_CW1, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5836_CW1, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5836_CW1, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5836_CW1, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5836_CW1, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5836_CW1, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5836_CW1, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5836_CW1, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(3, 5, 22, 17))) <- "NK/T Cells"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(6))) <- "Dendritic Cells"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(10, 13, 24))) <- "Monocytes"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(20, 21, 2))) <- "Macrophage"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(25, 8))) <- "Venous"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(11, 16, 9, 15))) <- "Vascular Smooth Muscle"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(19, 14, 12, 18))) <- "Capillary"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(28))) <- "Mature Adipocytes"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(23, 0, 26, 4, 1, 7))) <- "Stromal cells"
Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(27))) <- "B Cells"
#Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(27))) <- "Transcriptional Doublets"
#Idents(P5836_CW1, cells = WhichCells(P5836_CW1, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5836_CW1, label = T) + NoLegend()
p2 <- DimPlot(P5836_CW1, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5836_CW1/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5836_CW1[['Gross_Annotation']] <- Idents(P5836_CW1)


saveRDS(P5836_CW1, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5836_CW1_singlet.rds")

###-------------------------------####
# P5836_CW2 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5836_CW2, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5836_CW2, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5836_CW2, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5836_CW2, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5836_CW2, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5836_CW2, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5836_CW2, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5836_CW2, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5836_CW2, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5836_CW2, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5836_CW2, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5836_CW2, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5836_CW2, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(21, 8, 4, 3))) <- "NK/T Cells"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(26))) <- "Dendritic Cells"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(19, 13, 25))) <- "Monocytes"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(22, 24, 15, 17))) <- "Macrophage"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(11))) <- "Venous"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(14, 5, 2, 18, 16))) <- "Vascular Smooth Muscle"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(10, 23, 29, 20, 9, 12, 30, 6))) <- "Capillary"
#Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(28))) <- "Mature Adipocytes"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(7, 0, 1))) <- "Stromal cells"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(27))) <- "B Cells"
Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(28))) <- "Transcriptional Doublets"
#Idents(P5836_CW2, cells = WhichCells(P5836_CW2, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5836_CW2, label = T) + NoLegend()
p2 <- DimPlot(P5836_CW2, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5836_CW2/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5836_CW2[['Gross_Annotation']] <- Idents(P5836_CW2)


saveRDS(P5836_CW2, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5836_CW2_singlet.rds")

###-------------------------------####
# P5836_CW3 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5836_CW3, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5836_CW3, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5836_CW3, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5836_CW3, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5836_CW3, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5836_CW3, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5836_CW3, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5836_CW3, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5836_CW3, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5836_CW3, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5836_CW3, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5836_CW3, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5836_CW3, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(10, 15, 22, 28, 21, 14))) <- "NK/T Cells"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(6))) <- "Dendritic Cells"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(25, 16, 27))) <- "Monocytes"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(24, 11, 26, 4))) <- "Macrophage"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(18))) <- "Venous"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(17, 2))) <- "Vascular Smooth Muscle"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(23, 19, 12, 9, 20, 30))) <- "Capillary"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(32))) <- "Mature Adipocytes"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(0, 5, 8, 1, 7, 3, 13, 29))) <- "Stromal cells"
Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(31))) <- "B Cells"
#Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(28))) <- "Transcriptional Doublets"
#Idents(P5836_CW3, cells = WhichCells(P5836_CW3, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5836_CW3, label = T) + NoLegend()
p2 <- DimPlot(P5836_CW3, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5836_CW3/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5836_CW3[['Gross_Annotation']] <- Idents(P5836_CW3)


saveRDS(P5836_CW3, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5836_CW3_singlet.rds")

###-------------------------------####
# P5877_CW1 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5877_CW1, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5877_CW1, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5877_CW1, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5877_CW1, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5877_CW1, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5877_CW1, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5877_CW1, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5877_CW1, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5877_CW1, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5877_CW1, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5877_CW1, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5877_CW1, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5877_CW1, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(16, 3, 17, 20, 13))) <- "NK/T Cells"
Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(25))) <- "Dendritic Cells"
Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(11, 9, 22))) <- "Monocytes"
Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(21))) <- "Macrophage"
Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(19, 6))) <- "Venous"
Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(18, 12, 4, 0))) <- "Vascular Smooth Muscle"
Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(8, 7, 14, 5))) <- "Capillary"
#Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(32))) <- "Mature Adipocytes"
Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(23, 10, 2, 1, 15))) <- "Stromal cells"
Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(24))) <- "B Cells"
#Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(28))) <- "Transcriptional Doublets"
#Idents(P5877_CW1, cells = WhichCells(P5877_CW1, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5877_CW1, label = T) + NoLegend()
p2 <- DimPlot(P5877_CW1, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5877_CW1/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5877_CW1[['Gross_Annotation']] <- Idents(P5877_CW1)


saveRDS(P5877_CW1, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5877_CW1_singlet.rds")

###-------------------------------####
# P5877_CW2 Manual Annotation
###------------------------------####
# A lot of transcriptional doublets in this lane

p1 <- DimPlot(P5877_CW2, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5877_CW2, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5877_CW2, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5877_CW2, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5877_CW2, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5877_CW2, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5877_CW2, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5877_CW2, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5877_CW2, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5877_CW2, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5877_CW2, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5877_CW2, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5877_CW2, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


#Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c())) <- "NK/T Cells"
Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(10))) <- "Dendritic Cells"
Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(24, 16, 26))) <- "Monocytes"
Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(15, 12, 4, 27))) <- "Macrophage"
Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(11))) <- "Venous"
Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(23, 6))) <- "Vascular Smooth Muscle"
Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(3, 18, 25, 9))) <- "Capillary"
#Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(32))) <- "Mature Adipocytes"
Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(8, 7, 1, 5, 0, 20, 13, 2, 14, 19, 17))) <- "Stromal cells"
#Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(24))) <- "B Cells"
Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(21, 22))) <- "Transcriptional Doublets"
#Idents(P5877_CW2, cells = WhichCells(P5877_CW2, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5877_CW2, label = T) + NoLegend()
p2 <- DimPlot(P5877_CW2, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5877_CW2/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5877_CW2[['Gross_Annotation']] <- Idents(P5877_CW2)


saveRDS(P5877_CW2, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5877_CW2_singlet.rds")

###-------------------------------####
# P5877_CW3 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5877_CW3, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5877_CW3, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5877_CW3, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5877_CW3, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5877_CW3, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5877_CW3, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5877_CW3, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5877_CW3, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5877_CW3, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5877_CW3, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5877_CW3, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5877_CW3, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5877_CW3, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(22, 17))) <- "NK/T Cells"
Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(12))) <- "Dendritic Cells"
Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(20, 28, 8))) <- "Monocytes"
Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(7, 29, 10, 23, 24))) <- "Macrophage"
Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(25))) <- "Venous"
Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(26, 3))) <- "Vascular Smooth Muscle"
Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(19, 15))) <- "Capillary"
#Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(32))) <- "Mature Adipocytes"
Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(21, 18, 4, 14, 2, 5, 11, 13, 16, 11, 6, 9, 0, 1))) <- "Stromal cells"
Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(27))) <- "B Cells"
#Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(21, 22))) <- "Transcriptional Doublets"
#Idents(P5877_CW3, cells = WhichCells(P5877_CW3, idents = c(33, 35))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5877_CW3, label = T) + NoLegend()
p2 <- DimPlot(P5877_CW3, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5877_CW3/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5877_CW3[['Gross_Annotation']] <- Idents(P5877_CW3)


saveRDS(P5877_CW3, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5877_CW3_singlet.rds")

###-------------------------------####
# P5903_CW1 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5903_CW1, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5903_CW1, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5903_CW1, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5903_CW1, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5903_CW1, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5903_CW1, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5903_CW1, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5903_CW1, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5903_CW1, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5903_CW1, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5903_CW1, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5903_CW1, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5903_CW1, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(23, 10, 18, 19, 4, 5))) <- "NK/T Cells"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(21))) <- "Dendritic Cells"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(13, 14, 24, 28, 17))) <- "Monocytes"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(20, 22))) <- "Macrophage"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(16))) <- "Venous"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(15, 1, 11))) <- "Vascular Smooth Muscle"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(2, 9, 8, 12, 7, 6))) <- "Capillary"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(27))) <- "Mature Adipocytes"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(0, 3))) <- "Stromal cells"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(26))) <- "B Cells"
#Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(21, 22))) <- "Transcriptional Doublets"
Idents(P5903_CW1, cells = WhichCells(P5903_CW1, idents = c(25))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5903_CW1, label = T) + NoLegend()
p2 <- DimPlot(P5903_CW1, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5903_CW1/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5903_CW1[['Gross_Annotation']] <- Idents(P5903_CW1)


saveRDS(P5903_CW1, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5903_CW1_singlet.rds")

###-------------------------------####
# P5903_CW2 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5903_CW2, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5903_CW2, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5903_CW2, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5903_CW2, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5903_CW2, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5903_CW2, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5903_CW2, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5903_CW2, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5903_CW2, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5903_CW2, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5903_CW2, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5903_CW2, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5903_CW2, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(30, 23, 25))) <- "NK/T Cells"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(10))) <- "Dendritic Cells"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(21, 26, 35))) <- "Monocytes"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(9, 2, 22, 14, 33))) <- "Macrophage"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(12, 31, 29))) <- "Venous"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(24, 19, 28, 5, 11))) <- "Vascular Smooth Muscle"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(0, 3, 16, 17, 15))) <- "Capillary"
#Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c())) <- "Mature Adipocytes"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(4, 8, 36, 13, 6, 7, 1, 20, 27, 18))) <- "Stromal cells"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(34))) <- "B Cells"
Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(32))) <- "Transcriptional Doublets"
#Idents(P5903_CW2, cells = WhichCells(P5903_CW2, idents = c(25))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5903_CW2, label = T) + NoLegend()
p2 <- DimPlot(P5903_CW2, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5903_CW2/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5903_CW2[['Gross_Annotation']] <- Idents(P5903_CW2)


saveRDS(P5903_CW2, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5903_CW2_singlet.rds")

###-------------------------------####
# P5903_CW3 Manual Annotation
###------------------------------####
# A lot of transcriptional doublets

p1 <- DimPlot(P5903_CW3, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5903_CW3, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5903_CW3, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5903_CW3, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5903_CW3, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5903_CW3, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5903_CW3, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5903_CW3, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5903_CW3, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5903_CW3, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5903_CW3, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5903_CW3, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5903_CW3, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(25, 14, 35))) <- "NK/T Cells"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(21, 32))) <- "Dendritic Cells"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(24, 31, 7))) <- "Monocytes"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(9, 3, 8, 4, 16, 12, 10, 1, 23, 34, 22))) <- "Macrophage"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(26))) <- "Venous"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(27, 13))) <- "Vascular Smooth Muscle"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(28, 15, 20, 19))) <- "Capillary"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(30))) <- "Mature Adipocytes"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(2, 5, 17, 5, 6, 11, 0, 29))) <- "Stromal cells"
#Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(34))) <- "B Cells"
Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(33, 18))) <- "Transcriptional Doublets"
#Idents(P5903_CW3, cells = WhichCells(P5903_CW3, idents = c(25))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5903_CW3, label = T) + NoLegend()
p2 <- DimPlot(P5903_CW3, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5903_CW3/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5903_CW3[['Gross_Annotation']] <- Idents(P5903_CW3)


saveRDS(P5903_CW3, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5903_CW3_singlet.rds")

###-------------------------------####
# P5963_CW1 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5963_CW1, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5963_CW1, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5963_CW1, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5963_CW1, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5963_CW1, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5963_CW1, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5963_CW1, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5963_CW1, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5963_CW1, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5963_CW1, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5963_CW1, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5963_CW1, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5963_CW1, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(26, 18, 6, 7, 31))) <- "NK/T Cells"
Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(24))) <- "Dendritic Cells"
Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(21, 23, 27))) <- "Monocytes"
Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(22, 19, 30, 11, 10))) <- "Macrophage"
Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(16))) <- "Venous"
Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(15, 25, 0, 29))) <- "Vascular Smooth Muscle"
Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(12, 1, 4, 20, 14, 13))) <- "Capillary"
Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(32))) <- "Mature Adipocytes"
Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(5, 28, 17, 8, 3, 2, 9))) <- "Stromal cells"
#Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(34))) <- "B Cells"
#Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(33, 18))) <- "Transcriptional Doublets"
#Idents(P5963_CW1, cells = WhichCells(P5963_CW1, idents = c(25))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5963_CW1, label = T) + NoLegend()
p2 <- DimPlot(P5963_CW1, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5963_CW1/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5963_CW1[['Gross_Annotation']] <- Idents(P5963_CW1)


saveRDS(P5963_CW1, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5963_CW1_singlet.rds")

###-------------------------------####
# P5963_CW2 Manual Annotation
###------------------------------####

p1 <- DimPlot(P5963_CW2, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5963_CW2, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5963_CW2, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5963_CW2, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5963_CW2, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5963_CW2, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5963_CW2, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5963_CW2, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5963_CW2, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5963_CW2, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5963_CW2, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5963_CW2, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5963_CW2, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(21, 25))) <- "NK/T Cells"
Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(22))) <- "Dendritic Cells"
Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(29, 10, 30))) <- "Monocytes"
Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(31, 6, 24, 13, 4, 27))) <- "Macrophage"
Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(11, 17))) <- "Venous"
Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(23, 2, 1))) <- "Vascular Smooth Muscle"
Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(33, 8, 18, 20, 19, 12))) <- "Capillary"
#Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(32))) <- "Mature Adipocytes"
Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(15, 14, 16, 5, 32, 0, 7, 9, 3))) <- "Stromal cells"
#Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(34))) <- "B Cells"
Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(34, 28, 26))) <- "Transcriptional Doublets"
#Idents(P5963_CW2, cells = WhichCells(P5963_CW2, idents = c(25))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5963_CW2, label = T) + NoLegend()
p2 <- DimPlot(P5963_CW2, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5963_CW2/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5963_CW2[['Gross_Annotation']] <- Idents(P5963_CW2)


saveRDS(P5963_CW2, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5963_CW2_singlet.rds")

###-------------------------------####
# P5963_CW3 Manual Annotation
###------------------------------####
# A lot of doublets, but data looks cleaner with doublets removed so will further process after merging
# in hopes to salvage some of the T cells

p1 <- DimPlot(P5963_CW3, label = T) + NoLegend()

# T/NK Cells
p2 <- FeaturePlot(P5963_CW3, features = c("CD3E", "CD8A", "NKG7", "IL7R"))
p1 + p2

# Dendritic
p2 <- FeaturePlot(P5963_CW3, features = c("CCR7", "CLEC9A", "CLEC10A", "LILRA4", "CD1C"))
p1 + p2

# Monocytes
p2 <- FeaturePlot(P5963_CW3, features = c("CD14", "VCAN", "LYZ", "FCGR3A"))
p1 + p2

# Macrophages
p2 <- FeaturePlot(P5963_CW3, features = c("CD68", "MARCO", "RNASE1", "APOE", "C1QB"))
p1 + p2

# Venous
p2 <- FeaturePlot(P5963_CW3, features = c("ACKR1", "NR2F2", "VCAM1"))
p1 + p2

# VSM
p2 <- FeaturePlot(P5963_CW3, features = c("ACTA2", "TAGLN", "RGS5"))
p1 + p2

# Capillary
p2 <- FeaturePlot(P5963_CW3, features = c("CA4", "PRX", "SPARC", "CLDN5"))
p1 + p2

# Mature Adipocytes
p2 <- FeaturePlot(P5963_CW3, features = c("CIDEC", "ADIPOQ", "LPL", "PPARG"))
p1 + p2

# Stromal cells
p2 <- FeaturePlot(P5963_CW3, features = c("CFD", "CCDC80", "LUM", "MGP"))
p1 + p2

# B Cells
p2 <- FeaturePlot(P5963_CW3, features = c("MS4A1", "JCHAIN", "CD79A", "BANK1"))
p1 + p2

# Hemoglobin
p2 <- FeaturePlot(P5963_CW3, features = c("HBB"))
p1 + p2

# Doublets
p2 <- FeaturePlot(P5963_CW3, features = c("CLDN5", "PTPRC", "COL1A2"))
p1 + p2


Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(25))) <- "NK/T Cells"
Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(26))) <- "Dendritic Cells"
Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(23, 16))) <- "Monocytes"
Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(3, 14, 29, 32, 28, 33))) <- "Macrophage"
Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(8, 20, 35))) <- "Venous"
Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(18, 34, 24, 7, 6, 9, 36))) <- "Vascular Smooth Muscle"
Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(10, 30, 12, 1, 22))) <- "Capillary"
#Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(32))) <- "Mature Adipocytes"
Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(15, 0, 2, 11, 17, 13, 5, 19, 27))) <- "Stromal cells"
#Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(34))) <- "B Cells"
Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(21, 4, 31))) <- "Transcriptional Doublets"
#Idents(P5963_CW3, cells = WhichCells(P5963_CW3, idents = c(25))) <- "Erythrocytes"

# Visualize Umap and Doublet Finder
p1 <-DimPlot(P5963_CW3, label = T) + NoLegend()
p2 <- DimPlot(P5963_CW3, group.by = "DoubletFinder")
pdf(file = "/data/p_koethe_lab/Atlas_AT/Analysis/P5963_CW3/Doublet_UMAP.pdf")
p1 + p2
dev.off()

# Set Idents to gross celltype
P5963_CW3[['Gross_Annotation']] <- Idents(P5963_CW3)

saveRDS(P5963_CW3, "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/P5963_CW3_singlet.rds")

# 3. LINK SOUPORCELL GENETIC DEMULTIPLEXING WITH CITE-SEQ HASHING-------------------------------------

###-------------------------------####
# Load Processed Singlet Lanes
###------------------------------####
path <- "/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Singlets_8_18/"
file_list <- list.files(path = path, pattern = ".rds")

names(file_list) <- c('P5344_CW1', 'P5344_CW2', 'P5344_CW3', 'P5544_CW1', 'P5544_CW2', 'P5544_CW3', 'P5573_CW1',
                      'P5573_CW2', 'P5573_CW3', 'P5657_CW1', 'P5657_CW3', 'P5836_CW1', 'P5836_CW2',
                      'P5836_CW3', 'P5877_CW1', 'P5877_CW2', 'P5877_CW3', 'P5903_CW1', 'P5903_CW2', 'P5903_CW3',
                      'P5963_CW1', 'P5963_CW2', 'P5963_CW3')

for (i in 1:length(file_list)){
  print(file_list[[i]])
  print(names(file_list[i]))
  seurat_object <- readRDS(paste(path, file_list[[i]], sep = ""))
  assign(names(file_list[i]), seurat_object)
}

###-------------------------------####
# Load Processed Singlet Lanes
###------------------------------####

###-------------------------------####
# Link SoupOrCell to Hashtag
###------------------------------####

# P5344_CW1
table(P5344_CW1$HTO_classification, P5344_CW1$assignment)
# assignment: 251 (1) 252 (0), 253 (3), 254 (2)

P5344_CW1[['SoupHashLink']] <- ifelse(P5344_CW1[['assignment']] == 1, 251, ifelse(P5344_CW1[['assignment']] == 0, 252,
                                                                                  ifelse(P5344_CW1[['assignment']] == 3, 253, 254)))
table(P5344_CW1$HTO_classification, P5344_CW1$SoupHashLink) # Check correct assignment

# P5344_CW2
table(P5344_CW2$HTO_classification, P5344_CW2$assignment)
# assignment: 251 (1) 252 (3), 253 (2), 254 (0)
P5344_CW2[['SoupHashLink']] <- ifelse(P5344_CW2[['assignment']] == 1, 251, ifelse(P5344_CW2[['assignment']] == 3, 252,
                                                                                  ifelse(P5344_CW2[['assignment']] == 2, 253, 254)))
table(P5344_CW2$HTO_classification, P5344_CW2$SoupHashLink) # Check correct assignment
table(P5344_CW2$SoupHashLink, P5344_CW2$assignment)

# P5344_CW3
table(P5344_CW3$HTO_classification, P5344_CW3$assignment)
# assignment: 251 (0) 252 (3), 253 (2), 254 (1)
P5344_CW3[['SoupHashLink']] <- ifelse(P5344_CW3[['assignment']] == 0, 251, ifelse(P5344_CW3[['assignment']] == 3, 252,
                                                                                  ifelse(P5344_CW3[['assignment']] == 2, 253, 254)))
table(P5344_CW3$HTO_classification, P5344_CW3$SoupHashLink) # Check correct assignment
table(P5344_CW3$SoupHashLink, P5344_CW3$assignment)

# P5544_CW1
table(P5544_CW1$HTO_classification, P5544_CW1$assignment)
# assignment: 251 (3) 252 (1), 253 (0), 254 (2)
P5544_CW1[['SoupHashLink']] <- ifelse(P5544_CW1[['assignment']] == 3, 251, ifelse(P5544_CW1[['assignment']] == 1, 252,
                                                                                  ifelse(P5544_CW1[['assignment']] == 0, 253, 254)))
table(P5544_CW1$HTO_classification, P5544_CW1$SoupHashLink) # Check correct assignment
table(P5544_CW1$SoupHashLink, P5544_CW1$assignment)

# P5544_CW2
table(P5544_CW2$HTO_classification, P5544_CW2$assignment)
# assignment: 251 (1) 252 (2), 253 (3), 254 (0)
P5544_CW2[['SoupHashLink']] <- ifelse(P5544_CW2[['assignment']] == 1, 251, ifelse(P5544_CW2[['assignment']] == 2, 252,
                                                                                  ifelse(P5544_CW2[['assignment']] == 3, 253, 254)))
table(P5544_CW2$HTO_classification, P5544_CW2$SoupHashLink) # Check correct assignment
table(P5544_CW2$SoupHashLink, P5544_CW2$assignment)

# P5544_CW3
table(P5544_CW3$HTO_classification, P5544_CW3$assignment)
# assignment: 251 (1) 252 (0), 253 (2), 254 (3)
P5544_CW3[['SoupHashLink']] <- ifelse(P5544_CW3[['assignment']] == 1, 251, ifelse(P5544_CW3[['assignment']] == 0, 252,
                                                                                  ifelse(P5544_CW3[['assignment']] == 2, 253, 254)))
table(P5544_CW3$HTO_classification, P5544_CW3$SoupHashLink) # Check correct assignment
table(P5544_CW3$SoupHashLink, P5544_CW3$assignment)

# P5573_CW1
table(P5573_CW1$HTO_classification, P5573_CW1$assignment)
# assignment: 251 (2) 252 (1), 253 (3), 254 (0)
P5573_CW1[['SoupHashLink']] <- ifelse(P5573_CW1[['assignment']] == 2, 251, ifelse(P5573_CW1[['assignment']] == 1, 252,
                                                                                  ifelse(P5573_CW1[['assignment']] == 3, 253, 254)))
table(P5573_CW1$HTO_classification, P5573_CW1$SoupHashLink) # Check correct assignment
table(P5573_CW1$SoupHashLink, P5573_CW1$assignment)

# P5573_CW2
table(P5573_CW2$HTO_classification, P5573_CW2$assignment)
# assignment: 251 (3) 252 (2), 253 (0), 254 (1)
P5573_CW2[['SoupHashLink']] <- ifelse(P5573_CW2[['assignment']] == 3, 251, ifelse(P5573_CW2[['assignment']] == 2, 252,
                                                                                  ifelse(P5573_CW2[['assignment']] == 0, 253, 254)))
table(P5573_CW2$HTO_classification, P5573_CW2$SoupHashLink) # Check correct assignment
table(P5573_CW2$SoupHashLink, P5573_CW2$assignment)

# P5573_CW3
table(P5573_CW3$HTO_classification, P5573_CW3$assignment)
# assignment: 251 (3) 252 (1), 253 (2), 254 (0)
P5573_CW3[['SoupHashLink']] <- ifelse(P5573_CW3[['assignment']] == 3, 251, ifelse(P5573_CW3[['assignment']] == 1, 252,
                                                                                  ifelse(P5573_CW3[['assignment']] == 2, 253, 254)))
table(P5573_CW3$HTO_classification, P5573_CW3$SoupHashLink) # Check correct assignment
table(P5573_CW3$SoupHashLink, P5573_CW3$assignment)

# P5657_CW1
table(P5657_CW1$HTO_classification, P5657_CW1$assignment)
# assignment: 251 (3) 252 (2), 253 (0), 254 (1)
P5657_CW1[['SoupHashLink']] <- ifelse(P5657_CW1[['assignment']] == 3, 251, ifelse(P5657_CW1[['assignment']] == 2, 252,
                                                                                  ifelse(P5657_CW1[['assignment']] == 0, 253, 254)))
table(P5657_CW1$HTO_classification, P5657_CW1$SoupHashLink) # Check correct assignment
table(P5657_CW1$SoupHashLink, P5657_CW1$assignment)

# P5657_CW3
table(P5657_CW3$HTO_classification, P5657_CW3$assignment)
# assignment: 251 (0) 252 (3), 253 (2), 254 (1)
P5657_CW3[['SoupHashLink']] <- ifelse(P5657_CW3[['assignment']] == 0, 251, ifelse(P5657_CW3[['assignment']] == 3, 252,
                                                                                  ifelse(P5657_CW3[['assignment']] == 2, 253, 254)))
table(P5657_CW3$HTO_classification, P5657_CW3$SoupHashLink) # Check correct assignment
table(P5657_CW3$SoupHashLink, P5657_CW3$assignment)

# P5836_CW1
table(P5836_CW1$HTO_classification, P5836_CW1$assignment)
# assignment: 251 (3) 252 (NA), 253 (2), 254 (0)
P5836_CW1[['SoupHashLink']] <- ifelse(P5836_CW1[['assignment']] == 3, 251, ifelse(P5836_CW1[['assignment']] == 1, 252,
                                                                                  ifelse(P5836_CW1[['assignment']] == 2, 253, 254)))
table(P5836_CW1$HTO_classification, P5836_CW1$SoupHashLink) # Check correct assignment
table(P5836_CW1$SoupHashLink, P5836_CW1$assignment)

# P5836_CW2
table(P5836_CW2$HTO_classification, P5836_CW2$assignment)
# assignment: 251 (1) 252 (2), 253 (0), 254 (3)
P5836_CW2[['SoupHashLink']] <- ifelse(P5836_CW2[['assignment']] == 1, 251, ifelse(P5836_CW2[['assignment']] == 2, 252,
                                                                                  ifelse(P5836_CW2[['assignment']] == 0, 253, 254)))
table(P5836_CW2$HTO_classification, P5836_CW2$SoupHashLink) # Check correct assignment
table(P5836_CW2$SoupHashLink, P5836_CW2$assignment)

# P5836_CW3
table(P5836_CW3$HTO_classification, P5836_CW3$assignment)
# assignment: 251 (2) 252 (3), 253 (0), 254 (1)
P5836_CW3[['SoupHashLink']] <- ifelse(P5836_CW3[['assignment']] == 2, 251, ifelse(P5836_CW3[['assignment']] == 3, 252,
                                                                                  ifelse(P5836_CW3[['assignment']] == 0, 253, 254)))
table(P5836_CW3$HTO_classification, P5836_CW3$SoupHashLink) # Check correct assignment
table(P5836_CW3$SoupHashLink, P5836_CW3$assignment)

# P5877_CW1
table(P5877_CW1$HTO_classification, P5877_CW1$assignment)
# assignment: 251 (1) 252 (3), 253 (2), 254 (0)
P5877_CW1[['SoupHashLink']] <- ifelse(P5877_CW1[['assignment']] == 1, 251, ifelse(P5877_CW1[['assignment']] == 3, 252,
                                                                                  ifelse(P5877_CW1[['assignment']] == 2, 253, 254)))
table(P5877_CW1$HTO_classification, P5877_CW1$SoupHashLink) # Check correct assignment
table(P5877_CW1$SoupHashLink, P5877_CW1$assignment)

# P5877_CW2
table(P5877_CW2$HTO_classification, P5877_CW2$assignment)
# assignment: 251 (2) 252 (0), 253 (1), 254 (3)
P5877_CW2[['SoupHashLink']] <- ifelse(P5877_CW2[['assignment']] == 2, 251, ifelse(P5877_CW2[['assignment']] == 0, 252,
                                                                                  ifelse(P5877_CW2[['assignment']] == 1, 253, 254)))
table(P5877_CW2$HTO_classification, P5877_CW2$SoupHashLink) # Check correct assignment
table(P5877_CW2$SoupHashLink, P5877_CW2$assignment)

# P5877_CW3
table(P5877_CW3$HTO_classification, P5877_CW3$assignment)
# assignment: 251 (1) 252 (3), 253 (0), 254 (2)
P5877_CW3[['SoupHashLink']] <- ifelse(P5877_CW3[['assignment']] == 1, 251, ifelse(P5877_CW3[['assignment']] == 3, 252,
                                                                                  ifelse(P5877_CW3[['assignment']] == 0, 253, 254)))
table(P5877_CW3$HTO_classification, P5877_CW3$SoupHashLink) # Check correct assignment
table(P5877_CW3$SoupHashLink, P5877_CW3$assignment)

# P5903_CW1
table(P5903_CW1$HTO_classification, P5903_CW1$assignment)
# assignment: 251 (1) 252 (2), 253 (3), 254 (0)
P5903_CW1[['SoupHashLink']] <- ifelse(P5903_CW1[['assignment']] == 1, 251, ifelse(P5903_CW1[['assignment']] == 2, 252,
                                                                                  ifelse(P5903_CW1[['assignment']] == 3, 253, 254)))
table(P5903_CW1$HTO_classification, P5903_CW1$SoupHashLink) # Check correct assignment
table(P5903_CW1$SoupHashLink, P5903_CW1$assignment)

# P5903_CW2
table(P5903_CW2$HTO_classification, P5903_CW2$assignment)
# assignment: 251 (1) 252 (0), 253 (3), 254 (2)
P5903_CW2[['SoupHashLink']] <- ifelse(P5903_CW2[['assignment']] == 1, 251, ifelse(P5903_CW2[['assignment']] == 0, 252,
                                                                                  ifelse(P5903_CW2[['assignment']] == 3, 253, 254)))
table(P5903_CW2$HTO_classification, P5903_CW2$SoupHashLink) # Check correct assignment
table(P5903_CW2$SoupHashLink, P5903_CW2$assignment)

# P5903_CW3
table(P5903_CW3$HTO_classification, P5903_CW3$assignment)
# assignment: 251 (0) 252 (3), 253 (1), 254 (2)
P5903_CW3[['SoupHashLink']] <- ifelse(P5903_CW3[['assignment']] == 0, 251, ifelse(P5903_CW3[['assignment']] == 3, 252,
                                                                                  ifelse(P5903_CW3[['assignment']] == 1, 253, 254)))
table(P5903_CW3$HTO_classification, P5903_CW3$SoupHashLink) # Check correct assignment
table(P5903_CW3$SoupHashLink, P5903_CW3$assignment)

# P5963_CW1
table(P5963_CW1$HTO_classification, P5963_CW1$assignment)
# assignment: 251 (2) 252 (3), 253 (0), 254 (1)
P5963_CW1[['SoupHashLink']] <- ifelse(P5963_CW1[['assignment']] == 2, 251, ifelse(P5963_CW1[['assignment']] == 3, 252,
                                                                                  ifelse(P5963_CW1[['assignment']] == 0, 253, 254)))
table(P5963_CW1$HTO_classification, P5963_CW1$SoupHashLink) # Check correct assignment
table(P5963_CW1$SoupHashLink, P5963_CW1$assignment)

# P5963_CW2
table(P5963_CW2$HTO_classification, P5963_CW2$assignment)
# assignment: 251 (0) 252 (2), 253 (1), 254 (3)
P5963_CW2[['SoupHashLink']] <- ifelse(P5963_CW2[['assignment']] == 0, 251, ifelse(P5963_CW2[['assignment']] == 2, 252,
                                                                                  ifelse(P5963_CW2[['assignment']] == 1, 253, 254)))
table(P5963_CW2$HTO_classification, P5963_CW2$SoupHashLink) # Check correct assignment
table(P5963_CW2$SoupHashLink, P5963_CW2$assignment)

# P5963_CW3
table(P5963_CW3$HTO_classification, P5963_CW3$assignment)
# assignment: 251 (1) 252 (3), 253 (0), 254 (2)
P5963_CW3[['SoupHashLink']] <- ifelse(P5963_CW3[['assignment']] == 1, 251, ifelse(P5963_CW3[['assignment']] == 3, 252,
                                                                                  ifelse(P5963_CW3[['assignment']] == 0, 253, 254)))
table(P5963_CW3$HTO_classification, P5963_CW3$SoupHashLink) # Check correct assignment
table(P5963_CW3$SoupHashLink, P5963_CW3$assignment)

###-------------------------------####
# Link HATIMID with Hashtag/SouporCell
###------------------------------####

# Hashtag link by Project
# P5344_CW1: 1141 - 251; 1147-252; 1157-253; 1158-254
# P5344_CW2: 3018 - 251; 3019-252; 3020-253; 3021-254
# P5344_CW3: 5003-251; 5022-252; 5032-253;5033-254

# P5544_CW1: 1117-251; 1137-252; 1138-253; 1139-254
# P5544_CW2: 3012-251; 3015-252; 3023-253; 3024-254
# P5544_CW3: 5024-251; 5027-252; 5029-253; 5030-254

# P5573_CW1:1103-251; 1126-252; 1172-253; 1175-254
# P5573_CW2:3007-251; 3014-252; 3026-253; 3028-254
# P5573_CW3:5011-251; 5012-252; 5013-253; 5015-254

# P5344_CW1
P5344_CW1[['HATIMID']] <- ifelse(P5344_CW1[['SoupHashLink']] == 251, 1141,
                                 ifelse(P5344_CW1[['SoupHashLink']] == 252, 1147,
                                        ifelse(P5344_CW1[["SoupHashLink"]] == 253, 1157, 1158)))
table(P5344_CW1$HATIMID, P5344_CW1$SoupHashLink)

# P5344_CW2
P5344_CW2[['HATIMID']] <- ifelse(P5344_CW2[['SoupHashLink']] == 251, 3018,
                                 ifelse(P5344_CW2[['SoupHashLink']] == 252, 3019,
                                        ifelse(P5344_CW2[["SoupHashLink"]] == 253, 3020, 3021)))
table(P5344_CW2$HATIMID, P5344_CW2$SoupHashLink)

# P5344_CW3
P5344_CW3[['HATIMID']] <- ifelse(P5344_CW3[['SoupHashLink']] == 251, 5003,
                                 ifelse(P5344_CW3[['SoupHashLink']] == 252, 5022,
                                        ifelse(P5344_CW3[["SoupHashLink"]] == 253, 5032, 5033)))
table(P5344_CW3$HATIMID, P5344_CW3$SoupHashLink)

# P5544_CW1
P5544_CW1[['HATIMID']] <- ifelse(P5544_CW1[['SoupHashLink']] == 251, 1117,
                                 ifelse(P5544_CW1[['SoupHashLink']] == 252, 1137,
                                        ifelse(P5544_CW1[["SoupHashLink"]] == 253, 1138, 1139)))
table(P5544_CW1$HATIMID, P5544_CW1$SoupHashLink)

# P5544_CW2
P5544_CW2[['HATIMID']] <- ifelse(P5544_CW2[['SoupHashLink']] == 251, 3012,
                                 ifelse(P5544_CW2[['SoupHashLink']] == 252, 3015,
                                        ifelse(P5544_CW2[["SoupHashLink"]] == 253, 3023, 3024)))
table(P5544_CW2$HATIMID, P5544_CW2$SoupHashLink)

# P5544_CW3
P5544_CW3[['HATIMID']] <- ifelse(P5544_CW3[['SoupHashLink']] == 251, 5024,
                                 ifelse(P5544_CW3[['SoupHashLink']] == 252, 5027,
                                        ifelse(P5544_CW3[["SoupHashLink"]] == 253, 5029, 5030)))
table(P5544_CW3$HATIMID, P5544_CW3$SoupHashLink)

# P5573_CW1
P5573_CW1[['HATIMID']] <- ifelse(P5573_CW1[['SoupHashLink']] == 251, 1103,
                                 ifelse(P5573_CW1[['SoupHashLink']] == 252, 1126,
                                        ifelse(P5573_CW1[["SoupHashLink"]] == 253, 1172, 1175)))
table(P5573_CW1$HATIMID, P5573_CW1$SoupHashLink)

# P5573_CW2
P5573_CW2[['HATIMID']] <- ifelse(P5573_CW2[['SoupHashLink']] == 251, 3007,
                                 ifelse(P5573_CW2[['SoupHashLink']] == 252, 3014,
                                        ifelse(P5573_CW2[["SoupHashLink"]] == 253, 3026, 3028)))
table(P5573_CW2$HATIMID, P5573_CW2$SoupHashLink)

# P5573_CW3
# P5573_CW3:5011-251; 5012-252; 5013-253; 5015-254
P5573_CW3[['HATIMID']] <- ifelse(P5573_CW3[['SoupHashLink']] == 251, 5011,
                                 ifelse(P5573_CW3[['SoupHashLink']] == 252, 5012,
                                        ifelse(P5573_CW3[["SoupHashLink"]] == 253, 5013, 5015)))
table(P5573_CW3$HATIMID, P5573_CW3$SoupHashLink)

# P5657_CW1: 1112-251; 1113-252; 1136-253; 254-1161
P5657_CW1[['HATIMID']] <- ifelse(P5657_CW1[['SoupHashLink']] == 251, 1112,
                                 ifelse(P5657_CW1[['SoupHashLink']] == 252, 1113,
                                        ifelse(P5657_CW1[["SoupHashLink"]] == 253, 1136, 1161)))
table(P5657_CW1$HATIMID, P5657_CW1$SoupHashLink)

# P5657_CW3: 5005-251; 5009-252; 5028-253; 5031-254
P5657_CW3[['HATIMID']] <- ifelse(P5657_CW3[['SoupHashLink']] == 251, 5005,
                                 ifelse(P5657_CW3[['SoupHashLink']] == 252, 5009,
                                        ifelse(P5657_CW3[["SoupHashLink"]] == 253, 5028, 5031)))
table(P5657_CW3$HATIMID, P5657_CW3$SoupHashLink)

# P5836_CW1: 1135-251; 1143-252; 1149-253; 1155-254
P5836_CW1[['HATIMID']] <- ifelse(P5836_CW1[['SoupHashLink']] == 251, 1135,
                                 ifelse(P5836_CW1[['SoupHashLink']] == 252, 1143,
                                        ifelse(P5836_CW1[["SoupHashLink"]] == 253, 1149, 1155)))
table(P5836_CW1$HATIMID, P5836_CW1$SoupHashLink)

# P5836_CW2: 3013-251; 3022-252; 3032-253; 3035-254
P5836_CW2[['HATIMID']] <- ifelse(P5836_CW2[['SoupHashLink']] == 251, 3013,
                                 ifelse(P5836_CW2[['SoupHashLink']] == 252, 3022,
                                        ifelse(P5836_CW2[["SoupHashLink"]] == 253, 3032, 3035)))
table(P5836_CW2$HATIMID, P5836_CW2$SoupHashLink)

# P5836_CW3: 5004-251; 5008-252; 5010-253; 5016-254
P5836_CW3[['HATIMID']] <- ifelse(P5836_CW3[['SoupHashLink']] == 251, 5004,
                                 ifelse(P5836_CW3[['SoupHashLink']] == 252, 5008,
                                        ifelse(P5836_CW3[["SoupHashLink"]] == 253, 5010, 5016)))
table(P5836_CW3$HATIMID, P5836_CW3$SoupHashLink)

# P5877_CW1: 1109-251; 1132-252; 1151-253; 1167-254
P5877_CW1[['HATIMID']] <- ifelse(P5877_CW1[['SoupHashLink']] == 251, 1109,
                                 ifelse(P5877_CW1[['SoupHashLink']] == 252, 1132,
                                        ifelse(P5877_CW1[["SoupHashLink"]] == 253, 1151, 1167)))
table(P5877_CW1$HATIMID, P5877_CW1$SoupHashLink)

# P5877_CW2: 3011-251; 3029-252; 3031-253; 3037-254
P5877_CW2[['HATIMID']] <- ifelse(P5877_CW2[['SoupHashLink']] == 251, 3011,
                                 ifelse(P5877_CW2[['SoupHashLink']] == 252, 3029,
                                        ifelse(P5877_CW2[["SoupHashLink"]] == 253, 3031, 3037)))
table(P5877_CW2$HATIMID, P5877_CW2$SoupHashLink)

# P5877_CW3: 5001-251; 5002-252; 5017-253; 5019-254
P5877_CW3[['HATIMID']] <- ifelse(P5877_CW3[['SoupHashLink']] == 251, 5001,
                                 ifelse(P5877_CW3[['SoupHashLink']] == 252, 5002,
                                        ifelse(P5877_CW3[["SoupHashLink"]] == 253, 5017, 5019)))
table(P5877_CW3$HATIMID, P5877_CW3$SoupHashLink)

# P5903_CW1: 1108-251; 1180-252; 1181-253; 1182-254
P5903_CW1[['HATIMID']] <- ifelse(P5903_CW1[['SoupHashLink']] == 251, 1108,
                                 ifelse(P5903_CW1[['SoupHashLink']] == 252, 1180,
                                        ifelse(P5903_CW1[["SoupHashLink"]] == 253, 1181, 1182)))
table(P5903_CW1$HATIMID, P5903_CW1$SoupHashLink)

# P5903_CW2: 2108-251; 2111-252; 2113-253; 2116-254
P5903_CW2[['HATIMID']] <- ifelse(P5903_CW2[['SoupHashLink']] == 251, 2108,
                                 ifelse(P5903_CW2[['SoupHashLink']] == 252, 2111,
                                        ifelse(P5903_CW2[["SoupHashLink"]] == 253, 2113, 2116)))
table(P5903_CW2$HATIMID, P5903_CW2$SoupHashLink)

# P5903_CW3: 5018-251; 5020-252; 5021-253; 5026-254
P5903_CW3[['HATIMID']] <- ifelse(P5903_CW3[['SoupHashLink']] == 251, 5018,
                                 ifelse(P5903_CW3[['SoupHashLink']] == 252, 5020,
                                        ifelse(P5903_CW3[["SoupHashLink"]] == 253, 5021, 5026)))
table(P5903_CW3$HATIMID, P5903_CW3$SoupHashLink)

# P5963_CW1: 1107-251; 1159-252; 1165-253; 1176-254
P5963_CW1[['HATIMID']] <- ifelse(P5963_CW1[['SoupHashLink']] == 251, 1107,
                                 ifelse(P5963_CW1[['SoupHashLink']] == 252, 1159,
                                        ifelse(P5963_CW1[["SoupHashLink"]] == 253, 1165, 1176)))
table(P5963_CW1$HATIMID, P5963_CW1$SoupHashLink)

# P5963_CW2: 2107-251; 2110-252; 2112-253; 2114-254
P5963_CW2[['HATIMID']] <- ifelse(P5963_CW2[['SoupHashLink']] == 251, 2107,
                                 ifelse(P5963_CW2[['SoupHashLink']] == 252, 2110,
                                        ifelse(P5963_CW2[["SoupHashLink"]] == 253, 2112, 2114)))
table(P5963_CW2$HATIMID, P5963_CW2$SoupHashLink)

# P5963_CW3: 5034-251; 5036-252; 5037-253; 5038-254
P5963_CW3[['HATIMID']] <- ifelse(P5963_CW3[['SoupHashLink']] == 251, 5034,
                                 ifelse(P5963_CW3[['SoupHashLink']] == 252, 5036,
                                        ifelse(P5963_CW3[["SoupHashLink"]] == 253, 5037, 5038)))
table(P5963_CW3$HATIMID, P5963_CW3$SoupHashLink)

# 4. ADD METADATA TO SEURAT OBJECTS-------------------------------------
###-------------------------------####
# Link Metadata
###------------------------------####
meta.data <- read.csv("/data/p_koethe_lab/Atlas_AT/MetaData/Reference/HATIMStudyVisit_DATA_2021-08-09_1359.csv")

# Rename metabolic groups
meta.data <- meta.data %>%
  dplyr::mutate(StudyGroup = factor(hatim_final_arm, levels = c(1, 2, 3, 4), labels = c("HIV+ non-diabetic", "HIV+ prediabetic",
                                                                                        "HIV+ diabetic", "HIV- diabetic")))
# Add HATIM
meta.data <- meta.data %>%
  dplyr::mutate(HATIMID = as.numeric(as.character(hatim_clin_visit_pid)))

# Pull metformin data from metadata
meta.data$metformin <- ifelse(stringr::str_detect(meta.data$other_meds, regex('metformin', ignore_case = T)), 1, 0)

# Pull integrase data from metadata
meta.data$INSTI <- ifelse(meta.data$meta_current_art == 1, 1, ifelse(meta.data$meta_current_art == 4, 1, 0))

# Create a loop to add metadata
seurat_list <- c(P5344_CW1, P5344_CW2, P5344_CW3, P5544_CW1, P5544_CW2, P5544_CW3, P5573_CW1, P5573_CW2, P5573_CW3, P5657_CW1, P5657_CW3,
                 P5836_CW1, P5836_CW2, P5836_CW3, P5877_CW1, P5877_CW2, P5877_CW3, P5903_CW1, P5903_CW2, P5903_CW3, P5963_CW1, P5963_CW2, P5963_CW3)
names(seurat_list) <- c('P5344_CW1', 'P5344_CW2', 'P5344_CW3', 'P5544_CW1', 'P5544_CW2', 'P5544_CW3', 'P5573_CW1',
                        'P5573_CW2', 'P5573_CW3', 'P5657_CW1', 'P5657_CW3', 'P5836_CW1', 'P5836_CW2',
                        'P5836_CW3', 'P5877_CW1', 'P5877_CW2', 'P5877_CW3', 'P5903_CW1', 'P5903_CW2', 'P5903_CW3',
                        'P5963_CW1', 'P5963_CW2', 'P5963_CW3')

for(i in 1:length(seurat_list)) {
  print(names(seurat_list[i]))
  seurat_list[[i]][["StudyGroup"]] <- meta.data$StudyGroup[match(seurat_list[[i]]$HATIMID, meta.data$HATIMID)]
  seurat_list[[i]][["Age"]] <- meta.data$meta_age[match(seurat_list[[i]]$HATIMID, meta.data$HATIMID)]
  seurat_list[[i]][["Sex"]] <- meta.data$meta_sex[match(seurat_list[[i]]$HATIMID, meta.data$HATIMID)]
  seurat_list[[i]][["BMI"]] <- meta.data$meta_bmi[match(seurat_list[[i]]$HATIMID, meta.data$HATIMID)]
  seurat_list[[i]][["Metformin"]] <- meta.data$metformin[match(seurat_list[[i]]$HATIMID, meta.data$HATIMID)]
  seurat_list[[i]][["INSTI"]] <- meta.data$INSTI[match(seurat_list[[i]]$HATIMID, meta.data$HATIMID)]
  seurat_list[[i]][["Project"]] <- names(seurat_list[i])
  seurat_list[[i]][["ATVL"]] <- meta.data$upitt_at_hiv_per_mil_cells[match(seurat_list[[i]]$HATIMID, meta.data$HATIMID)]
}

###-------------------------------####
# Import TCR data and link to seurat object
###------------------------------####
path <- "/data/p_koethe_lab/Atlas_AT/CellRanger/Enclone/"
TCR_files <- list.files(path = path, pattern = "enclone.csv")
names(TCR_files) <- c('P5344_CW1_TCR', 'P5344_CW2_TCR', 'P5344_CW3_TCR', 'P5544_CW1_TCR', 'P5544_CW2_TCR', 'P5544_CW3_TCR', 'P5573_CW1_TCR',
                      'P5573_CW2_TCR', 'P5573_CW3_TCR', 'P5657_CW1_TCR', 'P5657_CW3_TCR', 'P5836_CW1_TCR', 'P5836_CW2_TCR',
                      'P5836_CW3_TCR', 'P5877_CW1_TCR', 'P5877_CW2_TCR', 'P5877_CW3_TCR', 'P5903_CW1_TCR', 'P5903_CW2_TCR', 'P5903_CW3_TCR',
                      'P5963_CW1_TCR', 'P5963_CW2_TCR', 'P5963_CW3_TCR')

for (i in 1:length(TCR_files)) {
  TCR <- read.csv(paste(path, TCR_files[[i]], sep = ""))
  rownames(TCR) <- TCR$barcode
  TCR <- TCR %>% dplyr::select(n) %>% dplyr::mutate(TCR_clonotype_size = n) %>% dplyr::select(TCR_clonotype_size)
  seurat_list[[i]] <- AddMetaData(seurat_list[[i]], TCR)
}

###-------------------------------####
# Import BCR data and link to seurat object
###------------------------------####
path <- "/data/p_koethe_lab/Atlas_AT/CellRanger/Enclone_BCR/"
BCR_files <- list.files(path = path, pattern = "enclone_bcr.csv")
names(BCR_files) <- c('P5344_CW1_BCR', 'P5344_CW2_BCR', 'P5344_CW3_BCR', 'P5544_CW1_BCR', 'P5544_CW2_BCR', 'P5544_CW3_BCR', 'P5573_CW1_BCR',
                      'P5573_CW2_BCR', 'P5573_CW3_BCR', 'P5657_CW1_BCR', 'P5657_CW3_BCR', 'P5836_CW2_BCR',
                      'P5836_CW3_BCR', 'P5877_CW1_BCR', 'P5877_CW2_BCR', 'P5877_CW3_BCR', 'P5903_CW1_BCR', 'P5903_CW2_BCR', 'P5903_CW3_BCR',
                      'P5963_CW1_BCR', 'P5963_CW2_BCR', 'P5963_CW3_BCR')

for (i in 1:length(BCR_files)) {
  if(names(seurat_list[i]) != "P5836_CW1") {
    BCR <- read.csv(paste(path, BCR_files[[i]], sep = ""))
    rownames(BCR) <- BCR$barcode
    BCR <- BCR %>% dplyr::select(n, cdr3_aa1, const1, const2) %>% dplyr::mutate(BCR_clonotype_size = n) %>% dplyr::select(BCR_clonotype_size, cdr3_aa1, const1, const2)
    seurat_list[[i]] <- AddMetaData(seurat_list[[i]], BCR)
  }
}

###-------------------------------####
# Save Processed Lanes
###------------------------------####
for (i in 1:length(seurat_list)) {
  print(paste("Saving file: ", names(seurat_list[i]), " in the path ", " /data/p_koethe_lab/Atlas_AT/Analysis/Processed_Final_objects/", ".",  sep = ""))
  seurat_object <- seurat_list[[i]]
  name <- names(seurat_list[i])
  saveRDS(seurat_object, file = paste("/data/p_koethe_lab/Atlas_AT/Analysis/Processed_Final_objects/", name, "_final", ".rds", sep = ""))
}

# 5. INTEGRATION--------------------------------------
# The next step will integrate the 23 lanes to generate one seurat
# object.
