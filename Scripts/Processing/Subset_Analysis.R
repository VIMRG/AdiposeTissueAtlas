##############################################################################################
##---------------- HATIM SINGLE CELL SUBCUTANEOUS ADIPOSE TISSUE ATLAS----------------------##
##------------------------- DATE: 10/29/2021 AUTHOR: -------------------------------##
## DESCRIPTION: The following code will analyze the cell subsets generated from merging all lanes
#  from the project. Each major subset (Lymphoid, Vascular, Endothelial, Myeloid) will be subset
#  and analyzed separately to remove heterotypic doublets and cells that cluster by mitochondrial
#  genes. Each finalized subset will then be merged to generate the final single cell atlas. This
#  script takes advantage of the future package to parallelize certain steps for efficiency.
##############################################################################################

# 1. SET UP_____________________________________________________
###-------------------####
# Set seed
###-------------------####
set.seed(7612) # Reproducibility

###-------------------####
# Set working directory
###-------------------####
setwd('/data/p_koethe_lab/Atlas_AT/Analysis')

###-------------------####
# Load Libraries
###-------------------####
library(Seurat) # Seurat 4.0 beta version includes uwot (note that several packages will also be attached)
library(tidyverse) # collection of packages for manipulating data (ggplot2, dplyr, tidyr, purr, stringr, tibble)
library(patchwork) # tool to assist with plotting
library(Matrix) # dense and sparse matrix manipulations
library(cowplot) # tool to assist with plotting
library(SeuratDisk) # Stores objects for reference assignments
library(stringr) # grab strings
library(data.table) # convert table to dataframe
library(MAST) # DGE
library(harmony) # Integration
library(future) # Paralellize Process

plan("multiprocess", workers = 5) # 5 cpus
options(future.globals.maxSize = 3000 * 1024^2)
options(future.seed = TRUE)

###-------------------####
# Load harmony integrated object
###-------------------####
integrated.data <- readRDS("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Integrated_Harmony_10.30.rds")

# 2. ANALYZE INTEGRATED OBJECT AND SUBSET_____________________________________________________
# Dataset (with doublets) contains 355,392 cells. The strategy will be to subset on major cell
# lineages, remove doublets and dying cells, and annotate the cell types, then merge all subsets back together

pdf("/data/p_koethe_lab/Atlas_AT/Analysis/Integrated_Objects/Main_Annotations/UMAP.pdf")
DimPlot(integrated.data, label = T, raster = FALSE) + NoLegend()
dev.off()

pdf("/data/p_koethe_lab/Atlas_AT/Analysis/Integrated_Objects/Main_Annotations/UMAP_Split.pdf")
DimPlot(integrated.data, label = F, split.by = "Project", ncol = 4, raster = FALSE) + NoLegend()
dev.off()

pdf("/data/p_koethe_lab/Atlas_AT/Analysis/Integrated_Objects/Main_Annotations/Features.pdf")
FeaturePlot(integrated.data, features = c("CLDN5", "COL1A2", "PTPRC"), raster = FALSE)
dev.off()

###-------------------####
# Myeloid Subset
###-------------------####
p1 <- DimPlot(integrated.data, label = T, raster = FALSE) + NoLegend()

# Macrophages
p2 <- FeaturePlot(integrated.data, features = c("TREM2", "RNASE1", "C1QB", "LYVE1"), raster = F) # Clusters 47, 35, 9, 27, 7, 8, 29, 33
p1 + p2

# Monocytes
p2 <- FeaturePlot(integrated.data, features = c("CD14", "VCAN", "FCGR3A", "S100A8"), raster = F) # Clusters 30, 6, 44, 19, 51
p1 + p2

# Dendritic Cells
p2 <- FeaturePlot(integrated.data, features = c("CCR7", "CD1C", "CLEC9A", "LILRA4"), raster = F) # Clusters 22, 39, 45
p1 + p2

# Cycling Myeloid
p2 <- FeaturePlot(integrated.data, features = c("MKI67", "CD1C", "TREM2"), raster = F) # Clusters 9, 29
p1 + p2

# Myeloid Subset
Myeloid <- subset(integrated.data, idents = c(6,7,8,9,19,22,27,29,30,33,35,39,44,45,47,51))

###-------------------####
# Stromal Subset
###-------------------####
p2 <- FeaturePlot(integrated.data, features = c("CCDC80", "LUM", "CFD", "ADIPOQ"), raster = F) # Clusters 1, 16, 3, 0, 14, 18, 20, 25, 41, 52, 36, 43, 33, 48
p1 + p2

Stromal <- subset(integrated.data, idents = c(0,1,3,14,16,18,20,25,33,36,41,43,48,52))

###-------------------####
# T/NK Cell Subset
###-------------------####
p2 <- FeaturePlot(integrated.data, features = c("CD3E", "CD8A", "IL7R", "NKG7"), raster = F) # Clusters 11, 12, 21, 15, 49, 40, 29
p1 + p2

Tcells <- subset(integrated.data, idents = c(11,12,15,21,29,40,49))

###-------------------####
# Vascular Subset
###-------------------####
p2 <- FeaturePlot(integrated.data, features = c("CLDN5", "ACKR1", "CA4", "TAGLN"), raster = F) # Clusters 13,23,24,2,32,26,28,37,53,50,31,4,5,17,38,42,10
p1 + p2

Endo <- subset(integrated.data, idents = c(2,4,5,10,13,17,23,24,26,28,31,32,37,38,42,50,51,53))


###-------------------####
# B Cell Subset
###-------------------####
p2 <- FeaturePlot(integrated.data, features = c("CD79A", "JCHAIN", "MS4A1"), raster = F) # Clusters 34,46,51,52
p1 + p2

Bcells <- subset(integrated.data, idents = c(34,46,51,52))

###-------------------####
# Erythropoietic Cluster
###-------------------####
p2 <- FeaturePlot(integrated.data, features = c("HBB"), raster = F) # Cluster 47
p1 + p2

# 2. T & NK Cell SUBSET_____________________________________________________
#Initial: 43,365 cells (including doublets and other cell types)
#==========================================
# Find variable features, rescale, and run PCA
#==========================================
DefaultAssay(Tcells) <- "RNA"
Tcells <- Tcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

#==========================================
# Run Harmony Correction on RNA PCA matrix (Note:Stochastic Process)
#==========================================
Tcells <- RunHarmony(Tcells, "Project", assay.use = "RNA", max.iter.harmony = 20)

#==========================================
# Run dimensional reduction and clustering on harmony corrected PC
# Note: This is stochastic process so results may not be exactly the same
#==========================================
ElbowPlot(Tcells, ndims = 40)
DimHeatmap(Tcells, cells = 500, balanced = T, dims = 15:20) # 20 clusters

Tcells <- Tcells %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 2.0)

#==========================================
# Evaluate for Doublets (Stromal doublet: 19; Vascular Doublet: 12,26,15,20; Myeloid doublet: 27,13,28,23,16,17,25)
#==========================================
p1 <- DimPlot(Tcells, label = T) + NoLegend()
p2 <- FeaturePlot(Tcells, features = c("CLDN5", "COL1A1", "LYZ", "CCDC80"))
p1 + p2

DimPlot(Tcells, group.by = "Gross_Annotation") # Same clusters
DimPlot(Tcells, group.by = "DoubletFinder") # Similar clusters

#==========================================
# Round 2: Remove doublet clusters 12,13,15,16,17,19,20,23,25,26,27,28 (29 probably doublet but will keep for now)
#==========================================
Tcells <- subset(Tcells, idents = c(12,13,15,16,17,19,20,23,25,26,27,28), invert = T)
DefaultAssay(Tcells) <- "RNA"
Tcells <- Tcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction (Note:Stochastic process)
Tcells <- Tcells %>% RunHarmony("Project", assay.use = "RNA", max.iter.harmony = 20)

# Downstream Analysis
Tcells <- Tcells %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:25) %>%
        FindClusters(resolution = 2.0)

# Evaluate For Doublets/Low Quality Cells: Myeloid/Stromal- 26,; HBB - 23; Junk - 5; Endo Doublet: 17 (SPARCL1, CD36, GNG11, CAV1, TYROBP, GSN, ADIRF) (Note: 22 suspicious for Myeloid but not clear)
p1 <- DimPlot(Tcells, label = T) + NoLegend()
p2 <- FeaturePlot(Tcells, features = c("CLDN5", "COL1A1", "LYZ", "CCDC80"))
p1 + p2

#==========================================
# Round 3 (34,684): Remove doublet clusters 5,17,23,26
#==========================================
Tcells <- subset(Tcells, idents = c(5,17,23,26), invert = T)
DefaultAssay(Tcells) <- "RNA"
Tcells <- Tcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction (Note:Stochastic process)
Tcells <- Tcells %>% RunHarmony("Project", assay.use = "RNA", max.iter.harmony = 20)

# Downstream Analysis
Tcells <- Tcells %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:25) %>%
        FindClusters(resolution = 2.0)

# Evaluate For Doublets/Low Quality Cells: Junk - 26; Myeloid - 24 (this was cluster 22)
p1 <- DimPlot(Tcells, label = T) + NoLegend()
p2 <- FeaturePlot(Tcells, features = c("CLDN5", "COL1A1", "LYZ", "CCDC80"))
p1 + p2

#==========================================
# Round 4 (31,502 cells): Remove cluster 24 & 26
#==========================================
# Would like to remove some TCR genes that are driving clustering of clonotypes (gamma delta TCR, MAIT TRAV1-2 okay since it will differentiate those cells and 
# not doing simultaneous TCR:GEX correlations here that would be problematic)
# Clustering driven by metallothionein genes. Will regress out since I want to characterize percentage of CD4 and CD8 subsets and 
# these are likely from multiple subsets (naive, TCM, TEM, TEMRA). May be interesting to look at in another context (some data
# shows association of CD8 T cells with diabetes e.g. Vijay et. al. 2020) but these occur in all cell lineages (myeloid, stromal) and not different by disease states.

Tcells <- subset(Tcells, idents = c(24,26), invert = T)
DefaultAssay(Tcells) <- "RNA"
TCR <- list("TRBV5-6","TRGV2", "TRAV14DV4", "TRAV13-1")
MTs <- list("MT1G", "MT1E","MT1X", "MT1F", "MT2A")
Tcells <- AddModuleScore(Tcells, features = MTs, name = "MTs")
Tcells <- AddModuleScore(Tcells, features = TCR, name = "TRAV")

Tcells <- Tcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "TRAV1", "MTs1")) %>%
        RunPCA()

# Run Harmony Correction (Note:Stochastic process)
Tcells <- RunHarmony(Tcells, "Project", assay.use = "RNA", max.iter.harmony = 20)

# Downstream Analysis (Take 23 because after that, more mitochondrial genes)
Tcells <- Tcells %>% RunUMAP(reduction = 'harmony', dims = 1:23) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:23) %>%
        FindClusters(resolution = 2.0)
        
# Evaluate For Doublets/Low Quality Cells
p1 <- DimPlot(Tcells, label = T) + NoLegend()
p2 <- FeaturePlot(Tcells, features = c("CLDN5", "COL1A1", "LYZ", "CCDC80"))
p1 + p2

#==========================================
# Final T Cells (31,140 cells): No further processing needed!
#==========================================
# Note: Some ambient RNA contamination driving clustering. Final dataset excludes T and NK cells with significant ambient contamination
# though they likely represent true lymphoid cells. Contamination is mainly from myeloid cells.
# Now will annotate CD4 & CD8 naive, TCM, TEM, TEMRA & Senescent; mature NK, CD16+ NK, and immature NK, ILCs, and CD4 Treg

# Get T cell markers
Tcells.markers <- FindAllMarkers(Tcells, only.pos = T, min.pct = 0.25, assay = "RNA")
write.csv(Tcells.markers, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcells_RNA_markers_11.1.csv')

# RNA Annotation markers
# CD4 naive: CCR7, SELL, LDHB, ribosomal proteins
# CD8 naive: CD8A, CCR7, SELL, NELL2, ribosomal proteins
# CD4 TCM: CCR7, CD27, GPR183, LTB, ribosomal proteins
# CD8 TCM: CD8A, CCR7, CD27, GZMA
# CD4 TEM/TR: IL32, VIM, ALOXA5, COTL1, CCL5
# CD8 TEM/TR: CD8A, VIM, ALOXA5, COTL1, CCL5, IFNG
# CD4 TEMRA & Senescent: CCL5, NKG7, GNLY, PRF1
# CD8 TEMRA & Senescent: CD8A, CCL5, NKG7, GNLY, PRF1
# CD4 Regulatory: CTLA, FOXP3
# MAIT: IL7R, CCR6
# Gamma Delta: KLRC2, TRDV1, TRDC, CCL5
# Cycling T & NK Cells: MKI67
# ILC: KIT, SPINK2, IL7R, SELL, XCL1
# CD57+ Mature NK: FCEGR3A, CCL5, PRF1, GNLY
# CD16+ mature NK: FCER1G, FCGR3A, CCL5, PRF1
# CD56+ immature NK: XCL1, XCL2, SELL
# No iNKT cells in data, unable to classify exhausted/activated clearly from transcriptome

# Cluster Assignment based on RNA
# 0: CD57+ CD16+ Mature NK Cells
# 1: CD4 TCM
# 2: CD16+ Mature NK Cells
# 3: CD8 TCM
# 4: Cytotoxic CD8 TCM/TEMRA
# 5: CD8 TEM
# 6: CD8 Senescent & TEMRA
# 7: CD16+ Mature NK Cells
# 8: Gamma Delta
# 9: CD4 TCM
# 10: CD8 TEM
# 11: CD56+ Immature NK
# 12: CD4 TEM
# 13: CD4 Naive
# 14: CD8 Naive
# 15: CD4 TEM/TR stressed or activated
# 16: CD8 TEM
# 17: Unclear
# 18: Regulatory CD4
# 19: Cycling T & NK Cells
# 20: MAIT
# 21: CD8 TEMRA & Senescent
# 22: ILC
# 23: CD4 or 8 TEM


# ADT Annotation markers
DefaultAssay(Tcells) <- "ADT"
p1 <- DimPlot(Tcells, label = T) + NoLegend()
p2 <- FeaturePlot(Tcells, features = c("CD4", "CD8", "CD45RA", "CD27"), min.cutoff = 'q10', max.cutoff = 'q90')
p3 <- VlnPlot(Tcells, features = c("CD4", "CD8", "CD45RA", "CD27"), pt.size = 0)
p4 <- FeaturePlot(Tcells, features = c("CD4", "CD8", "CD45RA", "CD57"), min.cutoff = 'q10', max.cutoff = 'q90')
p5 <- FeaturePlot(Tcells, features = c("CD16", "CD56-recom", "CD45RA", "CD57"), min.cutoff = 'q10', max.cutoff = 'q90')
p6 <- VlnPlot(Tcells, features = c("CD16", "CD56-recom", "CD45RA", "CD57"), pt.size = 0)
p1 + p2

# Cluster Assignment based on ADT (agreement with RNA)
# 0: CD57+ CD16+ Mature NK Cells (Y)
# 1: CD4+CD45RO+CD27+: CD4 TCM (Y)
# 2: CD16+CD57- Mature NK Cells (Y)
# 3: CD8+CD45RO+CD27+: CD8 TCM (Y)
# 4: CD4+CD8+CD27-CD45RA+: Mixed CD4 and CD8 TEMRA & Senescent (N)
# 5: CD8+CD45RO+CD27+/-: CD8 TEM (Y)
# 6: CD8+CD45RA+CD27-: CD8 Senescent & TEMRA (Y)
# 7: CD16+CD57+/-:Mixed CD57+/CD57- Mature NK Cells (N)
# 8: Gamma Delta (Y)
# 9: CD4+CD45RO+CD27+: CD4 TCM (Y)
# 10: CD8+CD45RA+CD27-: CD8 TEMRA & Senescent (N)
# 11: CD56+CD16-: Immature NK (Y)
# 12: CD4+CD45RO+CD27-: CD4 TEM (Y)
# 13: CD4+CD45RA+CD27+: CD4 Naive (Y)
# 14: CD8+CD45RA+CD27+: CD8 Naive (Y)
# 15: CD4+CD8+CD45RO+CD27-: CD4 TEM/TR stressed or activated (Y)
# 16: CD8+CD45RO+CD27-: CD8 TEM (Y)
# 17: CD8+/-CD45RA+CD27-: Mixed CD8 TEMRA & GD (N)
# 18: CD4+CD45RO+CD27+CTLA4+: Regulatory CD4 (Y)
# 19: Cycling T & NK Cells (Y)
# 20: MAIT (Y)
# 21: CD8+CD45RA+CD27-:CD8 TEMRA & Senescent (Y)
# 22: ILC (Y)
# 23: CD4+CD45RO+CD27-: CD4 TEM (N)

# Resolve Disagreements
# Cluster 4 (contains smaller subset of CD4 TEMRA & Senescent)
clust4 <- subset(Tcells, idents = 4)
RidgePlot(clust4, features = c("CD4", "CD8")) # Use plot to identify demarkation of CD8 antibody (around 1.1)
Idents(Tcells, cells = WhichCells(clust4, expression = CD8 >=1.2)) <- "CD8 TEMRA & Senescent"
Idents(Tcells, cells = WhichCells(clust4, expression = CD8 < 1.2)) <- "CD4 TEMRA & Senescent"

# Cluster 8 and 17 (This is composed of both CD8 TEMRA & Gamma Delta)
clust8 <- subset(Tcells, idents = c(8,17))
clust8[['TCR']] <- ifelse(is.na(clust8$TCR_clonotype_size), "no", "yes") # Exploiting TCR data, gamma delta should not have TCRA or TCRB, these are likely CD8 TEMRA
Idents(clust8) <- "TCR"
Idents(Tcells, cells = WhichCells(clust8, idents = 'yes')) <- "CD8 TEMRA & Senescent"
Idents(Tcells, cells = WhichCells(clust8, idents = 'no')) <- "Gamma Delta"

Idents(Tcells, cells = WhichCells(Tcells, idents = c(1,9))) <- "CD4 TCM"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(12,23,15))) <- "CD4 TEM"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(13))) <- "CD4 Naive"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(14))) <- "CD8 Naive"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(20))) <- "MAIT"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(18))) <- "CD4 Regulatory"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(3))) <- "CD8 TCM"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(5,16))) <- "CD8 TEM"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(6,10,21))) <- "CD8 TEMRA & Senescent"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(19))) <- "Cycling T & NK"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(22))) <- "ILC"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(11))) <- "Immature NK"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(0,7))) <- "CD57+ mNK"
Idents(Tcells, cells = WhichCells(Tcells, idents = c(2))) <- "CD16+ mNK"

# DotPlot (in no particular order)
DefaultAssay(Tcells) <- "RNA"
DotPlot(Tcells, features = c("IL7R", "SELL", "LEF1", "CCR7", "LTB", "LDHB", "GPR183", "COTL1", "ALOX5AP", "CCL5", "NKG7", "GNLY", "PRF1", "FOXP3", "CD8A", "CCR6", "TRDV1", "KLRC2", 
"FCER1G", "FCGR3A", "XCL1", "KIT", "GATA2", "SPINK2", "MKI67")) + RotatedAxis()

# Save Annotations
Tcells[['Manual_Annotation']] <- Idents(Tcells)

# Proportion
write.table(prop.table(table(Tcells$HATIMID, Tcells$Manual_Annotation), margin = 1) * 100, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcells.Prop.by.HATIMID.11.1.txt')

# Save T Cells Subcluster
saveRDS(Tcells, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_10.30.rds')


# 3. CD4 SUBCLUSTERING_____________________________________________________
# Need to Separate CD4 and CD8 T cells using CD8 CITE-seq antibody. 
#==========================================
# Round 1: CD4 subclustering
#==========================================
Tcells <- readRDS("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_10.30.rds")

# Examine ADT Expression of CD4 and CD8
DefaultAssay(Tcells) <- "ADT"
RidgePlot(Tcells, features = c('CD4', 'CD8')) # CD4 expression essentially non-existent on NK and GD cells as expected. 

# Subset on T cells (20,493 cells)
CD4 <- subset(Tcells, idents = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Regulatory", "CD4 TEMRA & Senescent", "CD8 TEMRA & Senescent", "CD8 TCM", "CD8 Naive", "CD8 TEM"))
RidgePlot(CD4, features = c("CD8", "CD4")) # Will keep all cells expressing < 1.1 CD8
CD4 <- subset(CD4, cells = WhichCells(CD4, expression = CD8 < 1.1))
RidgePlot(CD4, features = c("CD8", "CD4")) # Will keep all cells expressing < 1.1 CD8, however likely some NK cells included

# Rerun Harmony Integration
DefaultAssay(CD4) <- "RNA"

CD4 <- CD4 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "TRAV1")) %>%
        RunPCA()

# Run Harmony Correction
CD4 <- RunHarmony(CD4, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(CD4, ndims = 30)
DimHeatmap(CD4, cells = 500, balanced = T, dims = 15:20)

# Run DownStream clustering
CD4 <- RunUMAP(CD4, reduction = 'harmony', dims = 1:20)
CD4 <- FindNeighbors(CD4, reduction = 'harmony', dims = 1:20)
CD4 <- FindClusters(CD4, resolution = 1.0)

# View Plots
DefaultAssay(CD4) <- "ADT"
p1 <- DimPlot(CD4, label = T) + NoLegend()
p2 <- FeaturePlot(CD4, features = c("CD45RA", "CD27", "CD57", "CD4"), min.cutoff = 'q15', max.cutoff = 'q90')
p3 <- VlnPlot(CD4, features = c("CD4", "CD45RA", "CD27", "CD57", "CD8"), pt.size = 0)
p1 + p3
CD4.markers <- FindAllMarkers(CD4, only.pos = T, min.pct = 0.25, assay = "RNA")

# Cluster Identity
# 0: CD4 TCM (CD45RO+CD27+)
# 1: CD4 TEMRA & Senescent (CD45RA+/-CD27-CD57+)
# 2: CD4 TEM/Transitional (CD45RO+CD27+/-)
# 3: CD4 Naive (CD45RA+CD27+CCR7+)
# 4: CD4 TEM (CD45RO+CD27-)
# 5: CD4 TEM/Transitional (CD45RO+CD27+/-)
# 6: Contaminating NK/CD8 cluster (no CD4 expression, reduced TCRs)
# 7: CD4 TEM (CD45RO+CD27-)
# 8: CD4 TCM (CD45RO+CD27+)
# 9: CD4 Regulatory (CD45RO+CD27+)
# 10: CD4 TEM (CD45RO+CD27-)

#==========================================
# Round 2: CD4 subclustering - Remove contaminating NK cell cluster
#==========================================
DefaultAssay(CD4) <- "RNA"
CD4 <- subset(CD4, idents = 6, invert = T)

# Scale and PCA
CD4 <- CD4 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "TRAV1", "MTs1")) %>%
        RunPCA()

# Run Harmony Correction
CD4 <- RunHarmony(CD4, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(CD4, ndims = 30)
DimHeatmap(CD4, cells = 500, balanced = T, dims = 15:20)

# Run DownStream clustering
CD4 <- RunUMAP(CD4, reduction = 'harmony', dims = 1:20)
CD4 <- FindNeighbors(CD4, reduction = 'harmony', dims = 1:20)
CD4 <- FindClusters(CD4, resolution = 1.0)

# View Plots
DefaultAssay(CD4) <- "ADT"
p1 <- DimPlot(CD4, label = T) + NoLegend()
p2 <- FeaturePlot(CD4, features = c("CD45RA", "CD27", "CD57", "CD4"), min.cutoff = 'q15', max.cutoff = 'q90')
p3 <- VlnPlot(CD4, features = c("CD4", "CD45RA", "CD27", "CD57","CD8"), pt.size = 0)
p1 + p3
CD4.markers <- FindAllMarkers(CD4, only.pos = T, min.pct = 0.25, assay = "RNA")

# Cluster Identity
# 0: CD4 TCM (CD45RO+CD27+), CCR7, LTB, IL7R, AQP3,
# 1: CD4 Naive (CD45RA+CD27+), SELL, LEF1, CCR7
# 2: CD4 TEM/Transitional (CD45RO+CD27+/-), VIM, FOSB, ZFP36, ALOX5AP
# 3: CD4 TEMRA & Senescent (CD45RA+CD27-CD57+), GNLY, NKG7, GZMH, GZMB, PRF1, CCL5
# 4: CD4 TCM (CD45RO+CD27+/-), GPR183, IL7R,
# 5: CD4 TEM (CD45RO+CD27-)
# 6: CD4 TEM (CD45RO+CD27-)
# 7: CD4 Regualtory (CD45RO+CD27+), FOXP3
# 8: CD4 TCM (CD45RO+CD27+), CCR7
# 9: CD4 Regulatory (CD45RO+CD27+)
# 10: Contaminating CD8
# 11: CD4 TEM (CD45RO+CD27-)
# 12: Contaminating CD8

#==========================================
# Round 3: CD4 subclustering - Remove contaminating CD8 clusters
#==========================================
DefaultAssay(CD4) <- "RNA"
CD4 <- subset(CD4, idents = c(10,12), invert = T)

# Scale and PCA
CD4 <- CD4 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "TRAV1", "MTs1")) %>%
        RunPCA()

# Run Harmony Correction
CD4 <- RunHarmony(CD4, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(CD4, ndims = 30)
DimHeatmap(CD4, cells = 500, balanced = T, dims = 15:20)

# Run DownStream clustering
CD4 <- RunUMAP(CD4, reduction = 'harmony', dims = 1:15)
CD4 <- FindNeighbors(CD4, reduction = 'harmony', dims = 1:15)
CD4 <- FindClusters(CD4, resolution = 1.5)

# View Plots
DefaultAssay(CD4) <- "ADT"
p1 <- DimPlot(CD4, label = T) + NoLegend()
p2 <- FeaturePlot(CD4, features = c("CD45RA", "CD27", "CD57", "CD4"), min.cutoff = 'q15', max.cutoff = 'q90')
p3 <- VlnPlot(CD4, features = c("CD4", "CD45RA", "CD27", "CD57","CD8"), pt.size = 0)
p1 + p3
CD4.markers <- FindAllMarkers(CD4, only.pos = T, min.pct = 0.25, assay = "RNA")

#==========================================
# Round 4: CD4 subclustering - Remove cluster 14 (Myeloid genes expressed, probably ambient RNA contamination rather than true doublet)
#==========================================
CD4 <- subset(CD4, idents = 14, invert = T)

# Recluster without rerunning PCA or UMAP
DefaultAssay(CD4) <- "RNA"
CD4 <- RunUMAP(CD4, reduction = 'harmony', dims = 1:15)
CD4 <- FindNeighbors(CD4, reduction = 'harmony', dims = 1:15)
CD4 <- FindClusters(CD4, resolution = 1.5)

DefaultAssay(CD4) <- "ADT"
p1 <- DimPlot(CD4, label = T) + NoLegend()
p2 <- FeaturePlot(CD4, features = c("CD45RA", "CD27", "CD57", "CD4"), min.cutoff = 'q15', max.cutoff = 'q90')
p3 <- VlnPlot(CD4, features = c("CD4", "CD45RA", "CD27", "CD57","CD8"), pt.size = 0)
p1 + p3
CD4.markers <- FindAllMarkers(CD4, only.pos = T, min.pct = 0.25, assay = "RNA")

# Add Labels
# CD45RA+: 0,3,10,15
# CD45RO+: 1,2,4,5,6,7,8,9,11,12,13,14
# CD27+: 0,1,2,4,8,13
# CD27-: 3,5,7,9,10,12,15
# CD27 intermediate: 14,11,6

# Assignment High-Confidence
# Naive: Cluster 0
# TCM: Cluster 1, 4
# TEM: 5,7,9,12
# Treg: Cluster 8
# CTL: 3,10,15

# Lower Confidence (14)
# TCM:13, 2 (some cells in cluster 2 may actually be TEM)
# TEM: 6, 11

Idents(CD4) <- "seurat_clusters"
Idents(CD4, cells = WhichCells(CD4, idents = c(0))) <- "CD4 Naive"
Idents(CD4, cells = WhichCells(CD4, idents = c(1,2,4,13))) <- "CD4 TCM"
Idents(CD4, cells = WhichCells(CD4, idents = c(5,7,9,12,6,11))) <- "CD4 TEM"
Idents(CD4, cells = WhichCells(CD4, idents = c(3,10,15))) <- "CD4 TEMRA & Senescent"
Idents(CD4, cells = WhichCells(CD4, idents = c(8))) <- "CD4 Regulatory"
Idents(CD4, cells = WhichCells(CD4, idents = c(14))) <- "CD4 Transitional"

CD4[['Tcell_subcluster']] <- "CD4"
CD4[['Fine_Annotation']] <- Idents(CD4)

# DotPlot
VlnPlot(CD4, features = c("CD45RA", "CD27", "CD57"), pt.size = 0)
DefaultAssay(CD4) <- "RNA"
DotPlot(CD4, features = c("SELL", "LEF1", "CCR7", "GPR183", "LTB", "LDHB", "ALOX5AP", "ANXA1", "CCL5", "NKG7", "PRF1", "FOXP3", "CD8A")) + RotatedAxis()

saveRDS(CD4, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/CD4_RNA_11.1.rds')

# Proportion Cells
write.table(prop.table(table(CD4$HATIMID, CD4$Fine_Annotation), margin = 1) * 100, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/CD4.Prop.by.HATIMID.11.1.txt')
write.table(table(CD4$HATIMID), file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/CD4.Num.by.HATIMID.11.1.txt')

# 4. CD8 SUBCLUSTERING_____________________________________________________
# Need to Separate CD8 T cells using CD8 CITE-seq antibody. 
#==========================================
# Round 1: CD8 subclustering
#==========================================
Tcells <- readRDS("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_10.30.rds")

# Examine ADT Expression of CD4 and CD8
DefaultAssay(Tcells) <- "ADT"
RidgePlot(Tcells, features = c('CD4', 'CD8')) # Some CD8 expression on NK and GD cells as expected. 

# Subset on T cells (28,720 cells: Will include gamma delta and mature NK as these have some cells expressing TCRA & TCRB, may actually be CD8 T cells given similarity of TEMRA and gamma delta/mNK
# gene expression profile)
DefaultAssay(CD8) <- "ADT"
CD8 <- subset(Tcells, idents = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 TEMRA & Senescent", "CD8 TEMRA & Senescent", "CD8 TCM", "CD8 Naive", "CD8 TEM", "Gamma Delta", "CD57+ mNK", "CD16+ mNK"))
RidgePlot(CD8, features = c("CD8", "CD4")) # Will keep all cells expressing >= 1.1 CD8, however likely some NK and gamma delta cells included
CD8 <- subset(CD8, cells = WhichCells(CD8, expression = CD8 >= 1.1))
RidgePlot(CD8, features = c("CD8", "CD4")) # Very few CD4 cell crossover except TEM compartment

# Rerun Harmony Integration
DefaultAssay(CD8) <- "RNA"

CD8 <- CD8 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "TRAV1")) %>%
        RunPCA()

# Run Harmony Correction
CD8 <- RunHarmony(CD8, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(CD8, ndims = 30)
DimHeatmap(CD8, cells = 500, balanced = T, dims = 15:20)

# Run DownStream clustering
CD8 <- RunUMAP(CD8, reduction = 'harmony', dims = 1:20)
CD8 <- FindNeighbors(CD8, reduction = 'harmony', dims = 1:20)
CD8 <- FindClusters(CD8, resolution = 1.0)

# View Plots
DefaultAssay(CD8) <- "ADT"
p1 <- DimPlot(CD8, label = T) + NoLegend()
p2 <- FeaturePlot(CD8, features = c("CD45RA", "CD27", "CD57", "CD8"), min.cutoff = 'q15', max.cutoff = 'q90')
p3 <- VlnPlot(CD8, features = c("CD4", "CD45RA", "CD27", "CD57", "CD8"), pt.size = 0)
p1 + p3
CD8.markers <- FindAllMarkers(CD8, only.pos = T, min.pct = 0.25, assay = "RNA")

# Cluster Identity
# 0: CD8 TEM (CD45RO+CD27+/-), GZMK, COTL1, DUSP2, ALOX5AP, VIM, IL32, CCL5
# 1: CD8 TCM (CD45RO+CD27+), IL7R, GZMK, CXCR3, LDHA, ribosomal
# 2: CD8 TEMRA & Senescent (CD45RA+CD27-CD57+), NKG7, GNLY, GZMH, GZMB
# 3: CD8 TCM (CD45RO+CD27+), GPR183, LTB, ILR7R, CCR7, AQP3, ribosomal
# 4: CD8 TEMRA & Senescent (CD45RA+CD27-CD57+), GZMH, GZMB, NKG7, B2M, GNLY, PRF1, ZNF683
# 5: CD8 Naive (CD45RA+CD27+), CCR7, SEEL, LEF1, MAL
# 6: CD8 TEMRA & Senescent + NK (CD45RA+CD27-CD57+),KLRC2, KLRC3, TYROBP, FCGR3A, NKG7 -> mNK cells
# 7: Mostly NK (CD8lowCD45RA+CD27-CD57+), FCER1G, TYROBP, FCGR3A -> mNK cells
# 8: CD8 TEMRA & Senescent (CD45RA+CD27-CD57+), TXNIP, GIMAP7, KLRG1, GZMA, CCL5
# 9: CD8 TEMRA & Senescent (CD45RA+CD27-CD57+), CCL4, FOS, GZMH, 
# 10: CD8 TEM (CD45RO+CD27-), 
# 11: Mixed CD8 TEMRA & Senescent (CD45RO+CD27-CD57+) -> gamma delta
# 12: CD8 TEMRA & Senescent (CD45RA+CD27-CD57+)
#13: Transcriptionally like TEM, but markers show CD45RA expression

# Subcluster on NK and gamma delta (see if CD8 TEMRA separate)
other <- subset(CD8, idents = c(7,6,11))
p1 <- DimPlot(other, group.by = "TCR_clonotype_size") + NoLegend() # highly clonal cells mixed in mostly gamma delta, less NK
other[['TCR_present']] <- ifelse(is.na(other$TCR_clonotype_size), 'no', 'yes')
DefaultAssay(other) <- "RNA"

other <- other %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "TRAV1")) %>%
        RunPCA()
        
other <- RunHarmony(other, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(other, ndims = 30)
other <- other %>% RunUMAP(reduction = 'harmony', dims = 1:10) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:10) %>%
        FindClusters(resolution = 1.0)
# View Plots
DefaultAssay(other) <- "ADT"
p1 <- DimPlot(other, label = T) + NoLegend()
p2 <- DimPlot(other, group.by = 'TCR_present')
p3 <- VlnPlot(other, features = c("CD4", "CD45RA", "CD27", "CD57", "CD8"), pt.size = 0)
p1 + p3 # Very little CD8 expression and TCR sequences from clusters 1,3,8 (comprised of mNK)
p1 + p2

#==========================================
# Round 2: CD8 subclustering - exclude mature NK (very little CD8 expression and TCR clones)
#==========================================
Idents(CD8) <- "Manual_Annotation"
CD8 <- subset(CD8, idents = c("CD57+ mNK", "CD16+ mNK"), invert = T)

# Rerun Harmony Integration
DefaultAssay(CD8) <- "RNA"

CD8 <- CD8 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "TRAV1")) %>%
        RunPCA()

# Run Harmony Correction
CD8 <- RunHarmony(CD8, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(CD8, ndims = 30)
DimHeatmap(CD8, cells = 500, balanced = T, dims = 15:20)

# Run DownStream clustering
CD8 <- RunUMAP(CD8, reduction = 'harmony', dims = 1:20)
CD8 <- FindNeighbors(CD8, reduction = 'harmony', dims = 1:20)
CD8 <- FindClusters(CD8, resolution = 1.0)

# View Plots
DefaultAssay(CD8) <- "ADT"
p1 <- DimPlot(CD8, label = T) + NoLegend()
p2 <- FeaturePlot(CD8, features = c("CD45RA", "CD27", "CD57", "CD8"), min.cutoff = 'q15', max.cutoff = 'q90')
p3 <- VlnPlot(CD8, features = c("CD4", "CD45RA", "CD27", "CD57", "CD8"), pt.size = 0)
p1 + p3
CD8.markers <- FindAllMarkers(CD8, only.pos = T, min.pct = 0.25, assay = "RNA")

# Cluster Identity
# 0: CD8 TEM
# 1: CD8 TCM
# 2: CD8 TEMRA & Senescent
# 3: CD8 TEMRA & Senescent
# 4: Mixed
# 5: CD8 TCM
# 6: CD8 Naive
# 7: CD8 TEMRA & Senescent
# 8: CD8 TCM
# 9: CD8 TEMRA & Senescent
# 10: CD8 TEM 
# 11: Odd Transcriptional TEM but somewhat higher expression CD45RA
# 12: CD8 TEMRA & Senescent

# Subcluster 4 (Mixed cluster)
clust4 <- subset(CD8, idents = 4)
clust4[['TCR_present']] <- ifelse(is.na(clust4$TCR_clonotype_size), 'no', 'yes')
table(clust4$Manual_Annotation, clust4$TCR_present) # Gamma Delta contain 0 TCRA and TCRB sequences. Therefore will remove them.

#==========================================
# Round 3: CD8 subclustering - exclude Gamma Delta cells (labeled from T cell data set)
#==========================================
Idents(CD8) <- "Manual_Annotation"
CD8 <- subset(CD8, idents = c("Gamma Delta"), invert = T)

# Rerun Harmony Integration
DefaultAssay(CD8) <- "RNA"

CD8 <- CD8 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "TRAV1")) %>%
        RunPCA()

# Run Harmony Correction
CD8 <- RunHarmony(CD8, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(CD8, ndims = 30)
DimHeatmap(CD8, cells = 500, balanced = T, dims = 15:20)

# Run DownStream clustering
CD8 <- RunUMAP(CD8, reduction = 'harmony', dims = 1:20)
CD8 <- FindNeighbors(CD8, reduction = 'harmony', dims = 1:20)
CD8 <- FindClusters(CD8, resolution = 1.0)

# View Plots
DefaultAssay(CD8) <- "ADT"
p1 <- DimPlot(CD8, label = T) + NoLegend()
p2 <- FeaturePlot(CD8, features = c("CD45RA", "CD27", "CD57", "CD8"), min.cutoff = 'q15', max.cutoff = 'q90')
p3 <- VlnPlot(CD8, features = c("CD4", "CD45RA", "CD27", "CD57", "CD8"), pt.size = 0)
p1 + p3
CD8.markers <- FindAllMarkers(CD8, only.pos = T, min.pct = 0.25, assay = "RNA")

# Add Labels
# CD45RA+: 1,4,5,6,7,8,11,12
# CD45RO+: 0,2,3,9,10
# CD27+: 2,3,5,9
# CD27-: 1,4,6,7,8,12
# CD27 intermediate: 0,10,11

# Assignment High-Confidence
# Naive: Cluster 5
# TCM: Cluster 2,3,9
# TEM: 0,10,11
# CTL: 1,4,6,7,8,12

Idents(CD8) <- "seurat_clusters"
Idents(CD8, cells = WhichCells(CD8, idents = c(5))) <- "CD8 Naive"
Idents(CD8, cells = WhichCells(CD8, idents = c(2,3,9))) <- "CD8 TCM"
Idents(CD8, cells = WhichCells(CD8, idents = c(0,10,11))) <- "CD8 TEM"
Idents(CD8, cells = WhichCells(CD8, idents = c(1,4,6,7,8,12))) <- "CD8 TEMRA & Senescent"

CD8[['Tcell_subcluster']] <- "CD8"
CD8[['Fine_Annotation']] <- Idents(CD8)

# DotPlot
VlnPlot(CD8, features = c("CD45RA", "CD27", "CD57"), pt.size = 0)
DefaultAssay(CD8) <- "RNA"
DotPlot(CD8, features = c("SELL", "LEF1", "CCR7", "GPR183", "LTB", "LDHB", "ALOX5AP", "ANXA1", "CCL5", "NKG7", "PRF1", "FOXP3", "CD8A")) + RotatedAxis()

saveRDS(CD8, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/CD8_RNA_11.1.rds')

# Proportion Cells
write.table(prop.table(table(CD8$HATIMID, CD8$Fine_Annotation), margin = 1) * 100, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/CD8.Prop.by.HATIMID.11.1.txt')
write.table(table(CD8$HATIMID), file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/CD8.Num.by.HATIMID.11.1.txt')

# 5. MYELOID SUBCLUSTERING_____________________________________________________

#==========================================
# Round 1: Myeloid Clustering - 87,105 cells
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Myeloid) <- "RNA"
Myeloid <- Myeloid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction on RNA PCA matrix (Note:Stochastic Process)
Myeloid <- RunHarmony(Myeloid, "Project", assay.use = "RNA", max.iter.harmony = 20)

# Run dimensional reduction and clustering on harmony corrected PC
ElbowPlot(Myeloid, ndims = 40)
DimHeatmap(Myeloid, cells = 500, balanced = T, dims = 20:25)

Myeloid <- Myeloid %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:25) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Stromal doublet: 38,9,18,22,24; Vascular Doublet: 13,23,36,27; Lymphoid doublet: 34,37,15; Junk - 29,validated by DGE)
p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- FeaturePlot(Myeloid, features = c("MS4A1", "JCHAIN", "CLDN5", "CCDC80", "CD3E", "NKG7"))
p1 + p2

DimPlot(Myeloid, group.by = "Gross_Annotation") # not very helpful
DimPlot(Myeloid, group.by = "DoubletFinder") # not very helpful
FeaturePlot(Myeloid, features = c("CLEC9A", "CD1C", "CCR7", "MKI67")) # Make sure we are keeping rarer cell types

#==========================================
# Round 2 (87,105 cells): Myeloid Clustering excluding clusters 9,13,15,18,22,23,24,27,29,34,36,37,38 
#==========================================
Myeloid <- subset(Myeloid, idents = c(9,13,15,18,22,23,24,27,29,34,36,37,38), invert = T)

# Find variable features, rescale, and run PCA
DefaultAssay(Myeloid) <- "RNA"
Myeloid <- Myeloid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction on RNA PCA matrix (Note:Stochastic Process)
Myeloid <- RunHarmony(Myeloid, "Project", assay.use = "RNA", max.iter.harmony = 20)

# Run dimensional reduction and clustering on harmony corrected PC
ElbowPlot(Myeloid, ndims = 40)
DimHeatmap(Myeloid, cells = 500, balanced = T, dims = 20:25)

Myeloid <- Myeloid %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:25) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Stromal doublet: ; Vascular Doublet: ; Lymphoid doublet: ; Junk - )
p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- FeaturePlot(Myeloid, features = c("MS4A1", "JCHAIN", "CLDN5", "CCDC80", "CD3E", "NKG7"))
p1 + p2

# View Plots
p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- DotPlot(Myeloid, features = c("TREM2", "CD9", "RNASE1", "LYVE1", "C1QB", "LILRA4", "CLEC9A", "CD1C", "CCR7", "CD14", "VCAN", "S100A8", "FCGR3A", "WARS", "CXCL2", "CXCL3", "TXNIP", "CLEC10A", "MKI67")) + RotatedAxis()
p2

# LAMS (1,3,5,8,11,12,15,17,19,20,23,24,27)
p3 <- VlnPlot(Myeloid, features = c("TREM2", "CD9", "FABP4"), pt.size = 0, ncol = 2)
p3

# PVMs (2,7,9,32)
p4 <- VlnPlot(Myeloid, features = c("C1QB", "RNASE1", "LYVE1", "CD68"), pt.size = 0, ncol = 2)
p4

# Intermediate Macrophages (30)
p5 <- VlnPlot(Myeloid, features = c("CD86", "CXCL2", "CXCL3", "TXNIP"), pt.size = 0, ncol = 2)
p5

# cDC2B (13,25,26)
p6 <- VlnPlot(Myeloid, features = c("CD1C", "CLEC10A"), pt.size = 0, ncol = 2)
p6

# CCR7 (21)
p7 <- VlnPlot(Myeloid, features = c("CCR7", "LAMP3"), pt.size = 0, ncol = 2)
p7

# Classical Monocyte (0,6,18,31)
p8 <- VlnPlot(Myeloid, features = c("CD14", "VCAN", "S100A8", "S100A12"), pt.size = 0, ncol = 2)
p8

# Non-classical Monocyte (10)
p9 <- VlnPlot(Myeloid, features = c("FCGR3A", "SIGLEC10", "APOBEC3A"), pt.size = 0, ncol = 2)
p9

# LAM: 1,3,5,8,11,15,20,23,27 (cluster 17,12, and 27 have lower CD11C expression but do express TREM2)
# PVM: 2,7,9
# CCR7 DC: 21
# cDC1: 29
# cDC2B: 13,25,26
# pDC: 28
# Cycling: 15,24
# nMo: 10
# cMo: 0,6,18,31
# Other Mo: 14
# Other Macrophage:
# Mo-Mac: 22,19,16,30,4

# CD11C+: 1,3,5,8,11,15,20,23,27 (LAM)
# CD11C-: 2,7,9,28,29 
# CD11C intermediate: 4,6,10,12,13,14,16,17,18,19,21,22,25,26,30,31,32
# CD86+: 1,3,5,8,9,11,12,15,17,20,21,23,24,27,30,32, ?19
# CLUSTER 27 is doublet (very difficult to spot but wreaking havoc on clustering)

Idents(Myeloid) <- 'seurat_clusters'
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(21))) <- "CCR7+ DC" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(29))) <- "cDC1" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(13,25,26))) <- "cDC2B" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(28))) <- "pDC" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(0,6,18,31,32))) <- "cMo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(10))) <- "nMo" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(4,16,19,22,30))) <- "Mo-Mac"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(2,7,9))) <- "PVM" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(1,3,5,8,11,20,23,12,17))) <- "LAM"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(15,24))) <- "Cycling Myeloid" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(14))) <- "Other Mo" #

Myeloid[['Manual_Annotation']] <- Idents(Myeloid)

saveRDS(Myeloid, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_RNA_11.1.rds')

#==========================================
# Round 3 (69,330 cells): Myeloid Clustering removing cluster 27
#==========================================
Idents(Myeloid) <- 'seurat_clusters'
Myeloid <- subset(Myeloid, idents = 27, invert = T)
DefaultAssay(Myeloid) <- "RNA"
Myeloid <- Myeloid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Harmony Correction
Myeloid <- RunHarmony(Myeloid, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 300, kmeans_init_iter_max = 1000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Myeloid, ndims = 40)
DimHeatmap(Myeloid, cells = 500, balanced = T, dims = 20:25)

Myeloid <- Myeloid %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:25) %>%
        FindClusters(resolution = 2.0)

#==========================================
# Final Round (68,652 cells): Myeloid Clustering
#==========================================
# Evaluate for Doublets (Stromal doublet: ; Vascular Doublet: ; Lymphoid doublet: ; Junk - )
p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- FeaturePlot(Myeloid, features = c("MS4A1", "JCHAIN", "CLDN5", "CCDC80", "CD3E", "NKG7"))
p1 + p2

p3 <- DotPlot(Myeloid, features = c("TREM2", "CD9", "RNASE1", "LYVE1", "C1QB", "LILRA4", "CLEC9A", "CD1C", "CCR7",
"CD14", "VCAN", "S100A8", "FCGR3A", "WARS", "CXCL2", "CXCL3", "TXNIP", "CLEC10A", "MKI67")) + RotatedAxis()
p3

# LAM: 0,2,6,7,10,15,17,22,32,20
# PVM: 3,4,12
# Cycling: 16,19,29
# pDC: 27
# cDC1: 28
# cDC2B: 13,21,30, ?18
# CCR7 DC: 26
# nMo: 11
# cMo: 1,24,25,8
# Other Mo:9
# Mo-Mac: 31, 
# Unknown: 5, 14,23

# LAM (TREM2+: 0,2,6,7,10,15,16,17,19,20,22,29,32; LIPA: 0,2,6,7,10,15,16,17,19,20,22,29,32)
p1 <- VlnPlot(Myeloid, features = c("C1QB", "TREM2", "CD9", "LIPA"), pt.size = 0, ncol = 3)
p1
# LAM (CD11C+: 0,2,6,7,15,16,17,19,22,29,32; CD11C intermediate: 10,11,13,20,21,26,30)
DefaultAssay(Myeloid) <- "ADT"
p2 <- VlnPlot(Myeloid, features = c("CD11C", "CD86"), pt.size = 0, ncol = 2)
p2

# PVM (LYVE1:3,4,12)
DefaultAssay(Myeloid) <- "RNA"
p3 <- VlnPlot(Myeloid, features = c("LYVE1", "RNASE1", "C1QB"), pt.size = 0, ncol = 3)
p3

# Classical Monocyte (1,8,24,25)
p4 <- VlnPlot(Myeloid, features = c("CD14", "VCAN", "S100A8", "S100A12"), pt.size = 0, ncol = 2)
p4

# Other Mo-Mac: 14,23
Idents(Myeloid) <- 'seurat_clusters'
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(26))) <- "CCR7+ DC" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(28))) <- "cDC1" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(13,18,21,30))) <- "cDC2B" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(27))) <- "pDC" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(1,8,24,25))) <- "cMo" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(11))) <- "nMo" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(5,31))) <- "Mo-Mac"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(3,4,12))) <- "PVM" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(0,2,6,7,10,15,17,22,32,20))) <- "LAM"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(16,19,29))) <- "Cycling Myeloid" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(9))) <- "Other Mo" #
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(14,23))) <- "Other Mac" #

Myeloid[['Manual_Annotation']] <- Idents(Myeloid)

# DotPlot
DefaultAssay(Myeloid) <- "RNA"
DotPlot(Myeloid, features = c("S100A8", "VCAN", "CD14", "FCGR3A", "SERPINA1", "WARS", "APOBEC3A", "AREG", "TREM2", "FABP4", "CD68",
"RNASE1", "LYVE1", "CLEC9A", "CD1C", "CLEC10A", "LILRA4", "CCR7", "HLA-DPA1", "HLA-DQB1", "CD86", "CXCL2", "CXCL3", "TXNIP", "MKI67")) + 
RotatedAxis()

saveRDS(Myeloid, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_RNA_11.1.rds')
write.table(prop.table(table(Myeloid$HATIMID, Myeloid$Manual_Annotation), margin = 1) * 100, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid.Prop.by.HATIMID_11.1.txt')
write.table(table(Myeloid$HATIMID), file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid.Num.by.HATIMID_11.1.txt')

# 5. B Cell SUBCLUSTERING_____________________________________________________

#==========================================
# Round 1: B Cells
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Bcells) <- "RNA"
Bcells <- Bcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction
Bcells <- RunHarmony(Bcells, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(Bcells, ndims = 30)
DimHeatmap(Bcells, balanced = T, cells = 500, dims =5:10) # 9 dimensions

# Run DownStream clustering
Bcells <- Bcells %>% RunUMAP(reduction = 'harmony', dims = 1:9) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:9) %>%
        FindClusters(resolution = 1.0)
        
# Evaluate for Doublets: (T cells: 7, Myeloid:11,8; Stromal/Endo: 5,10)
p1 <- DimPlot(Bcells, label = T) + NoLegend()
p2 <- FeaturePlot(Bcells, features = c("TAGLN", "LYZ", "CLDN5", "CCDC80", "CD3E", "NKG7", "CD1C", "CCR7", "LILRA4"))
p1 + p2

#==========================================
# Round 2 (3,926 cells): B Cells recluster removing clusters 5,7,8,10,11
#==========================================
DefaultAssay(Bcells) <- "RNA"
Bcells <- subset(Bcells, idents = c(5,7,8,10,11), invert = T)

# Find variable features, rescale, and run PCA
Bcells <- Bcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction
Bcells <- RunHarmony(Bcells, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(Bcells, ndims = 30) 
DimHeatmap(Bcells, balanced = T, cells = 500, dims =5:10) # 9 dimensions

# Run DownStream clustering
Bcells <- Bcells %>% RunUMAP(reduction = 'harmony', dims = 1:9) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:9) %>%
        FindClusters(resolution = 1.0)
        
# Evaluate for Doublets: Vascular (11), HPSC (9)
p1 <- DimPlot(Bcells, label = T) + NoLegend()
p2 <- FeaturePlot(Bcells, features = c("TAGLN", "LYZ", "CLDN5", "CCDC80", "CD3E", "NKG7", "CD1C", "CCR7", "LILRA4"))
p1 + p2

#==========================================
# Round 3 (3018 cells): B Cells recluster removing clusters 11,9
#==========================================
DefaultAssay(Bcells) <- "RNA"
Bcells <- subset(Bcells, idents = c(9,11), invert = T)

# Find variable features, rescale, and run PCA
Bcells <- Bcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction
Bcells <- RunHarmony(Bcells, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(Bcells, ndims = 30) 
DimHeatmap(Bcells, balanced = T, cells = 500, dims =5:10) # 9 dimensions

# Run DownStream clustering
Bcells <- Bcells %>% RunUMAP(reduction = 'harmony', dims = 1:9) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:9) %>%
        FindClusters(resolution = 1.0)
        
# Evaluate for Doublets: T/NK cell (13)
p1 <- DimPlot(Bcells, label = T) + NoLegend()
p2 <- FeaturePlot(Bcells, features = c("TAGLN", "LYZ", "CLDN5", "CCDC80", "CD3E", "NKG7", "CD1C", "CCR7", "LILRA4"))
p1 + p2

#==========================================
# Round 4 (2893 cells): B Cells recluster removing clusters 13
#==========================================
DefaultAssay(Bcells) <- "RNA"
Bcells <- subset(Bcells, idents = c(13), invert = T)

# Find variable features, rescale, and run PCA
Bcells <- Bcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction
Bcells <- RunHarmony(Bcells, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(Bcells, ndims = 30) 
DimHeatmap(Bcells, balanced = T, cells = 500, dims =5:10) # 7 dimensions

# Run DownStream clustering
Bcells <- Bcells %>% RunUMAP(reduction = 'harmony', dims = 1:7) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:7) %>%
        FindClusters(resolution = 1.0)
        
p1 <- FeaturePlot(Bcells, features = c("MS4A1", "BANK1", "MZB1", "JCHAIN"))
p1

#==========================================
# Round 5 (2867 cells): B Cells recluster removing clusters 7 (Myeloid contaminant) and 9 (Junk) (Final: 2614)
#==========================================
DefaultAssay(Bcells) <- "RNA"
Bcells <- subset(Bcells, idents = c(7,9), invert = T)

# Find variable features, rescale, and run PCA
Bcells <- Bcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction
Bcells <- RunHarmony(Bcells, "Project", assay.use = "RNA", max.iter.harmony = 20)
ElbowPlot(Bcells, ndims = 30) 
DimHeatmap(Bcells, balanced = T, cells = 500, dims =5:10) # 7 dimensions

# Run DownStream clustering
Bcells <- Bcells %>% RunUMAP(reduction = 'harmony', dims = 1:8) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:8) %>%
        FindClusters(resolution = 1.0)

DefaultAssay(Bcells) <- "ADT"
p1 <- DimPlot(Bcells, label = T) + NoLegend()
p2 <- FeaturePlot(Bcells, features = c("CD19", "CD20", "CD27", "CD38"), min.cutoff = 'q10', max.cutoff = 'q90')

# Memory (CD19+CD27+): 2,3,9
# Naive (CD19+CD20+CD27-): 0,1,4,6,7
# PB (CD27+CD38+):5,8
Idents(Bcells) <- 'seurat_clusters'
Idents(Bcells, cells = WhichCells(Bcells, idents = c(2,3,9))) <- "Memory B cell" #
Idents(Bcells, cells = WhichCells(Bcells, idents = c(0,1,4,6,7))) <- "Naive B cell" #
Idents(Bcells, cells = WhichCells(Bcells, idents = c(5,8))) <- "Plasmablast" #

Bcells[['Manual_Annotation']] <- Idents(Bcells)

saveRDS(Bcells, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Bcells_RNA_11.1.rds') 

# 6. STROMAL CLUSTERING_____________________________________________________
#==========================================
# Round 1: Stromal Clustering (118,807 cells)
#==========================================

# Find variable features, rescale, and run PCA
DefaultAssay(Stromal) <- "RNA"
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()


# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 30:35) # 23 clusters

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:25) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Myeloid-21,31,25,26,30; Lymphoid - 14; Vascular: 12,29,22)
p1 <- DimPlot(Stromal, label = T, raster = F, repel = T) + NoLegend()
p2 <- FeaturePlot(Stromal, features = c("LYZ", "CLDN5", "CD3E", "NKG7"), raster = F)
p1 + p2

#==========================================
# Round 2 (118,807 cells): Stromal Clustering removing clusters 12,14,21,22,25,26,29,30,31
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Stromal) <- "RNA"
Stromal <- subset(Stromal, idents = c(12,14,21,22,25,26,29,30,31), invert = T)
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 25:30)

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Myeloid - 33, Dying - 26,31)
p1 <- DimPlot(Stromal, label = T, repel = T, raster = F) + NoLegend()
p2 <- FeaturePlot(Stromal, features = c("LYZ", "CLDN5", "CD3E", "NKG7"), raster = F)
p3 <- DotPlot(Stromal, features = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "CD3E", "NKG7", "LYZ")) + RotatedAxis()
p1 + p2
p3

#==========================================
# Round 3 (104,008 cells): Stromal Clustering removing clusters 33 (Myeloid) and 26, 31 (dying)
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Stromal) <- "RNA"
Stromal <- subset(Stromal, idents = c(26,31,33), invert = T)
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()


# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 25:30)

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (24 - Junk)
p1 <- DimPlot(Stromal, label = T, repel = T) + NoLegend()
p2 <- FeaturePlot(Stromal, features = c("LYZ", "CLDN5", "CD3E", "NKG7"))
p1 + p2

p3 <- DotPlot(Stromal, features = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "CD3E", "NKG7", "LYZ")) + RotatedAxis()
p3
VlnPlot(Stromal, features = c("MT-CO1", "MT-ATP6"), pt.size = 0)

#==========================================
# Round 4 (102,608): Recluster Stromal with removal of cluster 24 (junk)
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Stromal) <- "RNA"
Stromal <- subset(Stromal, idents = c(24), invert = T)
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()


# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 25:30)

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Junk 30)
p1 <- DimPlot(Stromal, label = T, repel = T) + NoLegend()
p2 <- FeaturePlot(Stromal, features = c("LYZ", "CLDN5", "CD3E", "NKG7"))
p1 + p2

p3 <- DotPlot(Stromal, features = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "CD3E", "NKG7", "LYZ", "CLDN5")) + RotatedAxis()
p3
VlnPlot(Stromal, features = c("MT-CO1", "MT-ATP6"), pt.size = 0)

#==========================================
# Round 5 (101,700): Recluster Stromal removing clusters 30 (Junk)
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Stromal) <- "RNA"
Stromal <- subset(Stromal, idents = c(30), invert = T)
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 25:30)

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Junk - 31)
p1 <- DimPlot(Stromal, label = T, repel = T) + NoLegend()
p2 <- FeaturePlot(Stromal, features = c("LYZ", "CLDN5", "CD3E", "NKG7"))
p1 + p2

p3 <- DotPlot(Stromal, features = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "CD3E", "NKG7", "LYZ")) + RotatedAxis()
p3
VlnPlot(Stromal, features = c("MT-CO1", "MT-ATP6"), pt.size = 0)

# Dot Plot
DotPlot(Stromal, features = c("ADIPOQ", "LPL", "C7", "CD55", "CFD", "PI16", "LUM", "MYOC", "ITM2A", "CXCL14", "APOD", "CXCL12", "GPC3", "CLU", 
"MT1A", "MT2A", "PTGDS", "JUN", "TIMP1", "MMP2", "COL1A1", "IFI6", "PLA2G2A", "POSTN", "MKI67")) + RotatedAxis()

# DGE
markers <- FindAllMarkers(Stromal, only.pos = T, min.pct = 0.25, assay = "RNA")
write.csv(markers, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_Markers.11.1.csv')

# Top Genes (JUN, FOS, IER3, ZFP36, GREB1L, FKBP1A, STK17B, PDK2, SAMHD1, TGFBI, ECHDC2, EBF1, DUSP1, PLXND1, CCL2, BCL6, PRRX1, ZFP36L2, CREB3L1, EGR1, CEBPA)
# JUN,FOS, MYC (early), CEBPB/PPARG (intermediate), CEBPA and lipid (late)

Idents(Stromal) <- "seurat_clusters"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(0,2))) <- "Preadipocyte 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(1))) <- "Preadipocyte 2"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(3))) <- "Preadipocyte 3"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(29))) <- "Transitional Adipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(7,24))) <- "MYOC+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(25))) <- "DKK3+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(16))) <- "Progenitor"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(14,20))) <- "Metallothionein+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(28))) <- "Proliferating Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(22))) <- "Lipomyofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(23))) <- "IFN+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(10,21))) <- "Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(17))) <- "CD55+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(9,13))) <- "PI16+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(5,15))) <- "Junk" # Mitochondrial genes (15) and ribosomal genes (5)
Idents(Stromal, cells = WhichCells(Stromal, idents = c(4))) <- "Preadipocyte 4"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(19))) <- "Preadipocyte 5"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(8,12,18))) <- "Preadipocyte 6"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(6,11))) <- "Preadipocyte 7"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(26,27))) <- "Mature Adipocyte"

Stromal[['Manual_Annotation']] <- Idents(Stromal)

#==========================================
# Final Round (101,512): Stromal removing Junk (92,292 cells)
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Stromal) <- "RNA"
Idents(Stromal) <- "Manual_Annotation"
Stromal <- subset(Stromal, idents = "Junk", invert = T)
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 25:30)

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Junk - 31)
p1 <- DimPlot(Stromal, label = T, repel = T) + NoLegend()
p2 <- FeaturePlot(Stromal, features = c("LYZ", "CLDN5", "CD3E", "NKG7"))
p1 + p2

p3 <- DotPlot(Stromal, features = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "CD3E", "NKG7", "LYZ")) + RotatedAxis()
p3
VlnPlot(Stromal, features = c("MT-CO1", "MT-ATP6"), pt.size = 0)

# Cell types are more for capturing heterogenity rather than identifying true cell types. In addition, not always clear
# where the line between "progenitor" and early "committed" preadipocyte is. In general, progenitors expressing for Stem 
# cell like genes (LUM, CFD, DCN, GSN) while more committed cells express increasing levels of CD36, CXCL12, and lipid
# accumulation genes. In addition, clear that there are adipocytes expressing classic adipogenesis genes (ZFP36, FOS, JUN, MYC, EGR1)
# while others don't express these to much extent. Clustering also shows that while CD55+ fibroblast seem to differentate into myofibroblast
# based on clustering, there are more classic preadipocyte cells expressing fibroblastic markers that also cluster closely with these cells.
# It is possible this could represent de-differentation or lineage from a common progenitor cell.

Idents(Stromal) <- "seurat_clusters"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(5))) <- "Preadipocyte 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(11,19))) <- "Metallothionein+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(24))) <- "IFN+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(30))) <- "Proliferating Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(23,26,21,13))) <- "Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(14))) <- "CD55+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(28,29))) <- "Mature Adipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(25))) <- "Transitional Adipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(0))) <- "Preadipocyte 2"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(2))) <- "Preadipocyte 3"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(27))) <- "DKK3+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(12,20))) <- "MYOC+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(31))) <- "Junk"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(17))) <- "Progenitor 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(22))) <- "Lipomyofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(8))) <- "Progenitor 2"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(10))) <- "Progenitor 3"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(1,3,9))) <- "Preadipocyte 4"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(4))) <- "Preadipocyte 5"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(6,7,18))) <- "PI16+ Fibroblast" # 1,4,7,8,10,11,12,18,19,20,24
Idents(Stromal, cells = WhichCells(Stromal, idents = c(15,16))) <- "CD9hi Progenitor"

Stromal <- subset(Stromal, idents = "Junk", invert = T)

Stromal[['Manual_Annotation']] <- Idents(Stromal)

#==========================================
# Final Round (92,292): Stromal removing Mature Adipocytes (42, 44, and 46), Junk (4, 32, )
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Stromal) <- "RNA"
Stromal <- subset(Stromal, idents = c(4, 32, 42, 44, 46), invert = T)
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 25:30)

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
        FindClusters(resolution = 3.5)

# Evaluate for Doublets (Junk - 31)
p1 <- DimPlot(Stromal, label = T, repel = T) + NoLegend()
p2 <- DimPlot(Stromal, group.by = "Manual_Annotation", label = T) + NoLegend()
p1 + p2

# Save seurat markers
Stromal[['cluster']] <- Idents(Stromal)

#==========================================
# Final Round (87,503): Removing Junk (43)
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Stromal) <- "RNA"
Stromal <- subset(Stromal, idents = c(43), invert = T)
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 25:30)

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
        FindClusters(resolution = 3.5)

# Plot 
p1 <- DimPlot(Stromal, label = T, repel = T) + NoLegend()
p2 <- DimPlot(Stromal, group.by = "Manual_Annotation", label = T) + NoLegend()
p1 + p2

write.csv(Markers, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal.seurat_markers.12.17.csv')

# 0 - Preadipocyte (EGR1, JUN, FOS)
# 1 - Preadipocyte (JUN, FOS, EGR1)
# 2 - Undifferentiated ACS (MYOC, CXCL14, MGP, APOD, IGFBP7)
# 3 - Preadipocyte Fibrobotic (CXCL14, APOE, APOD, MMP2)
# 4 - Mature Preadipocyte (C7)
# 5 - Undifferentiated ACS (PI16, CCN5, DCN)
# 6 - Undifferentiated ACS (PI16, CCN, BGN)
# 7 - Preadipocyte (SOX4, JUN, FOSB)
# 8 - MYOC+ Fibroblast (IGFBP7, MYOC, APOD)
# 9 - Preadipocyte (ZFP36, JUN, ADM)
# 10 - Undifferentiated ACS (MFAP5, IGFBP6, CFD, DCN)
# 11 - CD9hi Progenitor (CD9, ZFP36, IGFBP6)
# 12 - Preadipocyte PTGDS (PTGDS, COL15A1)
# 13 - Preadipocyte (ADM, DEPP, GPC3)
# 14 - CD9hi Progenitor (CD9, CLU, S100A10)
# 15 - Myofibroblast (POSTN, FN1, COL3A1)
# 16 - Undifferentiated ACS (PI16, FBLN, CCN5, DCN)
# 17 - Early Preadipocyte (CXCL12, GPC3, FMO2)
# 18 - Early Preadipocyte Fibrotic (SOX4, MMP2, CCN5)
# 19 - Metallothionein+ Preadipocyte (MT1G, MT1X)
# 20 - Metallothionein+ Preadipocyte (MT1X, PLA2G2A)
# 21 - MYOC + Fibroblast (MYOC, IGFBP7)
# 22 - Preadipocyte (APOD, GPX3, CXCL14, GPC3)
# 23 - CD55+ Fibroblast (CD55, MFAP5)
# 24 - CD55+ Fibroblast (CD55, MFAP5)
# 25 - Myofibroblast (POSTN, FN1, COL3A1)
# 26 - Preadipocyte (FOSB, EGR1, MAFF)
# 27 - Fibroblastic Progenitor (GSN, MMP2, CD36, CFD, THY1)
# 28 - Undifferentiated ACS (MFAP5, CFD, DCN, CCN5)
# 29 - Myofibroblast (TIMP1, FN1)
# 30 - MYOC+ Fibroblast (MYOC, IGFBP7)
# 31 - Lipomyofibroblast (APOE, BGN, APOC)
# 32 - Mature Preadipocyte (C7, APOD, COL15A1, SRPX)
# 33 - Undifferentiated ACS (CFD, GSN, MYOC, DCN)
# 34 - Transitional Adipocyte (C7, CIDEC, FABP4)
# 35 - Metallothionein+ Preadipocyte (CCL2)
# 36 - CD55+ Fibroblast (CD55, PRG4)
# 37 - ISG+ Preadipocyte
# 38 - Myofibroblast (POSTN1, FN1)
# 39 - DKK3+ Fibroblast (DKK3, CD9)
# 40 - ISG+ Preadipocyte (ISG15)
# 41 - Proliferating Myofibroblast (MKI67)
# 42 - Contaminating Cluster

Idents(Stromal) <- "seurat_clusters"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(41))) <- "Proliferating Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(37,40))) <- "ISG+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(39))) <- "DKK3+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(15,29,38))) <- "Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(23,24,36))) <- "CD55+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(19,20,35))) <- "Metallothionein+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(34))) <- "Transitional Adipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(10,33))) <- "Undifferentiated ACS 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(32,3))) <- "Mature Preadipocyte 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(25,31))) <- "Lipomyofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(8,21,30))) <- "MYOC+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(2,6))) <- "Undifferentiated ACS 3"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(5))) <- "Undifferentiated ACS 2"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(4))) <- "Mature Preadipocyte 2"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(14,27))) <- "CD9hi Progenitor"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(0,1,26))) <- "Early Preadipocyte 6"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(13,16))) <- "Early Preadipocyte 5"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(42))) <- "Contamination"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(7,9))) <- "Early Preadipocyte 4"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(11,22))) <- "Early Preadipocyte 3"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(18,17))) <- "Early Preadipocyte 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(12,28))) <- "Early Preadipocyte 2"

Stromal[['Manual_Annotation']] <- Idents(Stromal)

#==========================================
# Final Round (87,128): Removing Contaminating Cluster
#==========================================
DefaultAssay(Stromal) <- "RNA"
Stromal <- subset(Stromal, idents = c("Contamination"), invert = T)
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()

# Run Harmony Correction on RNA PCA matrix
Stromal <- RunHarmony(Stromal, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Stromal, ndims = 40)
DimHeatmap(Stromal, cells = 500, balanced = T, dims = 25:30)

Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
        FindClusters(resolution = 3.5)

# Plot 
p1 <- DimPlot(Stromal, label = T, repel = T) + NoLegend()
p2 <- DimPlot(Stromal, group.by = "Manual_Annotation", label = T) + NoLegend()
p1 + p2


# Top Genes (JUN, FOS, IER3, ZFP36, GREB1L, FKBP1A, STK17B, PDK2, SAMHD1, TGFBI, ECHDC2, EBF1, DUSP1, PLXND1, CCL2, BCL6, PRRX1, ZFP36L2, CREB3L1, EGR1, CEBPA)
# JUN,FOS, MYC (early), CEBPB/PPARG (intermediate), CEBPA and lipid (late)
DotPlot(Stromal, features = c("CD55", "CLU", "DCN", "LUM", "PI16", "JUN", "FOS", "ZFP36", "EGR1", "PTGDS", "SOX4", "MMP2", "MYOC", "ITM2A", "APOE", "POSTN", "TIMP1", "ISG15", "MT1A", "MT2A", "CXCL14", "CXCL12", "GPC3", "CD36",
"C7", "LPL", "CIDEC", "ADIPOQ", "FABP4", "MKI67", "CD9"), group.by = "Manual_Annotation") + RotatedAxis()

saveRDS(Stromal, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_RNA_12.17.rds')
write.table(prop.table(table(Stromal$HATIMID, Stromal$Manual_Annotation), margin = 1) * 100, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal.Prop.by.HATIMID_11.1.txt')
write.table(table(Stromal$HATIMID), file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal.Num.by.HATIMID_11.1.txt')

# 7. VASCULAR CLUSTERING_____________________________________________________
#==========================================
# Round 1 (110,110): 
#==========================================
# Find variable features, rescale, and run PCA
DefaultAssay(Endo) <- "RNA"
Endo <- Endo %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction on RNA PCA matrix
Endo <- RunHarmony(Endo, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)


# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Endo, ndims = 40)
DimHeatmap(Endo, cells = 500, balanced = T, dims = 20:26)


Endo <- Endo %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:25) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Lymphoid - 12,30,33; Myeloid - 22,23 ; Stromal - 10,27,28,32; Junk - 31,26)
p1 <- DimPlot(Endo, label = T, repel = T) + NoLegend()
p2 <- FeaturePlot(Endo, features = c("LYZ", "CCDC80", "CD3E", "NKG7"))
p1 + p2

p3 <- DotPlot(Endo, features = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "CD3E", "NKG7", "LYZ", "CCDC80")) + RotatedAxis()
p3

#==========================================
# Round 2 (110,110 cells): Recluster removing clusters 12, 30, 33 (Lymphoid), 22, 23 (Myeloid), 10, 27, 28, 32 (Stromal), and 26, 31 (Junk/dying)
#==========================================
# Find variable features, rescale, and run PCA
Endo <- subset(Endo, idents = c(10,12,22,23,26,27,28,30,31,32,33), invert = T)
DefaultAssay(Endo) <- "RNA"
Endo <- Endo %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction on RNA PCA matrix
Endo <- RunHarmony(Endo, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)


# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Endo, ndims = 40)
DimHeatmap(Endo, cells = 500, balanced = T, dims = 15:20)


Endo <- Endo %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets (Doublet - 24,25, Junk - 23)
p1 <- DimPlot(Endo, label = T, repel = T) + NoLegend()
p2 <- FeaturePlot(Endo, features = c("LYZ", "CLDN5", "CD3E", "NKG7"))
p1 + p2

p3 <- DotPlot(Endo, features = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "CD3E", "NKG7", "LYZ", "CLDN5", "ACKR1")) + RotatedAxis()
p3

#==========================================
# Round 3 (95,962 cells): Recluster removing clusters 23 (Junk) ,24,25 (immune cell contaminant): Final (94,040 cells)
#==========================================
# Find variable features, rescale, and run PCA
Endo <- subset(Endo, idents = c(23,24,25), invert = T)
DefaultAssay(Endo) <- "RNA"
Endo <- Endo %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction on RNA PCA matrix
Endo <- RunHarmony(Endo, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)


# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Endo, ndims = 40)
DimHeatmap(Endo, cells = 500, balanced = T, dims = 15:20)


Endo <- Endo %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 4.5)

# Evaluate for Doublets
p1 <- DimPlot(Endo, label = T, repel = T) + NoLegend()
p2 <- FeaturePlot(Endo, features = c("LYZ", "CD3E", "NKG7", "CCDC80", "CLDN5"))
p1 + p2

p3 <- DotPlot(Endo, features = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "CD3E", "NKG7", "LYZ", "CD1C", "CCDC80", "CLDN5", 
"GJA5", "GJA4", "HEY1", "CA4", "VWF", "ACKR1", "ACTA2", "TAGLN", "MYH11",
"RGS5", "COX4I2", "PDGFRB", "CSPGR", "NOTCH3", "SEMA3G")) + RotatedAxis()
p3

#==========================================
# Round 4 (94,040 cells): Recluster removing 55 (Myeloid), 49 (VSM/Pericyte)
#==========================================

# Find variable features, rescale, and run PCA
Endo <- subset(Endo, idents = c(49,55), invert = T)
DefaultAssay(Endo) <- "RNA"
Endo <- Endo %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>%
        RunPCA()
        
# Run Harmony Correction on RNA PCA matrix
Endo <- RunHarmony(Endo, "Project", assay.use = "RNA", max.iter.harmony=20, kmeans_init_nstart = 1000, kmeans_init_iter_max = 5000)

# Run dimensional reduction and clustering on harmony corrected PCA
ElbowPlot(Endo, ndims = 40)
DimHeatmap(Endo, cells = 500, balanced = T, dims = 15:20)


Endo <- Endo %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
        FindClusters(resolution = 2.0)

# Evaluate for Doublets
p1 <- DimPlot(Endo, label = T, repel = T) + NoLegend()
p2 <- FeaturePlot(Endo, features = c("LYZ", "CD3E", "NKG7", "CCDC80", "CLDN5"))
p1 + p2

#==========================================
# Final (94,040 cells)
#==========================================
DotPlot(Endo, features = c("GJA5", "GJA4", "HEY1", "CA4", "VWF", "ACKR1", "ACTA2", "TAGLN", "MYH11",
"RGS5", "COX4I2", "PDGFRB", "CSPGR", "NOTCH3", "SEMA3G")) + RotatedAxis()


pVein <- FeaturePlot(Endo, features = "ACKR1", min.cutoff = 'q10', max.cutoff = 'q90', raster = F) # clusters 6,10,13,21,22,27
pART <- FeaturePlot(Endo, features = c("GJA5", "GJA4", "ACTA2", "HEY1", "SEMA3G", "BMX"), min.cutoff = 'q10', max.cutoff = 'q90', raster = F) # 19,12,14
pCap <- FeaturePlot(Endo, features = c("CA4", "GNG11", "VWF"), min.cutoff = 'q10', max.cutoff = 'q90', raster = F) # 1,3,15,5,24,7,4,
pVSM <- FeaturePlot(Endo, features = c("NOTCH3", "MFAP5", "PDGFRB", "MYH11", "ACTA2", "TAGLN2"), min.cutoff = 'q10', max.cutoff = 'q90', raster = F) # 9,18,8
pPeri <- FeaturePlot(Endo, features = c("RSG5", "COX4I2", "PDGFRB", "CSPG4"), min.cutoff = 'q10', max.cutoff = 'q90', raster = F) # 16,8,26


Idents(Endo) <- "seurat_clusters"
Idents(Endo, cells = WhichCells(Endo, idents = c(7,13,19,20,26))) <- "Venous"
Idents(Endo, cells = WhichCells(Endo, idents = c(15))) <- "Arterial"
Idents(Endo, cells = WhichCells(Endo, idents = c(11,14))) <- "Intermediate Capillary"
Idents(Endo, cells = WhichCells(Endo, idents = c(2,3,4,5,6,21))) <- "Capillary"
Idents(Endo, cells = WhichCells(Endo, idents = c(10,18,24))) <- "Pericyte"
Idents(Endo, cells = WhichCells(Endo, idents = c(9,16))) <- "VSMC"
Idents(Endo, cells = WhichCells(Endo, idents = c(12,25))) <- "Venous-EndoMT-like"
Idents(Endo, cells = WhichCells(Endo, idents = c(23))) <- "Arterial-EndoMT-like"
Idents(Endo, cells = WhichCells(Endo, idents = c(0,1,8,17,22))) <- "Capillary-EndoMT-like"

Endo[['Manual_Annotation']] <- Idents(Endo)

saveRDS(Endo, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Endo_RNA_1.3.rds')
write.table(prop.table(table(Endo$HATIMID, Endo$Manual_Annotation), margin = 1) * 100, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Endo.Prop.by.HATIMID.1.3.txt')
write.table(table(Endo$HATIMID), file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Endo.Num.by.HATIMID.1.3.txt')
