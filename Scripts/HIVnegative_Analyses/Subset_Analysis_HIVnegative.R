##############################################################################################
##------------------------- SUBSET ANALYSIS ---------------------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will clean the Lymphoid (minus B cell) dataset and then pass the 
# cleaned object to the scib pipeline to determine which integration technique preserves
# biological variation most effectively. For efficiency and speed, I will use harmony to 
# clean the dataset and annotate cell populations. I will save the unintegrated version and 
# import this into the scib pipeline. I will also convert this into h5ad format.
##############################################################################################

# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(future)
library(tidyverse)
library(harmony)

plan("multiprocess", workers = 5)
options(future.globals.maxSize = 5000 * 1024^2)
options(future.seed = TRUE)

set.seed(7412)

date = "6.27"

filtered_dir <- paste0("../HATIM_Analysis/", date, "/Integrated_Filtered_HIVneg")
subset_dir <- paste0("../HATIM_Analysis/", date, "/SubsetAnalysis_HIVneg")
marker_dir <- "../Markers" 

data_hatim <- readr::read_csv('../HATIMStudyVisit_DATA_2021-08-09_1359.csv')

data_hatim <- data_hatim %>% mutate(HATIMID = as.factor(hatim_clin_visit_pid),
                                    StudyGroup = factor(hatim_final_arm, levels = c(1,2,3,4),
                                                        labels = c("HIV+ non-diabetic", "HIV+ prediabetic", 
                                                        "HIV+ diabetic", "HIV- diabetic")),
                                    HIV = case_when(hatim_final_arm == 4 ~ "HIVneg",
                                                    TRUE ~ "HIVpos"),
                                    Glucose = case_when(hatim_final_arm == 1 ~ "Glucose Tolerant",
                                                        TRUE ~ "Glucose Intolerant"),
                                    Sex = factor(meta_sex, levels = c(0,1), labels= c("Male", "Female")),
                                    homa2_ir = as.numeric(homa2_ir), 
                                    Age = meta_age,
                                    BMI = meta_bmi,
                                    FBG = meta_fbg,
                                    hba1c = meta_hba1c)

Hash_link <- read.csv("../Hash_Link.csv")
Hash_link <- Hash_link %>% mutate(Hash = as.character(Hash))
colnames(Hash_link) <- c("Lane", "HATIMID", "Soup_Hash")

#2. MYELOID ANALYSIS_____________________________________________________
###############################
# Import Data
###############################
Myeloid <- readRDS(paste(filtered_dir, "Myeloid_HIVneg.rds", sep = "/"))

################################
# Plot Doublets (clusters 95,77,44)
###############################
DotPlot(Myeloid, features = c("CD3E", "NKG7", "CCDC80", "CLDN5", 'FCGR3A', "VCAN", "CD1C", "CLEC9A", "CCR7", "MKI67", "TREM2", "LYVE1")) + RotatedAxis() + coord_flip()

p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- FeaturePlot(Myeloid, features = c("CD3E", "NKG7", "CCDC80", "CLDN5"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

# Remove Doublet Clusters (44,76,81,94), nearly all from HIV-negative
Myeloid <- subset(Myeloid, idents = c(44,76,81,94), invert = T)

# Junk (clusters 48), but will keep for now, remove later
Idents(Myeloid) <- "seurat_clusters"
VlnPlot(Myeloid, features = "percent.mt", pt.size = 0) + NoLegend()

################################
# Round 1: n = 67,791 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Myeloid) <- "RNA"
Myeloid <- Myeloid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Myeloid <- RunHarmony(Myeloid, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Myeloid <- Myeloid %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:25) %>% 
  FindClusters(resolution = 1.5)
  
###############################
# Visualize
###############################
# Junk: None
p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- VlnPlot(Myeloid, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: cluster 20, 15, 16 (4 [DC],13 [Mo-Mac TIMP+],22 [PVM],24 [ISG+],26 [Mo-Mac])
DotPlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "CCDC80", "TAGLN", "ACTA2", "MS4A1", "MKI67", 'FCGR3A', "VCAN", "CD1C", "CLEC9A", "CCR7", "TREM2", "LYVE1", "LILRA4")) + RotatedAxis() + coord_flip()
p2 <- FeaturePlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "CCDC80", "MKI67"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: cluster 27
p2 <- VlnPlot(Myeloid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Remove Doublets
###############################
Myeloid <- subset(Myeloid, idents = c(15,16,20,27), invert = T)

################################
# Round 2: n = 63,798 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Myeloid) <- "RNA"
Myeloid <- Myeloid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Myeloid <- RunHarmony(Myeloid, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Myeloid <- Myeloid %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:25) %>% 
  FindClusters(resolution = 1.5)
  
###############################
# Visualize
###############################
# Junk: Cluster 27
p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- VlnPlot(Myeloid, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 26
DotPlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "CCDC80", "TAGLN", "ACTA2", "MS4A1", "MKI67", 'FCGR3A', "VCAN", "CD1C", "CLEC9A", "CCR7", "TREM2", "LYVE1", "LILRA4", "CXCL3", "TIMP1", "ISG15", "HLA-DRA")) + RotatedAxis()
p2 <- FeaturePlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "CCDC80"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: ?15, Macrophage so will keep for now, remove later if needed
p2 <- VlnPlot(Myeloid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Remove Doublets/Junk
###############################
Myeloid <- subset(Myeloid, idents = c(26,27), invert = T)

################################
# Round 3: n = 63,309 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Myeloid) <- "RNA"
Myeloid <- Myeloid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Myeloid <- RunHarmony(Myeloid, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Myeloid <- Myeloid %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:25) %>% 
  FindClusters(resolution = 1.0)

p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- DimPlot(Myeloid, group.by = "Annotation", label = T) + NoLegend()
p1 + p2

# Remove cluster 23 (2 cells and doublet)
Myeloid <- subset(Myeloid, idents = 23, invert = T)

Idents(Myeloid) <- "seurat_clusters"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(3,2,13))) <- "IM"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(4,6,16))) <- "cMo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(11))) <- "Other Mo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(7))) <- "nMo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(1,8))) <- "PVM"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(0,12))) <- "LAM"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(14,15))) <- "Cycling Myeloid"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(18))) <- "ISG+ Mo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(5))) <- "Mo-Mac 2"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(20))) <- "pDC"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(19))) <- "Migratory DC"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(10))) <- "Mo-Mac 1"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(9,17,21))) <- "cDC2B"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(22))) <- "cDC1"
Myeloid[['Annotation']] <- Idents(Myeloid)

###############################
# Macrophage subset
###############################
Macrophage <- subset(Myeloid, idents = c("PVM", "IM", "LAM", "Mo-Mac 1", "Mo-Mac 2"))

# Recluster
DefaultAssay(Macrophage) <- "RNA"
Macrophage <- Macrophage %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

# Run Harmony
Macrophage <- RunHarmony(Macrophage, "Lane", assay.use = "RNA", max.iter.harmony = 30)

# Run Donwstream Analysis
Macrophage <- Macrophage %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(resolution = 1.0)

# Visualize
p1 <- DimPlot(Macrophage, label = T) + NoLegend()
p2 <- DimPlot(Macrophage, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

p2 <- VlnPlot(Macrophage, features = 'percent.mt', pt.size = 0) + NoLegend()
p1 + p2

# Cluster 15
p2 <- VlnPlot(Macrophage, features = 'nCount_RNA', pt.size = 0) + NoLegend()
p1 + p2

Doublet <- WhichCells(Macrophage, idents = c(15))
Macrophage <- subset(Macrophage, idents = c(15), invert = T)

# Recluster
DefaultAssay(Macrophage) <- "RNA"
Macrophage <- Macrophage %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

# Run Harmony
Macrophage <- RunHarmony(Macrophage, "Lane", assay.use = "RNA", max.iter.harmony = 30)

# Run Donwstream Analysis
Macrophage <- Macrophage %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(resolution = 0.4)

###############################
# Annotation
###############################
DotPlot(Macrophage, features = c("CD9", "TREM2", "LYVE1", "LYZ", "VCAN", "CXCL2", "APOE", "SPP1", "APOC1", "FCER1A", "S100A8")) + RotatedAxis()

Idents(Macrophage) <- "seurat_clusters"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(0))) <- "IM"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(2))) <- "PVM"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(1))) <- "LAM"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(4))) <- "Mo-Mac 1"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(3))) <- "Mo-Mac 2"
Macrophage[["Annotation"]] <- Idents(Macrophage)

# Markers
markers <- FindAllMarkers(Macrophage, only.pos = T, min.pct = 0.1, assay = "RNA")
write.csv(markers, file = paste(marker_dir, "Macrophage.Markers.csv", sep = "/"))

################################
# Add Metadata
###############################
Macrophage@meta.data[c(38:46)] <- NULL
metadata <- Macrophage@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG", "hba1c")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG, hba1c)
rownames(metadata) <- colnames(Macrophage)

Macrophage <- AddMetaData(Macrophage, metadata)

###############################
# Save RDS Macrophage
###############################
saveRDS(Macrophage, file = paste(subset_dir, "Macrophage_HIVneg.rds", sep = "/"))

markers <- FindAllMarkers(Myeloid, only.pos = T, min.pct = 0.1, assay = "RNA")
marker_dir <- "../Markers" 
dir.create(marker_dir)
write.csv(markers, paste(marker_dir, "Myeloid.Markers.csv", sep = "/"))

###############################
# Add Macrophage Annotation
###############################
# Subset out Junk clusters
Myeloid <- subset(Myeloid, cells = Doublet, invert = T)

# Now relabel Macrophages
Idents(Myeloid, cells = WhichCells(Macrophage, idents = c("IM"))) <- "IM"
Idents(Myeloid, cells = WhichCells(Macrophage, idents = c("LAM"))) <- "LAM"
Idents(Myeloid, cells = WhichCells(Macrophage, idents = c("PVM"))) <- "PVM"
Idents(Myeloid, cells = WhichCells(Macrophage, idents = c("Mo-Mac 1"))) <- "Mo-Mac 1"
Idents(Myeloid, cells = WhichCells(Macrophage, idents = c("Mo-Mac 2"))) <- "Mo-Mac 2"
Myeloid[['Annotation']] <- Idents(Myeloid)


################################
# Add Metadata
###############################
Myeloid@meta.data[c(38:46)] <- NULL
metadata <- Myeloid@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG", "hba1c")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG, hba1c)
rownames(metadata) <- colnames(Myeloid)

Myeloid <- AddMetaData(Myeloid, metadata)

################################
# Save Object
###############################
saveRDS(Myeloid, paste(subset_dir, "Myeloid_HIVneg.rds", sep = "/"))

#2. LYMPHOID ANALYSIS_____________________________________________________
###############################
# Import Data
###############################
Lymphoid <- readRDS(paste(filtered_dir, "Lymphoid_HIVneg.rds", sep = "/"))

###############################
# PLOT Doublets
###############################
DotPlot(Lymphoid, features = c("LYZ", "CD68", "CLDN5", "GNG11", "CCDC80", "LUM")) + RotatedAxis()

# Remove Cluster 75
Lymphoid <- subset(Lymphoid, idents = c(75), invert = T)

###############################
# ROUND 1: n = 34,955 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Lymphoid) <- "RNA"
Lymphoid <- Lymphoid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Lymphoid <- RunHarmony(Lymphoid, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Lymphoid <- Lymphoid %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(resolution = 1.5)
  
###############################
# Visualize
###############################
# Junk: None
p1 <- DimPlot(Lymphoid, label = T) + NoLegend()
p2 <- VlnPlot(Lymphoid, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 11
DotPlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67", "ACTA2", "TAGLN", "CD3D", "CD8A", "IL7R", "LEF1", "CCR7", "GPR183", "LDHB", "IL32", "COTL1",
"CD69", "TNF", "IFNG", "PRF1", "CCR6", "KIT", "FCGR3A", "TRDC", "TRDV1", "CTLA4")) + RotatedAxis()
p2 <- FeaturePlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: None
p2 <- VlnPlot(Lymphoid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Remove Doublets
###############################
Lymphoid <- subset(Lymphoid, idents = c(11), invert = T)

###############################
# ROUND 2: n = 33,776 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Lymphoid) <- "RNA"
Lymphoid <- Lymphoid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Lymphoid <- RunHarmony(Lymphoid, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Lymphoid <- Lymphoid %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(resolution = 1.5)
  
###############################
# Visualize
###############################
p1 <- DimPlot(Lymphoid, label = T) + NoLegend()
p2 <- DimPlot(Lymphoid, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

# Doublets: Cluster 23
DotPlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67", "ACTA2", "TAGLN", "CD3D", "CD8A", "IL7R", "LEF1", "CCR7", "GPR183", "LDHB", "IL32", "COTL1",
"CD69", "TNF", "IFNG", "PRF1", "CCR6", "KIT", "FCGR3A", "TRDC", "TRDV1", "CTLA4")) + RotatedAxis()

###############################
# Remove Doublets
###############################
Lymphoid <- subset(Lymphoid, idents = c(23), invert = T)

###############################
# ROUND 3: n = 33,676 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Lymphoid) <- "RNA"
Lymphoid <- Lymphoid %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Lymphoid <- RunHarmony(Lymphoid, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Lymphoid <- Lymphoid %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(resolution = 1.5)
  
###############################
# Visualize
###############################
# Doublets: None
DotPlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "ACTA2", "TAGLN", "MKI67")) + RotatedAxis()

###############################
# Add Annotation Labels
###############################
p1 <- DimPlot(Lymphoid, group.by = "Annotation", label = T) + NoLegend()
p2 <- DimPlot(Lymphoid, label = T) + NoLegend()
p1 + p2

DefaultAssay(Lymphoid) <- "ADT"
p2 <- VlnPlot(Lymphoid, features = c("CD8", "CD4", "CD27", "CD57", "CD45RA", "CD3"), pt.size = 0, ncol = 2) + NoLegend()
p1 + p2

DotPlot(Lymphoid, features = c("CD3D", "CD8A", "SELL", "LEF1", "GPR183", "LDHB", "COTL1", "IL32", 'CD69', "VIM", "CCL5", "PRF1", "CTLA4", "CCR6", "TRDC", "TRDV1", "FCGR3A", "XCL1", "KIT")) + RotatedAxis()

clust8 <- subset(Lymphoid, idents = 8)
DefaultAssay(clust8) <- "RNA"
clust8 <- clust8 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()
clust8 <- RunHarmony(clust8, "Lane", assay.use = "RNA", max.iter.harmony = 30)
clust8 <- clust8 %>% RunUMAP(reduction = 'harmony', dims = 1:6) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:6) %>% 
  FindClusters(resolution = 0.1)

cluster11 <- subset(Lymphoid, idents = 11)
DefaultAssay(cluster11) <- "ADT"
Cytotoxic <- subset(Lymphoid, idents = c(0,12,18))
DefaultAssay(Cytotoxic) <- "ADT"

clust14 <- subset(Lymphoid, idents = 14)
DefaultAssay(clust14) <- "RNA"
clust14 <- clust14 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()
clust14 <- RunHarmony(clust14, "Lane", assay.use = "RNA", max.iter.harmony = 30)
clust14 <- clust14 %>% RunUMAP(reduction = 'harmony', dims = 1:6) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:6) %>% 
  FindClusters(resolution = 0.1)

DefaultAssay(Lymphoid) <- "RNA"
Idents(Lymphoid) <- "seurat_clusters"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(10))) <- "CD4 Naive"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(13))) <- "CD8 Naive"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(4))) <- "CD4 TCM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(3))) <- "CD4 TEM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(16))) <- "CD4 Regulatory"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(17))) <- "MAIT"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(5))) <- "CD8 TCM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(2))) <- "CD8 TEM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(0,12,18))) <- "CD8 Cytotoxic"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(15))) <- "Cycling T & NK"
Idents(Lymphoid, cells = WhichCells(clust8, idents = c(2))) <- "ILC"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(8))) <- "Immature NK"
Idents(Lymphoid, cells = WhichCells(clust14, idents = c(0))) <- "Gamma Delta"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(14))) <- "CD8 TEM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(6))) <- "CD16 mNK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(1,9))) <- "CD57 mNK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(7))) <- "Gamma Delta"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(19))) <- "MT+ T cells"
Idents(Lymphoid, cells = WhichCells(cluster11, expression = CD8 > 1.0)) <- "CD8 TEM"
Idents(Lymphoid, cells = WhichCells(cluster11, expression = CD8 <= 1.0)) <- "CD4 TEM"
Idents(Lymphoid, cells = WhichCells(Cytotoxic, expression = CD8 > 1.0)) <- "CD8 Cytotoxic"
Idents(Lymphoid, cells = WhichCells(Cytotoxic, expression = CD8 <= 1.0)) <- "CD4 Cytotoxic"
Lymphoid[['Annotation']] <- Idents(Lymphoid)

#Mixed: 14
###############################
# Verify Annotations
###############################
DotPlot(Lymphoid, features = c("CD3D", "CD8A", "SELL", "LEF1", "GPR183", "LDHB", "COTL1", "IL32", 'CD69', "VIM", "CCL5", "PRF1", "CTLA4", "CCR6", "TRDC", "TRDV1", "FCGR3A", "XCL1", "KIT")) + RotatedAxis()

DefaultAssay(Lymphoid) <- "ADT"
VlnPlot(Lymphoid, features = c("CD4", "CD8", "CD27", "CD57", "CD45RA", "CD16"), stack = T, flip = T)

markers <- FindAllMarkers(Lymphoid, only.pos = T, min.pct = 0.1, assay = "RNA")
write.csv(markers, file = paste(marker_dir, "Lymphoid.Markers.csv", sep = "/"))

################################
# Add Metadata
###############################
Lymphoid@meta.data[c(38:46)] <- NULL
metadata <- Lymphoid@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG", "hba1c")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG, hba1c)
rownames(metadata) <- colnames(Lymphoid)

Lymphoid <- AddMetaData(Lymphoid, metadata)
###############################
# Save Harmony Integration
###############################
saveRDS(Lymphoid, file = paste0(subset_dir, "/", "Lymphoid_HIVneg.rds"))

#3. B Cell Analysis_____________________________________________________
###############################
# Import Data
###############################
Bcells <- readRDS(paste(filtered_dir, "Bcells_HIVneg.rds", sep = "/"))

###############################
# PLOT Doublets
###############################
DotPlot(Bcells, features = c("LYZ", "CD68", "CLDN5", "GNG11", "CCDC80", "LUM", "CD3E", "NKG7", "CD79A", "JCHAIN")) + RotatedAxis()

###############################
# ROUND 1: n = 3,034 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Bcells) <- "RNA"
Bcells <- Bcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Bcells <- RunHarmony(Bcells, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Bcells <- Bcells %>% RunUMAP(reduction = 'harmony', dims = 1:10) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:10) %>% 
  FindClusters(resolution = 1.0)
  
###############################
# Visualize
###############################
# Junk: 9
p1 <- DimPlot(Bcells, label = T) + NoLegend()
p2 <- VlnPlot(Bcells, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 7, 10,13
DotPlot(Bcells, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "CD3E", "NKG7", "TAGLN")) + RotatedAxis()

#  RNA Count: None
p2 <- VlnPlot(Bcells, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk and Doublets
###############################
Bcells <- subset(Bcells, idents = c(7,10,13), invert = T)

###############################
# ROUND 2: n = 2713 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Bcells) <- "RNA"
Bcells <- Bcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Bcells <- RunHarmony(Bcells, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Bcells <- Bcells %>% RunUMAP(reduction = 'harmony', dims = 1:6) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:6) %>% 
  FindClusters(resolution = 0.2)
  
###############################
# Visualize
###############################
# Doublet: None
p1 <- DimPlot(Bcells, label = T) + NoLegend()
p2 <- FeaturePlot(Bcells, features = c("CLDN5", "LYZ", "C1QB", "CCDC80", "NKG7", "CD3E"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

###############################
# Annotation
###############################
DefaultAssay(Bcells) <- "ADT"
p3 <- VlnPlot(Bcells, features = c("CD19", "CD20", "CD38"), stack = T, flip = T)
p1 + p3

DefaultAssay(Bcells) <- "RNA"
Idents(Bcells) <- "seurat_clusters"
Idents(Bcells, cells = WhichCells(Bcells, idents = c(2,3))) <- "Plasmablast"
Idents(Bcells, cells = WhichCells(Bcells, idents = c(0,1))) <- "B Cell"
Bcells[['Annotation']] <- Idents(Bcells)

markers <- FindAllMarkers(Bcells, only.pos = T, min.pct = 0.1, assay = "RNA")
write.csv(markers, file = paste(marker_dir, "Bcells.Markers.csv", sep = "/"))

###############################
# Add Metadata
###############################
Bcells@meta.data[c(38:46)] <- NULL
metadata <- Bcells@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG", "hba1c")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG, hba1c)
rownames(metadata) <- colnames(Bcells)

Bcells <- AddMetaData(Bcells, metadata)

###############################
# Save Harmony Integration
###############################
saveRDS(Bcells, file = paste0(subset_dir, "/", "Bcells_HIVneg.rds"))

#4. Stromal Analysis_____________________________________________________
###############################
# Import Data
###############################
Stromal <- readRDS(paste(filtered_dir, "Stromal_HIVneg.rds", sep = "/"))

###############################
# PLOT Doublets
###############################
# Doublets: Clusters 32,78
DotPlot(Stromal, features = c("LYZ", "CD68", "CLDN5", "GNG11", "CD3E", "NKG7", "TAGLN", "MKI67")) + RotatedAxis()
p1 <- DimPlot(Stromal, label = T) + NoLegend()
p2 <- FeaturePlot(Stromal, features = c("CD3E", "NKG7", "CLDN5", "TAGLN", "CD68", "LYZ"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

# Remove Doublets
Stromal <- subset(Stromal, idents = c(32,78), invert = T)

###############################
# ROUND 1: n = 96,892 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Stromal) <- "RNA"
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Stromal <- RunHarmony(Stromal, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>% 
  FindClusters(resolution = 2.5)
  
###############################
# Visualize
###############################
# Junk: None
p1 <- DimPlot(Stromal, label = T) + NoLegend()
p2 <- VlnPlot(Stromal, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 36,35,26,28,29
DotPlot(Stromal, features = c("CLDN5", "LYZ", "CD68", "NKG7", "CD3E", "TAGLN", "ISG15", "MKI67")) + RotatedAxis()
p2 <- FeaturePlot(Stromal, features = c("CLDN5", "LYZ", "NKG7", "CD3E", "TAGLN"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

# Ribosomal:30
p2 <- DotPlot(Stromal, features = c("RPL11", "RPS8", "RPS27", "RPL12", "RPL30")) + RotatedAxis()
p1 + p2

# RNA Count (low: 30 high: 34 [cycling fibroblast]
p2 <- VlnPlot(Stromal, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk, Doublets
###############################
Stromal <- subset(Stromal, idents = c(26,28,29,30,35,36), invert = T)

###############################
# ROUND 2: n = 92,550 cells
###############################
DefaultAssay(Stromal) <- "RNA"
Stromal <- Stromal %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Stromal <- RunHarmony(Stromal, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Stromal <- Stromal %>% RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>% 
  FindClusters(resolution = 1.0)
  
###############################
# Visualize
###############################
# Junk: None
p1 <- DimPlot(Stromal, label = T) + NoLegend()
p2 <- VlnPlot(Stromal, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: None
DotPlot(Stromal, features = c("CLDN5", "LYZ", "NKG7", "CD3E", "TAGLN")) + RotatedAxis()
p2 <- FeaturePlot(Stromal, features = c("CLDN5", "LYZ", "NKG7", "CD3E"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

# RNA Count: None
p2 <- VlnPlot(Stromal, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

p1 <- DimPlot(Stromal, label = T) + NoLegend()
p2 <- DimPlot(Stromal, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

###############################
# Add Metadata
###############################
Stromal@meta.data[c(38:46)] <- NULL
metadata <- Stromal@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG", "hba1c")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG, hba1c)
rownames(metadata) <- colnames(Stromal)

Stromal <- AddMetaData(Stromal, metadata)

###############################
# Annotation
###############################
DotPlot(Stromal, features = c("CD55", "DPP4", "MFAP5", "MYOC", "POSTN", "MKI67", "TIMP3", "MT1X", "FABP5", "FABP4", "C7", "LPL", "EGR1", "JUN", "FOS", "DCN", "LUM", "PI16", "GSN", "CFD", "APOD", "THY1",
"TMEM176B", "COL1A1", "COL1A2", "TIMP1", "ACTA2", 'MYH11')) + RotatedAxis()

Idents(Stromal) <- "seurat_clusters"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(9))) <- "PCOLCE+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(7))) <- "MYOC+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(6,11))) <- "Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(12))) <- "Cycling Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(8))) <- "Metallothionein+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(2))) <- "Mature Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(4))) <- "Adipose Progenitor Cell 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(1,3))) <- "Early Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(0))) <- "ECM-Producing Early Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(10,5))) <- "Adipose Progenitor Cell 2"
Stromal[['Annotation']] <- Idents(Stromal)

###############################
# Save RDS
###############################
saveRDS(Stromal, file = paste0(subset_dir, "/", "Stromal_HIVneg.rds"))

#5. Vascular Analysis_____________________________________________________
###############################
# Import Data
###############################
Vascular <- readRDS(paste(filtered_dir, "Vascular_HIVneg.rds", sep = "/"))

###############################
# PLOT Doublets
###############################
# Clusters 74,90, 91 (but cycling), 52, 21, 31, 41,61, 60 (dual smooth muscle-endothelial expression)
p1 <- DimPlot(Vascular, label = T) + NoLegend()
DotPlot(Vascular, features = c("LYZ", "CD68", "CCDC80", "CD3E", "NKG7", "MKI67", "GJA5", "GJA4", "HEY1", "SEMA3G", "CA4", "VWF", "GNG11", "ACKR1", "ACTA2", "TAGLN", "MYH11", "NOTCH3", "PDGFRB", "MFAP5",
"RGS5", "COX4I2", "CSPG4", "SNAI1", "SNAI2")) + RotatedAxis()

# Possible Junk - none
VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()

# RNA count High: ?41
VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()

# Removing Doublets
Vascular <- subset(Vascular, idents = c(21,31,41,52,60,61,74,90), invert = T)

###############################
# ROUND 1: n = 63,125 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Vascular) <- "RNA"
Vascular <- Vascular %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Vascular <- RunHarmony(Vascular, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Vascular <- Vascular %>% RunUMAP(reduction = 'harmony', dims = 1:25) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:25) %>% 
  FindClusters(resolution = 2.5)
  
###############################
# Visualize
###############################
# Junk: Clusters 26
p1 <- DimPlot(Vascular, label = T) + NoLegend()
p2 <- VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 35, 34, 32, 30, 33, part of 29
DotPlot(Vascular, features = c("LYZ", "CD68", "CCDC80", "CD3E", "NKG7", "MKI67", "GJA5", "GJA4", "HEY1", "SEMA3G", "CA4", "VWF", "GNG11", "ACKR1", "ACTA2", "TAGLN", "MYH11", "NOTCH3", "PDGFRB", "MFAP5",
"RGS5", "COX4I2", "CSPG4", "SNAI1", "SNAI2")) + RotatedAxis()

# RNA Count High: 14, 31
p2 <- VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk, Doublets
###############################
clust29 <- subset(Vascular, idents = 29)
cycling = WhichCells(clust29, expression = RGS5 > 2 | VWF > 2)
Idents(Vascular, cells = cycling) <- "Cycling"
Vascular <- subset(Vascular, idents = c(14,26,29,30,31,32,33,34,35), invert = T)

###############################
# ROUND 2: n = 58,199 cells
###############################
###############################
# Variable Features, scaling, regress, PCA
###############################
DefaultAssay(Vascular) <- "RNA"
Vascular <- Vascular %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()

###############################
# Run Harmony
###############################
Vascular <- RunHarmony(Vascular, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
Vascular <- Vascular %>% RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(resolution = 0.5)
  
###############################
# Visualize
###############################
# Junk: Pass QC
p1 <- DimPlot(Vascular, label = T) + NoLegend()
p2 <- VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Cluster
p2 <- DotPlot(Vascular, features = c("CCDC80", "LUM", "LYZ", "NKG7", "CD3E", "MS4A1", "MKI67")) + RotatedAxis()
p1 + p2

# RNA Count High:
p2 <- VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Add Metadata
###############################
Vascular@meta.data[c(38:46)] <- NULL
metadata <- Vascular@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG", "hba1c")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG, hba1c)
rownames(metadata) <- colnames(Vascular)

Vascular <- AddMetaData(Vascular, metadata)

###############################
# Add Annotations
###############################
p2 <- DotPlot(Vascular, features = c("GJA5", "GJA4", "HEY1", "SEMA3G", "CA4", "VWF", "GNG11", "ACKR1", "ACTA2", "TAGLN", "MYH11", "NOTCH3", "PDGFRB", "MFAP5",
"RGS5", "COX4I2", "CSPG4", "SNAI1", "SNAI2", "LPL", "PPARG", "CD36")) + RotatedAxis()
p1 + p2

Idents(Vascular) <- "seurat_clusters"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(4))) <- "Arterial EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(5))) <- "VSMC 1"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(8))) <- "VSMC 2"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(3))) <- "Pericyte"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(9))) <- "Cycling Vascular Cells"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(2))) <- "Venous EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(0,1,6,7))) <- "Capillary EC"
Vascular[['Annotation']] <- Idents(Vascular)

###############################
# Save
###############################
saveRDS(Vascular, file = paste0(subset_dir, "/", "Vascular_HIVneg.rds"))

# 6. T Cells_____________________________________________________
Lymphoid <- readRDS(paste0(subset_dir, "/", "Lymphoid_HIVneg.rds"))

# Isolate T Cells: To do this, I will grab clusters that are identified as T cells and are expressing CD3 ADT. This will include MAIT and Gamma Delta initially
Tcells <- subset(Lymphoid, idents = c("CD4 Cytotoxic", "CD8 Cytotoxic", "CD4 Regulatory", "MAIT", "CD8 TEM", "CD4 TEM", "CD8 TCM", "CD8 Naive", "CD4 Naive", "Gamma Delta", "CD4 TCM", "MT+ T cells"))
DefaultAssay(Tcells) <- "RNA"
###############################
# Variable Selection, PCA
###############################
# Regress out TCRs to prevent clustering driven by TCR
tcr_genes <- row.names(Tcells)[grep(pattern="^TR(A|B|G|D)V", row.names(Tcells))]
Tcells <- AddModuleScore(object = Tcells, features = list(tcr_genes), name = 'tcr')

DefaultAssay(Tcells) <- "RNA"
Tcells <- Tcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "tcr1"), assay = "RNA") %>%
        RunPCA()

# Run Harmony
Tcells <- RunHarmony(Tcells, "Lane", assay.use = "RNA", max.iter.harmony = 30)

# Run Donwstream Analysis
Tcells <- Tcells %>% RunUMAP(reduction = 'harmony', dims = 1:15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:15) %>% 
  FindClusters(resolution = 1.5)

###############################
# Visualize
###############################
p1 <- DimPlot(Tcells, label = T) + NoLegend()
p2 <- DimPlot(Tcells, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

DefaultAssay(Tcells) <- "ADT"
p2 <- VlnPlot(Tcells, features = c("CD3", "CD4", "CD8"), pt.size = 0, ncol = 2)
p1 + p2

###############################
# Exclude MAIT and Gamma Delta
###############################
clust13 <- FindMarkers(Tcells, ident.1 = 13, only.pos = T, min.pct = 0.25, assay = "RNA") # MAIT confirmed
clust2 <- FindMarkers(Tcells, ident.1 = c(2,15), only.pos = T, min.pct = 0.25, assay = "RNA") # Gamma Delta confirmed

Tcells <- subset(Tcells, idents = c(2,13,15), invert = T)

###############################
# Recluster after exclusion
###############################
DefaultAssay(Tcells) <- "RNA"
Tcells <- Tcells %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "tcr1"), assay = "RNA") %>%
        RunPCA()

# Run Harmony
Tcells <- RunHarmony(Tcells, "Lane", assay.use = "RNA", max.iter.harmony = 30)

# Run Donwstream Analysis
Tcells <- Tcells %>% RunUMAP(reduction = 'harmony', dims = 1:15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:15) %>% 
  FindClusters(resolution = 1.5)

# Visualize
p1 <- DimPlot(Tcells, label = T) + NoLegend()
p2 <- DimPlot(Tcells, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

DefaultAssay(Tcells) <- "ADT"
p2 <- VlnPlot(Tcells, features = c("CD3", "CD4", "CD8"), stack = T, flip = T)
p1 + p2

###############################
# Separate CD4 and CD8 T cells
###############################
# CD3: All
# CD4: clusters 4, 6, 9, 11, 12, 14
# CD8 expression > 1.2 seems to differentiate cells quite well based on ridge plot

DefaultAssay(Tcells) <- "ADT"
Idents(Tcells) <- "seurat_clusters"
Idents(Tcells, WhichCells(Tcells, expression = CD8 > 1.2)) <- "CD8"
Idents(Tcells, WhichCells(Tcells, expression = CD8 <= 1.2)) <- "CD4"
Tcells[['Tcell_subset']] <- Idents(Tcells)

# 7. CD4 T Cells_____________________________________________________
###############################
# CD4 Analysis Round 1
###############################
Idents(Tcells) <- "Tcell_subset"
CD4 <- subset(Tcells, idents = "CD4")

# Recluster
DefaultAssay(CD4) <- "RNA"
CD4 <- CD4 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "tcr1"), assay = "RNA") %>%
        RunPCA()

# Run Harmony
CD4 <- RunHarmony(CD4, "Lane", assay.use = "RNA", max.iter.harmony = 30)

# Run Donwstream Analysis
CD4 <- CD4 %>% RunUMAP(reduction = 'harmony', dims = 1:15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:15) %>% 
  FindClusters(resolution = 1.0)

# Visualize
p1 <- DimPlot(CD4, label = T) + NoLegend()
p2 <- DimPlot(CD4, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

DefaultAssay(CD4) <- "ADT"
p2 <- VlnPlot(CD4, features = c("CD4", "CD27", "CD8", "CD45RA"), pt.size = 0, ncol = 2)
p1 + p2

# Exclusion: Cluster 7 doesn't have much CD4 expression so will exclude. Cluster 1 & 9 also have disproportionate negative CD4. 
CD4 <- subset(CD4, idents = 7, invert = T)
cytotoxic <- subset(CD4, idents = c(0))
DefaultAssay(cytotoxic) <- "ADT"
Cells <- WhichCells(cytotoxic, expression = CD4 <= 1.2)

CD4 <- subset(CD4, cells = Cells, invert = T)

###############################
# CD4 Analysis Round 2
###############################
# Recluster
DefaultAssay(CD4) <- "RNA"
CD4 <- CD4 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "tcr1"), assay = "RNA") %>%
        RunPCA()

# Run Harmony
CD4 <- RunHarmony(CD4, "Lane", assay.use = "RNA", max.iter.harmony = 30)

# Run Donwstream Analysis
CD4 <- CD4 %>% RunUMAP(reduction = 'harmony', dims = 1:10) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:10) %>% 
  FindClusters(resolution = 0.3)

# Visualize
p1 <- DimPlot(CD4, label = T) + NoLegend()
p2 <- DimPlot(CD4, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

DefaultAssay(CD4) <- "ADT"
p2 <- VlnPlot(CD4, features = c("CD4", "CD27", "CD57", "CD45RA", "CD69"), stack = T, flip = T)
p1 + p2

###############################
# Annotation
###############################

Idents(CD4) <- "seurat_clusters"
Idents(CD4, cells = WhichCells(CD4, idents = c(1))) <- "CD4 TEM"
Idents(CD4, cells = WhichCells(CD4, idents = c(2))) <- "CD4 Cytotoxic"
Idents(CD4, cells = WhichCells(CD4, idents = c(4))) <- "CD4 Regulatory"
Idents(CD4, cells = WhichCells(CD4, idents = c(3))) <- "CD4 Naive"
Idents(CD4, cells = WhichCells(CD4, idents = c(0))) <- "CD4 TCM"
CD4[['Annotation']] <- Idents(CD4)

markers <- FindAllMarkers(CD4, only.pos = T, min.pct = 0.1, assay = "RNA")
write.csv(markers, file = paste(marker_dir, "CD4.Markers.csv", sep = "/"))
###############################
# Save RDS
###############################
saveRDS(CD4, file = paste(subset_dir, "CD4_HIVneg.rds", sep = "/"))

#7. CD8 T Cells_____________________________________________________
###############################
# CD8 Analysis Round 1
###############################
Idents(Tcells) <- "Tcell_subset"
CD8 <- subset(Tcells, idents = "CD8")

# Recluster
DefaultAssay(CD8) <- "RNA"
CD8 <- CD8 %>% FindVariableFeatures(verbose = T, nfeatures = 3000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "tcr1"), assay = "RNA") %>%
        RunPCA()

# Run Harmony
CD8 <- RunHarmony(CD8, "Lane", assay.use = "RNA", max.iter.harmony = 30)

# Run Donwstream Analysis
CD8 <- CD8 %>% RunUMAP(reduction = 'harmony', dims = 1:15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:15) %>% 
  FindClusters(resolution = 0.3)

# Visualize
p1 <- DimPlot(CD8, label = T) + NoLegend()
p2 <- DimPlot(CD8, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

DefaultAssay(CD8) <- "ADT"
p2 <- VlnPlot(CD8, features = c("CD8", "CD27", "CD57", "CD45RA"), pt.size = 0, ncol = 2)
p1 + p2

###############################
# Annotation
###############################
Idents(CD8) <- "seurat_clusters"
Idents(CD8, cells = WhichCells(CD8, idents = c(1))) <- "CD8 TCM"
Idents(CD8, cells = WhichCells(CD8, idents = c(2))) <- "CD8 TEM"
Idents(CD8, cells = WhichCells(CD8, idents = c(0,4,5))) <- "CD8 Cytotoxic"
Idents(CD8, cells = WhichCells(CD8, idents = c(3))) <- "CD8 Naive"
CD8[['Annotation']] <- Idents(CD8)

markers <- FindAllMarkers(CD8, only.pos = T, min.pct = 0.1, assay = "RNA")
write.csv(markers, file = paste(marker_dir, "CD8.Markers.csv", sep = "/"))

DefaultAssay(CD8) <- "ADT"
p1 <- DimPlot(CD8, label = T) + NoLegend()
p2 <- VlnPlot(CD8, features = c("CD8", "CD27", "CD57", "CD45RA"), pt.size = 0, ncol = 2)
p1 + p2

###############################
# Save RDS
###############################
saveRDS(CD8, file = paste(subset_dir, "CD8_HIVneg.rds", sep = "/"))

