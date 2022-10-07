##############################################################################################
##------------------------- SUBSET ANALYSIS ---------------------------------##
##------------------------- DATE: 6/27/2022 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will clean the subsets (major cell lineages) and annotate cells
## based on identified gene markers. 
##############################################################################################

# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(future)
library(tidyverse)
library(harmony)
#library(anndata)
#library(leiden)

plan("multiprocess", workers = 5)
options(future.globals.maxSize = 5000 * 1024^2)
options(future.seed = TRUE)

set.seed(7412)

date = "6.27"
tmp <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Integrated_Filtered")

data_hatim <- readr::read_csv('/data/p_koethe_lab/Atlas_AT/MetaData/Reference/HATIMStudyVisit_DATA_2021-08-09_1359.csv')

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
                                    FBG = meta_fbg) %>%
                            filter(HIV == "HIVpos")

Hash_link <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis_1.22/Metadata/Hash_Link.csv")
Hash_link <- Hash_link %>% mutate(Hash = as.character(Hash))
colnames(Hash_link) <- c("Lane", "HATIMID", "Soup_Hash")

#2. MYELOID ANALYSIS_____________________________________________________
###############################
# Import Data
###############################
Myeloid <- readRDS(paste(tmp, "Myeloid.rds", sep = "/"))

################################
# Plot Doublets (clusters 84, 82, 66, 67)
###############################
DotPlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "GNG11", "CCDC80", "LUM")) + RotatedAxis()
p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- FeaturePlot(Myeloid, features = c("CD3E", "NKG7"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

# Remove Clusters 84,82,66,67
Myeloid <- subset(Myeloid, idents = c(66,67,82,84), invert = T)

# All clusters labeled heterotypic were removed above
DotPlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "CCDC80"), group.by = "Heterotypic") + RotatedAxis()

# Junk (clusters 17 and 74)
Idents(Myeloid) <- "seurat_clusters"
VlnPlot(Myeloid, features = "percent.mt", pt.size = 0) + NoLegend()

################################
# Round 1: n = 45,903 cells
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

# Doublets: cluster 18,19,22,13
p2 <- FeaturePlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "CCDC80", "MKI67"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: cluster 17, unreasonably high compared with other cell types
p2 <- VlnPlot(Myeloid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Remove Doublets
###############################
Myeloid <- subset(Myeloid, idents = c(13,17,18,19,22), invert = T)

################################
# Round 2: n = 41,293 cells
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
# Junk: cluster 23
p1 <- DimPlot(Myeloid, label = T) + NoLegend()
p2 <- VlnPlot(Myeloid, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: None (some RNA contamination of nMo)
p2 <- FeaturePlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "CCDC80"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: Cluster 18
p2 <- VlnPlot(Myeloid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Remove Doublets
###############################
Myeloid <- subset(Myeloid, idents = c(18,23), invert = T)

################################
# Round 3: n = 40,585 cells
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

# Doublets: None
p2 <- FeaturePlot(Myeloid, features = c("CD3E", "NKG7", "CLDN5", "CCDC80"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: None
p2 <- VlnPlot(Myeloid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2
################################
# Cell Annotation
###############################
markers <- FindAllMarkers(Myeloid, only.pos = T, min.pct = 0.25, assay = "RNA")
write.csv(markers, paste(tmp, "Myeloid.Markers.csv", sep = "/"))

Idents(Myeloid) <- "seurat_clusters"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(1,3,4))) <- "IM"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(0,5))) <- "cMo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(6))) <- "Other Mo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(7))) <- "nMo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(2,10))) <- "PVM"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(8,11))) <- "LAM"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(14))) <- "Cycling Myeloid"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(15))) <- "ISG+ Mo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(17))) <- "cMo"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(18,22))) <- "Other Mac"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(19))) <- "pDC"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(20))) <- "Migratory DC"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(12,13,21))) <- "Mo-Mac"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(9,16))) <- "cDC2B"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c(23))) <- "cDC1"
Myeloid[['Annotation']] <- Idents(Myeloid)

DotPlot(Myeloid, features = c("TREM2", "APOE", "CXCL2", "CXCL3", "LYVE1", "RNASE1", "MKI67", "TIMP3", "CD1C", "CLEC9A", "CCR7", "LILRA4", "VCAN", "S100A8", "FCGR3A", "ISG15")) + 
RotatedAxis()

################################
# Add Metadata
###############################
metadata <- Myeloid@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG)
rownames(metadata) <- colnames(Myeloid)

Myeloid <- AddMetaData(Myeloid, metadata)

################################
# Save Object
###############################
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")
dir.create(tmp_dir)

saveRDS(Myeloid, paste(tmp_dir, "Myeloid.rds", sep = "/"))

#2. LYMPHOID ANALYSIS_____________________________________________________
###############################
# Import Data
###############################
tmp <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Integrated_Filtered")
Lymphoid <- readRDS(paste(tmp, "Lymphoid.rds", sep = "/"))

###############################
# PLOT Doublets
###############################
DotPlot(Lymphoid, features = c("LYZ", "CD68", "CLDN5", "GNG11", "CCDC80", "LUM")) + RotatedAxis()
# Cluster 71

# Remove Cluster 71
Lymphoid <- subset(Lymphoid, idents = 71, invert = T)

# No cells labeled heterotypic after removal of cluster 71
DotPlot(Lymphoid, features = c("LYZ", "CLDN5", "CCDC80"), group.by = "Heterotypic") + RotatedAxis() # Left over heterotypic are Vascular and will remove

# Hold off on removing clusters for now
VlnPlot(Lymphoid, features = "percent.mt", pt.size = 0) + NoLegend()

###############################
# ROUND 1: n = 31,563 cells
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
# Junk: Cluster 17
p1 <- DimPlot(Lymphoid, label = T) + NoLegend()
p2 <- VlnPlot(Lymphoid, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 9, 14 and 18 (cluster 22 is likely HPSC so will remove as well)
DotPlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67")) + RotatedAxis()
p2 <- FeaturePlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: None
p2 <- VlnPlot(Lymphoid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Remove Doublets
###############################
Lymphoid <- subset(Lymphoid, idents = c(9,14,17,18,22), invert = T)

###############################
# ROUND 2: n = 28,433 cells
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

# Doublets: Cluster 16 (HBB)
DotPlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67", "HBB")) + RotatedAxis()
p2 <- FeaturePlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: None
p2 <- VlnPlot(Lymphoid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Remove Doublets
###############################
Lymphoid <- subset(Lymphoid, idents = c(16), invert = T)

###############################
# ROUND 3: n = 28,061 cells
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
  FindClusters(resolution = 2.0)
  
###############################
# Visualize
###############################
# Junk: None
p1 <- DimPlot(Lymphoid, label = T) + NoLegend()
p2 <- VlnPlot(Lymphoid, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: None
DotPlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67")) + RotatedAxis()
p2 <- FeaturePlot(Lymphoid, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "MKI67"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

#  RNA Count: None
p2 <- VlnPlot(Lymphoid, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Add Annotation Labels
###############################
markers <- FindAllMarkers(Lymphoid, only.pos = T, min.pct = 0.25, assay = "RNA")
write.csv(markers, file = paste(tmp, "Lymphoid.Markers.csv", sep = "/"))

DefaultAssay(Lymphoid) <- "ADT"
p2 <- VlnPlot(Lymphoid, features = c("CD8", "CD4", "CD27", "CD57", "CD45RA", "CD3"), pt.size = 0, ncol = 2) + NoLegend()
p1 + p2


cluster0 <- subset(Lymphoid, idents = 0)
DefaultAssay(cluster0) <- "ADT"
cluster8 <- subset(Lymphoid, idents = 8)
DefaultAssay(cluster8) <- "ADT"
cluster17 <- subset(Lymphoid, idents = 17)
DefaultAssay(cluster17) <- "ADT"

DefaultAssay(Lymphoid) <- "RNA"
Idents(Lymphoid) <- "seurat_clusters"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(1))) <- "CD4 TCM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(2))) <- "CD57 mNK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(3))) <- "CD4 TEM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(4))) <- "CD8 TCM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(5))) <- "CD16 mNK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(6))) <- "CD8 TEM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(7))) <- "CD57 mNK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(8))) <- "CD8 Cytotoxic"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(9))) <- "CD57 mNK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(10))) <- "Gamma Delta"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(11))) <- "CD8 Cytotoxic"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(12))) <- "Immature NK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(13))) <- "CD8 Cytotoxic"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(14))) <- "CD4 Naive"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(15))) <- "CD8 Naive"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(16))) <- "CD8 TCM"
Idents(Lymphoid, cells = WhichCells(cluster17, expression = CD8 > 1.0)) <- "CD8 TEM"
Idents(Lymphoid, cells = WhichCells(cluster17, expression = CD8 <= 1.0)) <- "CD4 TEM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(18))) <- "CD8 TEM"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(19))) <- "Cycling T/NK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(20))) <- "MAIT"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(21))) <- "CD4 Regulatory"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(22))) <- "CD16 mNK"
Idents(Lymphoid, cells = WhichCells(Lymphoid, idents = c(23))) <- "ILC"
Idents(Lymphoid, cells = WhichCells(cluster0, expression = CD8 > 1.0)) <- "CD8 Cytotoxic"
Idents(Lymphoid, cells = WhichCells(cluster0, expression = CD8 <= 1.0)) <- "CD4 Cytotoxic"
Idents(Lymphoid, cells = WhichCells(cluster8, expression = CD8 > 1.0)) <- "CD8 Cytotoxic"
Idents(Lymphoid, cells = WhichCells(cluster8, expression = CD8 <= 1.0)) <- "CD4 Cytotoxic"
Lymphoid[['Annotation']] <- Idents(Lymphoid)

###############################
# Verify Annotations
###############################
DefaultAssay(Lymphoid) <- "ADT"
VlnPlot(Lymphoid, features = c("CD4", "CD8", "CD27", "CD57", "CD45RA", "CD16"), pt.size = 0, ncol = 2)

################################
# Add Metadata
###############################
metadata <- Lymphoid@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG)
rownames(metadata) <- colnames(Lymphoid)

Lymphoid <- AddMetaData(Lymphoid, metadata)

###############################
# Save Harmony Integration
###############################
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")
dir.create(tmp_dir)

saveRDS(Lymphoid, file = paste0(tmp_dir, "/", "Lymphoid.rds"))

#3. B Cell Analysis_____________________________________________________
###############################
# Import Data
###############################
tmp <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Integrated_Filtered")
Bcells <- readRDS(paste(tmp, "Bcells.rds", sep = "/"))

###############################
# PLOT Doublets
###############################
DotPlot(Bcells, features = c("LYZ", "CD68", "CLDN5", "GNG11", "CCDC80", "LUM", "CD3E", "NKG7")) + RotatedAxis()

# Remove Cluster 89, 95, 96
Bcells <- subset(Bcells, idents = c(89,95,96), invert = T)

VlnPlot(Bcells, features = "percent.mt", pt.size = 0) + NoLegend()

###############################
# ROUND 1: n = 2,993 cells
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
# Junk: None
p1 <- DimPlot(Bcells, label = T) + NoLegend()
p2 <- VlnPlot(Bcells, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 9, 11, 13, 14
DotPlot(Bcells, features = c("CD68", "LYZ", "CLDN5", "CCDC80", "CD3E", "NGK7", "TAGLN")) + RotatedAxis()

#  RNA Count: None
p2 <- VlnPlot(Bcells, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk and Doublets
###############################
Bcells <- subset(Bcells, idents = c(9,11,13,14), invert = T)

###############################
# ROUND 2: n = 2679 cells
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
# Junk: None
p1 <- DimPlot(Bcells, label = T) + NoLegend()
p2 <- VlnPlot(Bcells, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublet: None
p2 <- FeaturePlot(Bcells, features = c("CLDN5", "LYZ", "C1QB", "CCDC80", "NKG7", "CD3E"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

###############################
# Annotation
###############################
markers <- FindAllMarkers(Bcells, only.pos = T, min.pct = 0.25, assay = "RNA")
write.csv(markers, file = paste(tmp, "Bcells.Markers.csv", sep = "/"))

DefaultAssay(Bcells) <- "ADT"
VlnPlot(Bcells, features = c("CD19", "CD20", "CD38"), pt.size = 0, ncol = 2)

DefaultAssay(Bcells) <- "RNA"
Idents(Bcells) <- "seurat_clusters"
Idents(Bcells, cells = WhichCells(Bcells, idents = c(2))) <- "Plasmablast"
Idents(Bcells, cells = WhichCells(Bcells, idents = c(0,1,3))) <- "B Cell"
Bcells[['Annotation']] <- Idents(Bcells)

###############################
# Add Metadata
###############################
metadata <- Bcells@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG)
rownames(metadata) <- colnames(Bcells)

Bcells <- AddMetaData(Bcells, metadata)

###############################
# Save Harmony Integration
###############################
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")
dir.create(tmp_dir)

saveRDS(Bcells, file = paste0(tmp_dir, "/", "Bcells.rds"))

#4. Stromal Analysis_____________________________________________________
###############################
# Import Data
###############################
tmp <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Integrated_Filtered")
Stromal <- readRDS(paste(tmp, "Stromal.rds", sep = "/"))

###############################
# PLOT Doublets
###############################
# Doublets: Clusters 59, 72, 86, 97 81 (Mature adipocyte)
DotPlot(Stromal, features = c("LYZ", "CD68", "CLDN5", "GNG11", "CD3E", "NKG7", "TAGLN", "MKI67")) + RotatedAxis()
p1 <- DimPlot(Stromal, label = T) + NoLegend()
p2 <- FeaturePlot(Stromal, features = c("CD3E", "NKG7", "CLDN5", "TAGLN", "CD68", "LYZ"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

# Remove Doublets
Stromal <- subset(Stromal, idents = c(59,72,81,86,97), invert = T)

# All clusters labeled heterotypic were removed above
DotPlot(Stromal, features = c("CD3E", "NKG7", "CLDN5", "CCDC80"), group.by = "Heterotypic") + RotatedAxis()

# Junk:Cluster 50 (hold off on removing for now)
Idents(Stromal) <- "seurat_clusters"
VlnPlot(Stromal, features = "percent.mt", pt.size = 0) + NoLegend()

###############################
# ROUND 1: n = 57,993 cells
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
# Junk: Cluster 23, 32
p1 <- DimPlot(Stromal, label = T) + NoLegend()
p2 <- VlnPlot(Stromal, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 20,28,30,33
DotPlot(Stromal, features = c("CLDN5", "LYZ", "NKG7", "CD3E", "TAGLN")) + RotatedAxis()
p2 <- FeaturePlot(Stromal, features = c("CLDN5", "LYZ", "NKG7", "CD3E", "TAGLN"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

# Ribosomal: Cluster 24
p2 <- DotPlot(Stromal, features = c("RPL11", "RPS8", "RPS27", "RPL12", "RPL30")) + RotatedAxis()
p1 + p2

# RNA Count (low: 23,32; high: 34,17 [cycling fibroblast]
p2 <- VlnPlot(Stromal, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk, Doublets
###############################
Stromal <- subset(Stromal, idents = c(20,23,24,28,30,32,33), invert = T)

###############################
# ROUND 2: n = 52,493 cells
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
# Junk: Cluster 26
p1 <- DimPlot(Stromal, label = T) + NoLegend()
p2 <- VlnPlot(Stromal, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: None (Clusters 16 and 17 express TAGLN but this seems like stem-cell expression)
DotPlot(Stromal, features = c("CLDN5", "LYZ", "NKG7", "CD3E", "TAGLN")) + RotatedAxis()
p2 <- FeaturePlot(Stromal, features = c("CLDN5", "LYZ", "NKG7", "CD3E"), min.cutoff = 'q10', max.cutoff = 'q90')
p1 + p2

# Ribosomal: None
p2 <- DotPlot(Stromal, features = c("RPL11", "RPS8", "RPS27", "RPL12", "RPL30")) + RotatedAxis()
p1 + p2

# RNA Count (low: 26, 29)
p2 <- VlnPlot(Stromal, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk, Doublets: Cluster 27 passes QC but does not form distinct cluster so will remove. May try leiden at some point to see if that resolves this issue
###############################
Stromal <- subset(Stromal, idents = c(26,29,27), invert = T)

###############################
# ROUND 3: n = 51,029 cells
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
  
# QC: PASS
p1 <- DimPlot(Stromal, label = T) + NoLegend()
p2 <- VlnPlot(Stromal, features = "percent.mt", pt.size = 0) + NoLegend()
p3 <- DotPlot(Stromal, features = c("CLDN5", "LYZ", "NKG7", "CD3E", "TAGLN")) + RotatedAxis()
p4 <- VlnPlot(Stromal, features = "nCount_RNA", pt.size = 0) + NoLegend()
  
###############################
# Add Metadata
###############################
metadata <- Stromal@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG)
rownames(metadata) <- colnames(Stromal)

Stromal <- AddMetaData(Stromal, metadata)

###############################
# Annotation
###############################
DotPlot(Stromal, features = c("CD55", "MFAP5", "MYOC", "POSTN", "MKI67", "TIMP3", "ISG15", "MT1X", "FABP5", "FABP4", "C7", "LPL", "EGR1", "JUN", "FOS", "DCN", "LUM", "PI16", "GSN", "CFD", "APOD")) + RotatedAxis()

Idents(Stromal) <- "seurat_clusters"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(12))) <- "PCOLCE+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(6))) <- "MYOC+ Fibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(7,14))) <- "Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(16))) <- "Cycling Myofibroblast"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(15))) <- "ISG+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(11))) <- "Metallothionein+ Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(8))) <- "Mature Preadipocyte 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(2))) <- "Mature Preadipocyte 2"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(5))) <- "Adipose Progenitor Cell 1"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(0,1,13))) <- "Early Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(3,4))) <- "ECM-Producing Early Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c(9,10))) <- "Adipose Progenitor Cell 2"
Stromal[['Annotation']] <- Idents(Stromal)

###############################
# Save RDS
###############################
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")
dir.create(tmp_dir)

saveRDS(Stromal, file = paste0(tmp_dir, "/", "Stromal.rds"))

#5. Vascular Analysis_____________________________________________________
###############################
# Import Data
###############################
tmp <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Integrated_Filtered")
Vascular <- readRDS(paste(tmp, "Vascular.rds", sep = "/"))

###############################
# PLOT Doublets
###############################
# Clusters 87
p1 <- DimPlot(Vascular, label = T) + NoLegend()
DotPlot(Vascular, features = c("LYZ", "CD68", "CLDN5", "GNG11", "CD3E", "NKG7", "MKI67")) + RotatedAxis()

# Cluster 78 is Junk
VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()

# RNA count High: 76, 42
VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()

# Removing Doublets (76, 87), Junk (78), and implausibly high nRNA_count (42,76)
Vascular <- subset(Vascular, idents = c(42,76,78,87), invert = T)

###############################
# ROUND 1: n = 61,886 cells
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
# Junk: Clusters 35, 36, and 31
p1 <- DimPlot(Vascular, label = T) + NoLegend()
p2 <- VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublets: Clusters 19, 29, 30
p2 <- DotPlot(Vascular, features = c("CCDC80", "LUM", "LYZ", "NKG7", "CD3E", "MS4A1", "MKI67")) + RotatedAxis()
p1 + p2

# RNA Count High: 27, 32 (28, 34 both have higher RNA count so will hold off on removing)
p2 <- VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk, Doublets
###############################
Vascular <- subset(Vascular, idents = c(19,27,29,30,31,35,36), invert = T)

###############################
# ROUND 2: n = 57,098 cells
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
  FindClusters(resolution = 2.5)
  
###############################
# Visualize
###############################
# Junk: Pass QC
p1 <- DimPlot(Vascular, label = T) + NoLegend()
p2 <- VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Clean
p2 <- DotPlot(Vascular, features = c("CCDC80", "LUM", "LYZ", "NKG7", "CD3E", "MS4A1", "MKI67")) + RotatedAxis()
p1 + p2

# RNA Count High: 30 implausibly high, 26 and 29 are high but so is 16 so may be biological. 
p2 <- VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk, Doublets
###############################
Vascular <- subset(Vascular, idents = c(30), invert = T)

###############################
# ROUND 3: n = 56,709 cells
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
Vascular <- Vascular %>% RunUMAP(reduction = 'harmony', dims = 1:15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:15) %>% 
  FindClusters(resolution = 2.5)
  
###############################
# Visualize
###############################
# Junk: None
p1 <- DimPlot(Vascular, label = T) + NoLegend()
p2 <- VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Doublet: None
p2 <- DotPlot(Vascular, features = c("CCDC80", "LUM", "LYZ", "NKG7", "CD3E", "MS4A1")) + RotatedAxis()
p1 + p2

# RNA Count High: 35 (29 and 30 also high but near other clusters with high counts)
p2 <- VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

p2 <- VlnPlot(Vascular, features = "nFeature_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk, Doublets
###############################
Vascular <- subset(Vascular, idents = c(35), invert = T)

###############################
# ROUND 4: n = 56,645 cells
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
Vascular <- Vascular %>% RunUMAP(reduction = 'harmony', dims = 1:15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:15) %>% 
  FindClusters(resolution = 1.0)
  
###############################
# Visualize
###############################
# Pass QC
p1 <- DimPlot(Vascular, label = T) + NoLegend()
p2 <- VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Clust 17
p2 <- DotPlot(Vascular, features = c("CCDC80", "LUM", "LYZ", "NKG7", "CD3E", "MS4A1"))
p1 + p2

# Pass QC
p2 <- VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Exclude Junk, Doublets
###############################
Vascular <- subset(Vascular, idents = c(17), invert = T)

###############################
# FINAL ROUND: n = 56,267 cells
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
Vascular <- Vascular %>% RunUMAP(reduction = 'harmony', dims = 1:15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:15) %>% 
  FindClusters(resolution = 1.0)
  
###############################
# Visualize
###############################
# Pass QC
p1 <- DimPlot(Vascular, label = T) + NoLegend()
p2 <- VlnPlot(Vascular, features = "percent.mt", pt.size = 0) + NoLegend()
p1 + p2

# Pass QC
p2 <- DotPlot(Vascular, features = c("CCDC80", "LUM", "LYZ", "NKG7", "CD3E", "MS4A1"))
p1 + p2

# Pass QC
p2 <- VlnPlot(Vascular, features = "nCount_RNA", pt.size = 0) + NoLegend()
p1 + p2

###############################
# Add Metadata
###############################
metadata <- Vascular@meta.data
metadata <- left_join(metadata, Hash_link, by = c("Soup_Hash", "Lane"))
metadata$HATIMID <- as.factor(metadata$HATIMID)
data_hatim <- droplevels(data_hatim)
metadata <- left_join(metadata, data_hatim[, c("HATIMID", "StudyGroup", "HIV", "Glucose", "Sex", "homa2_ir", "Age", "BMI", "FBG")], by = "HATIMID")
metadata <- metadata %>% select(HATIMID, StudyGroup, HIV, Glucose, Sex, homa2_ir, Age, BMI, FBG)
rownames(metadata) <- colnames(Vascular)

Vascular <- AddMetaData(Vascular, metadata)

###############################
# Add Annotations
###############################
p2 <- DotPlot(Vascular, features = c("GJA5", "GJA4", "HEY1", "SEMA3G", "CA4", "VWF", "GNG11", "ACKR1", "ACTA2", "TAGLN", "MYH11", "NOTCH3", "PDGFRB", "MFAP5",
"RGS5", "COX4I2", "CSPG4", "SNAI1", "SNAI2")) + RotatedAxis()
p1 + p2

Idents(Vascular) <- "seurat_clusters"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(8))) <- "Arterial EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(9,10))) <- "VSMC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(3))) <- "Pericyte"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(2,6))) <- "Capillary EndoMT-like"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(11))) <- "Arterial EndoMT-like"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(16))) <- "Cycling Vascular Cells"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(14))) <- "Venous EndoMT-like"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(4,15))) <- "Venous EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(12))) <- "Ven-Cap EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(0,7,13))) <- "Capillary EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c(1,5))) <- "Intermediate Capillary EC"
Vascular[['Annotation']] <- Idents(Vascular)

###############################
# Save
###############################
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")
dir.create(tmp_dir)

saveRDS(Vascular, file = paste0(tmp_dir, "/", "Vascular.rds"))

# 6. T Cells_____________________________________________________
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))


# Isolate T Cells: To do this, I will grab clusters that are identified as T cells and are expressing CD3 ADT. This will include MAIT and Gamma Delta initially
Tcells <- subset(Lymphoid, idents = c("CD4 Cytotoxic", "CD8 Cytotoxic", "CD4 Regulatory", "MAIT", "CD8 TEM", "CD4 TEM", "CD8 TCM", "CD8 Naive", "CD4 Naive", "Gamma Delta", "CD4 TCM"))

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
clust14 <- FindMarkers(Tcells, ident.1 = 14, only.pos = T, min.pct = 0.25) # MAIT confirmed
clust6 <- FindMarkers(Tcells, ident.1 = 6, only.pos = T, min.pct = 0.25) # Gamma Delta confirmed

Tcells <- subset(Tcells, idents = c(6,14), invert = T)

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
p2 <- VlnPlot(Tcells, features = c("CD3", "CD4", "CD8"), pt.size = 0, ncol = 2)
p1 + p2

###############################
# Separate CD4 and CD8 T cells
###############################
# CD4: clusters 10,13,8,6,5, and 15 pure
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
CD4 <- CD4 %>% RunUMAP(reduction = 'harmony', dims = 1:12) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:12) %>% 
  FindClusters(resolution = 1.0)

# Visualize
p1 <- DimPlot(CD4, label = T) + NoLegend()
p2 <- DimPlot(CD4, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

DefaultAssay(CD4) <- "ADT"
p2 <- VlnPlot(CD4, features = c("CD4", "CD27", "CD8", "CD45RA"), pt.size = 0, ncol = 2)
p1 + p2

# Exclusion: Cluster 6 doesn't have much CD4 expression so will exclude. Cluster 1 & 9 also have disproportionate negative CD4. 
CD4 <- subset(CD4, idents = 6, invert = T)
cytotoxic <- subset(CD4, idents = c(1,9))
DefaultAssay(cytotoxic) <- "ADT"
Cells <- WhichCells(cytotoxic, expression = CD4 <= 1.0)

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
  FindClusters(resolution = 0.4)

# Visualize
p1 <- DimPlot(CD4, label = T) + NoLegend()
p2 <- DimPlot(CD4, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

DefaultAssay(CD4) <- "ADT"
p2 <- VlnPlot(CD4, features = c("CD4", "CD27", "CD57", "CD45RA"), pt.size = 0, ncol = 2)
p1 + p2

###############################
# Annotation
###############################
markers <- FindAllMarkers(CD4, only.pos = T, min.pct = 0.25, assay = "RNA")

Idents(CD4) <- "seurat_clusters"
Idents(CD4, cells = WhichCells(CD4, idents = c(1,5))) <- "CD4 TEM"
Idents(CD4, cells = WhichCells(CD4, idents = c(2))) <- "CD4 Cytotoxic"
Idents(CD4, cells = WhichCells(CD4, idents = c(4))) <- "CD4 Regulatory"
Idents(CD4, cells = WhichCells(CD4, idents = c(3))) <- "CD4 Naive"
Idents(CD4, cells = WhichCells(CD4, idents = c(0))) <- "CD4 TCM"
CD4[['Annotation']] <- Idents(CD4)

###############################
# Save RDS
###############################
saveRDS(CD4, file = paste(tmp_dir, "CD4.rds", sep = "/"))

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
  FindClusters(resolution = 0.4)

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
Idents(CD8, cells = WhichCells(CD8, idents = c(0,4))) <- "CD8 Cytotoxic"
Idents(CD8, cells = WhichCells(CD8, idents = c(3))) <- "CD8 Naive"
CD8[['Annotation']] <- Idents(CD8)

###############################
# Save RDS
###############################
saveRDS(CD8, file = paste(tmp_dir, "CD8.rds", sep = "/"))

#8. Macrophage_____________________________________________________
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))

Macrophage <- subset(Myeloid, idents = c("PVM", "IM", "LAM", "Mo-Mac"))

###############################
# Macrophage Round 1
###############################
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

p2 <- VlnPlot(Macrophage, features = 'nCount_RNA', pt.size = 0) + NoLegend()
p1 + p2

# Markers
markers <- FindAllMarkers(Macrophage, only.pos = T, min.pct = 0.25, assay = "RNA")

# Remove cluster 12 (more dendritic), and cluster 13 (high RNA count)
Macrophage <- subset(Macrophage, idents = c(12,13), invert = T)

###############################
# Macrophage Round 2
###############################
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

# Visualize
p1 <- DimPlot(Macrophage, label = T) + NoLegend()
p2 <- DimPlot(Macrophage, label = T, group.by = "Annotation") + NoLegend()
p1 + p2

p2 <- VlnPlot(Macrophage, features = 'percent.mt', pt.size = 0) + NoLegend()
p1 + p2

p2 <- VlnPlot(Macrophage, features = 'nCount_RNA', pt.size = 0) + NoLegend()
p1 + p2

p2 <- VlnPlot(Macrophage, features = c("TREM2", "CXCL3", "C1QB", "LYVE1"), pt.size = 0, ncol = 2)
p1 + p2

###############################
# Annotation
###############################
Idents(Macrophage) <- "seurat_clusters"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(0))) <- "IM"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(1))) <- "PVM"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(2))) <- "LAM"
Idents(Macrophage, cells = WhichCells(Macrophage, idents = c(3))) <- "Mo-Mac"
Macrophage[["Annotation"]] <- Idents(Macrophage)

###############################
# Save RDS
###############################
saveRDS(Macrophage, file = paste(tmp_dir, "Macrophage.rds", sep = "/"))



