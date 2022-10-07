##############################################################################################
##------------------------- DOUBLET REMOVAL ---------------------------------##
##------------------------- DATE: 6/27/2022 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: After applying  DoubletFinder, the parameters that take into
## account the homotypic doublets performed the best with 98.6% specificity for ground truth doublets
## (defined as the final dataset). Additionally, the 4,000 cell discrepancy did not remove any one
## cell type. Therefore, doublets identified by DoubletFinder, Genetic demultiplexing, or clusters 
## with > 60% of cells identified as doublets will be remove prior to further processing.
##############################################################################################

# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(future)
library(tidyverse)

plan("multiprocess", workers = 5)
options(future.globals.maxSize = 5000 * 1024^2)
options(future.seed = TRUE)

set.seed(7412)

path = "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis"
date = "6.27"

#2. LOAD DATA_____________________________________________________
Doublet <- readRDS("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Merged_Initial/Integrated_HIVnegative.rds")

#3. REMOVE DOUBLETS_____________________________________________________
###############################
# Doublet Removal
###############################
Doublet[['DF_GT']] <- ifelse(Doublet$GT == "Doublet" | Doublet$DoubletFinder == "Doublet", "Doublet", "Singlet") # Combine GT and DF
metadata <- data.frame(Doublet@meta.data)
Doublet_clusters <- metadata %>% group_by(seurat_clusters, DF_GT) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>% filter(DF_GT == "Doublet" & freq > 0.6) %>% select(seurat_clusters) # Filter out clusters with > 60% doublets
Doublet[['Subset_DF']] <- ifelse(Doublet$DF_GT == "Doublet" | Doublet$seurat_clusters %in% Doublet_clusters, "Exclude", "Keep") # Remove clusters with > 60% doublet cells

# Plot (looks reasonable)
p1 <- DimPlot(Doublet, group.by = "Subset_DF", raster = F)
p2 <- FeaturePlot(Doublet, features = c("LYZ", "CD68", "CD3E", "NKG7", "CLDN5", "CCDC80", "LUM"), min.cutoff = 'q10', max.cutoff = 'q90', raster = F)
p3 <- DimPlot(Doublet, group.by = "seurat_clusters", raster = F, label = T) + NoLegend()
p2 + p1 + p3

DotPlot(Doublet, features = c("LYZ", "CD68", "CD3E", "NKG7", "CLDN5", "GNG11", "MS4A1", "JCHAIN", "CCDC80", "LUM", "LILRA4")) + RotatedAxis() + coord_flip()

# 26,31,34,36,40,43,48,49,50,51,52,55,69,68,67,62,57
Idents(Doublet) <- "Subset_DF"
Singlet <- subset(Doublet, idents = "Keep")
rm(Doublet)

#4. TAG JUNK_____________________________________________________
###############################
# Overcluster
###############################
# Overcluster to find MT clusters
Singlet <- FindNeighbors(Singlet, dims = 1:30)
Singlet <- FindClusters(Singlet, resolution = 5.0)
DimPlot(Singlet, group.by = "seurat_clusters", label = T, raster = F) + NoLegend()

metadata <- Singlet@meta.data
###############################
# Evaluate for MT Junk
###############################
VlnPlot(Singlet, features = c("percent.mt"), pt.size = 0) + NoLegend()
junk_clusters <- metadata %>% group_by(seurat_clusters) %>% summarise(percent.mt = median(percent.mt)) %>% filter(percent.mt > 8.0) %>% select(seurat_clusters) # clusters with > 8% doublets

###############################
# Evaluate for Heterotypic Doublets
###############################
DotPlot(Singlet, features = c("LYZ", "CD68", "CD3E", "NKG7", "CLDN5", "GNG11", "MS4A1", "JCHAIN", "CCDC80", "LUM", "LILRA4")) + RotatedAxis() + coord_flip()
# Clusters 72, 73, 75,76, 82,95
###############################
# Tag Junk and Doublets
###############################
Singlet[['Junk']] <- ifelse(Singlet$seurat_clusters %in% c(33,34,87,88), "Junk", "Not Junk")
Singlet[['Heterotypic']] <- ifelse(Singlet$seurat_clusters %in% c(72,73,75,76,82,95), "Heterotypic", "Singlet")

#5. SAVE_____________________________________________________
# Save this dataset then perform subset analysis
tmp <- paste(path, date, "Integrated_Filtered_HIVneg", sep = "/")
dir.create(tmp)

saveRDS(Singlet, file = paste(tmp, "Singlet_Integrated.RDS", sep = "/"))

#6. SUSBSET ANALYSIS_____________________________________________________
###############################
# Save Subsets (except 88[mature adipocyte] and 82)
###############################
p1 <- DimPlot(Singlet, group.by = "CellType", raster = F, label = T) + NoLegend()
p2 <- DimPlot(Singlet, raster = F, label = T) + NoLegend()
p1 + p2
p3 <- DimPlot(Singlet, raster = F, label = T, repel = T) + NoLegend()

Lymphoid <- subset(Singlet, idents = c(86,41,30,69,51,42,55,76,58,15,39,73,29)) #13
saveRDS(Lymphoid, file = paste(tmp, "Lymphoid_HIVneg.rds", sep = "/"))

Stromal <- subset(Singlet, idents = c(90,7,6,9,2,36,5,0,81,1,57,4,56,34,10,75,85,46,37,8,50)) #21
saveRDS(Stromal, file = paste(tmp, "Stromal_HIVneg.rds", sep = "/"))

Myeloid <- subset(Singlet, idents = c(60,84,94,13,87,31,18,44,47,25,77,43,95,26,33,71,49,54,45,48,16,23,79,40,67,70,80,83)) #28 clusters
saveRDS(Myeloid, file = paste(tmp, "Myeloid_HIVneg.rds", sep = "/"))

Bcells <- subset(Singlet, idents = c(53,89)) #2
saveRDS(Bcells, file = paste(tmp, "Bcells_HIVneg.rds", sep = "/"))

Vascular <- subset(Singlet, idents = c(92,66,38,21,91,61,35,93,14,78,28,20,11,59,24,52,65,12,19,72,3,74,32,22,17,64,63,62,27,68)) # 30
saveRDS(Vascular, file = paste(tmp, "Vascular_HIVneg.rds", sep = "/"))





