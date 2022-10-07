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
Doublet <- readRDS(paste(path, date, "Merged_Initial/Integrated_Doublets.rds", sep = "/"))

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

# Clusters 17,50,78,81 expressing high percent mitochondria

###############################
# Evaluate for Heterotypic Doublets
###############################
DotPlot(Singlet, features = c("LYZ", "CD68", "CD3E", "NKG7", "CLDN5", "GNG11", "MS4A1", "JCHAIN", "CCDC80", "LUM", "LILRA4")) + RotatedAxis() + coord_flip()
# Clusters 97,96,95,89,84,82,72,71,67,66,59
###############################
# Tag Junk and Doublets
###############################
Singlet[['Junk']] <- ifelse(Singlet$seurat_clusters %in% c(17,50,78,81), "Junk", "Not Junk")
Singlet[['Heterotypic']] <- ifelse(Singlet$seurat_clusters %in% c(59,66,67,71,72,82,84,89,95,96,97), "Heterotypic", "Singlet")

#5. SAVE_____________________________________________________
# Save this dataset then perform subset analysis
tmp <- paste(path, date, "Integrated_Filtered", sep = "/")
dir.create(tmp)

saveRDS(Singlet, file = paste(tmp, "Singlet_Integrated.RDS", sep = "/"))

#6. SUSBSET ANALYSIS_____________________________________________________
###############################
# Save Subsets
###############################
Lymphoid <- subset(Singlet, idents = c(6,14,25,65,35,46,36,28,55,27,91,94,85,31,71))
saveRDS(Lymphoid, file = paste(tmp, "Lymphoid.rds", sep = "/"))

Stromal <- subset(Singlet, idents = c(90,18,21,12,51,7,37,54,3,58,50,59,86,97,81,72,1,4,32,68,29,10,43,22))
saveRDS(Stromal, file = paste(tmp, "Stromal.rds", sep = "/"))

Myeloid <- subset(Singlet, idents = c(56,74,75,5,19,82,40,38,82,73,83,24,41,13,53,30,17,84,44,67,48,34,77,70,79,66,80))
saveRDS(Myeloid, file = paste(tmp, "Myeloid.rds", sep = "/"))

Bcells <- subset(Singlet, idents = c(33,88,89,95,96))
saveRDS(Bcells, file = paste(tmp, "Bcells.rds", sep = "/"))

Vascular <- subset(Singlet, idents = c(62,45,87,2,20,49,26,9,61,42,39,76,15,60,69,52,93,16,64,8,0,78,11,47,57,23,63,92))
saveRDS(Vascular, file = paste(tmp, "Vascular.rds", sep = "/"))





