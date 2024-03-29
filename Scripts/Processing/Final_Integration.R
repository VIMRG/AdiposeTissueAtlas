##############################################################################################
##------------------------- FINAL INTEGRATION ---------------------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will merge all of the subsets that have been clean.
##############################################################################################

# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(future)
library(tidyverse)
library(harmony)

plan("multiprocess", workers = 8)
options(future.globals.maxSize = 17000 * 1024^2)
options(future.seed = TRUE)

set.seed(7412)

date = "6.27"

###############################
# Load seurat objects
###############################
tmp_dir <- paste0("../HATIM_Analysis/", date, "/SubsetAnalysis")

Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Bcells <- readRDS(paste0(tmp_dir, "/", "Bcells.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))

# 2. MODIFY OBJECTS_____________________________________________________

###############################
# RNA Assay
###############################
DefaultAssay(Vascular) <- "RNA"
DefaultAssay(Stromal) <- "RNA"
DefaultAssay(Bcells) <- "RNA"
DefaultAssay(Myeloid) <- "RNA"
DefaultAssay(Lymphoid) <- "RNA"

###############################
# Create Celltype Label
###############################
Vascular[['CellType']] <- ifelse(Vascular$Annotation %in% c("Venous EC", "Capillary EC", "Arterial EC"), "Endothelial", "Vascular Smooth Muscle")
Stromal[['CellType']] <- "Stromal"
Bcells[['CellType']] <- "B Cells"
Myeloid[['CellType']] <- "Myeloid"
Lymphoid[['CellType']] <- "T/NK Cells"

###############################
# Create Major Category
###############################
Vascular[['Major_CellType']] <- "Vascular"
Stromal[['Major_CellType']] <- "Stromal"
Bcells[['Major_CellType']] <- "Lymphoid"
Myeloid[['Major_CellType']] <- "Myeloid"
Lymphoid[['Major_CellType']] <- "Lymphoid"

# 3. INTEGRATE OBJECTS_____________________________________________________

###############################
# Merge Datasets
###############################
integrated.data <- merge(Stromal, y = c(Vascular, Bcells, Myeloid, Lymphoid))

# Clear Space
rm(Stromal)
rm(Vascular)
rm(Bcells)
rm(Myeloid)
rm(Lymphoid)

###############################
# Run variable features, normalization, and scale
###############################
integrated.data <- integrated.data %>% FindVariableFeatures(verbose = T, nfeatures = 5000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()
        
###############################
# Run Harmony
###############################
integrated.data <- RunHarmony(integrated.data, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
integrated.data <- integrated.data %>% RunUMAP(reduction = 'harmony', dims = 1:40) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:40) %>% 
  FindClusters(resolution = 2.0)

###############################
# QC Check
###############################
VlnPlot(integrated.data, features = c("RGS5", "VWF"), pt.size = 0)

# Cluster 35 is a endothelial/vascular stromal doublet. It does not appear in the subset version but regardless, will remove here since that is not informative.
integrated.data <- subset(integrated.data, idents = c("31"), invert = T)

# 4. SAVE RDS_____________________________________________________
tmp_dir <- paste0("../HATIM_Analysis/", date, "/Final_Integration")
dir.create(tmp_dir)

saveRDS(integrated.data, file = paste0(tmp_dir, "/", "Integrated.rds"))
