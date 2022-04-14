##############################################################################################
##------------------------- HARMONY INTEGRATION ---------------------------------##
##------------------------- DATE: 8/18/2021 AUTHOR: -------------------------------##
## DESCRIPTION: The following code will integrate the final seurat objects
##############################################################################################


# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(harmony)
library(future)

plan("multiprocess", workers = 8)
options(future.globals.maxSize = 15000 * 1024^2)
options(future.seed = TRUE)

###############################
# Set Seed
##############################
# Set Seed for reproducibility
set.seed(7612) # Reproducibility

###############################
# Load seurat objects
###############################
Endo <- readRDS("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo_RNA.rds")
Stromal <- readRDS("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal_RNA.rds")
Bcells <- readRDS("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Bcells_RNA.rds")
Myeloid <- readRDS("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid_RNA.rds")
Tcells <- readRDS("/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_Final.rds")

# 2. MODIFY OBJECTS_____________________________________________________

###############################
# RNA Assay
###############################
DefaultAssay(Endo) <- "RNA"
DefaultAssay(Stromal) <- "RNA"
DefaultAssay(Bcells) <- "RNA"
DefaultAssay(Myeloid) <- "RNA"
DefaultAssay(Tcells) <- "RNA"

###############################
# Create Celltype Label
###############################
Endo[['CellType']] <- "Endothelium"
Stromal[['CellType']] <- "Progenitors"
Bcells[['CellType']] <- "B Cells"
Myeloid[['CellType']] <- "Myeloid"
Tcells[['CellType']] <- "T Cells"

###############################
# Create Major Category
###############################
Endo[['Major_CellType']] <- "Endothelium"
Stromal[['Major_CellType']] <- "Stromal"
Bcells[['Major_CellType']] <- "Lymphoid"
Myeloid[['Major_CellType']] <- "Myeloid"
Tcells[['Major_CellType']] <- "Lymphoid"

###############################
# Check Duplicated Cells
###############################
Myeloid_dup <- data.frame(Myeloid[['CellType']])
Endo_dup <- data.frame(Endo[['CellType']])
stromal_dup <- data.frame(Stromal[['CellType']])
bcells_dup <- data.frame(Bcells[['CellType']])
tcells_dup <- data.frame(Tcells[['CellType']])

Myeloid_dup$Cell <- rownames(Myeloid_dup)
Endo_dup$Cell <- rownames(Endo_dup)
stromal_dup$Cell <- rownames(stromal_dup)
bcells_dup$Cell <- rownames(bcells_dup)
tcells_dup$Cell <- rownames(tcells_dup)

 # Merge
 Merged_dup <- rbind(Myeloid_dup, Endo_dup, stromal_dup, bcells_dup, tcells_dup)

# Find Duplicates
dup_cells <- duplicated(Merged_dup$Cell)
cells_dup <- Merged_dup[dup_cells,]

###############################
# Remove Duplicated B Cells
###############################
bcell_duplicates <- cells_dup[cells_dup == "B Cells",]
#bcell_singlets <- colnames(Bcells)[colnames(Bcells) %notin% bcell_duplicates$Cell]
Bcells <- subset(Bcells, cells = bcell_duplicates$Cell, invert = T)

# 2. INTEGRATE OBJECTS_____________________________________________________

###############################
# Merge Datasets
###############################
integrated.data <- merge(Stromal, y = c(Endo, Bcells, Myeloid, Tcells))

# Clear Space
rm(Stromal)
rm(Endo)
rm(Bcells)
rm(Myeloid)
rm(Tcells)

###############################
# Run variable features, normalization, and scale
###############################
integrated.data <- NormalizeData(integrated.data, assay = "RNA")
integrated.data <- NormalizeData(integrated.data, assay = "ADT", normalization.method = "CLR", margin = 2)

integrated.data <- integrated.data %>% FindVariableFeatures(verbose = T, nfeatures = 5000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()
        
###############################
# Run Harmony
###############################
integrated.data <- RunHarmony(integrated.data, "Project", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
integrated.data <- integrated.data %>% RunUMAP(reduction = 'harmony', dims = 1:40) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:40) %>% 
  FindClusters(resolution = 1.0)
  
# 3. SAVE OBJECT_____________________________________________________

saveRDS(integrated.data, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Final_Integrated_8.29.rds')

