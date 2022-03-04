##############################################################################################
##---------------- HATIM SINGLE CELL SUBCUTANEOUS ADIPOSE TISSUE ATLAS----------------------##
##------------------------- DATE: 10/29/2021 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: The following code will merge all processed subsets (Myeloid, Lymphoid (including 
# B cells), Stromal, and Vascular) and integrate using the harmony algorithm. This will serve
# as Figure 1 in the manuscript.
##############################################################################################

# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(harmony)
library(future)

plan("multiprocess", workers = 8)
options(future.globals.maxSize = 17000 * 1024^2)
options(future.seed = TRUE)

###############################
# Set Seed
##############################
# Set Seed for reproducibility
set.seed(7612) # Reproducibility

###############################
# Load seurat objects
###############################
Endo <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Endo/Endo_RNA_1.3.rds')
Stromal <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_1.3.rds')
Bcells <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Bcells_RNA_11.1.rds')
Myeloid <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_RNA_11.1.rds')
Tcells <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_10.30.rds')


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
Stromal[['CellType']] <- "Stromal"
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

 # Merged
MergedT_dup <- rbind(Myeloid_dup, Endo_dup, stromal_dup, bcells_dup, tcells_dup)
MergedM_dup <- rbind(tcells_dup, Endo_dup, stromal_dup, bcells_dup, Myeloid_dup)

# Find Duplicates in T cells (21 cells)
dup_cells <- duplicated(MergedT_dup$Cell)
Tcells_dup <- MergedT_dup[dup_cells,]

# Find Duplicates in Myeloid (21 cells)
dup_cells <- duplicated(MergedM_dup$Cell)
Mcells_dup <- MergedM_dup[dup_cells,]

###############################
# Remove Duplicated Cells (Affected Cycling T and Myeloid Cells): 21 cells
###############################
Tcells <- subset(Tcells, cells = Tcells_dup$Cell, invert = T)
Myeloid <- subset(Myeloid, cells = Mcells_dup$Cell, invert = T)

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
# Run variable features and scale (data already normalized)
###############################
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
  FindClusters(resolution = 2.0)
  
# 3. SAVE OBJECT_____________________________________________________

saveRDS(integrated.data, file = '/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Final_Integrated.1.3.rds')

