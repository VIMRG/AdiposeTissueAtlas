##############################################################################################
##---------------- HATIM SINGLE CELL SUBCUTANEOUS ADIPOSE TISSUE ATLAS----------------------##
##------------------------- DATE: 8/18/2021 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: The following code will process and analyze 10x data for the HATIM cohort with
## the goal to characterize the cell populations in persons with HIV and in HIV- diabetic persons. 
## Runs were multiplexed in groups of 4 and 96 individuals were sequenced in a total of 24   
## runs. CITE-seq antibodies were used to Hashtag individuals to link the samples with the 
## metadata. SouporCell was used to genetically demultiplex the cells and link them back to 
## the hashtag. DoubletFinder along with visualizating were used to identify doublets that 
## were subsequently removed from the dataset. Lanes were integrated using Seurat integration 
## function with reciprocal PCA and reference-based integration. Cells were manually identifed 
## using cannonical genes. Later, I will write code to use Harmony algorithm and compare results.
##############################################################################################

# 1. SETTING UP ENVIRONMENT------------------------------------------------------------------

###-------------------####
# Set seed
###-------------------####
set.seed(7612) # Reproducibility
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Functions.R') # Helper functions used throughout script

###-------------------####
# Set working directory
###-------------------####
setwd('/data/p_koethe_lab/Atlas_AT/Analysis')

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
library(SoupX) # Decontaminate count matrix

# 2. SOUPX BACKGROUND DECONTAMINATION------------------------------------------------------------------
# Most runs based on Souporcell had low ambient RNA (< 6%) but a few had ambient RNA between 10-20%,
# increasing noise. Therefore, SoupX was run on the 10x count matrices and the corrected matrices will
# be used for further downstream processing and analysis.

###-------------------####
# Get 10x Directories in a list
###-------------------####
dataset_local <- '/data/p_koethe_lab/Atlas_AT/CellRanger/cellranger_count'
inputfile_local <- 'outs/'

# List of Cellranger output directories (not stored in the same place)
list <- c('sampleP5344', 'P5344/P5344_CW_2', 'P5344/P5344_CW3', 'P5544/P5544_CW_1', 'P5544/P5544_CW_2', 'P5544/P5544_CW_3',
          'P5573/P5573_CW_1', 'P5573/P5573_CW_2', 'P5573/P5573_CW_3', 'P5657/P5657_CW_1', 'P5657/P5657_CW_3',
          'P5836/P5836_CW_1', 'P5836/P5836_CW_2', 'P5836/P5836_CW_3', 'P5877/P5877_CW_1', 'P5877/P5877_CW_2', 'P5877/P5877_CW_3',
          'P5903/P5903_CW_1', 'P5903/P5903_CW_2', 'P5903/P5903_CW_3', 'P5963/P5963_CW_1', 'P5963/P5963_CW_2', 'P5963/P5963_CW_3')

names(list) <- c('P5344_CW1', 'P5344_CW2', 'P5344_CW3', 'P5544_CW1', 'P5544_CW2', 'P5544_CW3',
                 'P5573_CW1', 'P5573_CW2', 'P5573_CW3', 'P5657_CW1', 'P5657_CW3', 
                 'P5836_CW1', 'P5836_CW2', 'P5836_CW3', 'P5877_CW1', 
                 'P5877_CW2', 'P5877_CW3', 'P5903_CW1', 'P5903_CW2', 
                 'P5903_CW3', 'P5963_CW1', 'P5963_CW2', 'P5963_CW3')

###-------------------####
# Call SoupX in a function to
# Write new count corrected matrices
# Wall Time approximation: 3 hrs (can easily parallelize to reduce time)
###-------------------####
lapply(list, FUN = SoupX_Fun)

rm(dataset_local)
rm(inputfile_local)

# 3. LOAD 10X DATA------------------------------------------------------------------
# The corrected count matrices generated from SoupX will be used for all downstream analyses
# from this point on.

###-------------------####
# Get SoupX matrix, barcodes, and features directory,
# will need to separately load antibody matrices
###-------------------####
matrix_location <- '/data/p_koethe_lab/Atlas_AT/SoupX'

list <- c('sampleP5344', 'P5344/P5344_CW_2', 'P5344/P5344_CW3', 'P5544/P5544_CW_1', 'P5544/P5544_CW_2', 'P5544/P5544_CW_3',
          'P5573/P5573_CW_1', 'P5573/P5573_CW_2', 'P5573/P5573_CW_3', 'P5657/P5657_CW_1', 'P5657/P5657_CW_3',
          'P5836/P5836_CW_1', 'P5836/P5836_CW_2', 'P5836/P5836_CW_3', 'P5877/P5877_CW_1', 'P5877/P5877_CW_2', 'P5877/P5877_CW_3',
          'P5903/P5903_CW_1', 'P5903/P5903_CW_2', 'P5903/P5903_CW_3', 'P5963/P5963_CW_1', 'P5963/P5963_CW_2', 'P5963/P5963_CW_3')

# Create names for Seurat sparse matrices
names(list) <- c('P5344_CW1', 'P5344_CW2', 'P5344_CW3', 'P5544_CW1', 'P5544_CW2', 'P5544_CW3',
                 'P5573_CW1', 'P5573_CW2', 'P5573_CW3', 'P5657_CW1', 'P5657_CW3', 
                 'P5836_CW1', 'P5836_CW2', 'P5836_CW3', 'P5877_CW1', 
                 'P5877_CW2', 'P5877_CW3', 'P5903_CW1', 'P5903_CW2', 
                 'P5903_CW3', 'P5963_CW1', 'P5963_CW2', 'P5963_CW3')

###-------------------####
# Call Read10X in function to generate 
# list of sparse matrices
###-------------------####
corrected.matrix_list <- lapply(list, FUN = Read_sparse, matrix_location = matrix_location)

rm(matrix_location)

###-------------------####
# Call Read10X in function to generate 
# list of sparse ADT matrices
###-------------------####
dataset_local <- '/data/p_koethe_lab/Atlas_AT/CellRanger/cellranger_count'
inputfile_local <- 'outs/filtered_feature_bc_matrix'

# List of Cellranger output directories (not stored in the same place)
list_ADT <- c('sampleP5344', 'P5344/P5344_CW_2', 'P5344/P5344_CW3', 'P5544/P5544_CW_1', 'P5544/P5544_CW_2', 'P5544/P5544_CW_3',
              'P5573/P5573_CW_1', 'P5573/P5573_CW_2', 'P5573/P5573_CW_3', 'P5657/P5657_CW_1', 'P5657/P5657_CW_3',
              'P5836/P5836_CW_1', 'P5836/P5836_CW_2', 'P5836/P5836_CW_3', 'P5877/P5877_CW_1', 'P5877/P5877_CW_2', 'P5877/P5877_CW_3',
              'P5903/P5903_CW_1', 'P5903/P5903_CW_2', 'P5903/P5903_CW_3', 'P5963/P5963_CW_1', 'P5963/P5963_CW_2', 'P5963/P5963_CW_3')

names(list_ADT) <- c('P5344_CW1.data', 'P5344_CW2.data', 'P5344_CW3.data', 'P5544_CW1.data', 'P5544_CW2.data', 'P5544_CW3.data',
                     'P5573_CW1.data', 'P5573_CW2.data', 'P5573_CW3.data', 'P5657_CW1.data', 'P5657_CW3.data', 
                     'P5836_CW1.data', 'P5836_CW2.data', 'P5836_CW3.data', 'P5877_CW1.data', 
                     'P5877_CW2.data', 'P5877_CW3.data', 'P5903_CW1.data', 'P5903_CW2.data', 
                     'P5903_CW3.data', 'P5963_CW1.data', 'P5963_CW2.data', 'P5963_CW3.data')

# Loop through to read 10x matrices
for (i in 1:length(list)){
  print(paste("Reading in 10x Matrix from Project ", names(list_ADT[i]), " located at ", dataset_local, "/", list_ADT[[i]], ".", sep = ""))
  seurat_data <- Read10X(data.dir = paste(dataset_local, list_ADT[[i]], inputfile_local, sep = "/"))
  assign(names(list_ADT[i]), seurat_data)
}


# 4. LOAD SOUPORCELL DATA------------------------------------------------------------------
# SouporCell is a genetic demultiplexing program that uses informative SNV to cluster cells
# into specific genotypic clusters belonging to individuals. This was run on the 10x CellRanger
# filtered matrices.

###-------------------####
# Load SouporCell
###-------------------####
# Path to SouporCell Results
dataset_local <- "/data/p_koethe_lab/Atlas_AT/SouporCell/Output"

# List directory of SouporCell
list_Soup <- c("SouporCell_P5344/CW1", "SouporCell_P5344/CW2", "SouporCell_P5344/CW3", "SouporCell_P5544/CW1", 
               "SouporCell_P5544/CW2", "SouporCell_P5544/CW3", "SouporCell_P5573/CW1", "SouporCell_P5573/CW2", "SouporCell_P5573/CW3",
               "SouporCell_P5657/CW1", "SouporCell_P5657/CW3", 
               "SouporCell_P5836/CW1", "SouporCell_P5836/CW2", "SouporCell_P5836/CW3", "SouporCell_P5877/CW1", "SouporCell_P5877/CW2", 
               "SouporCell_P5877/CW3", "SouporCell_P5903/CW1",
               "SouporCell_P5903/CW2", "SouporCell_P5903/CW3", "SouporCell_P5963/CW1", "SouporCell_P5963/CW2", "SouporCell_P5963/CW3")

names(list_Soup) <- c("Soup_P5344_CW1", "Soup_P5344_CW2", "Soup_P5344_CW3", "Soup_P5544_CW1", "Soup_P5544_CW2", "Soup_P5544_CW3",
                      "Soup_P5573_CW1", "Soup_P5573_CW2", "Soup_P5573_CW3", "Soup_P5657_CW1", "Soup_P5657_CW3",
                      "Soup_P5836_CW1", "Soup_P5836_CW2", "Soup_P5836_CW3", "Soup_P5877_CW1", "Soup_P5877_CW2", "Soup_P5877_CW3",
                      "Soup_P5903_CW1", "Soup_P5903_CW2", "Soup_P5903_CW3", "Soup_P5963_CW1", "Soup_P5963_CW2", "Soup_P5963_CW3")

###-------------------####
# Call function to load SouporCell
# data
###-------------------####
list_Soup <- lapply(list_Soup, FUN = Soup_tables, dataset_local = dataset_local)

rm(dataset_local)

# 5. PROCESS THE CITE-SEQ ANTIBODIES------------------------------------------------------------------
# Antibody names need to be the same between runs. A few names differ in the earlier runs compared to later
# runs so will fix this and use this data to create a seurat object.

###-------------------####
# Fix run P5344_CW1
###-------------------####
# rename CITE-SEQ Antibodies
rownames(P5344_CW1.data[['Antibody Capture']]) <- gsub("TotalSeq_C", "", rownames(P5344_CW1.data$`Antibody Capture`))
rownames(P5344_CW1.data[['Antibody Capture']]) <- mgsub(rownames(P5344_CW1.data[['Antibody Capture']]), c("Antibody", "anti_human_", "anti_human", "(^[0-9]+)"), c("", "", "", ""))

# rename Hashtag antibodies
rownames(P5344_CW1.data[['Antibody Capture']])[c(1:5)] <- c("251", "252", "253", "254", "255") # change naming of hashtag to conform with others

# rename antibodies that are not the same as the others
rownames(P5344_CW1.data[['Antibody Capture']])[c(8,16,19,20,23,25,27,30,31,33,34,35,36,41,42,44,45,47, 48, 50)] <- c("CD56_recom", "ICOS", "OX40", "HLA-DR", "CD40L", "BTLA", "CCR7", "Fas", "PD-1", 
                                                                                                                     "LAG-3", "CTLA-4", "TIM-3", "CD1C", "TCRYd", "TCRVa24-Ja18", "TCRVa7.2",
                                                                                                                     "CCR6", "CCR4", "CD11C", "CXCR3") # Change names to match others
# Remove Hashtag 255 since that was not used
P5344_CW1.data[['Antibody Capture']] <- P5344_CW1.data[['Antibody Capture']][-5,] # Remove hasthag 255 as this was not used

###-------------------####
# Fix run P5344_CW2
###-------------------####
# Remove extra labels
rownames(P5344_CW2.data[['Antibody Capture']]) <- mgsub(rownames(P5344_CW2.data$`Antibody Capture`), c("TotalSeq_C_", "TotalSeq_C"), c("", ""))

# Change D11C to CD11C (error in typing) then remove hashtag 255
rownames(P5344_CW2.data[['Antibody Capture']])[48] <- "CD11C" # Obnoxious error of misspelling
P5344_CW2.data[['Antibody Capture']] <- P5344_CW2.data[['Antibody Capture']][-5,] # Remove extra hashtag

###-------------------####
# Fix run P5344_CW3
###-------------------####
# Remove extra labels and hashtag 255
rownames(P5344_CW3.data[['Antibody Capture']]) <- mgsub(rownames(P5344_CW3.data$`Antibody Capture`), c("TotalSeq_C_"), c(""))
P5344_CW3.data[['Antibody Capture']] <- P5344_CW3.data[['Antibody Capture']][-5,] # Remove extra hashtag

###-------------------####
# Fix run P5544_CW1
###-------------------####
# Remove hashtag 255 and extra Antibody labeles
rownames(P5544_CW1.data[['Antibody Capture']]) <- mgsub(rownames(P5544_CW1.data$`Antibody Capture`), c("TotalSeq_C_"), c(""))
P5544_CW1.data[['Antibody Capture']] <- P5544_CW1.data[['Antibody Capture']][-5,] # Remove extra hashtag

###-------------------####
# All other runs
###-------------------####
# Generate list of sparse CITE-Seq counts
list_ab <- c(P5544_CW2.data$`Antibody Capture`, 
             P5544_CW3.data$`Antibody Capture`, P5573_CW1.data$`Antibody Capture`,
             P5573_CW2.data$`Antibody Capture`, P5573_CW3.data$`Antibody Capture`, P5657_CW1.data$`Antibody Capture`,
             P5657_CW3.data$`Antibody Capture`, P5836_CW1.data$`Antibody Capture`, 
             P5836_CW2.data$`Antibody Capture`, P5836_CW3.data$`Antibody Capture`, P5877_CW1.data$`Antibody Capture`,
             P5877_CW2.data$`Antibody Capture`, P5877_CW3.data$`Antibody Capture`, P5903_CW1.data$`Antibody Capture`, 
             P5903_CW2.data$`Antibody Capture`, P5903_CW3.data$`Antibody Capture`, P5963_CW1.data$`Antibody Capture`,
             P5963_CW2.data$`Antibody Capture`, P5963_CW3.data$`Antibody Capture`)


# Loop through to change antibody label (remove 'total_seq' e.g.)
list_ab <- lapply(list_ab, FUN = Citeseq_name)

###-------------------####
# Check that all antibodies are
# names the same
###-------------------####
all.equal(rownames(P5344_CW1.data[['Antibody Capture']]), rownames(P5344_CW2.data[["Antibody Capture"]]))
all.equal(rownames(P5344_CW1.data[['Antibody Capture']]), rownames(P5344_CW3.data[["Antibody Capture"]]))
all.equal(rownames(P5344_CW1.data[['Antibody Capture']]), rownames(P5544_CW1.data[["Antibody Capture"]]))


for (files in list_ab) {
  t <- all.equal(rownames(P5344_CW1.data[['Antibody Capture']]), rownames(files))
  print(t)
}

# 6. GENERATE SEURAT OBJECTS------------------------------------------------------------------
# This section will create seurat objects with a gene expression, HTO, and ADT assay

###-------------------------------####
# Load Gene Counts from Corrected Matrix
###------------------------------####
# Load the corrected matrix from 'corrected.matrix_list' into seurat objects

list_seurat <- lapply(corrected.matrix_list, FUN = function(x) {
  print(paste("Creating Gene Expression Seurat Objects for Projects.", sep = ""))
  x <- CreateSeuratObject(counts = x, min.cells = 3)
})

###-------------------------------####
# Load HTO and ADT from original matrix
# since this data does not need correcting
###------------------------------####
# Merge antibody list with first 4 projects
list_ab_capture <- c(P5344_CW1.data$`Antibody Capture`, P5344_CW2.data$`Antibody Capture`, P5344_CW3.data$`Antibody Capture`, P5544_CW1.data$`Antibody Capture`, list_ab)

names(list_ab_capture) <-  c('P5344_CW1.data', 'P5344_CW2.data', 'P5344_CW3.data', 'P5544_CW1.data', 'P5544_CW2.data', 'P5544_CW3.data',
                             'P5573_CW1.data', 'P5573_CW2.data', 'P5573_CW3.data', 'P5657_CW1.data', 'P5657_CW3.data', 
                             'P5836_CW1.data', 'P5836_CW2.data', 'P5836_CW3.data', 'P5877_CW1.data', 
                             'P5877_CW2.data', 'P5877_CW3.data', 'P5903_CW1.data', 'P5903_CW2.data', 
                             'P5903_CW3.data', 'P5963_CW1.data', 'P5963_CW2.data', 'P5963_CW3.data')

# Add HTO and ADT assays to all seurat objects
for (i in 1:length(list_ab_capture)) {
  print(names(list_seurat[i]))
  print(names(list_ab_capture[i]))
  list_seurat[[i]][['HTO']] <- CreateAssayObject(list_ab_capture[[i]][1:4,][, colnames(list_seurat[[i]])])
  list_seurat[[i]][['ADT']] <- CreateAssayObject(list_ab_capture[[i]][5:49,][, colnames(list_seurat[[i]])])
}


# 7. CHECK THAT OBJECTS ARE CORRECT------------------------------------------------------------------
# Check to see that the above generated the correct seurat objects and match the 
list_seurat

# 8.  DEMULTIPLEX SAMPLES BY HASHTAG--------------------------------------------------------

###-------------------------------####
# Change Default Assay to HOT
###------------------------------####
for (i in 1:length(list_seurat)) {
  DefaultAssay(list_seurat[[i]]) <- "HTO"
}

###-------------------------------####
# Normalize and use demux function
###------------------------------####
list_seurat <- lapply(list_seurat, FUN = function(x) {
  x <- NormalizeData(x, assay = "HTO", normalization.method = "CLR")
  x <- HTODemux(x, assay = "HTO", positive.quantile = 0.99)
})

###-------------------------------####
# Visualize Heatmaps of HTO
###------------------------------####

lapply(list_seurat, FUN = function(x) {
  HTOHeatmap(x, assay = "HTO", ncells = 10000)
})

###-------------------------------####
# Print Global Classification of HTO
###------------------------------####

lapply(list_seurat, FUN = function(x) {
  print(table(x[['HTO_classification.global']]))
})

# 9. ADD METADATA (will add individual metadata later) -------------------------------------

###-------------------------------####
# Change Assay to RNA
###------------------------------####
for (i in 1:length(list_seurat)) {
  DefaultAssay(list_seurat[[i]]) <- "RNA"
}

###-------------------------------####
# Add Percent Mitochondria
###------------------------------####
list_seurat <- lapply(list_seurat, FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = '^MT-', col.name = 'percent.mt')
})

###-------------------------------####
# Add SouporCell
###------------------------------####
# Will use list of soup cluster table we obtained earlier
for (i in 1:length(list_seurat)) {
  print(names(list_seurat[i]))
  print(names(list_Soup[i]))
  list_seurat[[i]] <- AddMetaData(list_seurat[[i]], list_Soup[[i]][,2:3])
}

# 10. SAVE RAW SEURAT OBJECTS PRIOR TO APPLYING QC-------------------------------------
for (i in 1:length(list_seurat)) {
  print(list_seurat[[i]])
  seurat_object <- list_seurat[[i]]
  name <- names(list_seurat[i])
  saveRDS(seurat_object, file = paste("/data/p_koethe_lab/Atlas_AT/Analysis/Raw/", name, ".rds", sep = ""))
}

# 11. APPLY QC FILTER TO SEURAT OBJECTS-------------------------------------
# For QC, will remove cells with >= 25% mitochondrial genes (stressed/dying), =< 200 features, or
# gene count of <= 800. 800 was chosen as a conservative cutoff to ensure enough informative SNVs
# to classify cell genotypes. Unfortunately, this may miss granulocytic cells, but they likely do
# not exist in our data due to the tissue processing.

###-------------------------------####
# Violin Plots
###------------------------------####
lapply(X = list_seurat, FUN = function(x){
  VlnPlot(x, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))
})


###-------------------------------####
# Apply QC Parameters
###------------------------------####
# Apply QC (mt.percent < 25%, UMI > 800, feature > 200)
list_seurat <- lapply(X = list_seurat, FUN = function(x) {
  x <- subset(x, subset = 
                (percent.mt < 25) &
                (nFeature_RNA > 200) &
                (nCount_RNA > 800))
})

###-------------------------------####
# Select only Singlets from SouporCell
###------------------------------####
list_seurat <- lapply(X = list_seurat, FUN = function(x) {
  Idents(x) <- "status"
  x <- subset(x, idents = "singlet")
})

# 12. APPLY PRELIMINARY NORMALIZATION, SCALING, AND CLUSTERING OF SEURAT OBJECTS-------------------------------------
# Normalization, scaling, and dimensional reduction/clustering is performed on each
# lane separately to examine quality.

###-------------------------------####
# Change assay to RNA
###------------------------------####
for (i in 1:length(list_seurat)) {
  DefaultAssay(list_seurat[[i]]) <- "RNA"
}

###-------------------------------####
# Normalize and Scale Data
###------------------------------####
load("/data/p_koethe_lab/Atlas_AT/MetaData/Reference/cycle.rda")

list_seurat <- lapply(X = list_seurat, FUN = function(x) {
  x <- NormalizeData(x, assay = "RNA")
  x <- NormalizeData(x, assay = "ADT", normalization.method = "CLR", margin = 2)
  x <- CellCycleScoring(x, g2m.features = g2m_genes, s.features = s_genes)
  x <- SCTransform(x, vars.to.regress = "percent.mt")
})


###-------------------------------####
# PCA, UMAP, and clustering with 20 PCs
###------------------------------####
list_seurat <- lapply(X = list_seurat, FUN = function(x) {
  x <- RunPCA(x)
  x <- RunUMAP(x, dims = 1:20)
  x <- FindNeighbors(x, dims = 1:20)
  x <- FindClusters(x, resolution = 2.0)
})


# 13. SAVE PROCESSED SEURAT OBJECTS TO FOLDER-------------------------------------
# generate for loop and save output
for (i in 1:length(list_seurat)) {
  print(list_seurat[[i]])
  seurat_object <- list_seurat[[i]]
  name <- names(list_seurat[i])
  saveRDS(seurat_object, file = paste("/data/p_koethe_lab/Atlas_AT/Analysis/Pre_Processed_objects/", name, ".rds", sep = ""))
}

# 13. HETEROTYPIC DOUBLET DETECTION-------------------------------------
# After genetic demultiplexing, the majority of the doublets in the datasets have been removed.
# However, genetic demultiplexing will miss heterotypic doublets. DoubletFinder can identify
# cells that represent heterotypic doublets. Standard settings will be used and will not adjust for 
# homotypic doublets since these have been removed from the dataset already with SouporCell. Using 20 PCs
# is reasonable based on examination of the datasets. Note: this will take a long time to run if not done in 
# parallel. This script will be run in parallel on the cluster but is here for record.

###-------------------------------####
# PCA, UMAP, and clustering with 20 PCs
###------------------------------####
for (i in 1:length(list_seurat)) {
  DefaultAssay(list_seurat[[i]]) <- "SCT"
}

PCs <- 1:20

# generate for loop and save output
for (i in 1:length(list_seurat)) {
  print(list_seurat[[i]])
  sweep.res.list <- paramSweep_v3(list_seurat[[i]], PCs = PCs, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
  nExp_poi <- round(0.05*nrow(list_seurat[[i]]@meta.data))
  list_seurat[[i]] <- doubletFinder_v3(list_seurat[[i]], PCs = PCs, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  seurat_object <- list_seurat[[i]]
  name <- names(list_seurat[i])
  saveRDS(seurat_object, file = paste("/data/p_koethe_lab/Atlas_AT/Analysis/Processed_objects/", name, ".rds", sep = ""))
}


