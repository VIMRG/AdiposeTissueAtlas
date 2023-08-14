##############################################################################################
##---------------- HATIM SINGLE CELL SUBCUTANEOUS ADIPOSE TISSUE ATLAS----------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: The following code will process and analyze 10x data for the HATIM cohort with
## the goal to characterize the cell populations in persons with HIV and in HIV- diabetic persons. 
## Runs were multiplexed in groups of 4 and 91 individuals were sequenced in a total of 23   
## runs. The HIV-negative individuals will be processed here at the same time but then
## integrated later. CITE-seq antibodies were used to Hashtag individuals to link the samples with the 
## metadata. SouporCell was used to genetically demultiplex the cells and link them back to 
## the hashtag. DoubletFinder along with visualizating were used to identify doublets that 
## were subsequently removed from the dataset. Lanes were integrated using Harmony Algorithm
##############################################################################################

# 1. SETTING UP ENVIRONMENT------------------------------------------------------------------

###############################
# set seed
###############################
set.seed(7612) # Reproducibility

###############################
# set directory
###############################
Analysis_dir <- "../Atlas_AT"
Cellranger_dir <- "../CellRanger_count"
date <- "6.27"

###############################
# Load Libraries
###############################
library(Seurat) # Seurat 4.0 version(note that several packages will also be attached)
library(tidyverse) # collection of packages for manipulating data (ggplot2, dplyr, tidyr, purr, stringr, tibble)
library(Matrix) # dense and sparse matrix manipulations
library(SeuratDisk) # Stores objects for reference assignments
library(DoubletFinder) # find heterotypic doublets
library(SoupX) # Decontaminate count matrix
library(future) # Parallelize certain functions

plan("multiprocess", workers = 5)
options(future.globals.maxSize = 3000 * 1024^2)
options(future.seed = TRUE)

# 2. SOUPX BACKGROUND DECONTAMINATION------------------------------------------------------------------
# Most runs based on Souporcell had low ambient RNA (< 6%) but a few had ambient RNA between 10-20%,
# increasing noise. Therefore, SoupX was run on the 10x count matrices and the corrected matrices will
# be used for further downstream processing and analysis. Run time is ~ 3 hours but can easily parallelize.
# Note, I already ran this previously so I will skip re-running for this version but is here for reproducibility

project <- c('P5344_CW1', 'P5344_CW2', 'P5344_CW3', 'P5544_CW1', 'P5544_CW2', 'P5544_CW3', 'P5573_CW1', 'P5573_CW2', 'P5573_CW3', 'P5657_CW1', 'P5657_CW3',  'P5836_CW1', 'P5836_CW2', 
'P5836_CW3', 'P5877_CW1', 'P5877_CW2', 'P5877_CW3', 'P5903_CW1', 'P5903_CW2', 'P5903_CW3', 'P5963_CW1', 'P5963_CW2', 'P5963_CW3')

SoupX_Fun <- function(project) {
    print(paste("Reading in 10x Matrix from Project located at ", Cellranger_dir, "/", project, "/outs", ".", sep = ""))
    sc = load10X(paste(Cellranger_dir, project, "outs", sep = "/")) # Will use Cellranger supplied clustering
    print("Calculating rho value of each channel.")
    sc = autoEstCont(sc) # Auto estimates reasonable and rho in line with Souporcell estimates of ambient RNA
    print("Adjusting Counts")
    out = adjustCounts(sc) # Generate corrected count matrix
    print("Writing corrected matrix")
    tmp_path <- paste("../Analysis_", date, "/SoupX/", project, sep = "")
    dir.create(file.path(tmp_path), recursive = T) # Create directory to store corrected matrix
    DropletUtils:::write10xCounts(tmp_path, out, overwrite = TRUE) # Use DropletUtils write10xCounts function to write matrix, feature, and barcode to directory
    print(paste("Corrected Matrix has successfully saved to ", tmp_path, ".", sep = ""))
}

lapply(project, SoupX_Fun)

# 3. LOAD 10X DATA------------------------------------------------------------------
# The corrected count matrices generated from SoupX will be used for all downstream analyses
# from this point on. ADT matrices will be separately loaded

###############################
# Load Corrected Matrices
###############################
soupX_path <- paste("../Analysis_", "1.22", "/SoupX", sep = "")

Read_sparse <- function(project, soupX_path) {
    tmp_dir <- paste(soupX_path, project, sep = "/")
    print(paste("Reading in 10x Matrix from Project located at ", tmp_dir, ".", sep = ""))
    matrix <- Read10X(data.dir = tmp_dir)
    return(matrix)
}

corrected.matrix_list <- lapply(project, FUN = Read_sparse, soupX_path = soupX_path)
names(corrected.matrix_list) <- project

###############################
# Load Original ADT Matrices
###############################
ADT_list <- vector("list", length = length(project))

ADT_fun <- function(project) {
    tmp_dir <- paste0(Cellranger_dir, "/", project, "/outs/filtered_feature_bc_matrix")
    print(paste("Reading in 10x Matrix from Project ", project, " located at ", tmp_dir, ".", sep = ""))
    matrix <- Read10X(data.dir = tmp_dir)
    adt_matrix <- matrix[['Antibody Capture']]
    return(adt_matrix)
}

ADT_list <- lapply(project, ADT_fun)
names(ADT_list) <- project

# 4. PROCESS THE CITE-SEQ ANTIBODIES------------------------------------------------------------------
# Antibody names need to be the same between runs. A few names differ in the earlier runs compared to later
# runs so will clean this and use this data to create a seurat object.

###############################
# Remove extra hashtag on first 4 runs
###############################
for (i in 1:4) {
    ADT_list[[i]] <- ADT_list[[i]][-5,]
}

###############################
# Make all ADT same name
###############################
adt_names <- rownames(ADT_list[[23]])
adt_names <- stringr::str_remove(adt_names, "TotalSeq_C_")

ADT_list <- lapply(ADT_list, function(xx) {
    rownames(xx) <- adt_names
    return(xx)
})

# 5. GENERATE SEURAT OBJECTS------------------------------------------------------------------
# Create seurat objects with a gene expression, HTO, and ADT assay

###############################
# Load Corrected Count Matrices
###############################
seurat_list <- vector("list", length = length(project))
names(seurat_list) <- project

seurat_list <- lapply(corrected.matrix_list, FUN = function(xx) {
  print(paste("Creating Gene Expression Seurat Objects for ", names(xx), ".", sep = ""))
  CreateSeuratObject(counts = xx, min.cells = 3)
})

###############################
# Load ADT and HTO data
###############################
# Add HTO and ADT assays to all seurat objects

for (i in 1:length(ADT_list)) {
  print(paste0("Creating HTO and ADT assay for project ", names(seurat_list[i]), "from ADT matrix ", names(ADT_list[i]), "."))
  seurat_list[[i]][['HTO']] <- CreateAssayObject(ADT_list[[i]][1:4,][, colnames(seurat_list[[i]])])
  seurat_list[[i]][['ADT']] <- CreateAssayObject(ADT_list[[i]][5:49,][, colnames(seurat_list[[i]])])
}

# 6.  DEMULTIPLEX SAMPLES BY HASHTAG--------------------------------------------------------
###############################
# Change to HTO assay
###############################
for (i in 1:length(seurat_list)) {
  DefaultAssay(seurat_list[[i]]) <- "HTO"
}

###############################
# Normalize and use Demux
###############################
seurat_list <- lapply(seurat_list, FUN = function(xx) {
  xx <- NormalizeData(xx, assay = "HTO", normalization.method = "CLR")
  xx <- HTODemux(xx, assay = "HTO", positive.quantile = 0.99)
})

# 7. ADD METADATA -------------------------------------

###############################
# Change to RNA Assay
###############################
for (i in 1:length(seurat_list)) {
  DefaultAssay(seurat_list[[i]]) <- "RNA"
}

###############################
# Add Percent Mitochondria
###############################
seurat_list <- lapply(seurat_list, FUN = function(xx) {
  xx <- PercentageFeatureSet(xx, pattern = '^MT-', col.name = 'percent.mt')
})

###############################
# Add SouporCell Data
###############################
soup_path <- paste("../Analysis_", "1.22", "/Soup", sep = "")

Soup_fun <- function(project) {
    tmp_dir <- paste(soup_path, project, sep = "/")
    print(paste("Loading Genetic Demultiplexing Cluster Assignment for ", project, " located at ", tmp_dir, ".", sep = ""))
    cluster <- read.table(paste(tmp_dir, "clusters.tsv", sep = "/"), header = T) # read in tsv file
    rownames(cluster) <- cluster$barcode # Need to set rownames to barcode to merge with seurat object later on
    return(cluster)
}

Soup_list <- lapply(project, FUN = Soup_fun)
names(Soup_list) <- project

for (i in 1:length(seurat_list)) {
  print(paste0("Adding Genetic Clustering Data to ", names(seurat_list[i]), " from the Souporcell cluster ", names(Soup_list[i])))
  seurat_list[[i]] <- AddMetaData(seurat_list[[i]], Soup_list[[i]][,1:3])
}

# 8. SAVE Unprocessed SEURAT OBJECTS PRIOR TO APPLYING QC-------------------------------------
tmp_dir <- paste("../HATIM_Analysis", date, "unprocessed", sep = "/")
dir.create(tmp_dir)

for (i in 1:length(seurat_list)) {
  print(paste("Saving Seurat object ", names(seurat_list[i]), ".", sep = ""))
  name <- names(seurat_list[i])
  seurat_object <- seurat_list[[i]]
  saveRDS(seurat_object, file = paste(tmp_dir, "/", name, ".rds", sep = ""))
}

# Read in seurat files
seurat_path = paste0("../HATIM_Analysis/", date, "/unprocessed")
seurat_list <- list.files(path = seurat_path, full.names = T)
names(seurat_list) <- project
seurat_list <- lapply(seurat_list, function(xx) {
    print(paste0("Reading saved seurat files for project ", "from ", xx, ".")) 
    readRDS(xx)
})

# 9. APPLY QC FILTER TO SEURAT OBJECTS-------------------------------------
# For QC, will remove cells with >= 25% mitochondrial genes (stressed/dying), =< 200 features, or
# gene count of <= 800. 800 was chosen as a conservative cutoff to ensure enough informative SNVs
# to classify cell genotypes. Unfortunately, this may miss granulocytic cells, but they likely do
# not exist in our data due to the tissue processing.


###############################
# Apply QC Parameters
###############################
# Apply QC (mt.percent < 25%, UMI > 800, feature > 200)
seurat_list <- lapply(X = seurat_list, FUN = function(xx) {
    subset(xx, subset = 
                (percent.mt < 25) &
                (nFeature_RNA > 200) &
                (nCount_RNA > 800))
})

# 10. APPLY PRELIMINARY NORMALIZATION, SCALING, AND CLUSTERING OF SEURAT OBJECTS-------------------------------------
# Normalization, scaling, and dimensional reduction/clustering is performed on each
# lane separately to examine quality. Holding off on removing genetic doublets as I
# will use these for ground truths in DoubletFinder. Running on SCTransform.

###-------------------------------####
# Change assay to RNA
###------------------------------####
for (i in 1:length(seurat_list)) {
  DefaultAssay(seurat_list[[i]]) <- "RNA"
}

###############################
# Normalize and Scale
###############################
#Will regress out cycling genes here for purposes of doublet finder

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_list <- lapply(X = seurat_list, FUN = function(xx) {
  xx <- NormalizeData(xx, assay = "RNA")
  xx <- NormalizeData(xx, assay = "ADT", normalization.method = "CLR", margin = 2)
  xx <- CellCycleScoring(xx, g2m.features = g2m.genes, s.features = s.genes)
  xx <- SCTransform(xx, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
})


###############################
# Apply PCA and DownStream Clustering
###############################
seurat_list <- lapply(X = seurat_list, FUN = function(xx) {
  xx <- RunPCA(xx)
  xx <- RunUMAP(xx, dims = 1:15)
  xx <- FindNeighbors(xx, dims = 1:15)
  xx <- FindClusters(xx, resolution = 2.0)
})

# 11. ADD METADATA TO SEURAT OBJECTS-------------------------------------
###############################
# Link Soup Cluster to Hashtag
###############################
Soup_Hash_fun <- function(seurat_object) {
    print(seurat_object)
    df <- data.frame(HTO_Global = seurat_object$HTO_classification.global, HTO_classification = seurat_object$HTO_classification, 
    Soup = seurat_object$assignment, Status = seurat_object$status) # Pull hashtag and soup cluster info from metadata
    df <- df%>% filter(HTO_Global == "Singlet" & Status == "singlet") # Filter so only singlet in both soup and hashtag (informative cells)
    df <- df %>% group_by(Soup, HTO_classification) %>% count(HTO_classification) %>% group_by(Soup) %>% mutate(Freq = n/sum(n)) %>% group_by(Soup) %>% slice_max(Freq) # Count cells in 251, 252, 253, 254 for each cluster and find frequency
    print(df)
    print
    seurat_object[["Soup_Hash"]] <- df$HTO_classification[match(seurat_object$assignment, df$Soup)]
    print(table(seurat_object$Soup_Hash))
    print(rowSums(table(seurat_object$assignment, seurat_object$HTO_classification)))
    return(seurat_object)
}

seurat_list <- lapply(seurat_list, Soup_Hash_fun)

# Manually fix a couple of lanes (20) P903_CW3 (error because 254 didn't work well at all) , P5903_CW3
seurat_list[[20]]$Soup_Hash <- ifelse(seurat_list[[20]]$assignment == 2, 254, seurat_list[[20]]$Soup_Hash) # Cite-seq didn't work that well in this run. 252 diffusely positive so by process of elimination,

# Of note, lane P5836_CW1 only had 3 genotypic clusters due to small sample size (technically 4 but only 2 cells in the genotypic cluster). Thus, only 3 individuals were meaningfully tagged back to hashtag.

###############################
# Save
###############################
tmp_dir <- paste("../HATIM_Analysis", date, "processed", sep = "/")
dir.create(tmp_dir)

for (i in 1:length(seurat_list)) {
  print(paste("Saving Seurat object ", names(seurat_list[i]), ".", sep = ""))
  name <- names(seurat_list[i])
  seurat_object <- seurat_list[[i]]
  saveRDS(seurat_object, file = paste(tmp_dir, "/", name, ".rds", sep = ""))
}
