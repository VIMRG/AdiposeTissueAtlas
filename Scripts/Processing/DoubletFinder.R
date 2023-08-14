##############################################################################################
##------------------------- DOUBLET FINDER PARALLELIZATION ---------------------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: The following code uses DoubletFinder to identify potential heterotypic doublets
# after genetic demultiplexing has removed homotypic doublets. Based on 10x Chromium expected
# doublets accounting for multiplexing of 4 unique individuals, I expect after genetic demultiplexing
# and removing homotypic doublets that ~ 5% of the total number of cells (25% of total doublets) likely
# remains and will use that as the estimate. 
##############################################################################################

###-------------------------------####
# Environment
###------------------------------####
set.seed(7612) # random seed

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID') # Get environmental variable from slurm array
i <- as.numeric(slurm_arrayid) # Convert to integer from 1-24

###-------------------------------####
# Import libraries
###------------------------------####
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(future) # Parallelize certain functions

plan("multiprocess", workers = 5)
options(future.globals.maxSize = 3000 * 1024^2)
options(future.seed = TRUE)
###-------------------------------####
# Import data
###------------------------------####
# Read in files
path <- "../processed"

file_list <- list.files(path = path, pattern = ".rds")

names(file_list) <- c('P5344_CW1', 'P5344_CW2', 'P5344_CW3', 'P5544_CW1', 'P5544_CW2', 'P5544_CW3', 'P5573_CW1', 'P5573_CW2', 'P5573_CW3', 'P5657_CW1', 'P5657_CW3',  'P5836_CW1', 'P5836_CW2', 
'P5836_CW3', 'P5877_CW1', 'P5877_CW2', 'P5877_CW3', 'P5903_CW1', 'P5903_CW2', 'P5903_CW3', 'P5963_CW1', 'P5963_CW2', 'P5963_CW3')

# grab only seurat rds for the specific array
file <- file_list[[i]]
names(file) <- names(file_list[i])

# Read in rds
seurat_object <- readRDS(paste(path, file, sep = ""))
Idents(seurat_object) <- "status"
seurat_object <- subset(seurat_object, idents = "unassigned", invert = T) # remove unassigned since not informative

seurat_object[['GT']] <- as.factor(ifelse(seurat_object$status == "singlet", "Singlet", "Doublet")) # Create ground truths

###-------------------------------####
# Call DoubletFinder
###------------------------------####
PCs <- 1:15 # Choose 15 PCs, reasonable based on examination
DoubletRate = (dim(seurat_object)[2]/0.57) * 4.6e-6 # Calculate doublet rate based on data from 10X
homotypic.prop <- 0.75

sweep.res.list <- paramSweep_v3(seurat_object, PCs = PCs, sct = TRUE)
gt.calls <- seurat_object@meta.data[rownames(sweep.res.list[[1]]), "GT"]
sweep.stats <- summarizeSweep(sweep.res.list, GT = TRUE, GT.calls = gt.calls)
bcmvn <- find.pK(sweep.stats)

pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK

nExp_poi <- round(DoubletRate*nrow(seurat_object@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seurat_object <- doubletFinder_v3(seurat_object, PCs = PCs, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)

DFname <- colnames(seurat_object@meta.data)[str_detect(string = colnames(seurat_object@meta.data), pattern = "DF.classifications*")] # pull out DF classification name
seurat_object[['DoubletFinder']] <- seurat_object[[DFname]] # Save DF to common name among all datasets
seurat_object[[DFname]] <- NULL # Delete DF classification

pANN <- colnames(seurat_object@meta.data)[str_detect(string = colnames(seurat_object@meta.data), pattern = "pANN")] # pull out unique pANN
seurat_object[['pANN']] <- seurat_object[[pANN]] # Save pANN to common name
seurat_object[[pANN]] <- NULL # Remove unique pANN

name <- names(file)
dir.create("../DoubletFinder")
saveRDS(seurat_object, file = paste("../DoubletFinder", name, ".rds", sep = ""))
