##############################################################################################
##------------------------- RNA MARKER PARALLELIZATION ---------------------------------##
##------------------------- DATE: 8/18/2021 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: The following code parallelizes the FindMarkers function from seurat for the
# large integrated dataset.
##############################################################################################


# 1. SET UP AND LOADING--------------------------------------------------------

###############################
# Set Seed for Reproducibility
###############################
set.seed(7612)

###############################
# Extract Slurm Array Data
###############################
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID') # Get environmental variable from slurm array
i <- as.numeric(slurm_arrayid) # Convert to integer from 1-55

###############################
# Import Libraries
###############################
library(Seurat)

###############################
# Function
###############################
RNAclust_Markers <- function(seurat_object, outpath, ident.1) {
    Idents(seurat_object) <- "Annotation"
    DGE <- FindMarkers(seurat_object, ident.1 = ident.1,
    assay = "RNA", min.pct = 0.1, only.pos = T) # Perform DGE
    ident.1 = gsub(" ", ".", ident.1)
    write.csv(DGE, file = paste(outpath, ident.1, ".csv", sep = ""))
    return(DGE)
}

# 2. DGE ANALYSIS------------------------------------------------

###############################
# Load Dataset
###############################
dir <- "../Final_Integration"
Integrated <- readRDS(paste(dir, "Integrated.rds", sep = "/"))

###############################
# DGE by Cell Type
###############################
clusters <- as.character(unique(Integrated$Annotation))
cluster <- clusters[i]

Integrated.markers <- RNAclust_Markers(Integrated, ident.1 = cluster, outpath = "../CellTypes")
