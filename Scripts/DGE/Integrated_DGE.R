##############################################################################################
##------------------------- MAST RNA MARKER PARALLELIZATION ---------------------------------##
##------------------------- DATE: 8/18/2021 AUTHOR: -------------------------------##
## DESCRIPTION: The following code parallelizes the FindMarkers function from seurat for the
# large integrated dataset using MAST.
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
i <- as.numeric(slurm_arrayid) # Convert to integer from 1-62

clusters <- c("ASC 3", "PCOLCE+ Fibroblast", "Myofibroblast", "ASC 2", "Preadipocyte", "ASC 1", "Mature Preadipocyte 2", "Extracellular Preadipocyte", "PTGDS+ Extracellular Preadipocyte",
"Early Preadipocyte 1", "Early Preadipocyte 2", "Metallothionein+ Preadipocyte/ASC", "PTGDS+ Preadipocyte", "Early Preadipocyte 3", "Mature Preadipocyte 1", "CD9hi Preadipocyte",
"Early Preadipocyte 4", "MYOC+ Fibroblast", "ISG+ Preadipocyte", "DKK3+ Fibroblast", "Lipomyofibroblast", "Proliferating Myofibroblast", "Intermediate Capillary", "Capillary-EndoMT-like",
"Arterial", "Capillary", "Venous", "VSMC", "Pericyte", "Venous-EndoMT-like", "Arterial-EndoMT-like", "Naive B cell", "Memory B cell", "Plasmablast", "nMo", "PVM", "LAM", "Other Mac",
"Cycling Myeloid", "cMo", "cDC1", "cDC2B", "CCR7+ DC", "Other Mo", "Mo-Mac", "pDC", "CD57+ mNK", "CD4 TEM", "CD16+ mNK", "CD4 TCM", "CD8 TEMRA & Senescent", "CD4 Naive", "CD8 TEM", "CD8 TCM",
"CD4 TEMRA & Senescent", "Gamma Delta", "CD4 Regulatory", "CD8 Naive", "Immature NK", "MAIT", "Cycling T & NK", "ILC")

cluster <- clusters[i]

###############################
# Import Libraries
###############################
library(Seurat)
library(future)

plan("multiprocess", workers = 8)

###############################
# Function
###############################
RNAclust_Markers <- function(seurat_object, outpath, ident.1) {
    Idents(seurat_object) <- "Manual_Annotation"
    DGE <- FindMarkers(seurat_object, ident.1 = ident.1,
    test.use = "MAST", assay = "RNA", min.pct = 0.25, only.pos = T) # Perform DGE
    write.csv(DGE, file = paste(outpath, ident.1, ".csv", sep = ""))
    return(DGE)
}

# 2. STROMAL SUBSET DGE ANALYSIS------------------------------------------------

###############################
# Load Dataset
###############################
Integrated <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Final_Integrated.1.3.rds')

###############################
# MAST DGE by Cell Type
###############################
Integrated.markers <- RNAclust_Markers(Integrated, ident.1 = cluster, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/DGE/")
