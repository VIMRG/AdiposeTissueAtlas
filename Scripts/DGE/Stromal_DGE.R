# 1. SET UP AND LOADING--------------------------------------------------------

###############################
# Set Seed for Reproducibility
###############################
set.seed(7612)

###############################
# Set Seed for Reproducibility
###############################
library(Seurat)
library(future)

source("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/DGE/DGE_Function.R")
plan("multiprocess", workers = 4)

# 2. STROMAL SUBSET DGE ANALYSIS------------------------------------------------

###############################
# Load Dataset
###############################
Stromal <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_1.3.rds')

###############################
# MAST DGE by Cell Type
###############################
Stromal.markers <- RNAclust_Markers(Stromal, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_cluster_RNAmarkers")
Stromal.ADT <- ADTclust_Markers(Stromal, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_cluster_ADTmarkers")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV+DM-)
###############################
# Combine celltype with study group
Stromal <- CellType_Fun(Stromal)

# Perform DGE
Stromal.HIVDMvHIVnoDM <- HIV_comparison(Stromal, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_DGE_by_DMstatus")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV-DM+)
###############################
Stromal.HIVDMvHIVnegDM <- DM_comparison(Stromal, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_DGE_by_HIVstatus")

# 3. STROMAL SEURAT CLUSTERS DGE ANALYSIS------------------------------------------------

###############################
# MAST DGE by Seurat Clusters
###############################
Stromal.markers <- RNAseurat_Markers(Stromal, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_seuratcluster_RNAmarkers")
Stromal.ADT <- ADTseurat_Markers(Stromal, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_seuratcluster_ADTmarkers")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV+DM-)
###############################
# Combine study group and seurat clusters
Stromal <- SeuratGroup_Fun(Stromal)

# Perform DGE
Stromal.HIVDMvHIVnoDM <- HIV_comp_seurat(Stromal, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_seurat.DGE_by_DMstatus")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV-DM+)
###############################
Stromal.HIVDMvHIVnegDM <- DM_comp_seurat(Stromal, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_seurat.DGE_by_HIVstatus")

