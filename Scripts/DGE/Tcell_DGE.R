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

# 2. T CELL SUBSET DGE ANALYSIS------------------------------------------------

###############################
# Load Dataset
###############################
Tcells <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells_RNA_10.30.rds')

###############################
# MAST DGE by Cell Type
###############################
Tcells.markers <- RNAclust_Markers(Tcells, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcell_cluster_RNAmarkers")
Tcells.ADT <- ADTclust_Markers(Tcells, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcell_cluster_ADTmarkers")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV+DM-)
###############################
# Combine celltype with study group
Tcells <- CellType_Fun(Tcells)

# Perform DGE
Tcells.HIVDMvHIVnoDM <- HIV_comparison(Tcells, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcells_DGE_by_DMstatus")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV-DM+)
###############################
Tcells.HIVDMvHIVnegDM <- DM_comparison(Tcells, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcells_DGE_by_HIVstatus")

# 3. T CELL SEURAT CLUSTERS DGE ANALYSIS------------------------------------------------

###############################
# MAST DGE by Seurat Clusters
###############################
Tcells.markers <- RNAseurat_Markers(Tcells, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcell_seuratcluster_RNAmarkers")
Tcells.ADT <- ADTseurat_Markers(Tcells, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcell_seuratcluster_ADTmarkers")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV+DM-)
###############################
# Combine study group and seurat clusters
Tcells <- SeuratGroup_Fun(Tcells)

# Perform DGE
Tcells.HIVDMvHIVnoDM <- HIV_comp_seurat(Tcells, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcells_seurat.DGE_by_DMstatus")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV-DM+)
###############################
Tcells.HIVDMvHIVnegDM <- DM_comp_seurat(Tcells, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Tcells/Tcells_seurat.DGE_by_HIVstatus")

