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

# 2. MONOS SUBSET DGE ANALYSIS------------------------------------------------

###############################
# Load Dataset
###############################
Myeloid <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_RNA_11.1.rds')

###############################
# MAST DGE by Cell Type
###############################
Myeloid.markers <- RNAclust_Markers(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_cluster_RNAmarkers")
Myeloid.ADT <- ADTclust_Markers(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_cluster_ADTmarkers")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV+DM-)
###############################
#Combine celltype with study group
Myeloid <- CellType_Fun(Myeloid)

# Perform DGE
Myeloid.HIVDMvHIVnoDM <- HIV_comparison(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_DGE_by_DMstatus")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV-DM+)
###############################
Myeloid.HIVDMvHIVnegDM <- DM_comparison(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_DGE_by_HIVstatus")

# 3. MONOS SEURAT CLUSTERS DGE ANALYSIS------------------------------------------------

###############################
# MAST DGE by Seurat Clusters
###############################
Myeloid.markers <- RNAseurat_Markers(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_seuratcluster_RNAmarkers")
Myeloid.ADT <- ADTseurat_Markers(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_seuratcluster_ADTmarkers")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV+DM-)
###############################
# Combine study group and seurat clusters
Myeloid <- SeuratGroup_Fun(Myeloid)

# Perform DGE
Myeloid.HIVDMvHIVnoDM <- HIV_comp_seurat(MMyeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_seurat.DGE_by_DMstatus")

###############################
# MAST DGE by Study Group (HIV+DM+ vs HIV-DM+)
###############################
Myeloid.HIVDMvHIVnegDM <- DM_comp_seurat(Myeloid, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Myeloid/Myeloid_seurat.DGE_by_HIVstatus")

