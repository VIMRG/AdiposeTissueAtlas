##############################################################################################
##------------------------- CELL CHAT ---------------------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will evaluate cell-cell communication pathways
##############################################################################################

# 1. SETTING UP ENVIRONMENT------------------------------------------------------------------

###############################
# set seed
###############################
set.seed(7612) # Reproducibility

###############################
# set directory
###############################
library(Seurat)
library(ggplot2)
library(dplyr)
library(gplots)
library(stats)
library(scater)
library(future)
library(scater)
library(tidyverse)
library(CellChat)
library(patchwork)

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 3000 * 1024^2)
options(future.seed = TRUE)

###############################
# Load Data
###############################
tmp_dir <- "../SubsetAnalysis"

Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))

# Collapse Stromal
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Early Preadipocyte", "ECM-Producing Early Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Adipose Progenitors"
Stromal[['Annotation']] <- Idents(Stromal)

# Subset out Monocytes
Idents(Myeloid) <- "Annotation"
Myeloid <- subset(Myeloid, idents = c("cMo", "nMo", "Other Mo", "ISG+ Mo", "pDC", "cDC1", "cDC2B", "Migratory DC"))

# T cells
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8.rds"))

# Lymphoid
cellNames <- c(colnames(CD8), colnames(CD4))
Lymphoid <- subset(Lymphoid, cells = cellNames, invert = T)
Lymphoid <- subset(Lymphoid, idents = c("CD57 mNK", "CD16 mNK", "Immature NK", "MAIT", "Gamma Delta", "ILC"))

# Vascular
Vascular <- subset(Vascular, idents = c("Venous EC", "Arterial EC", "Capillary EC"))

CellChat <- merge(x = Stromal, y = c(Macrophage, Myeloid, CD4, CD8, Vascular, Lymphoid))

rm(Stromal)
rm(Macrophage)
rm(Myeloid)
rm(Vascular)
rm(Lymphoid)
rm(CD4)
rm(CD8)

# 2. GENERATE CELL CHAT OBJECT_______________________________________________________
#########################
# Create CellChat Objects
#########################

Idents(CellChat) <- "StudyGroup"
HIV_nonDM <- subset(CellChat, idents = "HIV+ non-diabetic")
HIV_GI <- subset(CellChat, idents = c("HIV+ prediabetic", "HIV+ diabetic"))
#rm(CellChat)

cellchat_HIVnonDM <- createCellChat(object = HIV_nonDM, meta = HIV_nonDM@meta.data, group.by = "Annotation")
cellchat_HIVGI <- createCellChat(object = HIV_GI, meta = HIV_GI@meta.data, group.by = "Annotation")
rm(HIV_GI)
rm(HIV_nonDM)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

cellchat_HIVnonDM@DB <- CellChatDB.use
cellchat_HIVGI@DB <- CellChatDB.use

#########################
# Preprocessing
#########################
cellchat_HIVnonDM <- subsetData(cellchat_HIVnonDM)
cellchat_HIVGI <- subsetData(cellchat_HIVGI)

cellchat_HIVnonDM <- identifyOverExpressedGenes(cellchat_HIVnonDM)
cellchat_HIVGI <- identifyOverExpressedGenes(cellchat_HIVGI)

cellchat_HIVnonDM <- identifyOverExpressedInteractions(cellchat_HIVnonDM)
cellchat_HIVGI <- identifyOverExpressedInteractions(cellchat_HIVGI)

cellchat_HIVnonDM <- projectData(cellchat_HIVnonDM, PPI.human)
cellchat_HIVGI <- projectData(cellchat_HIVGI, PPI.human)

#########################
# Inference HIV+ non-DM
#########################
cellchat_HIVnonDM <- computeCommunProb(cellchat_HIVnonDM, population.size = TRUE, type = "truncatedMean", trim = 0.1)
cellchat_HIVnonDM <- filterCommunication(cellchat_HIVnonDM, min.cells = 10)
cellchat_HIVnonDM <- computeCommunProbPathway(cellchat_HIVnonDM)
cellchat_HIVnonDM <- aggregateNet(cellchat_HIVnonDM)
cellchat_HIVnonDM <- netAnalysis_computeCentrality(cellchat_HIVnonDM, slot.name = "netP")

tmp_dir <- "../CellChat"
dir.create(tmp_dir)

saveRDS(cellchat_HIVnonDM, file = paste0(tmp_dir, "/", "cellchat_HIVnonDM.rds"))
rm(cellchat_HIVnonDM)

#########################
# Inference HIV+ DM
#########################
cellchat_HIVGI <- computeCommunProb(cellchat_HIVGI, population.size = TRUE, , type = "truncatedMean", trim = 0.1)
cellchat_HIVGI <- filterCommunication(cellchat_HIVGI, min.cells = 10)
cellchat_HIVGI <- computeCommunProbPathway(cellchat_HIVGI)
cellchat_HIVGI <- aggregateNet(cellchat_HIVGI)
cellchat_HIVGI <- netAnalysis_computeCentrality(cellchat_HIVGI, slot.name = "netP")

saveRDS(cellchat_HIVGI, file = paste0(tmp_dir, "/", "cellchat_HIVGI.rds"))
