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
tmp_dir <- "../SubsetAnalysis_HIVneg"

Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal_HIVneg.rds")) 
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage_HIVneg.rds")) 
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4_HIVneg.rds")) 
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8_HIVneg.rds")) 
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid_HIVneg.rds")) 
Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular_HIVneg.rds")) 
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid_HIVneg.rds")) 

# Collapse Stromal
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Early Preadipocyte", "ECM-Producing Early Preadipocyte", "Mature Preadipocyte"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Adipose Progenitors"
Stromal[['Annotation']] <- Idents(Stromal)

# Subset out Monocytes
Idents(Myeloid) <- "Annotation"
Myeloid <- subset(Myeloid, idents = c("cMo", "nMo", "Other Mo", "ISG+ Mo", "pDC", "cDC1", "cDC2B", "Migratory DC"))

# Lymphoid
cellNames <- c(colnames(CD8), colnames(CD4))
Lymphoid <- subset(Lymphoid, cells = cellNames, invert = T)
Lymphoid <- subset(Lymphoid, idents = c("CD57 mNK", "CD16 mNK", "Immature NK", "MAIT", "Gamma Delta", "ILC"))

# Vascular
Vascular <- subset(Vascular, idents = c("Venous EC", "Arterial EC", "Capillary EC"))

Lymphoid[['SCT']] <- NULL # Get rid of SCT slot otherwise it won't merge
CD4[['SCT']] <- NULL
CD8[['SCT']] <- NULL
Vascular[['SCT']] <- NULL
Stromal[['SCT']] <- NULL
Myeloid[['SCT']] <- NULL
Macrophage[['SCT']] <- NULL

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
HIV_pos <- subset(CellChat, idents = "HIV+ diabetic")
HIV_neg <- subset(CellChat, idents = c("HIV- diabetic"))

cellchat_HIVpos <- createCellChat(object = HIV_pos, meta = HIV_pos@meta.data, group.by = "Annotation")
cellchat_HIVneg <- createCellChat(object = HIV_neg, meta = HIV_neg@meta.data, group.by = "Annotation")
rm(HIV_pos)
rm(HIV_neg)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

#cellchat@DB <- CellChatDB.use
cellchat_HIVpos@DB <- CellChatDB.use
cellchat_HIVneg@DB <- CellChatDB.use

#########################
# Preprocessing
#########################
cellchat_HIVpos <- subsetData(cellchat_HIVpos)
cellchat_HIVneg <- subsetData(cellchat_HIVneg)

cellchat_HIVpos <- identifyOverExpressedGenes(cellchat_HIVpos)
cellchat_HIVneg <- identifyOverExpressedGenes(cellchat_HIVneg)

cellchat_HIVpos <- identifyOverExpressedInteractions(cellchat_HIVpos)
cellchat_HIVneg <- identifyOverExpressedInteractions(cellchat_HIVneg)

cellchat_HIVpos <- projectData(cellchat_HIVpos, PPI.human)
cellchat_HIVneg <- projectData(cellchat_HIVneg, PPI.human)

#########################
# Inference HIV+ non-DM
#########################
cellchat_HIVpos <- computeCommunProb(cellchat_HIVpos, population.size = TRUE, type = "truncatedMean", trim = 0.1)
cellchat_HIVpos <- filterCommunication(cellchat_HIVpos, min.cells = 10)
cellchat_HIVpos <- computeCommunProbPathway(cellchat_HIVpos)
cellchat_HIVpos <- aggregateNet(cellchat_HIVpos)
cellchat_HIVpos <- netAnalysis_computeCentrality(cellchat_HIVpos, slot.name = "netP")

tmp_dir <- "../CellChat"
dir.create(tmp_dir)

saveRDS(cellchat_HIVpos, file = paste0(tmp_dir, "/", "cellchat_HIVpos.rds"))
rm(cellchat_HIVpos)

#########################
# Inference HIV+ DM
#########################
cellchat_HIVneg <- computeCommunProb(cellchat_HIVneg, population.size = TRUE, , type = "truncatedMean", trim = 0.1)
cellchat_HIVneg <- filterCommunication(cellchat_HIVneg, min.cells = 10)
cellchat_HIVneg <- computeCommunProbPathway(cellchat_HIVneg)
cellchat_HIVneg <- aggregateNet(cellchat_HIVneg)
cellchat_HIVneg <- netAnalysis_computeCentrality(cellchat_HIVneg, slot.name = "netP")

saveRDS(cellchat_HIVneg, file = paste0(tmp_dir, "/", "cellchat_HIVneg.rds"))
