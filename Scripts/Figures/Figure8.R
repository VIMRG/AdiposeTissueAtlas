###############################################################################
######----------------------------FIGURE 8--------------------------------#####
###############################################################################

###############################
# set directory
###############################
library(Seurat)
library(ggplot2)
library(CellChat)

###-------------------####
# Directories
###-------------------####
fig_dir <- "../Figures"
tmp_dir <- "../CellChat"

###############################
# Load Data
###############################
cellchat_HIVnonDM <- readRDS(paste0(tmp_dir, "/", "cellchat_HIVnonDM.rds")) # non-diabetic
cellchat_HIVGI <- readRDS(paste0(tmp_dir, "/", "cellchat_HIVGI.rds")) # Glucose Intolerant

object.list <- list(nonDM = cellchat_HIVnonDM, GI = cellchat_HIVGI)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/CodeRevised/Utils.R')

###############################
# Figure 8A: Overview of Cell-Cell Interactions
###############################
png(paste0(fig_dir, "/", "Figure8A.png"), units = 'in',height = 20, width = 25, res = 300)
gg1 <- netVisual_heatmap(cellchat, font.size = 12, cluster.rows = TRUE, cluster.cols = TRUE)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", font.size = 12, cluster.rows = TRUE, cluster.cols = TRUE)
gg1 + gg2
dev.off()

###############################
# Figure 8B: IM to Myofibroblast/Cycling Myofibroblast
###############################
png(paste0(fig_dir, "/", "Figure8B.png"), units = 'in',height = 5, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("IM"), targets.use = c(18,31), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("IM"), targets.use = c(18,31), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Figure 8C: IM to Preadipocyte
###############################
png(paste0(fig_dir, "/", "Figure8C.png"), units = 'in',height = 5, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("IM"), targets.use = c(1,36), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("IM"), targets.use = c(1,36), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Figure 8D: Myofibroblast/Cycling Myofibroblast to IM
###############################
png(paste0(fig_dir, "/", "Figure8D.png"), units = 'in',height = 6, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("Myofibroblast", "Cycling Myofibroblast"), targets.use = c(21), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("Myofibroblast", "Cycling Myofibroblast"), targets.use = c(21), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Figure 8E: Preadipocyte to IM
###############################
png(paste0(fig_dir, "/", "Figure8E.png"), units = 'in',height = 6, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("Preadipocyte", "Adipose Progenitors"), targets.use = c(21), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("Preadipocyte", "Adipose Progenitors"), targets.use = c(21), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Figure 8F: CD4 TEM to Myofibroblast/Cycling Myofibroblast
###############################
png(paste0(fig_dir, "/", "Figure8F.png"), units = 'in',height = 5, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("CD4 TEM"), targets.use = c(18,31), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("CD4 TEM"), targets.use = c(18,31), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Figure 8G: CD4 TEM to Preadipocyte
###############################
png(paste0(fig_dir, "/", "Figure8G.png"), units = 'in',height = 5, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("CD4 TEM"), targets.use = c(1,36), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("CD4 TEM"), targets.use = c(1,36), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()



