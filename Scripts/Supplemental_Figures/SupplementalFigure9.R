###############################################################################
######----------------------------Supplemental Figure 9--------------------------------#####
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

###############################
# Supplemental Figure 9A: LAM to Myofibroblast/Cycling Myofibroblast
###############################
png(paste0(fig_dir, "/", "SF9A.png"), units = 'in',height = 5, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("LAM"), targets.use = c(18,31), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("LAM"), targets.use = c(18,31), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Supplemental Figure 9B: LAM to Preadipocyte
###############################
png(paste0(fig_dir, "/", "SF9B.png"), units = 'in',height = 5, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("LAM"), targets.use = c(1,36), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("LAM"), targets.use = c(1,36), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Supplemental Figure 9C: Myofibroblast/Cycling Myofibroblast to LAM
###############################
png(paste0(fig_dir, "/", "SF9C.png"), units = 'in',height = 6, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("Myofibroblast", "Cycling Myofibroblast"), targets.use = c(24), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("Myofibroblast", "Cycling Myofibroblast"), targets.use = c(24), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Supplemental Figure 9D: Preadipocyte to LAM
###############################
png(paste0(fig_dir, "/", "SF9D.png"), units = 'in',height = 6, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("Preadipocyte", "Adipose Progenitors"), targets.use = c(24), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("Preadipocyte", "Adipose Progenitors"), targets.use = c(24), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Supplemental Figure 9E: TEM to IM/LAM
###############################
png(paste0(fig_dir, "/", "SF9E.png"), units = 'in',height = 6, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("CD4 TEM", "CD8 TEM"), targets.use = c(21,24), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("CD4 TEM", "CD8 TEM"), targets.use = c(21,24), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()

###############################
# Supplemental Figure 9F: IM/LAM to TEM
###############################
png(paste0(fig_dir, "/", "SF9F.png"), units = 'in',height = 6, width = 10, res = 300)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c("IM", "LAM"), targets.use = c(9,14), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c("IM", 'LAM'), targets.use = c(9,14), font.size = 12, color.use = c("#136bcf", "#cf4813"), x.angle = 45)
gg1 + gg2
dev.off()
