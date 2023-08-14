###############################################################################
######------------------SUPPLEMENTAL FIGURE 8-----------------------------#####
###############################################################################

###-------------------####
# Import Libraries
###-------------------####
library(ggplot2)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(PResiduals)
library(rms)
library(gridExtra)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(stringi)
library(circlize)

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/CodeRevised/Utils.R')

###-------------------####
# Directories
###-------------------####
Subset_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/SubsetAnalysis"
AllCells_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/Dialogue/AllCells" 
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/FiguresRevised"

###-------------------####
# Load Seurat Object
###-------------------####
Stromal <- readRDS(paste0(Subset_dir, "/", "Stromal_1.7.rds"))
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage_revised.rds"))
CD4 <- readRDS(paste0(Subset_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(Subset_dir, "/", "CD8.rds"))

# Process Seurat Objects to Match DIALOGUE
# CD4 T cells All (>30 cells)
# Get samples to use
metadata <- CD4@meta.data
HATIMID_CD4 <- metadata %>% dplyr::count(HATIMID) %>% filter(n > 30) %>% dplyr::select(HATIMID)
metadata <- CD8@meta.data
HATIMID_CD8 <- metadata %>% dplyr::count(HATIMID) %>% filter(n > 30) %>% dplyr::select(HATIMID)

CD4[["Annotation"]] <- "CD4"
CD4 <- subset(CD4, cells = colnames(CD4)[CD4$HATIMID %in% HATIMID_CD4$HATIMID])
CD4[['HATIMID']] <- droplevels(CD4$HATIMID)
CD4[['HATIMID']] <- as.character(CD4$HATIMID)

metadata <- CD4@meta.data
metadata <- left_join(metadata, data_hatim[,c("HATIMID", "hba1c")], by = "HATIMID")
hba1c <- metadata %>% dplyr::select(hba1c)
rownames(hba1c) <- colnames(CD4)
CD4 <- AddMetaData(CD4, hba1c)

# CD8 T cells
CD8[["Annotation"]] <- "CD8"
CD8 <- subset(CD8, cells = colnames(CD8)[CD8$HATIMID %in% HATIMID_CD8$HATIMID])
CD8[['HATIMID']] <- droplevels(CD8$HATIMID)
CD8[['HATIMID']] <- as.character(CD8$HATIMID)

metadata <- CD8@meta.data
metadata <- left_join(metadata, data_hatim[,c("HATIMID", "hba1c")], by = "HATIMID")
hba1c <- metadata %>% dplyr::select(hba1c)
rownames(hba1c) <- colnames(CD8)
CD8 <- AddMetaData(CD8, hba1c)

# Macrophage Collapse
metadata <- Macrophage@meta.data
HATIMID_Macrophage <- metadata %>% dplyr::count(HATIMID) %>% filter(n > 30) %>% dplyr::select(HATIMID)
Macrophage <- subset(Macrophage, cells = colnames(Macrophage)[Macrophage$HATIMID %in% HATIMID_Macrophage$HATIMID])
Macrophage[['HATIMID']] <- droplevels(Macrophage$HATIMID)
Macrophage[['HATIMID']] <- as.character(Macrophage$HATIMID)

metadata <- Macrophage@meta.data
metadata <- left_join(metadata, data_hatim[,c("HATIMID", "hba1c")], by = "HATIMID")
hba1c <- metadata %>% dplyr::select(hba1c)
rownames(hba1c) <- colnames(Macrophage)
Macrophage <- AddMetaData(Macrophage, hba1c)

# Stromal
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 1", "Mature Preadipocyte 2", "Early Preadipocyte", "ECM-Producing Early Preadipocyte"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Progenitor"
Stromal[['Annotation']] <- Idents(Stromal)
Stromal[['HATIMID']] <- droplevels(Stromal$HATIMID)
Stromal[['HATIMID']] <- as.character(Stromal$HATIMID)

metadata <- Stromal@meta.data
metadata <- left_join(metadata, data_hatim[,c("HATIMID", "hba1c")], by = "HATIMID")
hba1c <- metadata %>% dplyr::select(hba1c)
rownames(hba1c) <- colnames(Stromal)
Stromal <- AddMetaData(Stromal, hba1c)

###-------------------####
# Load Data
###-------------------####
R <- readRDS(paste(AllCells_dir, "Results.rds", sep = "/"))

###-------------------####
# Supplemental Figure 8A: PVM
###-------------------####
plot <- Heatmap_fun(seurat_object = Macrophage, subset = T, metadata = data_hatim, cluster = "PVM", genes = c("ATF4", "CEBPB", "CCL2", "CCL8", "EGR1", "EGR2", "FOSB", "KLF4", 
"A2M", "COX5B", "FABP4", "IFI27", "IFITM1", "IL1RN", "LPL", "TREM2", "TYROBP"), 
filename = "SF8A.png",tmp_dir = fig_dir, R = R, MCP_dim = 1, height = 10, width = 15)

###-------------------####
# Supplemental Figure 8B: Mo-Mac 1
###-------------------####
plot <- Heatmap_fun(seurat_object = Macrophage, subset = T, metadata = data_hatim, cluster = "Mo-Mac 1", genes = c("ATF4", "CCL2", "CEBPB", "EGR1", "EGR2", "KLF4", "SERPINE1", "TNF", "AREG", "B2M",
"CASP9", "CD74", "CSTA", "IFI27", "ISG15", "MT2A"), 
filename = "SF8B.png", tmp_dir = fig_dir, R = R, height = 10, width = 15)

###-------------------####
# Supplemental Figure 8C:IM
###-------------------####
plot <- Heatmap_fun(seurat_object = Macrophage, subset = T, metadata = data_hatim, cluster = "IM", genes = c("AREG", "C1QB", "IL18", "TYROBP", "KLF4", "CEBPB", "TNF","TREM2"), 
filename = "SF8C.png", tmp_dir = fig_dir, R = R, height = 10, width = 15)

###-------------------####
#  Supplemental Figure 8D: Myofibroblast
###-------------------####
plot <- Heatmap_fun(seurat_object = Stromal, subset = T, metadata = data_hatim, cluster = "Myofibroblast", genes = c("PRSS23", "TAGLN", "BGN", "TNC", "COL4A1", "CCN1", "CD34", "CD55", "MMP2"),tmp_dir = fig_dir, filename = "SF8D.png", R = R, MCP_dim = 1,height = 10, width = 15)

###-------------------####
# Supplemental Figure 8E: Preadipocyte
###-------------------####
##############################
# Generate Data for Genes Up in MCP1
##############################
MCP_dim = 1
clusters <- c("CD4", "CD8", "IM", "PVM", "LAM", "Mo-Mac 1", "Myofibroblast") # clusters to include
HATIM <- unique(Stromal$HATIMID) # Grab all HATIM IDs from Stromal
HATIM <- HATIM[HATIM %in% unique(CD4$HATIMID) & HATIM %in% unique(Macrophage$HATIMID) & HATIM %in% unique(CD8$HATIMID)] # Filter to make sure only include individuals in all cell types

# Empty Lists/datasets
MCP_up_list <- list()
cluster_up_list <- list()
gene_up_list <- list()
merged_up <- data.frame(matrix(data = NA, nrow = length(HATIM), ncol = 0))

for(cluster in clusters) { # Specify what to set seurat_object equal to
    if(cluster %in% "CD4") {
        seurat_object <- CD4
    } else if(cluster %in% "CD8") {
        seurat_object <- CD8
    } else if (cluster %in% c("LAM", "PVM", "IM", "Mo-Mac 1")) {
        seurat_object <- Macrophage
    } else {
        seurat_object <- Stromal
    }
    
    clustergene <- gsub(" ", ".", cluster) # Add "." to name to match Dialogue names
    cluster.up <- R$MCPs[[paste0("MCP",MCP_dim)]][[paste0(clustergene, ".up")]] # Grab all Up genes
    cluster_genes <- c(cluster.up)
    
    if(cluster %in% c("LAM", "PVM", "IM", "Mo-Mac 1", "Myofibroblast")) { # Subset on these cell
        Idents(seurat_object) <- "Annotation"
        seurat_object <- subset(seurat_object, idents = cluster)
    }
    
    Idents(seurat_object) <- "HATIMID"
    seurat_object <- seurat_object[,seurat_object$HATIMID %in% HATIM] # filter to only HATIM ID in all cell types
    
    Avg <- AverageExpression(seurat_object, group.by = "HATIMID") # Grab average expression by participant
    Avg <- Avg[[1]]
    
    # Pull out cluster-specific genes
    Avg <- Avg[rownames(Avg) %in% cluster_genes,] # only grab genes in MCP1 Up
    
    # Filter out lowly expressed genes (likely contaminating)
    Exp <- data.frame(Expression = rowMeans(seurat_object@assays$RNA@data[rownames(seurat_object@assays$RNA@data) %in% cluster_genes,]))
    Exp$Gene <- rownames(Exp)
    Exp <- Exp[order(Exp$Expression, decreasing = T),]
    Exp <- Exp[Exp$Expression > 0.15,]
    Avg <- Avg[rownames(Avg) %in% Exp$Gene,]
    names <- colnames(Avg)
    names(names) <- NULL
    
    # Scale Expression
    mat_scaled <- t(apply(Avg,1,scale))
    colnames(mat_scaled) <- names
    
    # Generate Metadata
    DF <- data.frame(mat_scaled)
    colnames(DF) <- names
    DF$MCP <- ifelse(rownames(DF) %in% cluster.up, "Up", "Down")
    DF$cluster <- cluster
    DF$gene <- rownames(DF)
    MCP_label <- DF$MCP
    cluster_label <- DF$cluster
    gene_label <- DF$gene
    DF$MCP <- NULL
    DF$cluster <- NULL
    DF$gene <- NULL
    tDF <- data.frame(t(DF))
    if(cluster == "Myofibroblast") {
        tDF$HATIMID <- as.character(rownames(tDF))
        DF <- left_join(tDF, data_hatim[,c("HATIMID", "StudyGroup", "BMI", "Age", "Sex", "FBG", "hba1c")], by = "HATIMID")
        merged_up <- cbind(merged_up, DF)
    } else {
        DF <- tDF
        merged_up <- cbind(merged_up, DF)
        print(paste0("Completed ", cluster))
    }
    
    # Append data to list
    MCP_up_list <- append(MCP_up_list, MCP_label)
    cluster_up_list <- append(cluster_up_list, cluster_label)
    gene_up_list <- append(gene_up_list, gene_label)
}

##############################
# Generate Data for Genes Down in MCP1
##############################
MCP_dim = 1
clusters <- c("CD4", "CD8", "IM", "PVM", "LAM", "Mo-Mac 1", "Myofibroblast")
HATIM <- unique(Stromal$HATIMID)
HATIM <- HATIM[HATIM %in% unique(CD4$HATIMID) & HATIM %in% unique(Macrophage$HATIMID) & HATIM %in% unique(CD8$HATIMID)]

MCP_down_list <- list()
cluster_down_list <- list()
gene_down_list <- list()
merged_down <- data.frame(matrix(data = NA, nrow = length(HATIM), ncol = 0))

for(cluster in clusters) {
    if(cluster %in% "CD4") {
        seurat_object <- CD4
    } else if(cluster %in% "CD8") {
        seurat_object <- CD8
    } else if (cluster %in% c("LAM", "PVM", "IM", "Mo-Mac 1")) {
        seurat_object <- Macrophage
    } else {
        seurat_object <- Stromal
    }
    
    clustergene <- gsub(" ", ".", cluster)
    #cluster.up <- R$MCPs[[paste0("MCP",MCP_dim)]][[paste0(clustergene, ".up")]]
    cluster.down <- R$MCPs[[paste0("MCP",MCP_dim)]][[paste0(clustergene, ".down")]]
    cluster_genes <- c(cluster.down)
    
    if(cluster %in% c("LAM", "PVM", "IM", "Mo-Mac 1", "Myofibroblast")) {
        Idents(seurat_object) <- "Annotation"
        seurat_object <- subset(seurat_object, idents = cluster)
    }
    
    Idents(seurat_object) <- "HATIMID"
    seurat_object <- seurat_object[,seurat_object$HATIMID %in% HATIM]
    
    Avg <- AverageExpression(seurat_object, group.by = "HATIMID")
    Avg <- Avg[[1]]
    
    # Pull out cluster-specific genes
    Avg <- Avg[rownames(Avg) %in% cluster_genes,]
    
    Exp <- data.frame(Expression = rowMeans(seurat_object@assays$RNA@data[rownames(seurat_object@assays$RNA@data) %in% cluster_genes,]))
    Exp$Gene <- rownames(Exp)
    Exp <- Exp[order(Exp$Expression, decreasing = T),]
    Exp <- Exp[Exp$Expression > 0.15,]
    Avg <- Avg[rownames(Avg) %in% Exp$Gene,]
    names <- colnames(Avg)
    names(names) <- NULL
    
    # Scale Expression
    mat_scaled <- t(apply(Avg,1,scale))
    colnames(mat_scaled) <- names
    
    # Generate Metadata
    DF <- data.frame(mat_scaled)
    colnames(DF) <- names
    DF$MCP <- ifelse(rownames(DF) %in% cluster.down, "Down", "Up")
    DF$cluster <- cluster
    DF$gene <- rownames(DF)
    MCP_label <- DF$MCP
    cluster_label <- DF$cluster
    gene_label <- DF$gene
    DF$MCP <- NULL
    DF$cluster <- NULL
    DF$gene <- NULL
    tDF <- data.frame(t(DF))
    if(cluster == "Myofibroblast") {
        tDF$HATIMID <- as.character(rownames(tDF))
        DF <- left_join(tDF, data_hatim[,c("HATIMID", "StudyGroup", "BMI", "Age", "Sex", "FBG", "hba1c")], by = "HATIMID")
        merged_down <- cbind(merged_down, DF)
    } else {
        DF <- tDF
        merged_down <- cbind(merged_down, DF)
    }
    
    MCP_down_list <- append(MCP_down_list, MCP_label)
    cluster_down_list <- append(cluster_down_list, cluster_label)
    gene_down_list <- append(gene_down_list, gene_label)
}

# Merge Up and Down lists
gene_list <- c(unlist(gene_up_list), unlist(gene_down_list))
MCP_list <- c(unlist(MCP_up_list), unlist(MCP_down_list))
cluster_list <- c(unlist(cluster_up_list), unlist(cluster_down_list))
merged <- cbind(merged_up[1:608], merged_down) # Only want metadata once so only include gene columns from Up genes

# Heatmap
    ha = HeatmapAnnotation(MP = unlist(MCP_list), Cluster = unlist(cluster_list), col = list(MP = c("Up" = "red", "Down" = "blue"), Cluster = c("CD4" = "blue", "CD8" = "green", "PVM" = "orange", "LAM" = "yellow", "IM" = "black", "Mo-Mac 1" = "purple", "Myofibroblast" = "darkorange")))
    col_fun = colorRamp2(c(25,35,45), c("blue", "white", "red"))
    col_fun2 = colorRamp2(c(25, 45, 65), c("#2D8A3D", "white", "#F49F2C"))
    col_fun3 = colorRamp2(c(80, 110, 150), c("blue", "white", "red"))
    col_fun4 = colorRamp2(c(5, 6, 8), c("#79FF33", "white", "#FF5733"))
    row_ha = rowAnnotation(Group = merged$StudyGroup, Sex = merged$Sex, BMI = merged$BMI, Age = merged$Age, FBG = merged$FBG, HbA1c = merged$hba1c, col = list(Group = c("HIV+ non-diabetic"="#91C46C", "HIV+ prediabetic" = "#7ABAF2", "HIV+ diabetic" = "#FFBE00"),
    Sex = c("Male" = "dodgerblue", "Female" = "orange"), BMI = col_fun, Age = col_fun2, FBG = col_fun3, HbA1c = col_fun4))

   index <- which(colnames(merged) %in% c("TREM2", "IFNG", "TNF", "AREG"))
    label <- colnames(merged)[index]
    column_ha <- columnAnnotation(foo = anno_mark(at = index, labels = label, which = "column", side = "bottom"))

png(file = paste(fig_dir, "SF8E.png", sep = "/"), res = 300, width = 20, height = 10, units = 'in') 
Heatmap(merged[,c(1:1054)], top_annotation = ha, right_annotation = row_ha, bottom_annotation = column_ha, show_column_names = F, cluster_columns = F) # take only gene columns (exclude metadata columns)
dev.off()

# Heatmap
    ha = HeatmapAnnotation(MP = unlist(MCP_list), Cluster = unlist(cluster_list), col = list(MP = c("Up" = "red", "Down" = "blue")), Cluster = c("CD4" = "blue", "CD8" = "green", "PVM" = "orange", "LAM" = "yellow", "IM" = "black", "Mo-Mac 1" = "purple", "Myofibroblast" = "darkorange"))
    col_fun = colorRamp2(c(25,35,45), c("blue", "white", "red"))
    col_fun2 = colorRamp2(c(25, 45, 65), c("#2D8A3D", "white", "#F49F2C"))
    col_fun3 = colorRamp2(c(80, 110, 150), c("blue", "white", "red"))
    col_fun4 = colorRamp2(c(5, 6, 8), c("#79FF33", "white", "#FF5733"))
    row_ha = rowAnnotation(Sex = merged$Sex, BMI = merged$BMI, Age = merged$Age, FBG = merged$FBG, HbA1c = merged$hba1c, Group = merged$StudyGroup, col = list(Group = c("HIV+ non-diabetic"="#91C46C", "HIV+ prediabetic" = "#7ABAF2", "HIV+ diabetic" = "#FFBE00"),
    Sex = c("Male" = "dodgerblue", "Female" = "orange"), BMI = col_fun, Age = col_fun2, FBG = col_fun3, HbA1c = col_fun4))

   index <- which(colnames(merged) %in% c("TREM2", "IFNG", "TNF", "AREG"))
    label <- colnames(merged)[index]
    column_ha <- columnAnnotation(foo = anno_mark(at = index, labels = label, which = "column", side = "bottom"))

png(file = paste(fig_dir, "SF8E_Legends.png", sep = "/"), res = 300, width = 20, height = 10, units = 'in') 
Heatmap(merged[,c(1:1054)], top_annotation = ha, right_annotation = row_ha, bottom_annotation = column_ha, show_column_names = F, cluster_columns = F) # take only gene columns (exclude metadata columns)
dev.off()

