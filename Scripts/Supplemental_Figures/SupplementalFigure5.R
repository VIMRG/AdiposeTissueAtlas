###############################################################################
######-----------------SUPPLEMENTAL FIGURE 5------------------------------#####
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

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Code/Utils.R')
fig_dir <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Figures/"

###-------------------####
# Load Data
###-------------------####
path = "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Pseudobulk/PreDMvsNonDM/"
file_list <- list.files(path = path, full.names = TRUE, pattern = ".csv")
names <- list.files(path= path, pattern = ".csv")
names <- str_remove(names, ".csv")
names(file_list) <- names

file_list <- lapply(file_list, function(x) {
  read.csv(x, row.names = 1)
})

###-------------------####
# Run GSEA on files
###-------------------####
# GO BP GSEA
pathways_GO <- lapply(file_list, GO_GSEA) # GO Pathway analysis
pathways_GO <- pathways_GO[names(pathways_GO) != "EC.HIVpreDM.HIVnonDM"] # No enriched terms
pathways_GO <- lapply(pathways_GO, GSEA_order) # Order by NES
saveRDS(pathways_GO, file = "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Pseudobulk/PreDMvsNonDM/GO_pathways.rds") # Save so can skip next time

# Kegg GSEA
pathways_Kegg <- lapply(file_list, KEGG_GSEA) # KEGG Pathway analysis
pathways_Kegg <- lapply(pathways_Kegg, GSEA_order)
saveRDS(pathways_Kegg, file = "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Pseudobulk/PreDMvsNonDM/Kegg_pathways.rds") # Save so can skip next time

###-------------------####
# Supplemental Figure 5A: Pseudotime Trajectory Slingshot Mo-Mac to PVM
###-------------------####
# see Pseudotime_PVM.R

###-------------------####
# Supplemental Figure 5B: Pseudotime Slingshot Mo-Mac to PVM
###-------------------####
# see Pseudotime_PVM.R

###-------------------####
# Supplemental Figure 5C: LAM GSEA Plot
###-------------------####
LAM_gsea_plot <- GSEA_plot(df = pathways_GO$LAM.HIVpreDM.HIVnonDM, filename = paste(fig_dir, "SF5C.png", sep = "/"), num_path = 6, legend.pos = "top", height = 6, width = 10)

###-------------------####
# Supplemental Figure 5D: Upset Plot Higher Expression Prediabetic vs Non-DM
###-------------------####
# Upregulated Genes
Gene_Up <- lapply(file_list, function(df) {
  df$gene <- rownames(df)
  df <- df %>% dplyr::filter(log2FoldChange > 0 & padj < 0.1)
})

# Remove Completely Null (EC and Preadipocyte)
Gene_Up <- Gene_Up[!names(Gene_Up) %in% c("EC.HIVpreDM.HIVnonDM", "Preadipocyte.HIVpreDM.HIVnonDM")]
genes <- bind_rows(Gene_Up)
genes <- unique(genes$gene)

mat <- data.frame(matrix(data = NA, nrow = length(genes), ncol = 0))
rownames(mat) <- genes

mat <- lapply(Gene_Up, upset_gene)
mat <- bind_cols(mat)

colnames(mat) <- gsub(".HIVpreDM.HIVnonDM", "", names(Gene_Up))

m2 = make_comb_mat(mat)
png(file = paste(fig_dir, "SF5D.png", sep = "/"), res = 300, units = "in", height = 10, width = 10)
UpSet(m2, comb_order = order(comb_size(m2), decreasing = T))
dev.off()

###-------------------####
# Supplemental Figure 5E: Upset Plot Down
###-------------------####
# Down Genes
Gene_Down <- lapply(file_list, function(df) {
  df$gene <- rownames(df)
  df <- df %>% dplyr::filter(log2FoldChange < 0 & padj < 0.1)
})

# Remove Completely Null (EC and Preadipocyte)
Gene_Down <- Gene_Down[!names(Gene_Down) %in% c("IM.HIVpreDM.HIVnonDM", "Preadipocyte.HIVpreDM.HIVnonDM")]
genes <- bind_rows(Gene_Down)
genes <- unique(genes$gene)

mat <- data.frame(matrix(data = NA, nrow = length(genes), ncol = 0))
rownames(mat) <- genes

mat <- lapply(Gene_Down, upset_gene)
mat <- bind_cols(mat)

colnames(mat) <- gsub(".HIVpreDM.HIVnonDM", "", names(Gene_Up))

m2 = make_comb_mat(mat)
png(file = paste(fig_dir, "SF5E.png", sep = "/"), res = 300, units = "in", height = 10, width = 10)
UpSet(m2, comb_order = order(comb_size(m2), decreasing = T))
dev.off()
