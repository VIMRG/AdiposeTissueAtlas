###############################################################################
######----------------------------FIGURE 6--------------------------------#####
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
source('../Utils.R')
fig_dir <- "../Figures"

###-------------------####
# Load Data: Note, this is Supplemental Table 8
###-------------------####
path = "../PreDMvsNonDM"
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
pathways_GO <- pathways_GO[names(pathways_GO) != "IM.HIVpreDM.HIVnonDM"] # No enriched terms
pathways_GO <- lapply(pathways_GO, GSEA_order) # Order by NES
saveRDS(pathways_GO, file = paste(path, "GO_pathways.rds", sep = "/")) # Save so can skip next time

# Kegg GSEA
pathways_Kegg <- lapply(file_list, KEGG_GSEA) # KEGG Pathway analysis
pathways_Kegg <- pathways_Kegg[names(pathways_Kegg) != "EndoMT.HIVpreDM.HIVnonDM"] # No enriched terms
pathways_Kegg <- lapply(pathways_Kegg, GSEA_order)
saveRDS(pathways_Kegg, file = paste(path, "Kegg_pathways.rds", sep = "/"))

###-------------------####
# Figure 6A: Macrophage Prediabetic vs non-diabetic Volcano Plot
###-------------------####
Mac_volcano <- VOLCANO(resOrderedDF = file_list$Macrophage.HIVpreDM.HIVnonDM, filename = paste(fig_dir, "Figure6A.png", sep = "/"), genes = c("HAMP", "FABP5", "C3", "DBI", "ETFB", "NDUFS5",
"CCL8", "CXCL1", "CCL2", "CCL4", "CCL3", "TNF", "TRIB1", "EGR1", "MRC1", "KLF6", "LYVE1", "ID3", "CAPG", "LALS1", "CSTB"), ylim = c(0,5))

###-------------------####
# Figure 6B: Macrophage Prediabetic vs diabetic GSEA
###-------------------####
Mac_gsea_plot <- GSEA_plot(df = pathways_Kegg$Macrophage.HIVpreDM.HIVnonDM, filename = paste(fig_dir, "Figure6B.png", sep = "/"), num_path = 5, width = 10)

###-------------------####
# Figure 6C: Pseudotime Mo-Mac to PVM
###-------------------####
# see Macrophage_Pseudotime.R

###-------------------####
# Figure 6D: Pseudotime Mo-Mac to PVM
###-------------------####
# see Macrophage_Pseudotime.R

###-------------------####
# Figure 6E: CD4 Volcano Plot
###-------------------####
CD4_volcano <- VOLCANO(resOrderedDF = file_list$CD4.HIVpreDM.HIVnonDM, filename = paste(fig_dir, "Figure6E.png", sep = "/"), genes = c("BATF", "CTSD", "S100A6", "PFN1", "TCF7", "CD63", "IL7R",
"IFITM3", "SELL", "LEF1", "ACTB", "MT2A", "IL32", "CSTB", "SH3BGRL3", "TCF7"), ylim = c(0,3.5), xlim = c(-1.5,1.5))

###-------------------####
# Figure 6F: CD8 Volcano Plot
###-------------------####
CD8_volcano <- VOLCANO(resOrderedDF = file_list$CD8.HIVpreDM.HIVnonDM, filename = paste(fig_dir, "Figure6F.png", sep = "/"), genes = c("CD7", "LDHB", "KLF2", "ALOX5AP", "CD63", "GZMA", 
"IL7R", "CLEC2B"), ylim = c(0,3.5), xlim = c(-1.5,1.5))
