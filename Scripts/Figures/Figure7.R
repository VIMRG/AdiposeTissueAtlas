###############################################################################
######----------------------------FIGURE 7--------------------------------#####
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
source('../Utils.R')

###-------------------####
# Directories
###-------------------####
Subset_dir <- "../SubsetAnalysis"
AllCells_dir <- "../AllCells" 
fig_dir <- "../Figures"

###-------------------####
# Load Seurat Object
###-------------------####
Stromal <- readRDS(paste0(Subset_dir, "/", "Stromal.rds"))
Macrophage <- readRDS(paste0(Subset_dir, "/", "Macrophage.rds"))
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
# Figure 7A: Dialogue Score
###-------------------####
# See output from Dialogue.R

###-------------------####
# Figure 7B: ORA Dialogue
###-------------------####
Gene_Up <- c(R$MCPs$MCP1[[1]], R$MCPs$MCP1[[2]], R$MCPs$MCP1[[3]], R$MCPs$MCP1[[4]], R$MCPs$MCP1[[5]], R$MCPs$MCP1[[6]], R$MCPs$MCP1[[7]])
Gene_Down <- c(R$MCPs$MCP1[[8]], R$MCPs$MCP1[[9]], R$MCPs$MCP1[[10]], R$MCPs$MCP1[[11]], R$MCPs$MCP1[[12]], R$MCPs$MCP1[[13]], R$MCPs$MCP1[[14]])

Gene_Up <- Gene_Up[!duplicated(Gene_Up)] # remove duplicates
Gene_Up <- Gene_Up[!str_detect(Gene_Up, "RP")] # Remove RP genes since ribosomal features will dominate pathway
Gene_Down <- Gene_Down[!duplicated(Gene_Down)] # remove duplicates
Gene_Down <- Gene_Down[!str_detect(Gene_Down, "RP")] # Remove RP genes since ribosomal features will dominate pathway

# GO Analysis
enrich_Up <- enrichGO(gene = Gene_Up,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  minGSSize = 10,
  maxGSSize = 500,
  keyType = "SYMBOL")
enrich_Up <- simplify(enrich_Up, cutoff=0.7, by="p.adjust", select_fun=min)

enrich_Down <- enrichGO(gene = Gene_Down,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  minGSSize = 10,
  maxGSSize = 500,
  keyType = "SYMBOL")
enrich_Down <- simplify(enrich_Down, cutoff=0.7, by="p.adjust", select_fun=min)

# Pull out Top pathways
Up_path <- enrich_Up[1:12, c("Description", "p.adjust", "GeneRatio")] # rna processes
Up_path$Type <- "Up"
Down_path <- enrich_Down[1:12, c("Description", "p.adjust", "GeneRatio")] # extracellular organization
Down_path$Type <- "Down"

pathways <- rbind(Up_path, Down_path)
pathways$Count <- as.numeric(stri_extract_first_regex(pathways$GeneRatio, "[0-9]+"))

plot <- pathways %>% group_by(Type) %>%
        mutate(Description = factor(Description, levels = Description)) %>%
        ggplot() + geom_point(aes(x = Type, y = Description, size = Count, fill = -log10(p.adjust)), shape = 21) +
        xlab("") + ylab("") + ggtitle("") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1, size = 20, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), legend.title = element_text(size = 14, face = "bold"), legend.position = "top") + scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
    
plot
ggsave(file = paste(fig_dir, "Figure7B.png", sep = "/"),  device = 'png', width = 15, height = 10, units = 'in', dpi = 600)

###-------------------####
# Figure 7C: CD4 T cells
###-------------------####
plot <- Heatmap_fun(seurat_object = CD4, subset = F, metadata = data_hatim, cluster = "CD4", genes = c("CCR7", "LEF1", "IL4R", "IL21R", "LDHB", "SELL", "ALOX5AP", "ANXA2", "CD36", "CCL5", "CCL4",
"CD40LG", "CD69", "CD99", "CTSD", "CXCR3", "CXCR4", "DUSP2", "GBP5", "FABP4", "FABP5", "GZMA", "GZMK", "HOPX", "HLA-DRA", "IFITM2", "IFNG", "IL32", "LGALS3", "KLRG1", "PRF1", "S100A4", "SH3BGRL3",
"TNF", "XCL1", "XCL2"), tmp_dir = fig_dir, filename = "Figure7C.png", R = R, MCP_dim = 1, height = 5, width = 12)

###-------------------####
# Figure 7D: CD8 T cells
###-------------------####
plot <- Heatmap_fun(seurat_object = CD8, subset = F, metadata = data_hatim, cluster = "CD8", genes = c("CCR7", "IL4R", "IL7R", "LEF1", "NELL2", "NFKB1", "SELL", "ACTB", "ALOX5AP", "ANXA2",
"B2M", "CCL3", "CCL4", "CCL5", "CD69", "CTSD", "FABP4", "GZMB", "HLA-DRA", "IFNG", "HOPX", "IFI16", "IL32", "KLRG1", "PRF1", "S100A4", "SH3BGRL3", "TNF", "XCL1", "XCL2"),
filename = "Figure7D.png", tmp_dir = fig_dir, R = R, MCP_dim = 1, height = 5, width = 12)
