###############################################################################
######--------------------SUPPLEMENTAL FIGURE 4---------------------------#####
###############################################################################
# 1. SET UP___________________________________________________________________
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
library(pheatmap)
library(clusterProfiler)
library(rWikiPathways)
library(org.Hs.eg.db)
library(DOSE)
library(ggrepel)

###-------------------####
# Source Code
###-------------------####
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Figure_Functions.R')
source('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/DGE/DGE_Function.R')

###-------------------####
# Import Data
###-------------------####
Stromal <- readRDS('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_1.3.rds')
Prop_Long <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Proportions_Long.csv') # Proportions Cell Types Long
Prop <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Merged_Proportions.csv')  # Proportions Cell Types Wide

Stromal.markers <- read.csv('/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Stromal_cluster_RNAmarkers.csv')

###-------------------####
# Environment
###-------------------####
set.seed(7612)

###-------------------####
# Set ggplot theme
###-------------------####
theme_set(theme_bw() +
            theme(axis.text = (element_text(size = 6)),
                  plot.caption = element_text(face = "italic", size = 9),
                  legend.position = "bottom",
                  legend.background = element_rect(fill = "transparent"), 
                  legend.text = element_text(size = 18),
                  panel.background = element_rect(fill = "transparent"),
                  panel.border = element_blank(),
                  axis.line = element_line(color="black"),
                  strip.text = element_text(size = 24), 
                  axis.title = element_text(size = 18),
                  plot.margin=unit(c(0.7,2.5,0,0),"cm"),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)
          ))

# 2. SUPPLEMENTAL FIGURE 4 ___________________________________________________________________
#########################
# Figure S4A: Stromal Heatmap
#########################
HEATMAP_FUN(Stromal, markers = Stromal.markers, path = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/SF4A.png', height = 16, Margins = c(0,6.5,0,0), plot_unit = 'cm')

#########################
# Figure S4B: VlnPlot Progenitor Markers
#########################
Idents(Stromal) <- "Manual_Annotation"
VlnPlot(Stromal, features = c("CFD", "LUM", "DCN", "MGP", "PI16", "GSN", "FBLN1"), pt.size = 0, stack = TRUE, flip = TRUE, fill.by = 'ident', col = col_shuffle) + 
NoLegend() + 
theme(axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold'), plot.margin = margin(c(0,0,0,2), unit = "cm")) +
xlab("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/SF4B.png", device = 'png', width = 20, height = 10, units = 'in', dpi = 600)

#########################
# Figure S4C: VlnPlot Preadipocyte Markers
#########################
VlnPlot(Stromal, features = c("CXCL14", "COL4A1", "ZFP36", "MYC", "JUN", "FOS", "CEBPD", "PPARG", "FABP4", "C7", "LPL", "CIDEC"),  pt.size = 0, stack = TRUE, flip = TRUE, fill.by = 'ident', col = col_shuffle) +
NoLegend() + 
theme(axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold'), plot.margin = margin(c(0,0,0,2), unit = "cm")) +
xlab("")
ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/SF4C.png', device = 'png', width = 20, height = 10, units = 'in', dpi = 600)

#########################
# Figure S4D: GO PATHWAY ANALYSIS
#########################
Preadipocytes <- Stromal.markers %>% filter(cluster %in% c("Early Preadipocyte 1", "Early Preadipocyte 2", "Early Preadipocyte 3", "Early Preadipocyte 4", "Preadipocyte", "CD9hi Preadipocyte", "ECM-Producing Preadipocyte",
"PTGDS+ Preadipocyte", "PTGDS+ ECM-Producing Preadipocyte", "ISG+ Preadipocyte", "Metallothionein+ Preadipocyte/ASC", "Mature Preadipocyte 1", "Mature Preadipocyte 2"))

Entrez <- clusterProfiler::bitr(Preadipocytes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Convert to ENTREZID
Preadipocytes <- left_join(Entrez, Preadipocytes, by = c("SYMBOL" = "gene")) # Add ENTREZID
Preadipocytes <- Preadipocytes %>% dplyr::filter(p_val_adj < 0.05) # Filter to genes < 0.05

Preadipocyte_list <- Preadipocytes %>% group_by(cluster) %>% group_split() # Make List
names(Preadipocyte_list) <- c("CD9hi Preadipocyte", "Early Preadipocyte 1", "Early Preadipocyte 2", "Early Preadipocyte 3", "Early Preadipocyte 4", "ECM-Producing Preadipocyte",
"ISG+ Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2", "Metallothionein+ Preadipocyte/ASC", "Preadipocyte", "PTGDS+ ECM-Producing Preadipocyte", "PTGDS+ Preadipocyte")

# Loop through list of DGE for each Preadipocyte cluster
DF <- data.frame()
DF <- lapply(Preadipocyte_list, function(x) {
    bp <- enrichGO(gene = x$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
    results <- data.frame(bp2)
    DF <- rbind(DF, results)
} )

# Plot Results
Pathways <- bind_rows(DF, .id = "Cluster") # Combine in data frame
Pathways$Cluster <- factor(Pathways$Cluster, levels = c("Early Preadipocyte 1", "Early Preadipocyte 2", "Early Preadipocyte 3", "Early Preadipocyte 4", "Preadipocyte", "CD9hi Preadipocyte",
"ECM-Producing Preadipocyte", "PTGDS+ Preadipocyte", "PTGDS+ ECM-Producing Preadipocyte", "Metallothionein+ Preadipocyte/ASC", "ISG+ Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2"))

Pathways %>% group_by(Cluster) %>%
        do(head(.,5)) %>% 
        mutate(Description = factor(Description, levels = Description)) %>%
        ggplot() + geom_point(aes(x = Cluster, y = Description, size = -log10(pvalue), fill = Count), shape = 21) +
        xlab("") + ylab("") + ggtitle("") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1, size = 14, face = "bold"), 
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 16, face = "bold")) + scale_y_discrete(labels = function(x) str_wrap(x, width = 70))

ggsave('/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/SF4D.png', device = 'png', width = 20, height = 20, units = 'in', dpi = 600)

#########################
# Figure S4E: Preadipocyte DGE HIV+ DM vs HIV+ non-DM
#########################
Idents(Stromal) <- "Manual_Annotation"
Preadipocytes <- subset(Stromal, idents = c("Early Preadipocyte 1", "Early Preadipocyte 2", "Early Preadipocyte 3", "Preadipocyte", "CD9hi Preadipocyte", "ECM-Producing Preadipocyte",
"PTGDS+ Preadipocyte", "PTGDS+ ECM-Producing Preadipocyte", "Mature Preadipocyte 1", "Mature Preadipocyte 2"))
Preadipocytes[['Manual_Annotation']] <- "Preadipocyte" # aggregate Preadipocytes
Preadipocytes <- CellType_Fun(Preadipocytes) # From DGE functions

Preadipocytes.DGE.HIVDMvHIVnonDM <- HIV_comparison.5(Preadipocytes, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Preadipocytes.HIVDM.HIVnonDM_1.3.csv")
volcano.HIVDMvHIVnonDM.plot <- VOLCANO(Preadipocytes.DGE.HIVDMvHIVnonDM, filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/SF4E.png', 
genes = c("GSN", "CD36", "POSTN", "KLF4", "TMEM176A", "EGR1", "MYOC", "ALPL", "DEPP1", "IGF2", "CEBPD", "MAFF", "IER2", "KLF2", "CTHRC1", "ZFP36", "TXNIP"), filter = "Preadipocyte", width = 5, height = 5)

#########################
# Figure S4F: Preadipocyte DGE HIV+ DM vs HIV- DM
#########################
Preadipocytes.DGE.HIVDMvHIVnegDM <- DM_comparison.5(Preadipocytes, "/data/p_koethe_lab/Atlas_AT/Analysis/Harmony_Integration/Stromal/Preadipocytes.HIVDM.HIVnegDM_1.3.csv")
volcano.HIVDMvHIVnegDM.plot <- VOLCANO(Preadipocytes.DGE.HIVDMvHIVnegDM, filename = '/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/SF4F.png', 
genes = c("DNAJB1", "HSPA8", "JUN", "ATF3", "CDKN1A", "JUNB", "EGR1", "FOS", "NFKBIZ", "ADIRF", "COL3A1", "COL1A2", "FABP5", "CCL2", "DEPP1"), filter = "Preadipocyte", width = 5, height = 5)

#########################
# Figure S4G: Comparison to Other Datasets (This code was adapted from Sarvari et al. 2021)
#########################
Sarvari_FAP1 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Sarvari_FAP1.csv")
Sarvari_FAP2 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Sarvari_FAP2.csv")
Sarvari_FAP3 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Savari_FAP3.csv")
Sarvari_FAP4 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Savari_FAP4.csv")

Helper_CAP <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Helper_CAP.csv")
Helper_FIP <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Helper_FIP.csv")
Helper_APC <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Helper_APC.csv")

Deplancke_P1 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Deplanke_P1.csv")
Deplancke_P2 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Deplanke_P2.csv")
Deplancke_P3 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Deplanke_P3.csv")


Burle_ASC1 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Burle_ASC1.csv")
Burle_ASC2 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Burle_ASC2.csv")

Seale_CD142 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Seale_CD142.csv")
Seale_ICAM <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Seale_ICAM.csv")
Seale_DPP4 <- read.csv("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Comparison/Seale_DPP4.csv")

GeneList <- list(as.character(Savari_FAP1[,5]), as.character(Savari_FAP2[,5]), as.character(Savari_FAP3[,5]), as.character(Savari_FAP4[,5]), as.character(Helper_CAP[,5]), as.character(Helper_FIP[,5]),
as.character(Helper_APC[,5]), as.character(Deplancke_P1[,5]), as.character(Deplancke_P2[,5]), as.character(Deplancke_P3[,5]), as.character(Burle_ASC1[,5]), as.character(Burle_ASC2[,5]), as.character(Seale_CD142[,5]),
as.character(Seale_ICAM[,5]), as.character(Seale_DPP4[,5]))

GeneList <- lapply(GeneList, function(x) {
    x[x != "N/A" & x != ""]
})

names(GeneList) <- c("Sarvari_FAP1", "Sarvari_FAP2", "Sarvari_FAP3", "Sarvari_FAP4", "Helper_CAP", "Helper_FIP", "Helper_APC", "Deplancke_P1", "Deplancke_P2", "Deplancke_P3", "Burle_ASC1", "Burle_ASC2",
"Seale_CD142", "Seale_ICAM", "Seale_DPP4")

Stromal <- AddModuleScore(Stromal, features = GeneList, name = "Comparison")

colnames(Stromal@meta.data)[95:109] <- names(GeneList)

# Scale the module scores
for (i in 95:109) {Stromal@meta.data[,i] <- scale(Stromal@meta.data[,i])}

# Sarvari
Stromal$Sarvari <- "None"
Stromal@meta.data[Stromal@meta.data$Sarvari_FAP1 > Stromal@meta.data$Sarvari_FAP2 & Stromal@meta.data$Sarvari_FAP1 > Stromal@meta.data$Sarvari_FAP3 & Stromal@meta.data$Sarvari_FAP1 > Stromal@meta.data$Sarvari_FAP4, "Sarvari"] <- "FAP1"
Stromal@meta.data[Stromal@meta.data$Sarvari_FAP2 > Stromal@meta.data$Sarvari_FAP1 & Stromal@meta.data$Sarvari_FAP2 > Stromal@meta.data$Sarvari_FAP3 & Stromal@meta.data$Sarvari_FAP2 > Stromal@meta.data$Sarvari_FAP4, "Sarvari"] <- "FAP2"
Stromal@meta.data[Stromal@meta.data$Sarvari_FAP3 > Stromal@meta.data$Sarvari_FAP1 & Stromal@meta.data$Sarvari_FAP3 > Stromal@meta.data$Sarvari_FAP2 & Stromal@meta.data$Sarvari_FAP3 > Stromal@meta.data$Sarvari_FAP4, "Sarvari"] <- "FAP3"
Stromal@meta.data[Stromal@meta.data$Sarvari_FAP4 > Stromal@meta.data$Sarvari_FAP1 & Stromal@meta.data$Sarvari_FAP4 > Stromal@meta.data$Sarvari_FAP2 & Stromal@meta.data$Sarvari_FAP4 > Stromal@meta.data$Sarvari_FAP3, "Sarvari"] <- "FAP4"

# Helper
Stromal$Helper <- "None"
Stromal@meta.data[Stromal@meta.data$Helper_CAP > Stromal@meta.data$Helper_FIP & Stromal@meta.data$Helper_CAP > Stromal@meta.data$Helper_APC, "Schwalie"] <- "CAP"
Stromal@meta.data[Stromal@meta.data$Helper_FIP > Stromal@meta.data$Helper_CAP & Stromal@meta.data$Helper_FIP > Stromal@meta.data$Helper_APC, "Schwalie"] <- "FIP"
Stromal@meta.data[Stromal@meta.data$Helper_APC > Stromal@meta.data$Helper_FIP & Stromal@meta.data$Helper_APC > Stromal@meta.data$Helper_CAP, "Schwalie"] <- "APC"

# Deplanke
Stromal$Deplancke <- "None"
Stromal@meta.data[Stromal@meta.data$Deplancke_P1 > Stromal@meta.data$Deplancke_P2 & Stromal@meta.data$Deplancke_P1 > Stromal@meta.data$Deplancke_P3, "Deplanke"] <- "P1"
Stromal@meta.data[Stromal@meta.data$Deplancke_P2 > Stromal@meta.data$Deplancke_P1 & Stromal@meta.data$Deplancke_P2 > Stromal@meta.data$Deplancke_P3, "Deplanke"] <- "P2"
Stromal@meta.data[Stromal@meta.data$Deplancke_P3 > Stromal@meta.data$Deplancke_P1 & Stromal@meta.data$Deplancke_P3 > Stromal@meta.data$Deplancke_P2, "Deplanke"] <- "P3"

# Burle
Stromal$Burle <- "None"
Stromal@meta.data[Stromal@meta.data$Burle_ASC1 > Stromal@meta.data$Burle_ASC2, "Burle"] <- "ASC1"
Stromal@meta.data[Stromal@meta.data$Burle_ASC2 > Stromal@meta.data$Burle_ASC1, "Burle"] <- "ASC2"

# Seale
Stromal$Seale <- "None"
Stromal@meta.data[Stromal@meta.data$Seale_CD142 > Stromal@meta.data$Seale_ICAM & Stromal@meta.data$Seale_CD142 > Stromal@meta.data$Seale_DPP4, "Seale"] <- "CD142"
Stromal@meta.data[Stromal@meta.data$Seale_ICAM > Stromal@meta.data$Seale_CD142 & Stromal@meta.data$Seale_ICAM > Stromal@meta.data$Seale_DPP4, "Seale"] <- "ICAM"
Stromal@meta.data[Stromal@meta.data$Seale_DPP4 > Stromal@meta.data$Seale_ICAM & Stromal@meta.data$Seale_DPP4 > Stromal@meta.data$Seale_CD142, "Seale"] <- "DPP4"

## DimPlots
Stromal <- SetIdent(Stromal, value = Stromal$Sarvari)
Stromal[['Sarvari']] <- factor(Stromal$Sarvari, levels = c("FAP1", "FAP2", "FAP3", "FAP4"))
color <- col_shuffle[c(1,25,8,29)]
DimPlot(Stromal, label=F, cols = color) + #scale_colour_manual(name = "Cell Type", breaks = c(1:4), labels = c("FAP1", "FAP2", "FAP3", "FAP4"), values = color) +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) + ggtitle("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Sarvari_Plot.png", dpi = 600, height = 5, width = 5, units = 'in', device = 'png')

Stromal <- SetIdent(Stromal, value = Stromal$Helper)
Stromal[['Helper']] <- factor(Stromal$Helper, levels = c("CAP", "FIP", "APC"))
color <- col_shuffle[c(1,25,8,29)]
DimPlot(Stromal, label=F, cols = color) + #scale_colour_manual(name = "Cell Type", breaks = c(1:4), labels = c("FAP1", "FAP2", "FAP3", "FAP4"), values = color) +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) + ggtitle("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Helper_Plot.png", dpi = 600, height = 5, width = 5, units = 'in', device = 'png')

Stromal <- SetIdent(Stromal, value = Stromal$Deplancke)
Stromal[['Deplancke']] <- factor(Stromal$Deplancke, levels = c("P1", "P2", "P3"))
color <- col_shuffle[c(1,25,8,29)]
DimPlot(Stromal, label=F, cols = color) + #scale_colour_manual(name = "Cell Type", breaks = c(1:4), labels = c("FAP1", "FAP2", "FAP3", "FAP4"), values = color) +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) + ggtitle("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Deplancke_Plot.png", dpi = 600, height = 5, width = 5, units = 'in', device = 'png')

Stromal <- SetIdent(Stromal, value = Stromal$Burle)
Stromal[['Burle']] <- factor(Stromal$Burle, levels = c("ASC1", "ASC2"))
color <- col_shuffle[c(1,25,8,29)]
DimPlot(Stromal, label=F, cols = color) + #scale_colour_manual(name = "Cell Type", breaks = c(1:4), labels = c("FAP1", "FAP2", "FAP3", "FAP4"), values = color) +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) + ggtitle("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Burle_Plot.png", dpi = 600, height = 5, width = 5, units = 'in', device = 'png')

Stromal <- SetIdent(Stromal, value = Stromal$Seale)
Stromal[['Seale']] <- factor(Stromal$Seale, levels = c("CD142", "ICAM", "DPP4"))
color <- col_shuffle[c(1,25,8,29)]
DimPlot(Stromal, label=F, cols = color) + #scale_colour_manual(name = "Cell Type", breaks = c(1:4), labels = c("FAP1", "FAP2", "FAP3", "FAP4"), values = color) +
theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14)) + ggtitle("")
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Seale_Plot.png", dpi = 600, height = 5, width = 5, units = 'in', device = 'png')

#########################
# Figure S4H: Feature Plots
#########################
# Sarvari
p1 <- FeaturePlot(Stromal, features = "Sarvari_FAP1", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p2 <- FeaturePlot(Stromal, features = "Sarvari_FAP2", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p3 <- FeaturePlot(Stromal, features = "Sarvari_FAP3", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
          
p4 <- FeaturePlot(Stromal, features = "Sarvari_FAP4", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

plot_grid(p1, p2, p3, p4)
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Sarvari_FeaturePlot.png", dpi = 600, height = 10, width = 10, units = 'in', device = 'png')

# Helper
p1 <- FeaturePlot(Stromal, features = "Helper_CAP", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p2 <- FeaturePlot(Stromal, features = "Helper_FIP", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p3 <- FeaturePlot(Stromal, features = "Helper_APC", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
          
plot_grid(p1, p2, p3)
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Helper_FeaturePlot.png", dpi = 600, height = 10, width = 10, units = 'in', device = 'png')

# Deplancke
p1 <- FeaturePlot(Stromal, features = "Deplancke_P1", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p2 <- FeaturePlot(Stromal, features = "Deplancke_P2", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p3 <- FeaturePlot(Stromal, features = "Deplancke_P3", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
          
plot_grid(p1, p2, p3)
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Deplancke_FeaturePlot.png", dpi = 600, height = 10, width = 10, units = 'in', device = 'png')

# Burle
p1 <- FeaturePlot(Stromal, features = "Burle_ASC1", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p2 <- FeaturePlot(Stromal, features = "Burle_ASC2", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
          
plot_grid(p1,p2)
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Burle_FeaturePlot.png", dpi = 600, height = 10, width = 10, units = 'in', device = 'png')

# Seale
p1 <- FeaturePlot(Stromal, features = "Seale_CD142", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p2 <- FeaturePlot(Stromal, features = "Seale_ICAM", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

p3 <- FeaturePlot(Stromal, features = "Seale_DPP4", min.cutoff = 'q10', max.cutoff = 'q90') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
          
plot_grid(p1, p2, p3)
ggsave("/data/p_koethe_lab/Atlas_AT/Analysis/Seurat_code/Harmony_Code/Supplementary_Figures/Supplemental_Figure4/Seale_FeaturePlot.png", dpi = 600, height = 10, width = 10, units = 'in', device = 'png')
