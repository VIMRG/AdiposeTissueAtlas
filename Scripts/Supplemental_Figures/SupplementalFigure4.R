###############################################################################
######-------------------SUPPLEMENTAL FIGURE 4----------------------------#####
###############################################################################

###-------------------####
# Import Libraries
###-------------------####
library(ggbeeswarm)
library(miloR)
library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(scProportionTest)
library(ggtext)
library(future)

plan("multiprocess", workers = 5)
options(future.globals.maxSize = 3000 * 1024^2)
options(future.seed = TRUE)

###-------------------####
# Source Code
###-------------------####
source('../Utils.R')

###-------------------####
# Directories
###-------------------####
tmp_dir <- "../SubsetAnalysis"
Prop_dir <- "../Proportion"
fig_dir <- "../Figures"

###-------------------####
# Load Seurat Object
###-------------------####
Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Lymphoid <- readRDS(paste0(tmp_dir, "/", "Lymphoid.rds"))
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage.rds"))
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8.rds"))

###-------------------####
# Load Cell Proportions
###-------------------####
Lymphoid.Prop <- read.table(paste(Prop_dir, "Lymphoid.Prop.txt", sep = "/"), check.names = FALSE)
Myeloid.Prop <- read.table(paste(Prop_dir, "Myeloid.Prop.txt", sep = "/"), check.names = F)
Stromal.Prop <- read.table(paste(Prop_dir, "Stromal.Prop.txt", sep = "/"), check.names = F)
Vascular.Prop <- read.table(paste(Prop_dir, "Vascular.Prop.txt", sep = "/"), check.names = F)
CD4.Prop <- read.table(paste(Prop_dir, "CD4.Prop.txt", sep = "/"), check.names = F)
CD8.Prop <- read.table(paste(Prop_dir, "CD8.Prop.txt", sep = "/"), check.names = F)
Macrophage.Prop <- read.table(paste(Prop_dir, "Macrophage.Prop.txt", sep = "/"), check.names = F)
CD4.Num <- read.table(paste(Prop_dir, "CD4.Num.txt", sep = "/"), check.names = F)
CD8.Num <- read.table(paste(Prop_dir, "CD8.Num.txt", sep = "/"), check.names = F)
Macrophage.Num <- read.table(paste(Prop_dir, "Macrophage.Num.txt", sep = "/"), check.names = F)

# Remove NAs
Lymphoid.Prop <- Lymphoid.Prop %>% na.omit()
Myeloid.Prop <- Myeloid.Prop %>% na.omit()
Stromal.Prop <- Stromal.Prop %>% na.omit()
Vascular.Prop <- Vascular.Prop %>% na.omit()
CD4.Prop <- CD4.Prop %>% na.omit()
CD8.Prop <- CD8.Prop %>% na.omit()
Macrophage.Prop <- Macrophage.Prop %>% na.omit()

###-------------------####
# Merge with Metadata
###-------------------####
# Main Subsets
Myeloid.prop <- Prop_Merge(Myeloid.Prop, data_hatim)
Lymphoid.prop <- Prop_Merge(Lymphoid.Prop, data_hatim)
Stromal.prop <- Prop_Merge(Stromal.Prop, data_hatim)
Vascular.prop <- Prop_Merge(Vascular.Prop, data_hatim)

# Macrophage Subset
Macrophage.prop <- Prop_Merge(Macrophage.Prop, data_hatim)
colnames(Macrophage.Num) <- c("HATIMID", "Macnum")
Macrophage.Num$HATIMID <- as.factor(as.character(Macrophage.Num$HATIMID))

Macrophage.prop <- left_join(Macrophage.prop, Macrophage.Num, by = "HATIMID")
Macrophage.prop <- Macrophage.prop %>% dplyr::filter(Macnum > 30) # Filter > 30
Macrophage.prop <- Macrophage.prop %>% dplyr::select(-c(Macnum))

# CD4 Subset
CD4.prop <- Prop_Merge(CD4.Prop, data_hatim)
colnames(CD4.Num) <- c("HATIMID", "CD4num")
CD4.Num$HATIMID <- as.factor(as.character(CD4.Num$HATIMID))

CD4.prop <- left_join(CD4.prop, CD4.Num, by = "HATIMID")
CD4.prop <- CD4.prop %>% dplyr::filter(CD4num > 30) # Filter > 30
CD4.prop <- CD4.prop %>% dplyr::select(-c(CD4num))

# CD8 Subset
CD8.prop <- Prop_Merge(CD8.Prop, data_hatim)
colnames(CD8.Num) <- c("HATIMID", "CD8num")
CD8.Num$HATIMID <- as.factor(as.character(CD8.Num$HATIMID))

CD8.prop <- left_join(CD8.prop, CD8.Num, by = "HATIMID")
CD8.prop <- CD8.prop %>% dplyr::filter(CD8num > 30) # Filter > 30
CD8.prop <- CD8.prop %>% dplyr::select(-c(CD8num))

###-------------------####
# Supplemental Figure 4A: Differential Abundance Myeloid Populations
###-------------------####
# Convert to SCE then Milo
Myeloid.sce <- as.SingleCellExperiment(Myeloid, assay = "RNA")
Myeloid.milo <- Milo(Myeloid.sce)

# Build KNN Graph
Myeloid.milo <- buildGraph(Myeloid.milo, k = 35, d = 25, reduced.dim = "HARMONY")

# Define Neighbourhoods
Myeloid.milo <- makeNhoods(Myeloid.milo, prop = 0.1, k = 35, d=25, refined = TRUE, reduced_dims = "HARMONY")
plotNhoodSizeHist(Myeloid.milo)

meta.data = as.data.frame(colData(Myeloid.milo))
meta.data = droplevels(meta.data)
Myeloid.milo <- countCells(Myeloid.milo, meta.data = meta.data, sample="HATIMID")
dim(nhoodCounts(Myeloid.milo))

# Prepare Metadata and Study Design
Myeloid.design <- data.frame(colData(Myeloid.milo))[,c("HATIMID", "Lane", "BMI", "StudyGroup")]
Myeloid.design <- distinct(Myeloid.design)
Myeloid.design <- droplevels(Myeloid.design)
rownames(Myeloid.design) <- Myeloid.design$HATIMID
Myeloid.design <- Myeloid.design[colnames(nhoodCounts(Myeloid.milo)), ]
Myeloid.design$Lane <- as.factor(Myeloid.design$Lane)
Myeloid.design <- Myeloid.design %>% mutate(StudyGroup = case_when(StudyGroup == "HIV+ non-diabetic"~"HIVnonDM",
    StudyGroup == "HIV+ prediabetic"~"HIVpreDM",
    TRUE~"HIVDM"))
table(Myeloid.design$StudyGroup)
rownames(Myeloid.design) <- Myeloid.design$HATIMID

# Set up Contrasts
model <- model.matrix(~ 0 + StudyGroup, data=Myeloid.design)
ave.contrast <- c("StudyGroupHIVnonDM - (StudyGroupHIVDM + StudyGroupHIVpreDM)/2")
mod.constrast <- makeContrasts(contrasts=ave.contrast, levels=model)
mod.constrast

Myeloid.milo <- calcNhoodDistance(Myeloid.milo, d=25, reduced.dim = "HARMONY")
da_results <- testNhoods(Myeloid.milo, design = ~ 0 + StudyGroup, design.df = Myeloid.design, model.contrasts = ave.contrast, fdr.weighting = "graph-overlap") # max, neighbor-distance, graph-overlap, or k-distance

# Plot Results
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

# Inspect distribution
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

# Volcano Plot
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)

# Plot on UMAP
Myeloid.milo <- buildNhoodGraph(Myeloid.milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(Myeloid.milo, dimred = "UMAP", colour_by="StudyGroup", text_by = "Annotation", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(Myeloid.milo, da_results, layout="UMAP",alpha=0.1) 
  
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
ggsave(filename = paste(fig_dir, "SF4B"), device = "png", height = 5, width = 10, units = 'in', dpi = 300)

# DA Annotation
da_results <- annotateNhoods(Myeloid.milo, da_results, coldata_col = "Annotation")

saveRDS(Myeloid.milo, file = paste(tmp_dir, "Myeloid.milo.R", sep = "/"))

plotDAbeeswarm(da_results, group.by = "Annotation")
ggsave(filename = paste(fig_dir, "SF4A.png", sep = "/"), device = "png", height = 10, width = 10, units = 'in', dpi = 300)

###-------------------####
# Supplemental Figure 4B: Scatter Plot IM and FBG
###-------------------####
plot <- scatter_plot(DF = Macrophage.prop, x = "IM", y = "meta_fbg", xlabel = "IM (% Macrophage)", ylabel = "Fasting Blood Glucose (mg/dl)", nonDM = T, dotsize = 4)
ggsave(filename = paste(fig_dir, "SF4B.png", sep = "/"), units = "in", height = 6, width = 6, dpi = 300, device = "png")

###-------------------####
# Supplemental Figure 4C: Scatter Plot LAM and HbA1c
###-------------------####
plot <- scatter_plot(DF = Macrophage.prop, x = "LAM", y = "meta_hba1c", xlabel = "LAM (% Macrophage)", ylabel = "HbA1c (%)", nonDM = T)
ggsave(filename = paste(fig_dir, "SF4C.png", sep = "/"), units = "in", height = 6, width = 6, dpi = 300, device = "png")

###-------------------####
# Supplemental Figure 4D: Scatter Plot CD8 TEM and FBG
###-------------------####
plot <- scatter_plot(DF = CD8.prop, x = "`CD8 TEM`", y = "meta_fbg", xlabel = "CD8 TEM (% CD8)", ylabel = "Fasting Blood Glucose (mg/dl)", nonDM = T)
ggsave(filename = paste(fig_dir, "SF4D.png", sep = "/"), units = "in", height = 6, width = 6, dpi = 300, device = "png")

###-------------------####
# Supplemental Figure 4E: CD4+CD69+ Cells by Group
###-------------------####
AT_CD4 <- data_hatim %>% dplyr::select(at_cd4_cd69_mem1, StudyGroup, HATIMID) %>% dplyr::filter(!is.na(at_cd4_cd69_mem1) & StudyGroup != "HIV- diabetic")
colnames(AT_CD4) <- c("CD4+CD69+", "StudyGroup", "HATIMID")
Prop_long <- AT_CD4 %>% pivot_longer(cols = -c(HATIMID, StudyGroup), names_to = "measure", values_to = "value")

plot <- box_plot_fun_paper(Prop_long, parent_cell_type = "CD4 Memory", 
connect = F, is_paired = FALSE, x_axis_text_size = 12, my_comparisons = list(c("HIV+ non-diabetic", "HIV+ prediabetic"), 
c("HIV+ non-diabetic", "HIV+ diabetic")), specify_x_order = NA, ncol = 1)
ggsave(filename = paste(fig_dir, "SF4E.png", sep = "/"), units = "in", height = 5, width = 8, dpi = 300, device = "png")

###-------------------####
# Supplemental Figure 4F: Scatter Plot CD4 TEM and FBG
###-------------------####
plot <- scatter_plot(DF = CD4.prop, x = "`CD4 TEM`", y = "meta_fbg", xlabel = "CD4 TEM (% CD4)", ylabel = "Fasting Blood Glucose (mg/dl)", nonDM = T)
ggsave(filename = paste(fig_dir, "SF4F.png", sep = "/"), units = "in", height = 6, width = 6, dpi = 300, device = "png")

###-------------------####
# Supplemental Figure 4G: Differential Abundance Lymphoid
###-------------------####
# Convert to SCE then Milo
Lymphoid.sce <- as.SingleCellExperiment(Lymphoid, assay = "RNA")
Lymphoid.milo <- Milo(Lymphoid.sce)

# Build KNN Graph
Lymphoid.milo <- buildGraph(Lymphoid.milo, k = 40, d = 20, reduced.dim = "HARMONY")

# Define Neighbourhoods
Lymphoid.milo <- makeNhoods(Lymphoid.milo, prop = 0.1, k = 40, d=20, refined = TRUE, reduced_dims = "HARMONY")
plotNhoodSizeHist(Lymphoid.milo)

meta.data = as.data.frame(colData(Lymphoid.milo))
meta.data = droplevels(meta.data)
Lymphoid.milo <- countCells(Lymphoid.milo, meta.data = meta.data, sample="HATIMID")
dim(nhoodCounts(Lymphoid.milo))

# Prepare Metadata and Study Design
Lymphoid.design <- data.frame(colData(Lymphoid.milo))[,c("HATIMID", "Lane", "BMI", "StudyGroup")]
Lymphoid.design <- distinct(Lymphoid.design)
Lymphoid.design <- droplevels(Lymphoid.design)
rownames(Lymphoid.design) <- Lymphoid.design$HATIMID
Lymphoid.design <- Lymphoid.design[colnames(nhoodCounts(Lymphoid.milo)), ]
Lymphoid.design$Lane <- as.factor(Lymphoid.design$Lane)
Lymphoid.design <- Lymphoid.design %>% mutate(StudyGroup = case_when(StudyGroup == "HIV+ non-diabetic"~"HIVnonDM",
    StudyGroup == "HIV+ prediabetic"~"HIVpreDM",
    TRUE~"HIVDM"))
table(Lymphoid.design$StudyGroup)
rownames(Lymphoid.design) <- Lymphoid.design$HATIMID

# Set up Contrasts
model <- model.matrix(~ 0 + StudyGroup, data=Lymphoid.design)
ave.contrast <- c("StudyGroupHIVnonDM - (StudyGroupHIVDM + StudyGroupHIVpreDM)/2")
mod.constrast <- makeContrasts(contrasts=ave.contrast, levels=model)
mod.constrast

Lymphoid.milo <- calcNhoodDistance(Lymphoid.milo, d=20, reduced.dim = "HARMONY")
#Myeloid.milo$StudyGroup <- ifelse(Myeloid.milo$StudyGroup == "HIV+ non-diabetic", "HIVnonDM", ifelse(Myeloid.milo$StudyGroup == "HIV+ prediabetic", "HIVpreDM", "HIVDM"))
da_results <- testNhoods(Lymphoid.milo, design = ~ 0 + StudyGroup, design.df = Lymphoid.design, model.contrasts = ave.contrast, fdr.weighting = "graph-overlap") # max, neighbor-distance, graph-overlap, or k-distance

# Plot on UMAP
Lymphoid.milo <- buildNhoodGraph(Lymphoid.milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(Lymphoid.milo, dimred = "UMAP", colour_by="StudyGroup", text_by = "Annotation", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(Lymphoid.milo, da_results, layout="UMAP",alpha=0.1) 
  
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

# DA Annotation
da_results <- annotateNhoods(Lymphoid.milo, da_results, coldata_col = "Annotation")
plotDAbeeswarm(da_results, group.by = "Annotation")

ggsave(filename = paste(fig_dir, "SF4G.png", sep = "/"), device = "png", height = 10, width = 10, units = 'in', dpi = 300)



