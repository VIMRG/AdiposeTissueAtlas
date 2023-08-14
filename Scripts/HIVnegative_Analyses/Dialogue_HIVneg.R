##############################################################################################
##------------------------- DIALOGUE ---------------------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will use Dialogue to discover gene regulatory programs across 
# cell types that vary with disease state. 
##############################################################################################

# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(future)
library(tidyverse)
library(harmony)
library(Matrix)
library(gridExtra)
library(DIALOGUE)
library(beanplot)
library(nnls)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)

plan("multiprocess", workers = 5)
options(future.globals.maxSize = 5000 * 1024^2)
options(future.seed = TRUE)

set.seed(7412)

date = "6.27"
tmp <- paste0("../", date)

data_hatim <- readr::read_csv('../HATIMStudyVisit_DATA_2021-08-09_1359.csv')
data_hatim <- data_hatim %>% mutate(HATIMID = as.factor(hatim_clin_visit_pid),
                                    StudyGroup = factor(hatim_final_arm, levels = c(1,2,3,4),
                                                        labels = c("HIV+ non-diabetic", "HIV+ prediabetic", 
                                                        "HIV+ diabetic", "HIV- diabetic")),
                                    HIV = case_when(hatim_final_arm == 4 ~ "HIVneg",
                                                    TRUE ~ "HIVpos"),
                                    Glucose = case_when(hatim_final_arm == 1 ~ "Normoglycemic",
                                                        TRUE ~ "Glucose Intolerant"),
                                    Sex = factor(meta_sex, levels = c(0,1), labels= c("Male", "Female")),
                                    homa2_ir = as.numeric(homa2_ir), 
                                    Age = meta_age,
                                    BMI = meta_bmi,
                                    metformin = case_when(stringr::str_detect(other_meds, regex('metformin', ignore_case = T))~ "Yes",stringr::str_detect(other_meds, regex('metphormin', ignore_case = T))~ "Yes",
                                    TRUE~"No"),
                                    FBG = meta_fbg,
                                    hba1c = meta_hba1c)

###############################
# Load Data
###############################
Macrophage <- readRDS(paste(tmp, "SubsetAnalysis_HIVneg", "Macrophage_HIVneg.rds", sep = "/"))
CD4 <- readRDS(paste(tmp, "SubsetAnalysis_HIVneg", "CD4_HIVneg.rds", sep = "/"))
CD8 <- readRDS(paste(tmp, "SubsetAnalysis_HIVneg", "CD8_HIVneg.rds", sep = "/"))
Stromal <- readRDS(paste(tmp, "SubsetAnalysis_HIVneg", "Stromal_HIVneg.rds", sep = "/"))

###############################
# Add Glucose Intolerance to Metadata
###############################
Macrophage[['GlucoseIntolerant']] <- ifelse(Macrophage$StudyGroup == "HIV+ non-diabetic", FALSE, TRUE)
CD4[['GlucoseIntolerant']] <- ifelse(CD4$StudyGroup == "HIV+ non-diabetic", FALSE, TRUE)
CD8[['GlucoseIntolerant']] <- ifelse(CD8$StudyGroup == "HIV+ non-diabetic", FALSE, TRUE)
Stromal[['GlucoseIntolerant']] <- ifelse(Stromal$StudyGroup == "HIV+ non-diabetic", FALSE, TRUE)

# 2. FUNCTIONS_____________________________________________________
# Dialogue requires cell-specific objects. Therefore, I will subset cells of interest from
# the seurat objects and extract the necessary information. To limit to many 
# interactions, I will subset on prespecified cells types that have either shown proportional
# differences or have biological relevance to T2DM. This will include the following:
# T cells: CD4, CD8
# Macrophage: IM, LAM, PVM
# Stromal: Myofibroblast, PCOLCE+ Fibroblast, MYOC+ Fibroblast, Preadipocyte
rA <- list()

###############################
# Scale and Log-Transform RNA Count (0-1)
###############################
Scale <- function(x) {
    (x-min(x))/(max(x)-min(x))
}

###############################
# Extract Data and make cell object
###############################
MakeCellObject <- function(subset, seurat_object) {
    print(paste0("Subsetting on ", subset, "."))
    Idents(seurat_object) <- "Annotation"
    DefaultAssay(seurat_object) <- "RNA"
    celltype <- subset(seurat_object, idents = subset)
    print(paste0("The object ", subset, " has ", dim(celltype)[1], " genes and ", dim(celltype)[2], " cells."))
    
    name <- gsub(" ", ".", subset) # name of cell object
    tpm <- celltype@assays$RNA@data # log normalized gene counts
    metadata <- data.frame(GlucoseIntolerant = celltype$GlucoseIntolerant, Sex = celltype$Sex, Age = celltype$Age, BMI = celltype$BMI, FBG = celltype$FBG, Hba1c = celltype$hba1c) # metadata
    samples <- as.character(celltype$HATIMID) # sample names
    
    print("Performing feature selection, data scaling, PCA.")
    celltype <- celltype %>% FindVariableFeatures() %>% 
                ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA")) %>% 
                RunPCA()
    X <- celltype@reductions$pca@cell.embeddings[,c(1:25)] # PCA embeddings
    cellQ <- Scale(log(celltype$nCount_RNA)) # log scaled (0-1) RNA counts per cell
    
    print(paste0("Making cell object ", name, "."))
    cell.object <- make.cell.type(name = name, tpm = as.matrix(tpm), samples = samples, X = X, metadata = metadata, cellQ = cellQ)
    return(cell.object)
}

# 3. GENERATE CELL SPECIFIC OBJECTS_____________________________________________________
###############################
# Generate List of Cell types HIV+ Only
###############################
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

rA_CD4 <- MakeCellObject(subset = "CD4", seurat_object = CD4)

# CD8 T cells
CD8[["Annotation"]] <- "CD8"
CD8 <- subset(CD8, cells = colnames(CD8)[CD8$HATIMID %in% HATIMID_CD8$HATIMID])
CD8[['HATIMID']] <- droplevels(CD8$HATIMID)
CD8[['HATIMID']] <- as.character(CD8$HATIMID)

rA_CD8 <- MakeCellObject(subset = "CD8", seurat_object = CD8)

# Macrophage Collapse
metadata <- Macrophage@meta.data
HATIMID_Macrophage <- metadata %>% dplyr::count(HATIMID) %>% filter(n > 30) %>% dplyr::select(HATIMID)
Macrophage <- subset(Macrophage, cells = colnames(Macrophage)[Macrophage$HATIMID %in% HATIMID_Macrophage$HATIMID])
Macrophage[['HATIMID']] <- droplevels(Macrophage$HATIMID)
Macrophage[['HATIMID']] <- as.character(Macrophage$HATIMID)

Macrophage_list <- c("LAM", "IM", "PVM", "Mo-Mac 1", "Mo-Mac 2")
rA_Macrophage <- lapply(Macrophage_list, MakeCellObject, seurat_object = Macrophage)

# Stromal
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte", "Early Preadipocyte", "ECM-Producing Early Preadipocyte"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Progenitor"
Stromal[['Annotation']] <- Idents(Stromal)
Stromal[['HATIMID']] <- droplevels(Stromal$HATIMID)
Stromal[['HATIMID']] <- as.character(Stromal$HATIMID)

Stromal_list <- c("Preadipocyte", "Progenitor", "PCOLCE+ Fibroblast", "Myofibroblast")
rA_Stromal <- lapply(Stromal_list, MakeCellObject, seurat_object = Stromal)

###############################
# Merge and save object
###############################
# Merge
rA <- c(rA_CD4, rA_CD8, rA_Macrophage, rA_Stromal)

tmp_dir <- paste0("../HIVneg")
dir.create(tmp_dir)
saveRDS(rA, file = paste(tmp_dir, "rA_HIVneg.rds", sep = "/"))

rm(Stromal)
rm(Macrophage)
rm(CD4)
rm(CD8)

# 4. Run Dialogue_____________________________________________________
###############################
# Run Dialogue All Cells
###############################
R <- DIALOGUE.run(rA = rA, main = "HIVneg", k = 5, results.dir = tmp_dir, pheno = "GlucoseIntolerant", plot.flag = T, conf = c('cellQ', 'Sex', 'BMI', 'Age'))
saveRDS(R, file = paste(tmp_dir, "Results_HIVneg.rds", sep = "/"))
