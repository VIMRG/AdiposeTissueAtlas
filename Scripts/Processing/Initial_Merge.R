##############################################################################################
##------------------------- INITIAL HIV+ INTEGRATION---------------------------------##
##------------------------- AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: The following code will perform integration with harmony of the HIV+ individuals.
## This includes doublets tagged by genetic demultiplexing and DoubletFinder.
##############################################################################################


# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(harmony)
library(SingleCellExperiment)


date = "6.27"

###############################
# Import Dataset
###############################
path <- "../DoubletFinder/"
file_list <- list.files(path = path, pattern = ".rds")

project <- c('P5344_CW1', 'P5344_CW2', 'P5544_CW1', 'P5544_CW2', 'P5573_CW1', 'P5573_CW2', 'P5657_CW1', 'P5836_CW1', 'P5836_CW2', 
'P5877_CW1', 'P5877_CW2', 'P5903_CW1', 'P5903_CW2', 'P5963_CW1', 'P5963_CW2')

seurat_path = paste0("../HATIM_Analysis/", date, "/DoubletFinder")
seurat_list <- list.files(path = seurat_path, full.names = T)
seurat_list <- seurat_list[grep(pattern = "(_CW1|_CW2)", seurat_list)]
names(seurat_list) <- project
seurat_list <- lapply(seurat_list, function(xx) {
    print(paste0("Reading saved seurat files for project ", "from ", xx, ".")) 
    readRDS(xx)
})

# Add Lane/Batch Information
for (i in 1:length(seurat_list)) {
    seurat_list[[i]]@meta.data['Lane'] <- names(seurat_list)[i]
}


# 2. MERGE AND INTEGRATE_____________________________________________________
###############################
# Merge Together
###############################
integrated.data <- merge(seurat_list[[1]], y = seurat_list[c(2:15)], 
add.cell.ids=c("5344CW1", "5344CW2", "5544CW1", "5544CW2", "5573CW1", "5573CW2", "5657CW1", "5836CW1", "5836CW2", "5877CW1", 
"5877CW2", "5903CW1", "5903CW2", "5963CW1", "5963CW2"))

rm(seurat_list)
###############################
# Run variable features, normalization, and scale
###############################
DefaultAssay(integrated.data) <- "RNA"
integrated.data <- NormalizeData(integrated.data, assay = "RNA")
integrated.data <- NormalizeData(integrated.data, assay = "ADT", normalization.method = "CLR", margin = 2)

integrated.data <- integrated.data %>% FindVariableFeatures(verbose = T, nfeatures = 5000) %>%
        ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA"), assay = "RNA") %>%
        RunPCA()
        
###############################
# Run Harmony
###############################
integrated.data <- RunHarmony(integrated.data, "Lane", assay.use = "RNA", max.iter.harmony = 30)

###############################
# Run Donwstream Analysis
###############################
integrated.data <- integrated.data %>% RunUMAP(reduction = 'harmony', dims = 1:40) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:40) %>% 
  FindClusters(resolution = 3.0)

# 3. SAVE_____________________________________________________
date = "6.27"
tmp_dir <- paste("../HATIM_Analysis/", date, "/Merged_Initial", sep = "")
dir.create(tmp_dir)

saveRDS(integrated.data, file = paste(tmp_dir, "Integrated_Doublets.rds", sep = "/"))

