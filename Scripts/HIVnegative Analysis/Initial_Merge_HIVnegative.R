##############################################################################################
##------------------------- INITIAL HIV+ and HIV- INTEGRATION---------------------------------##
##------------------------- DATE: 6/27/2022 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: The following code will perform integration with harmony of the HIV+ and HIV- individuals.
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
path <- "/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/DoubletFinder/"
file_list <- list.files(path = path, pattern = ".rds")

project <- c('P5344_CW3', 'P5544_CW3', 'P5573_CW3', 'P5657_CW3', 'P5836_CW3', 'P5877_CW3', 'P5903_CW3', 'P5963_CW3')

seurat_path = paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/DoubletFinder")
seurat_list <- list.files(path = seurat_path, full.names = T)
seurat_list <- seurat_list[grep(pattern = "(_CW3)", seurat_list)]
names(seurat_list) <- project
seurat_list <- lapply(seurat_list, function(xx) {
    print(paste0("Reading saved seurat files for project ", "from ", xx, ".")) 
    readRDS(xx)
})

# Add Lane/Batch Information
for (i in 1:length(seurat_list)) {
    seurat_list[[i]]@meta.data['Lane'] <- names(seurat_list)[i]
}

HIV <- readRDS("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Final_Integration/Integrated.rds")

# 2. MERGE AND INTEGRATE_____________________________________________________
###############################
# Merge Together
###############################
HIVneg <- merge(x = seurat_list[[1]], y = seurat_list[c(2:8)], add.cell.ids = c("5344CW3", "5544CW3", "5573CW3", "5657CW3", "5836CW3", "5877CW3", "5903CW3", "5963CW3"))
integrated.data <- merge(x = HIV, y = HIVneg)

rm(seurat_list)
rm(HIVneg)
rm(HIV)

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
tmp_dir <- paste("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Merged_Initial", sep = "")
dir.create(tmp_dir)

saveRDS(integrated.data, file = paste(tmp_dir, "Integrated_HIVnegative.rds", sep = "/"))
