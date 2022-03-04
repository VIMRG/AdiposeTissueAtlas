############-----------Helper Functions-------------###########

############-----------------------------------------###########
#SoupX
# Details: Ambient RNA contamination is common in single cell. Often
# genes of the most common cell types contaminate the background of
# droplets resulting of expression in cells that shouldn't express those
# genes. In adipose tissue data, endothelial cells are the most abundant
# and hence, markers for endothelial cells (e.g. CLDN5, GNG11, etc.) likely
# represent most of the background contamination. SoupX will identify these
# ambient RNA genes and remove from the count matrix.
###########-----------------------------------------#############
# input list of 10x matrices locations
SoupX_Fun <- function(matrix) {
    print(paste("Reading in 10x Matrix from Project located at ", dataset_local, "/", matrix, ".", sep = ""))
    sc = load10X(paste(dataset_local, matrix, inputfile_local, sep = "/")) # Will use Cellranger supplied clustering
    print("Calculating rho value of each channel.")
    sc = autoEstCont(sc) # Auto estimates reasonable and rho in line with Souporcell estimates of ambient RNA
    print("Adjusting Counts")
    out = adjustCounts(sc) # Generate corrected count matrix
    print("Writing corrected matrix")
    dir.create(file.path("/data/p_koethe_lab/Atlas_AT/SoupX", matrix), recursive = T) # Create directory to store corrected matrix
    tmpdir <- paste("/data/p_koethe_lab/Atlas_AT/SoupX", matrix, sep = "/")
    DropletUtils:::write10xCounts(tmpdir, out, overwrite = TRUE) # Use DropletUtils write10xCounts function to write matrix, feature, and barcode to directory
    print(paste("Corrected Matrix has successfully saved to ", tmpdir, ".", sep = ""))
}

############-----------------------------------------###########
#Read10X
# Details: Takes as input directory containing count matrix,
# to load as sparse matrix
###########-----------------------------------------#############
Read_sparse <- function(matrix, matrix_location) {
    print(paste("Reading in 10x Matrix from Project located at ", matrix_location, "/", matrix, ".", sep = ""))
    matrix <- Read10X(data.dir = paste(matrix_location, matrix, sep = "/"))
    return(matrix)
}

############-----------------------------------------###########
#SouporCell Tables
# Details: reads soup or cell .tsv cluster assigment files which will
# be merged with the seurat object to select only singlet cells later on
###########-----------------------------------------#############
Soup_tables <- function(table, dataset_local) {
    print(paste("Loading Genetic Demultiplexing Cluster Assignment for ", names(table), " located at ", dataset_local, "/", table, ".", sep = ""))
    cluster <- read.table(paste(dataset_local, table, "clusters.tsv", sep = "/"), header = T) # read in tsv file
    rownames(cluster) <- cluster$barcode # Need to set rownames to barcode to merge with seurat object later on
    return(cluster)
}

############-----------------------------------------###########
#SouporCell Tables
# Details: reads soup or cell .tsv cluster assigment files which will
# be merged with the seurat object to select only singlet cells later on
###########-----------------------------------------#############
Citeseq_name <- function(ADT_matrix) {
     rownames(ADT_matrix) <- mgsub(rownames(ADT_matrix), c("TotalSeq_C_", "TotalSeq_C"), c("", ""))
     return(ADT_matrix)
}

############-----------------------------------------###########
#Create Seurat Object
# Details: read in counts from corrected assay
###########-----------------------------------------#############
CreateSeurat <- function(seurat_data) {
    print(paste("Creating Gene Expression Seurat Object for Project ", names(seurat_data), ".", sep = ""))
    seurat_obj <- CreateSeuratObject(counts = seurat_data[['Gene Expression']], min.cells = 3)
    #seurat_obj[['HTO']] <- CreateAssayObject(seurat_data$`Antibody Capture`[1:4,][,colnames(seurat_obj)])
    #seurat_obj[['ADT']] <- CreateAssayObject(seurat_data$`Antibody Capture`[5:49,][,colnames(seurat_obj)])
    #assign(names(seurat_data), seurat_obj)
}


