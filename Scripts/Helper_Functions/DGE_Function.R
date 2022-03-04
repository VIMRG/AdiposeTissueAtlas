###############################
# MAST RNA DGE BY CELL TYPE
###############################
RNAclust_Markers <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "Manual_Annotation"
    DGE <- FindAllMarkers(seurat_object,
    test.use = "MAST", assay = "RNA", min.pct = 0.25, only.pos = T) # Perform DGE
    write.csv(DGE, file = paste(outpath, ".csv", sep = ""))
    return(DGE)
}

###############################
# MAST ADT DGE BY CELL TYPE
###############################
ADTclust_Markers <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "Manual_Annotation"
    DGE <- FindAllMarkers(seurat_object,
    test.use = "MAST", assay = "ADT", min.pct = 0.25, latent.vars = c("Project"), only.pos = T) # Perform DGE
    write.csv(DGE, file = paste(outpath, ".csv", sep = ""))
    return(DGE)
}

###############################
# Create New Column with Celltype and Study Group
###############################
CellType_Fun <- function(seurat_object) {
    print("Creating new metadta column combining cell type and Study Group.")
    seurat_object[['celltype.group']] <- paste(seurat_object$Manual_Annotation, seurat_object$StudyGroup, sep = "_")
    Idents(seurat_object) <- "celltype.group"
    return(seurat_object)
}

###############################
# MAST HIV+DM+ VS HIV+DM-
###############################
HIV_comparison <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "celltype.group"
    DF <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))
    for (i in unique(seurat_object$Manual_Annotation)) {
        ident.1 = paste0(i, "_HIV+ diabetic")
        ident.2 = paste0(i, "_HIV+ non-diabetic")
        print(paste("Finding differentially expressed genes between HIV+DM+ vs HIV+DM- in ", i, " using the seurat MAST wrapper implementation.", sep = ""))
        DGE <- FindMarkers(seurat_object, ident.1 = ident.1, ident.2 = ident.2, assay = "RNA", test.use = "MAST", min.pct = 0.25, latent.vars = c("Age", "Sex", "BMI"))
        DGE$gene <- rownames(DGE)
        DGE$cluster <- i
        DF <- rbind(DF, DGE)
    }
    print(paste("Saving CSV file of differentially expressed genes in the location ", outpath, ".csv", sep = ""))
    write.csv(DF, file = paste(outpath, ".csv", sep = ""))
    return(DF)
}

###############################
# MAST HIV+DM+ VS HIV+DM- No min percent
###############################
HIV_comparison.5 <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "celltype.group"
    DF <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))
    for (i in unique(seurat_object$Manual_Annotation)) {
        ident.1 = paste0(i, "_HIV+ diabetic")
        ident.2 = paste0(i, "_HIV+ non-diabetic")
        print(paste("Finding differentially expressed genes between HIV+DM+ vs HIV+DM- in ", i, " using the seurat MAST wrapper implementation.", sep = ""))
        DGE <- FindMarkers(seurat_object, ident.1 = ident.1, ident.2 = ident.2, assay = "RNA", test.use = "MAST", min.pct = 0.1, latent.vars = c("Age", "Sex", "BMI"))
        DGE$gene <- rownames(DGE)
        DGE$cluster <- i
        DF <- rbind(DF, DGE)
    }
    print(paste("Saving CSV file of differentially expressed genes in the location ", outpath, ".csv", sep = ""))
    write.csv(DF, file = paste(outpath, ".csv", sep = ""))
    return(DF)
}

###############################
# MAST HIV+DM+ VS HIV-DM+ No min percent
###############################
DM_comparison.5 <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "celltype.group"
    DF <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))
    for (i in unique(seurat_object$Manual_Annotation)) {
        ident.1 = paste0(i, "_HIV+ diabetic")
        ident.2 = paste0(i, "_HIV- diabetic")
        print(paste("Finding differentially expressed genes between HIV+DM+ vs HIV-DM+ in ", i, " using the seurat MAST wrapper implementation.", sep = ""))
        DGE <- FindMarkers(seurat_object, ident.1 = ident.1, ident.2 = ident.2, assay = "RNA", test.use = "MAST", min.pct = 0.1, latent.vars = c("Age", "Sex", "BMI"))
        DGE$gene <- rownames(DGE)
        DGE$cluster <- i
        DF <- rbind(DF, DGE)
    }
    print(paste("Saving CSV file of differentially expressed genes in the location ", outpath, ".csv", sep = ""))
    write.csv(DF, file = paste(outpath, ".csv", sep = ""))
    return(DF)
}

###############################
# MAST HIV+DM+ vs HIV-DM+
###############################
DM_comparison <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "celltype.group"
    DF <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))
    for (i in unique(seurat_object$Manual_Annotation)) {
        ident.1 = paste0(i, "_HIV+ diabetic")
        ident.2 = paste0(i, "_HIV- diabetic")
        print(paste("Finding differentially expressed genes between HIV+DM+ vs HIV-DM+ in ", i, " using the seurat MAST wrapper.", sep = ""))
        DGE <- FindMarkers(seurat_object, ident.1 = ident.1, ident.2 = ident.2, assay = "RNA", test.use = "MAST", min.pct = 0.25, latent.vars = c("Age", "Sex", "BMI"))
        DGE$gene <- rownames(DGE)
        DGE$cluster <- i
        DF <- rbind(DF, DGE)
    }
    
    # Save CSV file
    print(paste("Saving CSV file of differentially expressed genes in the location ", outpath, ".csv", sep = ""))
    write.csv(DF, file = paste(outpath, ".csv", sep = ""))
    
    return(DF)
}

###############################
# MAST HIV+DM+ vs HIV-DM+ no min percent
###############################
DM_comparison.5 <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "celltype.group"
    DF <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))
    for (i in unique(seurat_object$Manual_Annotation)) {
        ident.1 = paste0(i, "_HIV+ diabetic")
        ident.2 = paste0(i, "_HIV- diabetic")
        print(paste("Finding differentially expressed genes between HIV+DM+ vs HIV-DM+ in ", i, " using the seurat MAST wrapper.", sep = ""))
        DGE <- FindMarkers(seurat_object, ident.1 = ident.1, ident.2 = ident.2, assay = "RNA", test.use = "MAST", min.pct = 0.05, latent.vars = c("percent.mt", "Project"))
        DGE$gene <- rownames(DGE)
        DGE$cluster <- i
        DF <- rbind(DF, DGE)
    }
    
    # Save CSV file
    print(paste("Saving CSV file of differentially expressed genes in the location ", outpath, ".csv", sep = ""))
    write.csv(DF, file = paste(outpath, ".csv", sep = ""))
    
    return(DF)
}


###############################
# MAST RNA DGE BY SEURAT CLUSTERS
###############################
RNAseurat_Markers <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "seurat_clusters"
    DGE <- FindAllMarkers(seurat_object,
    test.use = "MAST", assay = "RNA", min.pct = 0.25, latent.vars = c("percent.mt", "Project"), only.pos = T) # Perform DGE
    write.csv(DGE, file = paste(outpath, ".csv", sep = ""))
    return(DGE)
}

###############################
# MAST ADT DGE BY Seurat Clusters
###############################
ADTseurat_Markers <- function(seurat_object, outpath) {
    Idents(seurat_object) <- "seurat_clusters"
    DGE <- FindAllMarkers(seurat_object,
    test.use = "MAST", assay = "ADT", min.pct = 0.25, latent.vars = c("Project"), only.pos = T) # Perform DGE
    write.csv(DGE, file = paste(outpath, ".csv", sep = ""))
    return(DGE)
}

###############################
# Create New Column with seurat clusters and Study Group
###############################
SeuratGroup_Fun <- function(seurat_object) {
    print("Creating new metadta column combining cell type and Study Group.")
    seurat_object[['seurat.group']] <- paste(seurat_object$seurat_clusters, seurat_object$StudyGroup, sep = "_")
    Idents(seurat_object) <- "seurat.group"
    return(seurat_object)
}

###############################
# MAST HIV+DM+ VS HIV+DM- by Seurat Clusters
###############################
HIV_comp_seurat <- function(seurat_object, outpath) {
    DF <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))
    for (i in unique(seurat_object$seurat_clusters)) {
        ident.1 = paste0(i, "_HIV_Diabetic")
        ident.2 = paste0(i, "_HIV_nonDiabetic")
        if (dim(subset(seurat_object, ident = ident.1))[2] > 50 & dim(subset(seurat_object, ident = ident.2))[2] > 50) {
            print(paste("Finding differentially expressed genes between HIV+DM+ vs HIV+DM- in ", i, " using the seurat MAST wrapper implementation.", sep = ""))
            DGE <- FindMarkers(seurat_object, ident.1 = ident.1, ident.2 = ident.2, assay = "RNA", test.use = "MAST", min.pct = 0.25, latent.vars = c("percent.mt", "Project"))
            DGE$gene <- rownames(DGE)
            DGE$cluster <- i
            DF <- rbind(DF, DGE)
        } else
        print(paste("One of the conditions has fewer than 50 cells in cluster ", i, " (skipping).", sep = ""))
        next
    }
    print(paste("Saving CSV file of differentially expressed genes in the location ", outpath, ".csv", sep = ""))
    write.csv(DF, file = paste(outpath, ".csv", sep = ""))
    return(DF)
}

###############################
# MAST HIV+DM+ vs HIV-DM+ by Seurat Clusters
###############################
DM_comp_seurat <- function(seurat_object, outpath) {
    DF <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))
    for (i in unique(seurat_object$seurat_clusters)) {
        ident.1 = paste0(i, "_HIV_Diabetic")
        ident.2 = paste0(i, "_HIVneg_Diabetic")
        if (dim(subset(seurat_object, ident = ident.1))[2] > 50 & dim(subset(seurat_object, ident = ident.2))[2] > 50) {
            print(paste("Finding differentially expressed genes between HIV+DM+ vs HIV-DM+ in ", i, " using the seurat MAST wrapper.", sep = ""))
            DGE <- FindMarkers(seurat_object, ident.1 = ident.1, ident.2 = ident.2, assay = "RNA", test.use = "MAST", min.pct = 0.25, latent.vars = c("percent.mt", "Project"))
            DGE$gene <- rownames(DGE)
            DGE$cluster <- i
            DF <- rbind(DF, DGE)
        } else
        print(paste("One of the conditions has fewer than 50 cells in cluster ", i, " (skipping).", sep = ""))
        next
    }
    # Save CSV file
    print(paste("Saving CSV file of differentially expressed genes in the location ", outpath, ".csv", sep = ""))
    write.csv(DF, file = paste(outpath, ".csv", sep = ""))
    return(DF)
}
