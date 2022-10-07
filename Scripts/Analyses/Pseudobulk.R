#############################################################################################
##------------------------- PSEUDOBULK ANALYSIS ---------------------------------##
##------------------------- DATE: 6/27/2022 AUTHOR: SAM BAILIN-------------------------------##
## DESCRIPTION: This script will evaluate DGE using Pseudobulk.
##############################################################################################

# 1. SET UP_____________________________________________________
###############################
# Import Libraries
###############################
library(Seurat)
library(tidyverse)
library(scuttle)
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(ggplot2)

set.seed(7412)

date = "6.27"

###############################
# Load seurat objects
###############################
tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/SubsetAnalysis")

Vascular <- readRDS(paste0(tmp_dir, "/", "Vascular.rds"))
Stromal <- readRDS(paste0(tmp_dir, "/", "Stromal.rds"))
Myeloid <- readRDS(paste0(tmp_dir, "/", "Myeloid.rds"))
Macrophage <- readRDS(paste0(tmp_dir, "/", "Macrophage.rds"))
CD4 <- readRDS(paste0(tmp_dir, "/", "CD4.rds"))
CD8 <- readRDS(paste0(tmp_dir, "/", "CD8.rds"))

tmp_dir <- paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Pseudobulk")
dir.create(tmp_dir)
DM_dir <- paste(tmp_dir, "DiabeticVsNonDM", sep = "/")
dir.create(DM_dir)
PreDM_dir <- paste(tmp_dir, "PreDMvsNonDM", sep = "/")
dir.create(PreDM_dir)

# 2. FUNCTIONS_____________________________________________________
###############################
# Aggregate Cell groups by individual
###############################
Aggregate <- function(seurat_object) {
    sce <- as.SingleCellExperiment(seurat_object, assay = "RNA") # Convert to SCE
    agg.sce <- aggregateAcrossCells(sce, ids = colData(sce)[,c("Pseudobulk", "HATIMID")]) # Aggregate by cluster and sample
    colData(agg.sce)$Pseudobulk <- as.character(colData(agg.sce)$Pseudobulk) # Change Manual Annotation to character for next step
    colData(agg.sce)$HATIMID <- as.character(colData(agg.sce)$HATIMID) # Change Manual Annotation to character for next step
    colnames(agg.sce) <- apply(colData(agg.sce)[,c("Pseudobulk", "HATIMID")], 1, function(x) paste(x, collapse="_")) # Add column names
    return(agg.sce)
}

###############################
# Evaluate Outliers
###############################
Outlier_Eval <- function(agg.sce, cluster, StudyGroup1 = "HIV+ diabetic", StudyGroup2 = "HIV+ non-diabetic") {
    agg.sce <- subset(agg.sce, , Pseudobulk == cluster) # subset on cluster of interest
    agg.sce <- subset(agg.sce, , StudyGroup %in% c(StudyGroup1, StudyGroup2)) # subset on studygroup of interest
    agg.sce <- subset(agg.sce, , ncells > 5) # take only clusters with > 5 cells
    
    if (StudyGroup2 == "HIV+ non-diabetic") {
        label2 = "HIVnonDM"
    } else if (StudyGroup2 == "HIV- diabetic") {
        label2 = "HIVnegDM"
    } else if (StudyGroup2 == "HIV+ prediabetic") {
        label2 = "HIVpreDM"
    }
    
        if (StudyGroup1 == "HIV+ diabetic") {
        label1 = "HIVDM"
    } else if (StudyGroup1 == "HIV+ prediabetic") {
        label1 = "HIVpreDM"
    }
    
    metadata <- colData(agg.sce)[,c("StudyGroup", "Age", "Sex", "BMI", "Lane")]
    metadata$StudyGroup <- factor(metadata$StudyGroup, levels = c(StudyGroup1, StudyGroup2), labels = c(label1, label2))
    metadata$Sex <- factor(metadata$Sex)
    metadata$Lane <- factor(metadata$Lane)
    
    if (all(rownames(metadata) != colnames(assay(agg.sce, 'counts')))) {
        stop("Metadata rownames do not equal assay column names.") # Using stop function
    }
    
    dds <- DESeqDataSetFromMatrix(round(assay(agg.sce, 'counts')), 
                    colData = metadata, 
                    design = ~ StudyGroup + Age + BMI + Sex)
    vst_dds <- vst(dds, blind = TRUE)
    return(vst_dds)
}

###############################
# Plot Outliers
###############################
Outlier_Plot <- function(agg.sce, vst, cluster, StudyGroup1 = "HIV+ diabetic", StudyGroup2 = "HIV+ non-diabetic") {
    agg.sce <- subset(agg.sce, , Pseudobulk == cluster) # subset on cluster of interest
    agg.sce <- subset(agg.sce, , StudyGroup %in% c(StudyGroup1, StudyGroup2)) # subset on studygroup of interest
    agg.sce <- subset(agg.sce, , ncells > 5) # take only clusters with > 5 cells
    
    if (StudyGroup2 == "HIV+ non-diabetic") {
        label2 = "HIVnonDM"
    } else if (StudyGroup2 == "HIV- diabetic") {
        label2 = "HIVnegDM"
    } else if (StudyGroup2 == "HIV+ prediabetic") {
        label2 = "HIVpreDM"
    }
    
        if (StudyGroup1 == "HIV+ diabetic") {
        label1 = "HIVDM"
    } else if (StudyGroup1 == "HIV+ prediabetic") {
        label1 = "HIVpreDM"
    }
    
    metadata <- colData(agg.sce)[,c("StudyGroup", "Age", "Sex", "BMI", "Lane")]
    metadata$StudyGroup <- factor(metadata$StudyGroup, levels = c(StudyGroup1, StudyGroup2), labels = c(label1, label2))
    metadata$Sex <- factor(metadata$Sex)
    metadata$Lane <- factor(metadata$Lane)
    vst_mat <- assay(vst)
    vst_cor <- cor(vst_mat)
    annotation <- data.frame(metadata[,c("StudyGroup", "Lane", "Sex", "Age", "BMI"), drop = F])
    plot <- pheatmap(vst_cor, annotation_col = annotation)
    return(plot)
}

###############################
# Perform Pseudobulk
###############################
Pseudobulk <- function(agg.sce, cluster, StudyGroup1 = "HIV+ diabetic", StudyGroup2 = "HIV+ non-diabetic", exclude, shrinkage = T) {
    agg.sce <- subset(agg.sce, , Pseudobulk == cluster) # subset on cluster of interest
    agg.sce <- subset(agg.sce, , StudyGroup %in% c(StudyGroup1, StudyGroup2)) # subset on studygroup of interest
    agg.sce <- subset(agg.sce, , ncells > 5) # take only clusters with > 5 cells
    agg.sce <- subset(agg.sce, , !HATIMID %in% exclude) # exclude outlier
    
    if (StudyGroup2 == "HIV+ non-diabetic") {
        label2 = "HIVnonDM"
    } else if (StudyGroup2 == "HIV- diabetic") {
        label2 = "HIVnegDM"
    } else if (StudyGroup2 == "HIV+ prediabetic") {
        label2 = "HIVpreDM"
    }
    
    if (StudyGroup1 == "HIV+ diabetic") {
        label1 = "HIVDM"
    } else if (StudyGroup1 == "HIV+ prediabetic") {
        label1 = "HIVpreDM"
    }
    
    metadata <- colData(agg.sce)[,c("StudyGroup", "Age", "Sex", "BMI", "Lane")]
    metadata$StudyGroup <- factor(metadata$StudyGroup, levels = c(StudyGroup2, StudyGroup1), labels = c(label2, label1))
    metadata$Sex <- factor(metadata$Sex)
    metadata$Lane <- factor(metadata$Lane)
    
    if (all(rownames(metadata) != colnames(assay(agg.sce, 'counts')))) {
        stop("Metadata rownames do not equal assay column names.") # Using stop function
    }
    
    dds <- DESeqDataSetFromMatrix(round(assay(agg.sce, 'counts')), 
                    colData = metadata, 
                    design = ~ StudyGroup + Age + BMI + Sex)
    
    keep <- rowSums(counts(dds)) >= 3
    dds <- dds[keep,]
    dds <- DESeq(dds)
    
    if(shrinkage) {
        res <- lfcShrink(dds, coef = paste0("StudyGroup_", label1, "_vs_", label2), type = "apeglm")
    } else {
        res <- results(dds, contrast=c("StudyGroup", label1, label2))
    }

    summary(res)
    resOrderedwald <- res[order(res$padj),] # order by pvalue
    resOrderedwald <- resOrderedwald[!is.na(resOrderedwald$padj),] # Remove NAs
    resOrderedDFwald <- as.data.frame(resOrderedwald)
    return(resOrderedDFwald)
}


# 3. MACROPHAGE PSEUDOBULK_____________________________________________________
###############################
# OVERALL MACROPHAGE
###############################
# HIV+ DM vs HIV+ non-diabetic
Macrophage[['Pseudobulk']] <- "Macrophage"
agg.macrophage <- Aggregate(Macrophage)

# HIV+ DM vs HIV+ nonDM
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Macrophage")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Macrophage")

Mac.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Macrophage", exclude = c("1109"))
write.csv(Mac.DM.NonDM, file = paste(DM_dir, "Macrophage.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Macrophage", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Macrophage", StudyGroup1 = "HIV+ prediabetic")

Mac.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Macrophage", exclude = c("1109"), StudyGroup1 = "HIV+ prediabetic")
write.csv(Mac.PreDM.NonDM, file = paste(PreDM_dir, "Macrophage.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# PVM
###############################
# Pseudobulk
Macrophage[['Pseudobulk']] <- Idents(Macrophage)
agg.macrophage <- Aggregate(Macrophage)

# HIV+ DM vs HIV+ non-DM
Outlier <- Outlier_Eval(agg.macrophage, cluster = "PVM")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "PVM")

PVM.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "PVM", exclude = c("1132")) # Remove substantial outliers (spearmans 0.6)
write.csv(PVM.DM.NonDM, file = paste(DM_dir, "PVM.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "PVM", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "PVM", StudyGroup1 = "HIV+ prediabetic")

PVM.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "PVM", exclude = c(""), StudyGroup1 = "HIV+ prediabetic")
write.csv(PVM.PreDM.NonDM, file = paste(PreDM_dir, "PVM.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# IM
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "IM")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "IM")

IM.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "IM", exclude = c("1132", "3013")) # Remove substantial outliers
write.csv(IM.DM.NonDM, file = paste(DM_dir, "IM.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "IM", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "IM", StudyGroup1 = "HIV+ prediabetic")

IM.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "IM", exclude = c("1132"), StudyGroup1 = "HIV+ prediabetic")
write.csv(IM.PreDM.NonDM, file = paste(PreDM_dir, "IM.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# LAM
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "LAM")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "LAM")

LAM.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "LAM", exclude = c("3014")) # Remove substantial outliers (spearmans 0.6)
write.csv(LAM.DM.NonDM, file = paste(DM_dir, "LAM.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "LAM", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "LAM", StudyGroup1 = "HIV+ prediabetic")

LAM.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "LAM", exclude = c("1137"), StudyGroup1 = "HIV+ prediabetic")
write.csv(LAM.PreDM.NonDM, file = paste(PreDM_dir, "LAM.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# Mo-Mac
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Mo-Mac")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Mo-Mac")

MoMac.DM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Mo-Mac", exclude = c("1109")) # Remove substantial outliers (spearmans 0.6)
write.csv(MoMac.DM.NonDM, file = paste(DM_dir, "MoMac.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.macrophage, cluster = "Mo-Mac", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.macrophage, vst = Outlier, cluster = "Mo-Mac", StudyGroup1 = "HIV+ prediabetic")

MoMac.PreDM.NonDM <- Pseudobulk(agg.macrophage, cluster = "Mo-Mac", exclude = c("1109", "1167"), StudyGroup1 = "HIV+ prediabetic")
write.csv(MoMac.PreDM.NonDM, file = paste(PreDM_dir, "MoMac.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 3. CD4 PSEUDOBULK_____________________________________________________
###############################
# OVERALL CD4
###############################
CD4[['Pseudobulk']] <- "CD4"
agg.CD4 <- Aggregate(CD4)

# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.CD4, cluster = "CD4")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
DESeq2::plotPCA(Outlier, intgroup = 'Lane')
Plot <- Outlier_Plot(agg.CD4, vst = Outlier, cluster = "CD4")

CD4.DM.NonDM <- Pseudobulk(agg.CD4, cluster = "CD4", exclude = c("1167"))
write.csv(CD4.DM.NonDM, file = paste(DM_dir, "CD4.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.CD4, cluster = "CD4", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.CD4, vst = Outlier, cluster = "CD4", StudyGroup1 = "HIV+ prediabetic")

CD4.PreDM.NonDM <- Pseudobulk(agg.CD4, cluster = "CD4", exclude = "", StudyGroup1 = "HIV+ prediabetic")
write.csv(CD4.PreDM.NonDM, file = paste(PreDM_dir, "CD4.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 4. MONOCYTE PSEUDOBULK_____________________________________________________
Idents(Myeloid) <- "Annotation"
Idents(Myeloid, cells = WhichCells(Myeloid, idents = c("cMo", "nMo", "Other Mo", "ISG+ Mo"))) <- "Monocyte"
Myeloid[['Pseudobulk']] <- Idents(Myeloid)
agg.Myeloid <- Aggregate(Myeloid)

# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Myeloid, cluster = "Monocyte")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Myeloid, vst = Outlier, cluster = "Monocyte")

Mo.DM.NonDM <- Pseudobulk(agg.Myeloid, cluster = "Monocyte", exclude = c(""))
write.csv(Mo.DM.NonDM, file = paste(DM_dir, "Monocyte.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Myeloid, cluster = "Monocyte", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Myeloid, vst = Outlier, cluster = "Monocyte", StudyGroup1 = "HIV+ prediabetic")

Mo.PreDM.NonDM <- Pseudobulk(agg.Myeloid, cluster = "Monocyte", exclude = "", StudyGroup1 = "HIV+ prediabetic")
write.csv(Mo.PreDM.NonDM, file = paste(PreDM_dir, "Monocyte.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 5. CD8 PSEUDOBULK_____________________________________________________
###############################
# OVERALL CD8
###############################
# HIV+ DM vs HIV+ non-diabetic
CD8[['Pseudobulk']] <- "CD8"
agg.CD8 <- Aggregate(CD8)

Outlier <- Outlier_Eval(agg.CD8, cluster = "CD8")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.CD8, vst = Outlier, cluster = "CD8")

CD8.DM.NonDM <- Pseudobulk(agg.CD8, cluster = "CD8", exclude = c("3015", "3037"))
write.csv(CD8.DM.NonDM, file = paste(DM_dir, "CD8.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.CD8, cluster = "CD8", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.CD8, vst = Outlier, cluster = "CD8", StudyGroup1 = "HIV+ prediabetic")

CD8.PreDM.NonDM <- Pseudobulk(agg.CD8, cluster = "CD8", exclude = "", StudyGroup1 = "HIV+ prediabetic")
write.csv(CD8.PreDM.NonDM, file = paste(PreDM_dir, "CD8.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 6. VASCULAR PSEUDOBULK_____________________________________________________
###############################
# OVERALL ENDOTHELIAL
###############################
Idents(Vascular) <- "Annotation"
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Venous EC", "Arterial EC", "Intermediate Capillary EC", "Capillary EC", "Ven-Cap EC"))) <- "EC"
Idents(Vascular, cells = WhichCells(Vascular, idents = c("Arterial EndoMT-like", "Venous EndoMT-like", "Capillary EndoMT-like"))) <- "EndoMT"
Vascular[['Pseudobulk']] <- Idents(Vascular)
agg.Vascular <- Aggregate(Vascular)

# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EC")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EC")

EC.DM.NonDM <- Pseudobulk(agg.Vascular, cluster = "EC", exclude = c("1113", "1182"))
write.csv(EC.DM.NonDM, file = paste(DM_dir, "EC.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EC", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EC", StudyGroup1 = "HIV+ prediabetic")

EC.PreDM.NonDM <- Pseudobulk(agg.Vascular, cluster = "EC", exclude = "", StudyGroup1 = "HIV+ prediabetic")
write.csv(EC.PreDM.NonDM, file = paste(PreDM_dir, "EC.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# OVERALL ENDOMT
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EndoMT")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EndoMT")

EndoMT.DM.NonDM <- Pseudobulk(agg.Vascular, cluster = "EndoMT", exclude = "")
write.csv(EndoMT.DM.NonDM, file = paste(DM_dir, "EndoMT.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Vascular, cluster = "EndoMT", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Vascular, vst = Outlier, cluster = "EndoMT", StudyGroup1 = "HIV+ prediabetic")

EndoMT.PreDM.NonDM <- Pseudobulk(agg.Vascular, cluster = "EndoMT", exclude = "", StudyGroup1 = "HIV+ prediabetic")
write.csv(EndoMT.PreDM.NonDM, file = paste(PreDM_dir, "EndoMT.HIVpreDM.HIVnonDM.csv", sep = "/"))

# 7. STROMAL PSEUDOBULK_____________________________________________________
###############################
# OVERALL STROMAL
###############################
Idents(Stromal) <- "Annotation"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Mature Preadipocyte 1", "Mature Preadipocyte 2", "Early Preadipocyte", "ECM-Producing Early Preadipocyte"))) <- "Preadipocyte"
Idents(Stromal, cells = WhichCells(Stromal, idents = c("Adipose Progenitor Cell 1", "Adipose Progenitor Cell 2"))) <- "Progenitor"
Stromal[['Pseudobulk']] <- Idents(Stromal)
agg.Stromal <- Aggregate(Stromal)

# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Preadipocyte")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Preadipocyte")

PreAd.DM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Preadipocyte", exclude = "")
write.csv(PreAd.DM.NonDM, file = paste(DM_dir, "Preadipocyte.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Preadipocyte", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Preadipocyte", StudyGroup1 = "HIV+ prediabetic")

Preadipocyte.PreDM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Preadipocyte", exclude = "", StudyGroup1 = "HIV+ prediabetic")
write.csv(Preadipocyte.PreDM.NonDM, file = paste(PreDM_dir, "Preadipocyte.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# Progenitor
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Progenitor")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Progenitor")

Progenitor.DM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Progenitor", exclude = "3013")
write.csv(Progenitor.DM.NonDM, file = paste(DM_dir, "Progenitor.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Progenitor", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Progenitor", StudyGroup1 = "HIV+ prediabetic")

Progenitor.PreDM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Progenitor", exclude = "1165", StudyGroup1 = "HIV+ prediabetic")
write.csv(Progenitor.PreDM.NonDM, file = paste(PreDM_dir, "Progenitor.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# PCOLCE+ Fibroblast
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "PCOLCE+ Fibroblast")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "PCOLCE+ Fibroblast")

PCOLCE.DM.NonDM <- Pseudobulk(agg.Stromal, cluster = "PCOLCE+ Fibroblast", exclude = "3013")
write.csv(PCOLCE.DM.NonDM, file = paste(DM_dir, "PCOLCE.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "PCOLCE+ Fibroblast", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "PCOLCE+ Fibroblast", StudyGroup1 = "HIV+ prediabetic")

PCOLCE.PreDM.NonDM <- Pseudobulk(agg.Stromal, cluster = "PCOLCE+ Fibroblast", exclude = "1109", StudyGroup1 = "HIV+ prediabetic")
write.csv(PCOLCE.PreDM.NonDM, file = paste(PreDM_dir, "PCOLCE.HIVpreDM.HIVnonDM.csv", sep = "/"))

###############################
# Myofibroblast
###############################
# HIV+ DM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Myofibroblast")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Myofibroblast")

Myofibroblast.DM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Myofibroblast", exclude = "1167")
write.csv(Myofibroblast.DM.NonDM, file = paste(DM_dir, "Myofibroblast.HIVDM.HIVnonDM.csv", sep = "/"))

# HIV+ PreDM vs HIV+ non-diabetic
Outlier <- Outlier_Eval(agg.Stromal, cluster = "Myofibroblast", StudyGroup1 = "HIV+ prediabetic")
DESeq2::plotPCA(Outlier, intgroup = 'StudyGroup')
Plot <- Outlier_Plot(agg.Stromal, vst = Outlier, cluster = "Myofibroblast", StudyGroup1 = "HIV+ prediabetic")

Myofibroblast.PreDM.NonDM <- Pseudobulk(agg.Stromal, cluster = "Myofibroblast", exclude = c("2110", "1167"), StudyGroup1 = "HIV+ prediabetic")
write.csv(Myofibroblast.PreDM.NonDM, file = paste(PreDM_dir, "Myofibroblast.HIVpreDM.HIVnonDM.csv", sep = "/"))

