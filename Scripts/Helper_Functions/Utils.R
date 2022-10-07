###############################################################################
######------------------------FUNCTIONS-------------------------------#####
###############################################################################
library(tidyverse)
library(PResiduals)
library(Seurat)
library(ggrepel)
library(ggforce)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(gridExtra)
library(scales)
library(reshape2)
library(cowplot)
library(rms)
library(viridis)

###############################
# Load Metadata
###############################
data_hatim <- readr::read_csv('/data/p_koethe_lab/Atlas_AT/MetaData/Reference/HATIMStudyVisit_DATA_2021-08-09_1359.csv')
data_hatim <- data_hatim %>% mutate(HATIMID = as.factor(hatim_clin_visit_pid),
                                    StudyGroup = factor(hatim_final_arm, levels = c(1,2,3,4),
                                                        labels = c("HIV+ non-diabetic", "HIV+ prediabetic", 
                                                        "HIV+ diabetic", "HIV- diabetic")),
                                    HIV = case_when(hatim_final_arm == 4 ~ "HIVneg",
                                                    TRUE ~ "HIVpos"),
                                    Glucose = case_when(hatim_final_arm == 1 ~ "Glucose Tolerant",
                                                        TRUE ~ "Glucose Intolerant"),
                                    Sex = factor(meta_sex, levels = c(0,1), labels= c("Male", "Female")),
                                    homa2_ir = as.numeric(homa2_ir), 
                                    Age = meta_age,
                                    BMI = meta_bmi,
                                    metformin = case_when(stringr::str_detect(other_meds, regex('metformin', ignore_case = T))~ "Yes",stringr::str_detect(other_meds, regex('metphormin', ignore_case = T))~ "Yes",
                                    TRUE~"No"),
                                    FBG = meta_fbg,
                                    hba1c = meta_hba1c,
                                    Older_NRTI = case_when(older_nrti == 1~"Yes",TRUE~"No"))

colors <- c("aliceblue", "aquamarine1", "azure", "bisque", "blue", "burlywood", "chartreuse", "coral", "cyan", "cornsilk", "chocolate", "darkolivegreen", "darkorange", "darksalmon", "deepskyblue", "dodgerblue",
"firebrick", "gold", "gray", "green", "greenyellow", "khaki", "lavenderblush", "lightblue", "lightcoral", "lightcyan", "lightgreen", "lightpink", "lightsalmon", "magenta", "maroon", "lightsteelblue",
"mediumturquoise", "navy", "orange", "orchid", "purple", "red", "rosybrown", "seagreen", "seashell", "sandybrown", "slategray", "sienna", "springgreen", "tan", "violet", "turquoise", "yellow", "salmon", 
"plum", "peru", "black", "yellow", "brown")

###-------------------####
# UMAP 
###-------------------####
UMAP_FUN <- function(seurat_object, path, color, idents = "Celltype_plot", breaks, label, xlab = "", ylab = "", plot_label = F, res = 600, height = 8, width = 10, units = 'in'){
    Idents(seurat_object) <- idents
    plot <- DimPlot(seurat_object, label = plot_label, repel = T, raster = FALSE, label.size = 6) + 
    scale_colour_manual(name = "Cell Type", breaks = breaks,
    labels = label, values = color) + xlab(xlab) + ylab(ylab) + 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14))
    ggsave(path, dpi = res, height = height, width = width, units = units, device = 'png')
    return(plot)
}

###-------------------####
# DotPlot
###-------------------####
DOTPLOT_FUN <- function(seurat_object, features, path, ident = "Annotation", width = 10, height = 10, res = 600, units = 'in') {
    DefaultAssay(seurat_object) <- "RNA"
    Idents(seurat_object) <- ident
    plot <- DotPlot(seurat_object, features = features, group.by = ident) + scale_color_gradient2(low = "steelblue", 
                        high = "red") + 
    RotatedAxis() +
    theme(axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
    # Save Plot
    ggsave(path, device = 'png', width = width, height = height, units = units, dpi = res)
    return(plot)
}

###-------------------####
# Calulate Average ADT per cluster
###-------------------####
CITESEQ_Avg <- function(seurat_object) {
    DefaultAssay(seurat_object) <- "ADT"
    Anno_list <- as.character(unique(seurat_object$Annotation))
    Results <- data.frame(matrix(ncol = length(rownames(seurat_object)) + 1, nrow = 0))
    colnames(Results) <- c("Annotation", rownames(seurat_object))
    for (i in 1:length(Anno_list)) {
        cluster <- seurat_object@assays$ADT@data[, seurat_object$Annotation == Anno_list[[i]]]
        Avg <- rowMeans(cluster)
        Avg <- c(Anno_list[[i]], Avg)
        df <- data.frame(t(Avg))
        colnames(df)[1] <- "Annotation"
        Results <- rbind(Results, df)
    }
    return(Results)
}

###-------------------####
# Prep Plotting
###-------------------####
CITESEQ_PlotPrep <- function(seurat_object, CITESeq_Averages, ADT) {
    metadata <- seurat_object@meta.data
    row_names <- rownames(metadata)
    CITESeq_Averages <- CITESeq_Averages[,c("Annotation",ADT)]
    i <- 2:dim(CITESeq_Averages)[2]
    CITESeq_Averages[, i] <- apply(CITESeq_Averages[, i], 2, function(x) as.numeric(as.character(x)))
    metadata <- left_join(metadata, CITESeq_Averages, by = "Annotation")
    metadata$barcode <- row_names
    UMAP <- data.frame(Embeddings(seurat_object, reduction = 'umap'))
    UMAP$barcode <- rownames(UMAP)
    joined <- inner_join(metadata, UMAP, by = "barcode")
    joined <- joined[,c(ADT, "UMAP_1", "UMAP_2")]
    return(joined)
}

###-------------------####
# Plot Average ADT
###-------------------####
CITESEQ_Plot <- function(DF) {  
    title <- DF$measure[1]
    plot <- DF %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = value)) + 
    geom_point(size = 0.2) +
    scale_color_viridis() +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    theme_classic() + 
    ggtitle(title) + 
    theme(plot.title = element_text(face = "bold",size=18, hjust = 0.5, color = "black"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line=element_blank())
    return(plot)
}

###-------------------####
# GO Pathway Analysis
###-------------------####
GO_Function <- function(markers, celltype, filename) {
    clusters <- markers %>% filter(cluster %in% celltype)
    Entrez <- clusterProfiler::bitr(clusters$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    cluster_markers <- left_join(Entrez, clusters, by = c("SYMBOL" = "gene"))
    cluster_markers <- cluster_markers %>% dplyr::filter(p_val_adj < 0.05) %>% arrange(cluster) # Filter to genes < 0.05
    cluster_list <- cluster_markers %>% group_by(cluster) %>% group_split()
    names(cluster_list) <- celltype

    # Loop through list of DGE for each Macrophage cluster
    DF <- data.frame()
    DF <- lapply(cluster_list, function(x) {
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
    })

    # Plot Results
    Pathways <- bind_rows(DF, .id = "Cluster") # Combine in data frame

    plot <- Pathways %>% group_by(Cluster) %>%
        do(head(.,5)) %>% 
        mutate(Description = factor(Description, levels = Description)) %>%
        ggplot() + geom_point(aes(x = Cluster, y = Description, size = -log10(pvalue), fill = Count), shape = 21) +
        xlab("") + ylab("") + ggtitle("") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1, size = 20, face = "bold"), 
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 20, face = "bold")) + scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
    
    plot
    ggsave(file = filename,  device = 'png', width = 10, height = 5, units = 'in', dpi = 600)
    return(plot)
}

###-------------------####
# KEGG Pathway Analysis
###-------------------####
KEGG_ORA <- function(markers, celltype, filename, flip = "No") {
    clusters <- markers %>% filter(cluster %in% celltype)
    uni <- clusterProfiler::bitr(clusters$gene, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = org.Hs.eg.db)
    uni <- uni[!duplicated(uni$SYMBOL),]
    cluster_markers <- left_join(uni, clusters, by = c("SYMBOL" = "gene"))
    cluster_markers <- cluster_markers %>% dplyr::filter(p_val_adj < 0.05) %>% arrange(cluster) # Filter to genes < 0.05
    cluster_list <- cluster_markers %>% group_by(cluster) %>% group_split()
    names(cluster_list) <- celltype

    # Loop through list of DGE for each Macrophage cluster
    DF <- data.frame()
    DF <- lapply(cluster_list, function(x) {
        KEGG <- enrichKEGG(gene = x$UNIPROT,
                    organism = 'hsa',
                    keyType = "uniprot",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)
        results <- data.frame(KEGG)
        DF <- rbind(DF, results)
    })

    # Plot Results
    Pathways <- bind_rows(DF, .id = "Cluster") # Combine in data frame

    plot <- Pathways %>% group_by(Cluster) %>%
        do(head(.,5)) %>% 
        mutate(Description = factor(Description, levels = Description)) %>%
        ggplot() + geom_point(aes(x = Cluster, y = Description, size = Count, fill = pvalue), shape = 21) +
        xlab("") + ylab("") + ggtitle("") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1, size = 14, face = "bold", color = "black"), 
        axis.text.y = element_text(size = 14, face = "bold", color = "black"), legend.title = element_text(size = 14, face = "bold")) + scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
    
    if(flip == "Yes") {
        plot <- plot + coord_flip()
    }
    plot
    ggsave(file = filename,  device = 'png', width = 15, height = 10, units = 'in', dpi = 600)
    return(plot)
}

###############################
# Box Plot: Plot Percent celltype by Study Group
############################### 
box_plot_fun_paper <- function(df, my_comparisons, 
                         parent_cell_type = "Lymphoid",
                         ncol = 4,
                         legend.pos = "right",
                         plot_cols = c("#91C46C", "#7ABAF2", "#FFBE00"), connect = F,
                         is_paired = FALSE, x_axis_text_size = 16,
                         specify_x_order = NA, specify_x_labs = NA,
                         level_p_val = F, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))){
  plot_size <- abs(diff(range(df$value)))
  if(level_p_val){
    stat_pos <- rep(max(df$value) + plot_size * 0.05,times = length(my_comparisons))
  } else if(!level_p_val){
    stat_pos <- c(max(df$value),
                  max(df$value) + plot_size * 0.10,
                  max(df$value) + plot_size * 0.20)
  }
  stat_pos <- stat_pos[1:length(my_comparisons)]
  
  plot_title <- parent_cell_type
  y_axis_title <- paste0("% of ",parent_cell_type)
  df <- df %>%
  dplyr::group_by(measure) %>%
  dplyr::mutate(Avg = mean(value)) %>%
  dplyr::arrange(desc(Avg))
  df$measure <- factor(df$measure, levels = df$measure[!duplicated(df$measure)])
  
  plot <- df %>%
  ggplot(aes(x = StudyGroup, y = value, fill = factor(StudyGroup))) +
  facet_wrap(~measure, ncol = ncol) +
  geom_boxplot(mapping = aes(fill = StudyGroup), outlier.shape = NA, alpha = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.75, lwd = 1) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.4, position = "dodge") +
  scale_fill_manual(values = plot_cols, name = "Study Group") + 
  ggtitle(plot_title) +
  ylab(y_axis_title) + 
  xlab("") +
  theme_classic() +
  ylim(0, stat_pos[[length(my_comparisons)]]*1.1) +
  stat_compare_means(comparisons = my_comparisons, label.y = stat_pos, paired = is_paired, method = "wilcox.test", symnum.args = symnum.args) +
   theme(plot.title = element_text(face = "bold",size=18, hjust = 0.5, color = "black"),
          axis.text.y = element_text(face = "bold", size = 13, color = "black"),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 16, color = "black"),
          axis.title.x = element_blank(),
          strip.text = element_text(size = 14, face = "bold"),
          legend.position = legend.pos,
          legend.title = element_text(colour="black", size=16, 
                                      face="bold"),
          legend.text = element_text(colour="black", size=13, 
                                     face="bold"))
    return(plot)
}

###############################
# Correlation function
###############################
Cor_rho <- function(df, adjusted, variable, nonDM = F) {
    Results <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(Results) <- c("rho", "pvalue", "Cluster")
    numb <- dim(df)[2]-10
    if (nonDM == T) {
        df <- df %>% filter(StudyGroup != "HIV+ diabetic")
    }
    
    if (adjusted == T) {
        if (variable == "BMI") {
            for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,variable] ~ Sex  + StudyGroup + Age, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i])
                Results <- rbind(Results, res)
            }
        } else if (variable == "Age") {
             for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,variable] ~ Sex  + StudyGroup + BMI, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i])
                Results <- rbind(Results, res)
            }
        } else if (variable %in% c("homa2_ir", "meta_hba1c", "meta_fbg")) {
            for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(df[,i] | df[,variable] ~ Sex + BMI + Age, data = df, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(df)[i])
                Results <- rbind(Results, res)
            }
        } 
    } else {
        for (i in 1:numb) {
            cor_res <- cor.test(df[,i], df[,variable], method = "spearman")
            res <- data.frame(rho = cor_res$estimate, pvalue = cor_res$p.value, colnames(df)[i])
            Results <- rbind(Results, res)
        }
    }
    return(Results)
}

###############################
# Correlation function: This will read in
# the cell proportions and perform spearman's (adjusted or unadjusted) 
# on the specified subset and return a dataframe of the results
###############################
Cor_rho_subset <- function(df1, df2, subset, HIV = "NA") {
    Results <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(Results) <- c("rho", "pvalue", "Cluster")
    numb <- dim(df1)[2]-10
    df2 <- df2[,c(1:(dim(df2)[2]-9))]
    merge <- inner_join(df1, df2, by = "HATIMID")
    HIVm <- merge %>% filter(HIV == "HIVpos")
    HIVneg <- merge %>% filter(HIV == "HIVneg")
    

        if (HIV == "T") {
             for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(HIVm[,i] | HIVm[,subset] ~ Sex + StudyGroup + BMI + Age, data = HIVm, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(HIVm)[i])
                Results <- rbind(Results, res)
            }
        } else if (HIV == "F") {
            for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(HIVneg[,i] | HIVneg[,subset] ~ Sex  + BMI + Age, data = HIVneg, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(HIVneg)[i])
                Results <- rbind(Results, res)
            }
        } else {
            for (i in 1:numb) {
                cor_res <- PResiduals::partial_Spearman(merge[,i] | merge[,subset] ~ Sex + StudyGroup + BMI + Age, data = merge, link.x = "logit", link.y = "logit")
                res <- data.frame(rho = cor_res$TS$TB$ts, pvalue = cor_res$TS$TB$pval, Cluster = colnames(merge)[i])
                Results <- rbind(Results, res)
        }
    }    
    return(Results)
}

###############################
# Merge Proportion and Metadata
############################### 
Prop_Merge <- function(Prop, metadata) {
    Prop$HATIMID <- as.factor(rownames(Prop))
    Prop <- left_join(Prop, metadata[,c("StudyGroup", "HATIMID", "HIV", "Glucose", "Sex", "Age", "BMI", "homa2_ir", "meta_fbg", "meta_hba1c")], by = "HATIMID")
    Prop <- Prop[!is.na(Prop[,1]),]
    return(Prop)
}

###############################
# Correlation Plotting: This will plot the rho per cluster
# super-imposed on the UMAP
###############################  
Cor_plot <- function(Results, seurat_object, title) {
    metadata <- seurat_object@meta.data
    New_met <- inner_join(metadata, Results, by = c("Annotation" = "Cluster"))
    New_met$barcode <- rownames(metadata)
    UMAP <- data.frame(Embeddings(seurat_object, reduction = 'umap'))
    UMAP$barcode <- rownames(UMAP)
    
    joined <- inner_join(New_met, UMAP, by = "barcode")
    plot <- joined %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = rho)) + 
    geom_point(size = 0.2) +
    scale_color_gradient2(low = "blue", high = "red", mid = "gray", midpoint = 0, limit = range(joined$rho), space = "lab", name = "Correlation") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    theme_classic() + 
    ggtitle(title) + 
    theme(plot.title = element_text(face = "bold",size=18, hjust = 0.5, color = "black"),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank())
    return(plot)
}

###############################
# ORM Model: This will use Ordinal
# regression to find the influence of Sex
# on proportions adjusted for other metadata. Ordinal
# was chosen due to the non-parametric distrubtion of the
# data.
############################### 
ORM_Model <- function(df) {
    df$StudyGroup <- droplevels(df$StudyGroup)
    dd <- rms::datadist(df)
    options(datadist = dd)
    
    Results <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(Results) <- c("Estimate", "Cluster")
    numb <- dim(df)[2]-10
    
    for (i in 1:numb) {
        f <- rms::orm(df[,i] ~ Sex + Age + BMI + StudyGroup, data = df)
        data_orm <- data.frame(summary(f))
        Sex_odds <- data.frame(Estimate = data_orm[rownames(data_orm) == "X.Odds.Ratio.2", 4], Cluster = colnames(df)[i])
        Results <- rbind(Results, Sex_odds)
    }
    return(Results)
}

###############################
# ORM Plot: This will use the Odds Ratio
# output from the ORM model and plot
# on each cluster for Sex
############################### 
ORM_Plot <- function(Results, seurat_object, title) {
    metadata <- seurat_object@meta.data
    New_met <- inner_join(metadata, Results, by = c("Annotation" = "Cluster"))
    New_met$barcode <- rownames(metadata)
    UMAP <- data.frame(Embeddings(seurat_object, reduction = 'umap'))
    UMAP$barcode <- rownames(UMAP)
    
    joined <- inner_join(New_met, UMAP, by = "barcode")
    plot <- joined %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = Estimate)) + 
    geom_point(size = 0.2) +
    scale_color_gradient2(low = "darkblue", mid = 'grey', high = "red", midpoint = 1, limit = range(joined$Estimate), space = "lab", name = "Odds Ratio F:M") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    theme_classic() + 
    ggtitle(title) + 
    theme(plot.title = element_text(face = "bold",size=18, hjust = 0.5, color = "black"),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank())
    return(plot)
}

###############################
# GO GSEA: This function will return
# a list of GSEA using GO BP
############################### 
GO_GSEA <- function(df, ont = "BP", minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05){
  df$gene <- rownames(df)
  df <- df %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::select(gene, log2FoldChange)
  geneList = df[,2]
  names(geneList) = as.character(df[,1])
  
  ego <- gseGO(geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  minGSSize = minGSSize,
  maxGSSize = maxGSSize,
  pvalueCutoff = pvalueCutoff,
  verbose = FALSE,
  keyType = "SYMBOL")
  ego <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  
  if(dim(ego)[1] > 0) { # Only for those with pathways, otherwise will error
    Pathways_gsea <- ego %>% dplyr::select(Description, setSize, NES, p.adjust, enrichmentScore, core_enrichment)
    return(Pathways_gsea)
  } 
}

###############################
# KEGG GSEA: This function will return
# a list of GSEA using KEGG
############################### 
KEGG_GSEA <- function(df, minGSSize = 50, pvalueCutoff = 0.05, keyType = "uniprot") {
  df$gene <- rownames(df)
  uni <- clusterProfiler::bitr(df$gene, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = org.Hs.eg.db)
  uni <- uni[!duplicated(uni$SYMBOL),] # Key returns > 1 mapping. Several don't actually map to gene. Only the first one does so therefore, will exclude duplicates
  df <- left_join(uni, df, by = c("SYMBOL" = "gene"))
  df <- df %>% arrange(desc(log2FoldChange)) %>% dplyr::select(UNIPROT, log2FoldChange)
  geneList = df[,2]
  names(geneList) = as.character(df[,1])

  kk2 <- gseKEGG(geneList = geneList,
  organism = "hsa",
  minGSSize = minGSSize,
  pvalueCutoff = pvalueCutoff,
  keyType = keyType)

    if(dim(kk2)[1] > 0) {
    Pathways_KEGG <- kk2 %>% dplyr::select(Description, setSize, NES, p.adjust, enrichmentScore, core_enrichment)
    return(Pathways_KEGG)
  } 
}

###############################
# GSEA order: Order dataframe 
# by NES
############################### 
GSEA_order <- function(df) {
  df <- as.data.frame(df)
  df <- df %>% arrange(desc(df$NES))
  return(df)
}

###############################
# Volcano Plot: Plot DGE
# of Pseudobulk
############################### 
VOLCANO <- function(resOrderedDF, units = "in", height = 8, width = 8, filename, dpi = 600, cols = c("Up" = "firebrick3", "Down" = "steelblue", "ns" = "grey"), genes, ylim = c("",""), 
xlim = c(-2.5,2.5)) {
    resOrderedDF$gene <- rownames(resOrderedDF)
    resOrderedDF <- resOrderedDF %>% dplyr::mutate(threshold = factor(case_when(log2FoldChange > 0.25 & padj < 0.05 ~ 'Up', log2FoldChange < -0.25 & padj < 0.05  ~ 'Down',
    TRUE ~ 'ns')),
    sizes = case_when(threshold %in% c("Up", "Down")~2, TRUE~1),
    alphas = case_when(threshold %in% c("Up", "Down")~1, TRUE~0.5))
    resOrderedDF <- resOrderedDF %>% mutate(padj = case_when(padj < 0.0000001~0.0000001, TRUE~padj))
    
    genes_interest <- resOrderedDF %>% dplyr::filter(gene %in% genes)
    
    plot <- ggplot(data = resOrderedDF, aes(x = log2FoldChange, y = -log10(padj), fill = threshold, size = sizes, alpha = alphas)) +
    geom_point(shape = 21, colour = "black") +
    xlim(xlim) +
    ylim(ylim) +
    geom_vline(xintercept = c(-0.25, 0.25), color = 'black', alpha = 1.0, linetype = 'longdash') +
    geom_hline(yintercept = 1.3, color = 'black', alpha = 1.0, linetype = 'longdash') +
    scale_fill_manual(values = cols) + 
    theme_bw() + 
    xlab(expression('Average Log'[2]*" "*'Fold Change')) +
    ylab(expression("-log"[10]*" "*"adjusted p-value")) + 
    theme(axis.text = (element_text(size = 16)),
                  plot.caption = element_text(face = "italic", size = 9),
                  legend.position = "none",
                  panel.background = element_rect(fill = "transparent"),
                  panel.border = element_blank(),
                  axis.line = element_line(color="black"),
                  strip.text = element_text(size = 24), 
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2),
                  axis.title.y = element_text(face = "bold", size = 16, color = "black"),
                  axis.title.x = element_text(face = "bold", size = 16, color = "black")) +
    geom_label_repel(data = genes_interest, aes(label = gene, color = "black"), force = 2, nudge_y = 1, color = "black", fill = "white")
    ggsave(filename, units = units, height = height, width = width, dpi = dpi)
    return(plot)
}

###############################
# GSEA Plot: Plot Pathways
# using dotplot
############################### 
GSEA_plot <- function(df, units = "in", height = 8, width = 8, filename, dpi = 600, device = 'png', num_path = 10, legend.pos = "right") {
  if (dim(df)[1] > 20) {
    df_max <- df %>% slice_max(NES, n = num_path)
    df_min <- df %>% slice_min(NES, n = num_path)
    df <- rbind(df_max, df_min)
  }
  list <- str_split(df$core_enrichment, "/")
  count <- data.frame(matrix(data = NA, nrow = 0, ncol = 1))
  for (i in 1:length(list)) {
    count_i <- length(list[[i]])
    count <- rbind(count, count_i)
  }
  colnames(count) <- "Count"
  df <- cbind(df, count)
  plot <- df %>% arrange(NES) %>% mutate(Description = factor(Description, levels = unique(Description))) %>%
    ggplot(aes(x = NES, y = Description, size = Count, colour = -log10(p.adjust)), shape = 21) + 
    geom_point() + 
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    scale_color_gradient(low = "darkorange", high = "#EB34D2") +
    theme(axis.text.x = element_text(size = 14, face = "bold"),
                  axis.text.y = element_text(size = 12, face = "bold"),
                  plot.caption = element_text(face = "italic", size = 9),
                  legend.position = legend.pos,
                  panel.background = element_rect(fill = "transparent"),
                  panel.border = element_blank(),
                  axis.line = element_line(color="black"),
                  strip.text = element_text(size = 24), 
                  axis.title = element_text(size = 18),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)) +
  ylab("")
  ggsave(filename, units = units, height = height, width = width, dpi = dpi, device = device)
  return(plot)
}


###############################
# Upset Plot: Generate Upset plot
# for genes by cell type
###############################
upset_gene <- function(df) {
  col <- ifelse(rownames(mat) %in% df$gene, 1, 0)
  name <- names(df)
  mat <- cbind(mat, data.frame(name = col))
  return(mat)
}

###############################
# Plot Dialogue MCP Heatmap
###############################
Heatmap_fun <- function(seurat_object, subset = F, R, metadata = data_hatim, cluster, genes, tmp_dir = paste0("/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/", date, "/Dialogue/", "All"),
    units = "in", height = 10, width = 15, res = 300, filename, MCP_dim = 1) {
    
    # Pull out genes
    clustergene <- gsub(" ", ".", cluster)
    cluster.up <- R$MCPs[[paste0("MCP",MCP_dim)]][[paste0(clustergene, ".up")]]
    cluster.down <- R$MCPs[[paste0("MCP",MCP_dim)]][[paste0(clustergene, ".down")]]
    cluster_genes <- c(cluster.up, cluster.down)
    
    # Get Expression
    if(subset) {
        Idents(seurat_object) <- "Annotation"
        seurat_object <- subset(seurat_object, idents = cluster)
    }
    Avg <- AverageExpression(seurat_object, group.by = "HATIMID")
    Avg <- Avg[[1]]
    
    # Pull out cluster-specific genes
    Avg <- Avg[rownames(Avg) %in% cluster_genes,]
    
    # Get rid of lowly expressed genes (likely ambient RNA contamination)
    Exp <- data.frame(Expression = rowMeans(seurat_object@assays$RNA@data[rownames(seurat_object@assays$RNA@data) %in% cluster_genes,]))
    Exp$Gene <- rownames(Exp)
    Exp <- Exp[order(Exp$Expression, decreasing = T),]
    Exp <- Exp[Exp$Expression > 0.15,]
    Avg <- Avg[rownames(Avg) %in% Exp$Gene,]
    names <- colnames(Avg)
    names(names) <- NULL
    
    # Scale Expression
    mat_scaled <- t(apply(Avg,1,scale))
    colnames(mat_scaled) <- names

    # Generate Metadata
    DF <- data.frame(mat_scaled)
    DF$MCP <- ifelse(rownames(DF) %in% cluster.up, "Up", "Down")
    MCP_label <- DF$MCP
    DF$MCP <- NULL
    colnames(DF) <- names
    tDF <- data.frame(t(DF))
    tDF$HATIMID <- as.character(rownames(tDF))
    merged <- left_join(tDF, metadata[,c("HATIMID", "StudyGroup", "BMI", "Age", "Sex", "FBG", "hba1c", "metformin")], by = "HATIMID")

    # Heatmap
    ha = HeatmapAnnotation(MP = MCP_label, col = list(MP = c("Up" = "red", "Down" = "blue")))
    col_fun = colorRamp2(c(25,35,45), c("blue", "white", "red"))
    col_fun2 = colorRamp2(c(25, 45, 65), c("#2D8A3D", "white", "#F49F2C"))
    col_fun3 = colorRamp2(c(80, 110, 150), c("blue", "white", "red"))
    col_fun4 = colorRamp2(c(5, 6, 8), c("#79FF33", "white", "#FF5733"))
    row_ha = rowAnnotation(Group = merged$StudyGroup, Sex = merged$Sex, BMI = merged$BMI, Age = merged$Age, FBG = merged$FBG, HbA1c = merged$hba1c, col = list(Group = c("HIV+ non-diabetic"="green", "HIV+ prediabetic" = "orange", "HIV+ diabetic" = "red"),
    Sex = c("Male" = "dodgerblue", "Female" = "orange"), BMI = col_fun, Age = col_fun2, FBG = col_fun3, HbA1c = col_fun4))

    index <- which(rownames(mat_scaled) %in% c(genes))
    label <- rownames(mat_scaled)[index]
    column_ha <- columnAnnotation(foo = anno_mark(at = index, labels = label, which = "column", side = "bottom"))
    
    png(file = paste(tmp_dir, filename, sep = "/"), res = res, width = width, height = height, units = units)
    print({
      Heatmap(t(mat_scaled), column_order = cluster_genes[cluster_genes %in% Exp$Gene], top_annotation = ha, right_annotation = row_ha, bottom_annotation = column_ha, show_column_names = F)
    })
    dev.off()
    
    plot <- Heatmap(t(mat_scaled), column_order = cluster_genes[cluster_genes %in% Exp$Gene], top_annotation = ha, right_annotation = row_ha, bottom_annotation = column_ha, show_column_names = F)
    return(plot)
}

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
Pseudobulk <- function(agg.sce, cluster, StudyGroup1 = "HIV+ diabetic", StudyGroup2 = "HIV+ non-diabetic", exclude) {
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
    resOrderedwaldt <- res[order(res$padj),] # order by pvalue
    resOrderedwaldt <- resOrderedwaldt[!is.na(resOrderedwaldt$padj),] # Remove NAs
    resOrderedDFwaldt <- as.data.frame(resOrderedwaldt)
    return(resOrderedDFwaldt)
}

######################
# Function for writing AnnData with normalized and log transformed counts, raw counts, and metadata (Seurat to H5ad conversion is clunky.)
######################
writeAnnData <- function(seurat_object, filename) {
    anndata::AnnData(X = t(as.matrix(GetAssayData(seurat_object, slot = "data", assay = "RNA"))),
    var = data.frame(gene = rownames(seurat_object),
    row.names = rownames(seurat_object)),
    layers = list(counts = t(as.matrix(GetAssayData(seurat_object, slot = "count", assay = "RNA")))),
    obs = data.frame(celltype = as.character(seurat_object$Annotation),
    louvain = as.character(seurat_object$seurat_clusters),
    StudyGroup = as.character(seurat_object$StudyGroup),
    Age = as.numeric(seurat_object$Age),
    BMI = as.numeric(seurat_object$BMI),
    HATIMID = as.character(seurat_object$HATIMID),
    Sex = as.character(seurat_object$Sex),
    batch = as.character(seurat_object$Lane),
    row.names = colnames(seurat_object)),
    obsm = list(X_umap = matrix(Embeddings(seurat_object, reduction = "umap"), ncol = 2)))$write_h5ad(filename, compression = "gzip")
}