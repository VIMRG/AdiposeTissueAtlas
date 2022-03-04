###############################################################################
######------------------------FIGURE FUNCTION-------------------------------#####
###############################################################################
###-------------------####
# Import Libraries
###-------------------####
library(ggplot2)
library(tidyverse)


###-------------------####
# UMAP Plot Function
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
# Proportion Plot Function
###-------------------####
Plot_prop <- function(Prop, filter, sort, y.axis = "all", filename, device = 'png', width = 14, height = 8, units = 'in', dpi = 600, Margins = c(5.5, 5.5, 5.5, 5.5), plot_unit = 'pt') {
    theme_set(theme_bw() +
            theme(axis.text = (element_text(size = 16)),
                  plot.caption = element_text(face = "italic", size = 9),
                  legend.position = "bottom",
                  legend.background = element_rect(fill = "transparent"), 
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 24, face = "bold"),
                  panel.background = element_rect(fill = "transparent"),
                  panel.border = element_blank(),
                  axis.line = element_line(color="black"),
                  strip.text = element_text(size = 24), 
                  axis.title = element_text(size = 20, face = "bold"),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_line(color = "grey90", linetype = 2)
          ))
          
    Prop <- Prop %>% filter(Origin == filter)
    
    if(sort){
        Prop <- Prop %>%
        dplyr::arrange(desc(Median)) %>%
        mutate(cell.f = factor(cell.f, levels = unique(cell.f))) %>%
        group_by(Major) %>%
        mutate(rank = mean(Median)) %>% # Order by Major Cell Type and Rank
        ungroup() %>%
        dplyr::arrange(desc(rank)) %>%
        mutate(Major = factor(Major, levels = unique(Major)),
        StudyGroup = factor(StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic", "HIV- diabetic")))
        
        
        ggplot(Prop, aes(x = cell.f, y = value, fill = StudyGroup)) +
        geom_boxplot(position = position_dodge(0.8), l = 45, outlier.shape = NA) + 
        facet_grid(~Major, scales = "free", space = "free") +
        scale_fill_manual(values = c("#F04832", "#F57636", "#EBAA3B", "#4287F5"), name = "Disease State") +
         theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1, size = 16, face = "bold"), plot.margin = margin(Margins, unit = plot_unit), axis.title.x = element_blank(), axis.text.y = element_text(size = 20)) +
        ylab(paste('% Total ', y.axis, ' cells'))
        
    } else {
        Prop <- Prop %>%
        dplyr::arrange(desc(Median)) %>%
        mutate(cell.f = factor(cell.f, levels = unique(cell.f)),
        StudyGroup = factor(StudyGroup, levels = c("HIV+ non-diabetic", "HIV+ prediabetic", "HIV+ diabetic", "HIV- diabetic")))
        
        ggplot(Prop, aes(x = cell.f, y = value, fill = StudyGroup)) +
        geom_boxplot(position = position_dodge(0.8), l = 45, outlier.shape = NA) +
        scale_fill_manual(values = c("#F04832", "#F57636", "#EBAA3B", "#4287F5"), name = "Disease State") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1, size = 20, face = "bold"), plot.margin = margin(Margins, unit = plot_unit), axis.title.x = element_blank()) +
        ylab(paste('% Total ', y.axis, ' cells'))

    }
    # Save Plot
    ggsave(filename, device = 'png', width = width, height = height, units = units, dpi = dpi)
}


###-------------------####
# DotPlot Function
###-------------------####
DOTPLOT_FUN <- function(seurat_object, features, path, ident = "Manual_Annotation", width = 12, height = 8, res = 600, units = 'in') {
    DefaultAssay(seurat_object) <- "RNA"
    Idents(seurat_object) <- ident
    Plot <- DotPlot(seurat_object, features = features, group.by = ident) + scale_color_gradient2(low = "dodgerblue", 
            high = "red") + 
            theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size = 16, face = "bold"), 
            legend.text = element_text(size = 14)) +
            RotatedAxis()
    # Save Plot
    ggsave(path, device = 'png', width = width, height = height, units = units, dpi = res)
    return(Plot)
}

###-------------------####
# Cell Marker Heatmap
###-------------------####
HEATMAP_FUN <- function(seurat_object, markers, path, res = 600, width = 10, height = 10, units = 'in', Margins = c(5.5, 5.5, 5.5, 5.5), plot_unit = "pt") {
    DefaultAssay(seurat_object) <- "RNA"
    Idents(seurat_object) <- "Manual_Annotation"
    markers.top <- markers %>%
    group_by(cluster) %>%
    do(head(.,10))
    seurat_object[['Manual_Annotation']] <- factor(seurat_object$Manual_Annotation, levels = unique(markers.top$cluster))
    
    maxcells <- min(table(seurat_object$Manual_Annotation))
    DoHeatmap(subset(seurat_object, downsample = maxcells), group.colors = col_shuffle, features = markers.top$gene, group.by = 'Manual_Annotation') + 
    theme(axis.text.y = element_text(size = 6), plot.margin = margin(Margins, unit = plot_unit)) + 
    scale_fill_gradientn(colors = c("dodgerblue", "black", "gold"))
    ggsave(path, dpi = res, width = width, height = height, units = units, device = 'png')
}

###-------------------####
# Volcano Plot
###-------------------####
VOLCANO <- function(DGE, filename, dpi = 600, width = 10, height = 10, units = 'in', genes, filter) {
    DGE <- DGE %>% dplyr::filter(cluster == filter)
    DGE <- DGE %>% dplyr::mutate(threshold = factor(case_when(avg_log2FC > 0.25 & p_val_adj < 0.05 ~ 'cond1', avg_log2FC < -0.25 & p_val_adj < 0.05  ~ 'cond2',
    TRUE ~ 'cond3')))
    DGE$label <- DGE$gene %in% genes
    
    plot <- ggplot(DGE, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = threshold, alpha = 0.5)) +
    xlim(-2.5, 2.5) +
    xlab('Average Log2 Fold Change') +
    ylab("-log10 adjusted p-value") + 
    geom_vline(xintercept = c(-0.25, 0.25), color = 'springgreen4', alpha = 1.0, linetype = 'longdash') +
    geom_hline(yintercept = 1.3, color = 'springgreen4', alpha = 1.0, linetype = 'longdash') +
    theme_bw() + scale_color_manual(name = "Condition", values = c("cond1" = 'firebrick3', 'cond2'= 'dodgerblue4', 'cond3' = 'grey')) + theme(legend.position="none") +
    geom_text_repel(aes(x = avg_log2FC, -log10(p_val_adj)),label = ifelse(DGE$label == TRUE, as.character(DGE$gene),""), box.padding = unit(0.45, "lines"),hjust=1) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 18, face = 'bold'),
    axis.title.x = element_text(size = 18, face = 'bold'))
    
    ggsave(filename, dpi = dpi, device = 'png', width = width, height = height, units = units)
    plot
    return(plot)
}
###-------------------####
# Colors
###-------------------####
col_shuffle <-c("#AAB7B8", "#C02FE0", "#9F9314", "#E3D559", "#A483C4", "#40E66F", "#52AB16",
"#F9C99B", "#F53722", "#14659F", "#607C3B","#868BCF", "#47AD02", "#E36659",
"#A2F390", "#FAE31B", "#37A44B", "#9F4E14", "#F6E89F", "#EB7C1C", "#1C9CD1",
"#FBC0C0", "#FFA047", "#0E97C9", "#9F1420", "#4FBDBC", "#FD5935", "#DE83BC",
"#2EA315", "#15C27A", "#0E6655", "#D98880", "#5DADE2", "#8EEDE3", "#C6A7F3",
"#EFF99E", "#8147D3", "#0E0C8A", "#C483A3", "#7AB056", "#9493F0", "#8B838F",
"#B2424C", "#0EC9C5", "#3432D4", "#AEEAD4", "#57BC42", "#50ABA9", "#13D2BE",
"#CADF15", "#FF5733", "#65AB8E", "#F9D198", "#EC95EE", "#699FB3", "#06C07C",
"#CB24CE", "#E4C627", "#8DD3F1", "#E5E206")