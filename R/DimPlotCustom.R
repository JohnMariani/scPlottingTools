#' @export
DimPlotCustom <- function(seurat, group.by = "orig.ident", pt.size = 1, plotLegend = T, labelFont = 6, titleFont = 8, nrow = NULL, ncol = NULL, split.by = NULL, plot = T){
  embeddings <- as.data.frame(seurat@reductions$umap@cell.embeddings)
  xLimits <- c(min(embeddings$UMAP_1), max(embeddings$UMAP_1))
  yLimits <- c(min(embeddings$UMAP_2), max(embeddings$UMAP_2))
  
  embeddings$group <- seurat@meta.data[,group.by]
  
  if(is.null(split.by)){
    p <- ggplot2::ggplot(data=embeddings, aes(x=UMAP_1, y=UMAP_2)) + ggplot2::geom_point(aes(fill= group), size = pt.size, colour = "black", stroke = .1, shape = 21) + ggplot2::theme_classic() +
      ggplot2::xlim(xLimits) + ggplot2::ylim(yLimits) + ggplot2::ylab("UMAP 2") + ggplot2::xlab("UMAP 1") + ggplot2::theme(plot.tag = element_text(size = 12), legend.position = "bottom", legend.direction = "horizontal", text = element_text(size = labelFont), legend.text = element_text(size = 6), plot.margin = unit(c(0,0,0,0), "cm"), plot.title = element_text(hjust = 0.5, size = titleFont))
  } else {
    splitDF <- seurat@meta.data[drop = F,,split.by]
    splits <- unique(splitDF[,1])
    p <- lapply(splits, function(y) {
      ggplot2::ggplot(data = embeddings[row.names(splitDF[splitDF[,split.by] %in% y,,drop = F]),], aes(x= UMAP_1, y=UMAP_2)) + ggplot2::geom_point(aes(fill= group), size = pt.size, colour = "black", stroke = .1, shape = 21) + ggplot2::theme_classic() +
        ggplot2::xlim(xLimits) + ggplot2::ylim(yLimits) + ggplot2::ylab("UMAP 2") + ggplot2::xlab("UMAP 1") + ggplot2::theme(plot.tag = element_text(size = 12), legend.position = "bottom", legend.direction = "horizontal", text = element_text(size = labelFont), legend.text = element_text(size = 6), plot.margin = unit(c(0,0,0,0), "cm"), plot.title = element_text(hjust = 0.5, size = titleFont)) + ggplot2::ggtitle(y) 
    })
  }
  if(plot == T){
    if(is.null(ncol) & is.null(nrow) & !is.null(split.by)){
      ncol = length(splits)
    }
    if(plotLegend == F){
      patchwork::wrap_plots(p, nrow = nrow, ncol = ncol) & theme(legend.position = "none") 
    } else {
      patchwork::wrap_plots(p, nrow = nrow, ncol = ncol) 
    }
  } else {
    return(p)
  }
}

