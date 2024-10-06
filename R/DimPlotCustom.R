#' @export
DimPlotCustom <- function(seurat, group.by = "orig.ident", 
                          pt.size = 1, 
                          plotLegend = T, 
                          label.size = 6, 
                          title.size = 8, 
                          nrow = NULL, 
                          ncol = NULL, 
                          split.by = NULL, 
                          plot = T, 
                          label = F,
                          rasterize = F,
                          raster.dpi = 600){
  embeddings <- as.data.frame(seurat@reductions$umap@cell.embeddings)
  xLimits <- c(min(embeddings$UMAP_1), max(embeddings$UMAP_1))
  yLimits <- c(min(embeddings$UMAP_2), max(embeddings$UMAP_2))
  
  embeddings$group <- seurat@meta.data[,group.by]
  
  if(is.null(split.by)){
    p <- ggplot2::ggplot(data=embeddings, ggplot2::aes(x=UMAP_1, y=UMAP_2)) + 
      ggplot2::geom_point(ggplot2::aes(fill= group), size = pt.size, colour = "black", stroke = .1, shape = 21) + 
      ggplot2::theme_classic() +
      ggplot2::xlim(xLimits) + 
      ggplot2::ylim(yLimits) + 
      ggplot2::ylab("UMAP 2") + 
      ggplot2::xlab("UMAP 1") + 
      ggplot2::theme(plot.tag = ggplot2::element_text(size = 12), legend.position = "bottom", legend.direction = "horizontal", text = ggplot2::element_text(size = label.size), legend.text = ggplot2::element_text(size = 6), plot.margin = ggplot2::unit(c(0,0,0,0), "cm"), plot.title = ggplot2::element_text(hjust = 0.5, size = title.size))
  } else {
    splitDF <- seurat@meta.data[drop = F,,split.by]
    if(!is.null(levels(splitDF[,1]))){
      splits <- levels(splitDF[,1])
    } else {
      splits <- unique(splitDF[,1])
    }
    p <- lapply(splits, function(y) {
      ggplot2::ggplot(data = embeddings[row.names(splitDF[splitDF[,split.by] %in% y,,drop = F]),], aes(x= UMAP_1, y=UMAP_2)) + 
        ggplot2::geom_point(aes(fill= group), size = pt.size, colour = "black", stroke = .1, shape = 21) +
        ggplot2::theme_classic() +
        ggplot2::xlim(xLimits) + 
        ggplot2::ylim(yLimits) + 
        ggplot2::ylab("UMAP 2") + 
        ggplot2::xlab("UMAP 1") + 
        ggplot2::theme(plot.tag = ggplot2::element_text(size = 12), legend.position = "bottom", legend.direction = "horizontal", text = ggplot2::element_text(size = label.size), legend.text = ggplot2::element_text(size = 6), plot.margin = unit(c(0,0,0,0), "cm"), plot.title = ggplot2::element_text(hjust = 0.5, size = title.size)) + ggplot2::ggtitle(y) 
    })
  }
  if(label ==T){
    if(is.null(split.by)){
      centroids <- aggregate(cbind(UMAP_1,UMAP_2) ~ group, data = embeddings, FUN=mean)
      p <- p + ggplot2::geom_text(data = centroids, size = label.size, mapping = aes(x=UMAP_1, y=UMAP_2, label=group))
    } else {
      embeddings$split <- seurat@meta.data[, split.by]
      centroids <- aggregate(cbind(UMAP_1,UMAP_2) ~ group + split, data= embeddings, FUN=mean)
      for(i in 1:length(splits)){
        p[[i]] <- p[[i]] + ggplot2::geom_text(data = centroids[centroids$split == splits[i],], size = label.size, mapping = aes(x=UMAP_1, y=UMAP_2, label=group))
      }
    }
  }
  if(plot == T){
    if(is.null(ncol) & is.null(nrow) & !is.null(split.by)){
      ncol = length(splits)
    }
    if(rasterize == T){
      p <- p + ggrastr::rasterise(ggplot2::geom_point(ggplot2::aes(fill= group), size = pt.size, colour = "black", stroke = .1, shape = 21), dpi = raster.dpi)
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

