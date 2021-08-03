#' @export
FeaturePlotCustom <- function(seurat, genes, plot = T, tag = element_blank(), plotLegend = T, cellNum = T, pt.size = .25, labelFont = 6, titleFont = 8, sharedScale = "All", nrow = NULL, ncol = NULL, split.by = NULL, sideBySide = T, color0 = "grey60", colPalette = c("dodgerblue2", "gold", "red2")){
  embeddings <- as.data.frame(seurat@reductions$umap@cell.embeddings)
  expr <- seurat@assays$RNA@data[genes,, drop = F]
  scatter_col = c("grey60",colorRampPalette(c("dodgerblue2", "gold", "red2"))(max(expr)*100))
  if(is.null(split.by)){
    p <- lapply(genes, function(x) {ggplot2::ggplot(data=embeddings[match(names(sort(expr[x,], decreasing = F)), row.names(embeddings)),], aes(x=UMAP_1, y=UMAP_2)) + ggplot2::geom_point(aes(color= sort(expr[x,], decreasing = F)), size = pt.size) +
        ggplot2::labs(col="Expression", title = x) + ggplot2::theme_classic() +
        ggplot2::ylab(element_blank()) + ggplot2::xlab(element_blank()) + ggplot2::theme(plot.tag = element_text(size = 12), text = element_text(size = labelFont), legend.text = element_text(size = 6), plot.margin = unit(c(0,0,0,0), "cm"), plot.title = element_text(hjust = 0.5, size = titleFont)) +
        if(sharedScale %in% c("All", "Gene")){
          ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)), limits = c(0, max(expr)))
        } else if(sharedScale == "None"){
          ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)), limits = c(0, max(expr[x,])))
        }})
  } else {
    xLimits <- c(min(embeddings$UMAP_1), max(embeddings$UMAP_1))
    yLimits <- c(min(embeddings$UMAP_2), max(embeddings$UMAP_2))
    splitDF <- seurat@meta.data[drop = F,,split.by]
    splits <- unique(splitDF[,1])
    q <- lapply(splits, function(y) {
      lapply(genes, function(x) {
        ggplot2::ggplot(data=embeddings[match(names(sort(expr[x,colnames(expr) %in% row.names(splitDF[splitDF[,split.by] %in% y,, drop = F])], decreasing = F)), row.names(embeddings)),], aes(x=UMAP_1, y=UMAP_2)) + ggplot2::geom_point(aes(color= sort(expr[x,colnames(expr) %in% row.names(splitDF[splitDF[,split.by] %in% y,, drop = F])], decreasing = F)), size = pt.size) +
          ggplot2::labs(col="Expression", title = paste0(x, " - ", y)) + ggplot2::theme_classic() +
          ggplot2::ylab(element_blank()) + ggplot2::xlab(element_blank()) + ggplot2::theme(plot.tag = element_text(size = 12), text = element_text(size = labelFont), legend.text = element_text(size = 6), plot.margin = unit(c(0,0,0,0), "cm"), plot.title = element_text(hjust = 0.5, size = titleFont)) +
          ggplot2::xlim(xLimits) + ggplot2::ylim(yLimits) +
          if(sharedScale == "All"){
            ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)), limits = c(0, max(expr)))
          } else if(sharedScale == "Gene"){
            ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)), limits = c(0, max(expr[x,])))
          } else if(sharedScale == "None"){
            ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)))
          }
      })})
    if(sideBySide == T & !is.null(split.by)){
      p <- list()
      loopCounter <- 1
      for(i in 1:length(genes)){
        for(j in 1:length(splits)) {
          p[[loopCounter]] <-q[[j]][i]
          loopCounter <- loopCounter + 1
        }
      }
    }
    p <- unlist(p, recursive = F)
  }
  if(plot == T){
    if(is.null(ncol) & is.null(nrow)){
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

