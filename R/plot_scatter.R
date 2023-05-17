#' Scatter plot of single cell data
#'
#' Creates a scatter plot of two features (typically feature expression), across a
#' set of single cells. Cells are colored by their identity class. Pearson
#' correlation between the two features is displayed above the plot.
#'
#' @param object Seurat object
#' @param feature1 First feature to plot. Typically feature expression but can also
#' be metrics, PC scores, etc. - anything that can be retreived with FetchData
#' @param feature2 Second feature to plot.
#' @param cells Cells to include on the scatter plot.
#' @param shuffle Whether to randomly shuffle the order of points. This can be
#' useful for crowded plots if points of interest are being buried. (default is FALSE)
#' @param seed Sets the seed if randomly shuffling the order of points.
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param cols Colors to use for identity class plotting.
#' @param pt.size Size of the points on the plot
#' @param shape.by Ignored for now
#' @param span Spline span in loess function call, if \code{NULL}, no spline added
#' @param smooth Smooth the graph (similar to smoothScatter)
#' @param slot Slot to pull data from, should be one of 'counts', 'data', or 'scale.data'
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' @param plot.cor Display correlation in plot title
#' @param raster Convert points to raster format, default is \code{NULL}
#' which will automatically use raster if the number of points plotted is greater than
#' 100,000
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param jitter Jitter for easier visualization of crowded points (default is FALSE)
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 geom_smooth aes_string
#' @importFrom patchwork wrap_plots
#'
#' @export
#' @concept visualization
#'
#' @aliases GenePlot
#'
#' @examples
#' data("pbmc_small")
#' FeatureScatter(object = pbmc_small, feature1 = 'CD9', feature2 = 'CD3E')
#'
FeatureScatter <-
  function(
    object,
    feature1,
    feature2,
    cells = NULL,
    shuffle = FALSE,
    seed = 1,
    group.by = NULL,
    cols = NULL,
    pt.size = 1,
    shape.by = NULL,
    span = NULL,
    smooth = FALSE,
    combine = TRUE,
    slot = 'data',
    plot.cor = TRUE,
    raster = NULL,
    raster.dpi = c(512, 512),
    jitter = FALSE
  ){
    # 1. Define cells to use: if unspecified, plot all cells
    cells <- cells %||% get_all_cells(object)

    ## 1.2. Process shuffle settings
    if (isTRUE(x = shuffle)) {
      set.seed(seed = seed)
      cells <- sample(x = cells)
    }

    object[['ident']] <- Idents(object = object)
    group.by <- group.by %||% 'ident'

    data <-  FetchData(
      object = object,
      vars = c(feature1, feature2, group.by),
      cells = cells,
      slot = slot
    )
    if (!grepl(pattern = feature1, x = colnames(x = data)[1])) {
      stop("Feature 1 (", feature1, ") not found.", call. = FALSE)
    }
    if (!grepl(pattern = feature2, x = colnames(x = data)[2])) {
      stop("Feature 2 (", feature2, ") not found.", call. = FALSE)
    }
    data <- as.data.frame(x = data)
    feature1 <-  colnames(x = data)[1]
    feature2 <-  colnames(x = data)[2]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
    plots <- lapply(
      X = group.by,
      FUN = function(x) {
        SingleCorPlot(
          data = data[,c(feature1, feature2)],
          col.by = data[, x],
          cols = cols,
          pt.size = pt.size,
          smooth = smooth,
          legend.title = 'Identity',
          span = span,
          plot.cor = plot.cor,
          raster = raster,
          raster.dpi = raster.dpi,
          jitter = jitter
        )
      }
    )
    if (isTRUE(x = length(x = plots) == 1)) {
      return(plots[[1]])
    }
    if (isTRUE(x = combine)) {
      plots <- wrap_plots(plots, ncol = length(x = group.by))
    }
    return(plots)
  }
