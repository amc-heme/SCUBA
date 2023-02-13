#------- Function copied and adapted from Seurat package ----------------------- 
# https://github.com/satijalab/seurat


#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a
#' cell and it's positioned based on the cell embeddings determined by the reduction technique. By
#' default, cells are colored by their identity class (can be changed with the group.by parameter).
#'
#' Code adapted from the Seurat package (https://github.com/satijalab/seurat/blob/master/R/visualization.R)
#' 
#' @param object Seurat object
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
#' or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
#' By default, ggplot2 assigns colors. We also include a number of palettes from the pals package.
#' See \code{\link{DiscretePalette}} for details.
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' see \code{\link{FetchData}} for more details
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells.  Only applicable if \code{raster = FALSE}.
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top)
#' @param shuffle Whether to randomly shuffle the order of points. This can be
#' useful for crowded plots if points of interest are being buried. (default is FALSE)
#' @param seed Sets the seed if randomly shuffling the order of points.
#' @param label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param label.color Sets the color of the label text
#' @param label.box Whether to put a box around the label text (geom_text vs
#' geom_label)
#' @param repel Repel labels
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#' will also resize to the size(s) passed to \code{sizes.highlight}
#' @param cols.highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells.highlight
#' @param sizes.highlight Size of highlighted cells; will repeat to the length
#' groups in cells.highlight
#' @param na.value Color value for NA points when using custom scale
#' @param ncol Number of columns for display when combining plots
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#' @param raster Convert points to raster format, default is \code{NULL} which
#' automatically rasterizes if plotting more than 100,000 cells
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom rlang !!
#' @importFrom ggplot2 facet_wrap vars sym labs
#' @importFrom patchwork wrap_plots
#'
#' @export
#' @concept visualization
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{CellSelector}, respectively.
#'
#' @aliases TSNEPlot PCAPlot ICAPlot
#' @seealso \code{\link{FeaturePlot}} \code{\link{HoverLocator}}
#' \code{\link{CellSelector}} \code{\link{FetchData}}
#'
#' @examples
#' data("pbmc_small")
#' DimPlot(object = pbmc_small)
#' DimPlot(object = pbmc_small, split.by = 'ident')
#'
DimPlot <- function(
    object,
    dims = c(1, 2),
    cells = NULL,
    cols = NULL,
    pt.size = NULL,
    reduction = NULL,
    group.by = NULL,
    split.by = NULL,
    shape.by = NULL,
    order = NULL,
    shuffle = FALSE,
    seed = 1,
    label = FALSE,
    label.size = 4,
    label.color = 'black',
    label.box = FALSE,
    repel = FALSE,
    cells.highlight = NULL,
    cols.highlight = '#DE2D26',
    sizes.highlight = 1,
    na.value = 'grey50',
    ncol = NULL,
    combine = TRUE,
    raster = NULL,
    raster.dpi = c(512, 512)
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  # Determine reduction to use (if NULL, use default))
  # %||% infix is from rlang
  reduction <- reduction %||% DefaultDimReduc(object = object)
  # Read parameter for cells to plot (defaults to all cells)
  cells <- cells %||% colnames(x = object)
  
  # Fetch data for chosen reduction from object 
  # subset for chosen cells and dims 
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  # Convert to data.frame (default is a matrix)
  data <- as.data.frame(x = data)
  # Form names for dim reduction coordinates on the x- and y- axis
  dims <- 
    paste0(
      # Use the key of the *reduction* chosen by the 
      # user (or the default reduction)
      Key(object = object[[reduction]]), 
      dims
      )
  
  # Store current ident class
  object[['ident']] <- Idents(object = object)
  # Store group by variable in `orig.groups`
  orig.groups <- group.by
  # Define the group by variable, using NULL if undefined
  group.by <- group.by %||% 'ident'
  
  # Add values of the group.by metadata variable(s) for 
  # each cell to the reduction data
  data <- 
    cbind(
      data, 
      object[[group.by]][cells, , drop = FALSE]
      )
  
  # Explicitly define group.by as the column names for each metadata variable
  # added to the reduction matrix above
  group.by <- colnames(x = data)[3:ncol(x = data)]

  # If any of the metadata columns are not factors, coerce them into factors.
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  
  # Add data for shape.by metadata variable if it exists
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  
  # Same for split.by metadata variable
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  
  # Randomly shuffle cells if specified by the user
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  
  # For each group by variable, create a DimPlot of cells grouped by that variable. 
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- Seurat::SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        # cols: a categorical pallete in this case
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        # Cells to highlight
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        raster = raster,
        raster.dpi = raster.dpi
      )
      if (label) {
        plot <- Seurat::LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by,
          box = label.box,
          color = label.color
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          ggplot2::facet_wrap(
            # rlang injection operator
            # ?`!!` for more information
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      plot <- if (is.null(x = orig.groups)) {
        plot + labs(title = NULL)
      } else {
        plot + CenterTitle()
      }
    }
  )
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    # %iff% infix (Seurat package)
    # Uses ncol if orig.groups is not NULL, otherwise uses ncol >:(
    plots <- patchwork::wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}