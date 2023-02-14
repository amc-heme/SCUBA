#------- Function copied and adapted from Seurat package -----------------------
# https://github.com/satijalab/seurat

#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a
#' cell and it's positioned based on the cell embeddings determined by the reduction technique. The function accepts both Seurat and SingleCellExperiment objects.
#' Seurat objects, cells are colored by their identity class by default, and for
#' SingleCellExperiment objects, cells are colored by the first metadata column in
#' colData(). The metadata variable used for coloring cells can be changed with the group_by parameter).
#'
#' The code for this function was from the {\link[Seurat Package]{https://github.com/satijalab/seurat/blob/master/R/visualization.R}} and adapted for use with
#' {\link[SingleCellExperiment]{https://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html}} objects.
#'
#' @param object  a Seurat object or a SingleCellExperiment object
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions i.e. c(1,2) to plot the first and second dimensions from the
#' reduction results.
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
#' or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
#' By default, ggplot2 assigns colors. We also include a number of palettes from the pals package.
#' See \code{\link{DiscretePalette}} for details.
#' @param pt_size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group_by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split_by Name of a metadata column to split plot by;
#' see \code{\link{FetchData}} for more details
#' @param shape_by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells.  Only applicable if \code{raster = FALSE}.
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top)
#' @param shuffle Whether to randomly shuffle the order of points. This can be
#' useful for crowded plots if points of interest are being buried. (default is FALSE)
#' @param seed Sets the seed if randomly shuffling the order of points.
#' @param label Whether to label the clusters
#' @param label_size Sets size of labels
#' @param label_color Sets the color of the label text
#' @param label_box Whether to put a box around the label text (geom_text vs
#' geom_label)
#' @param repel Repel labels
#' @param cells_highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols_highlight} and other cells black (white if dark.theme = TRUE);
#' will also resize to the size(s) passed to \code{sizes_highlight}
#' @param cols_highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells_highlight
#' @param sizes_highlight Size of highlighted cells; will repeat to the length
#' groups in cells_highlight
#' @param na_value Color value for NA points when using custom scale
#' @param ncol Number of columns for display when combining plots
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#' @param raster Convert points to raster format, default is \code{NULL} which
#' automatically rasterizes if plotting more than 100,000 cells
#' @param raster_dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom rlang !!
#' @importFrom ggplot2 facet_wrap vars sym labs
#' @importFrom patchwork wrap_plots
#' @importFrom SingleCellExperiment reducedDimNames reducedDims colData
#' @importFrom Seurat Embeddings SingleDimPlot LabelClusters
#'
#' @export
#'
DimPlot <- function(
    object,
    dims = c(1, 2),
    cells = NULL,
    cols = NULL,
    pt_size = NULL,
    reduction = NULL,
    group_by = NULL,
    split_by = NULL,
    shape_by = NULL,
    order = NULL,
    shuffle = FALSE,
    seed = 1,
    label = FALSE,
    label_size = 4,
    label_color = 'black',
    label_box = FALSE,
    repel = FALSE,
    cells_highlight = NULL,
    cols_highlight = '#DE2D26',
    sizes_highlight = 1,
    na_value = 'grey50',
    ncol = NULL,
    combine = TRUE,
    raster = NULL,
    raster_dpi = c(512, 512)
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }

  # Compile data for plotting: methods depend on object type
  if (is(object, "SingleCellExperiment")){
    # For SingleCellExperiment objects
    # 1. Define reduction (defaults to the first reduction stored for sce objects)
    # Uses rlang %||% infix
    reduction <- reduction %||% SingleCellExperiment::reducedDimNames(object)[1]
    # 2. Define cells to include in plot
    ## Same as for seurat object ##
    cells <- cells %||% BiocGenerics::colnames(object)

    # 3. Fetch dimensional reduction data from object
    data <- SingleCellExperiment::reducedDims(object)[[reduction]][cells, dims]
    data <- BiocGenerics::as.data.frame(data)

    # 4. Fetch names of dimensions to plot
    # For SCE objects, use the column names for the requested dim indices
    # (there is no `Key()` method for SingleCellExperiment objects)
    dim_names <- colnames(data)[dims]

    # 5. Process group by selection
    # There is no "ident" property for SingleCellExperiment objects, so defaults
    # will be chosen based on the first metadata column, rather than the current ident class
    # Store original entry for group_by
    orig_groups <- group_by

    if (is.null(group_by)){
      warning("group_by is NULL. Defaulting to the first column in colData.")
      group_by <- names(colData(object))[1]
    }

    # 6. Bind group by metadata to the table of reduction coordinates
    data <-
      cbind(
        data,
        # Subsets for selected cells, and the names of entered group_by columns
        colData(object)[cells, group_by, drop = FALSE]
      )

    # 7. Define group_by variables to iterate through
    group_by_cols <-
      colnames(data)[3:ncol(data)]

    # 8. Convert group by columns to factors if they are not already
    for (group in group_by_cols) {
      if (!is.factor(data[, group])) {
        data[, group] <- factor(data[, group])
      }
    }

    # 9. Add shape_by data if it exists
    if (!is.null(x = shape_by)) {
      data[, shape_by] <- colData(object)[cells, shape_by, drop = TRUE]
    }

    # 10. Same for split_by data
    if (!is.null(x = split_by)) {
      data[, split_by] <- colData(object)[cells, split_by, drop = TRUE]
    }

    # 11. If sufffle is TRUE, randomly shuffle cells
    if (isTRUE(shuffle)) {
      set.seed(seed = seed)
      data <- data[sample(x = 1:nrow(x = data)), ]
    }
  } else if (is(object, "Seurat")){
    # For Seurat Objects
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
    dim_names <-
      paste0(
        # Use the key of the *reduction* chosen by the
        # user (or the default reduction)
        Key(object = object[[reduction]]),
        dims
      )

    # Set group by to the current ident class if it is NULL
    # Store current ident class
    object[['ident']] <- Idents(object = object)
    # Store group by variable in `orig_groups`
    orig_groups <- group_by
    group_by <- group_by %||% 'ident'

    # Add values of the group_by metadata variable(s) for
    # each cell to the reduction data
    data <-
      cbind(
        data,
        object[[group_by]][cells, , drop = FALSE]
      )

    # Explicitly define group_by as the column names for each metadata variable
    # added to the reduction matrix above
    group_by_cols <- colnames(x = data)[3:ncol(x = data)]

    # If any of the metadata columns are not factors, coerce them into factors.
    for (group in group_by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }

    # Add data for shape_by metadata variable if it exists
    if (!is.null(x = shape_by)) {
      data[, shape_by] <- object[[shape_by, drop = TRUE]]
    }

    # Same for split_by metadata variable
    if (!is.null(x = split_by)) {
      data[, split_by] <- object[[split_by, drop = TRUE]]
    }

    # Randomly shuffle cells if specified by the user
    if (isTRUE(x = shuffle)) {
      set.seed(seed = seed)
      data <- data[sample(x = 1:nrow(x = data)), ]
    }
  } else {
    # Throw an error for unsupported data types
    stop("Object entered is not a SinglecellExperiment or a Seurat object.")
  }

  # For each group by variable, create a DimPlot of cells grouped by that variable.
  plots <- lapply(
    X = group_by_cols,
    FUN = function(group) {
      plot <- Seurat::SingleDimPlot(
        data = data[, c(dim_names, group, split_by, shape_by)],
        dims = dim_names,
        col.by = group,
        # cols: a categorical pallete in this case
        cols = cols,
        pt.size = pt_size,
        shape.by = shape_by,
        order = order,
        label = FALSE,
        # Cells to highlight
        cells.highlight = cells_highlight,
        cols.highlight = cols_highlight,
        sizes.highlight = sizes_highlight,
        na.value = na_value,
        raster = raster,
        raster.dpi = raster_dpi
      )
      if (label) {
        plot <-
          Seurat::LabelClusters(
            plot = plot,
            id = group,
            repel = repel,
            size = label_size,
            split.by = split_by,
            box = label_box,
            color = label_color
          )
      }
      if (!is.null(x = split_by)) {
        plot <-
          plot +
          Seurat:::FacetTheme() +
          ggplot2::facet_wrap(
            # rlang injection operator
            # ?`!!` for more information
            facets = vars(!!sym(x = split_by)),
            ncol = if (length(x = group_by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split_by]))
            } else {
              ncol
            }
          )
      }
      plot <- if (is.null(x = orig_groups)) {
        plot + labs(title = NULL)
      } else {
        plot + CenterTitle()
      }
    }
  )
  if (!is.null(x = split_by)) {
    ncol <- 1
  }
  if (combine) {
    plots <-
      patchwork::wrap_plots(
        plots,
        # %iff% infix (Seurat package)
        # Uses ncol if orig_groups is not NULL, otherwise uses NULL
        ncol = orig_groups %iff% ncol
        )
  }
  return(plots)
}
