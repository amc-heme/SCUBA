#' Plot feature expression by group
#'
#' Common codebase for violin, ridge plots.
#'
#' @param object Seurat object
#' @param type Plot type, choose from 'ridge', 'violin', or 'splitViolin'
#' @param features Features to plot (gene expression, metrics, PC scores,
#'  anything that can be retreived by FetchData)
#' @param idents Which classes to include in the plot (default is all)
#' @param ncol Number of columns if multiple plots are displayed
#' @param sort Sort identity classes (on the x-axis) by the average expression
#' of the attribute being potted, or, if stack is True, sort both identity
#' classes and features by hierarchical clustering
#' @param y_max Maximum y axis value
#' @param same_y_lims Set all the y-axis limits to the same values
#' @param adjust Adjust parameter for geom_violin
#' @param pt_size Point size for geom_violin
#' @param cols Colors to use for plotting
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param split_by A variable to split the plot by
#' @param log plot Y axis on log scale
#' @param slot Slot to pull expression data from (e.g. "counts" or "data")
#' @param stack Horizontally stack plots for multiple feature
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#' @param fill_by Color violins/ridges based on either 'feature' or 'ident'
#' @param flip flip plot orientation (identities on x-axis)
#' @param add_noise determine if adding a small noise for plotting
#' @param raster Convert points to raster format, default is \code{NULL} which
#' automatically rasterizes if plotting more than 100,000 cells
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom scales hue_pal
#' @importFrom ggplot2 xlab ylab
#' @importFrom patchwork wrap_plots
#'
ExPlot <-
  function(
    object,
    features,
    type = 'violin',
    idents = NULL,
    ncol = NULL,
    sort = FALSE,
    # assay = NULL,
    y_max = NULL,
    same_y_lims = FALSE,
    adjust = 1,
    cols = NULL,
    pt_size = 0,
    group_by = NULL,
    split_by = NULL,
    log = FALSE,
    slot = 'data',
    stack = FALSE,
    combine = TRUE,
    fill_by = NULL,
    flip = FALSE,
    add_noise = TRUE,
    raster = NULL
    ){
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    if (isTRUE(x = stack)) {
      if (!is.null(x = ncol)) {
        warning(
          "'ncol' is ignored with 'stack' is TRUE",
          call. = FALSE,
          immediate. = TRUE
        )
      }
      if (!is.null(x = y_max)) {
        warning(
          "'y_max' is ignored when 'stack' is TRUE",
          call. = FALSE,
          immediate. = TRUE
        )
      }
    } else {
      ncol <- ncol %||% ifelse(
        test = length(x = features) > 9,
        yes = 4,
        no = min(length(x = features), 3)
      )
    }
    data <- FetchData(object = object, vars = features, slot = slot)
    pt_size <- pt_size %||% AutoPointSize(data = object)
    # Set features equal to the returned features *with the assay keys added*
    features <- colnames(x = data)
    # Set cells based on the value of idents
    if (is.null(x = idents)) {
      cells <- colnames(x = object)
    } else {
      # ???
      #cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
      cells <- names(Idents(object)[Idents(object) %in% idents])
    }
    data <- data[cells, , drop = FALSE]
    idents <- if (is.null(x = group_by)) {
      Idents(object = object)[cells]
    } else {
      object[[group_by, drop = TRUE]][cells]
    }
    if (!is.factor(x = idents)) {
      idents <- factor(x = idents)
    }
    if (is.null(x = split_by)) {
      split <- NULL
    } else {
      # set `split` to the split by metadata column, subsetted for the cells in
      # the object, and convert to a factor if it is not already
      split <- object[[split_by, drop = TRUE]][cells]
      if (!is.factor(x = split)) {
        split <- factor(x = split)
      }
      # Set up colors for a split plot
      if (is.null(x = cols)) {
        cols <- hue_pal()(length(x = levels(x = idents)))
        cols <- Interleave(cols, InvertHex(hexadecimal = cols))
      } else if (length(x = cols) == 1 && cols == 'interaction') {
        # Splits each split_by group into secondary violins by ident class and
        # colors each ident-split_by combination separately
        # base interaction(): shows interaction of the idents and split factors
        split <- interaction(idents, split)
        cols <- hue_pal()(length(x = levels(x = idents)))
      } else {
        cols <- Col2Hex(cols)
      }
      if (length(x = cols) < length(x = levels(x = split))) {
        cols <- Interleave(cols, InvertHex(hexadecimal = cols))
      }
      cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
      names(x = cols) <- levels(x = split)
      if ((length(x = cols) > 2) & (type == "splitViolin")) {
        warning("Split violin is only supported for <3 groups, using multi-violin.")
        type <- "violin"
      }
    }
    if (same_y_lims && is.null(x = y_max)) {
      y_max <- max(data)
    }
    if (isTRUE(x = stack)) {
      return(MultiExIPlot(
        type = type,
        data = data,
        idents = idents,
        split = split,
        sort = sort,
        same.y.lims = same_y_lims,
        adjust = adjust,
        cols = cols,
        pt.size = pt_size,
        log = log,
        fill.by = fill_by,
        add.noise = add_noise,
        flip = flip
      ))
    }
    plots <- lapply(
      X = features,
      FUN = function(x) {
        return(SingleExIPlot(
          type = type,
          data = data[, x, drop = FALSE],
          idents = idents,
          split = split,
          sort = sort,
          y.max = y_max,
          adjust = adjust,
          cols = cols,
          pt.size = pt_size,
          log = log,
          add.noise = add_noise,
          raster = raster
        ))
      }
    )
    label.fxn <- switch(
      EXPR = type,
      'violin' = if (stack) {
        xlab
      } else {
        ylab
      },
      "splitViolin" = if (stack) {
        xlab
      } else {
        ylab
      },
      'ridge' = xlab,
      stop("Unknown ExIPlot type ", type, call. = FALSE)
    )
    for (i in 1:length(x = plots)) {
      key <- paste0(unlist(x = strsplit(x = features[i], split = '_'))[1], '_')
      obj <- names(x = which(x = Key(object = object) == key))
      if (length(x = obj) == 1) {
        if (inherits(x = object[[obj]], what = 'DimReduc')) {
          plots[[i]] <- plots[[i]] + label.fxn(label = 'Embeddings Value')
        } else if (inherits(x = object[[obj]], what = 'Assay')) {
          next
        } else {
          warning("Unknown object type ", class(x = object), immediate. = TRUE, call. = FALSE)
          plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
        }
      } else if (!features[i] %in% rownames(x = object)) {
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    }
    if (combine) {
      plots <- wrap_plots(plots, ncol = ncol)
      if (length(x = features) > 1) {
        plots <- plots & NoLegend()
      }
    }
    return(plots)
  }
