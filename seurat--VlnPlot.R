#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @inheritParams RidgePlot
#' @param pt.size Point size for geom_violin
#' @param split.by A variable to split the violin plots by,
#' @param split.plot  plot each group of the split violin plots by multiple or
#' single violin shapes.
#' @param adjust Adjust parameter for geom_violin
#' @param flip flip plot orientation (identities on x-axis)
#' @param add.noise determine if adding a small noise for plotting
#' @param raster Convert points to raster format. Requires 'ggrastr' to be installed.
# default is \code{NULL} which automatically rasterizes if ggrastr is installed and
# number of points exceed 100,000.
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @export
#' @concept visualization
#'
#' @seealso \code{\link{FetchData}}
#'
#' @examples
#' data("pbmc_small")
#' VlnPlot(object = pbmc_small, features = 'PC_1')
#' VlnPlot(object = pbmc_small, features = 'LYZ', split.by = 'groups')
#'
VlnPlot <- function(
    object,
    features,
    cols = NULL,
    pt.size = NULL,
    idents = NULL,
    sort = FALSE,
    assay = NULL,
    group.by = NULL,
    split.by = NULL,
    adjust = 1,
    y.max = NULL,
    same.y.lims = FALSE,
    log = FALSE,
    ncol = NULL,
    slot = 'data',
    split.plot = FALSE,
    stack = FALSE,
    combine = TRUE,
    fill.by = 'feature',
    flip = FALSE,
    add.noise = TRUE,
    raster = NULL
) {
  if (
    !is.null(x = split.by) &
    getOption(x = 'Seurat.warn.vlnplot.split', default = TRUE)
  ) {
    message(
      "The default behaviour of split.by has changed.\n",
      "Separate violin plots are now plotted side-by-side.\n",
      "To restore the old behaviour of a single split violin,\n",
      "set split.plot = TRUE.
      \nThis message will be shown once per session."
    )
    options(Seurat.warn.vlnplot.split = FALSE)
  }
  return(ExIPlot(
    object = object,
    type = ifelse(test = split.plot, yes = 'splitViolin', no = 'violin'),
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    adjust = adjust,
    pt.size = pt.size,
    cols = cols,
    group.by = group.by,
    split.by = split.by,
    log = log,
    slot = slot,
    stack = stack,
    combine = combine,
    fill.by = fill.by,
    flip = flip,
    add.noise = add.noise,
    raster = raster
  ))
}

# Plot feature expression by identity
#
# Basically combines the codebase for VlnPlot and RidgePlot
#
# @param object Seurat object
# @param type Plot type, choose from 'ridge', 'violin', or 'splitViolin'
# @param features Features to plot (gene expression, metrics, PC scores,
# anything that can be retreived by FetchData)
# @param idents Which classes to include in the plot (default is all)
# @param ncol Number of columns if multiple plots are displayed
# @param sort Sort identity classes (on the x-axis) by the average expression of the attribute being potted,
# or, if stack is True, sort both identity classes and features by hierarchical clustering
# @param y.max Maximum y axis value
# @param same.y.lims Set all the y-axis limits to the same values
# @param adjust Adjust parameter for geom_violin
# @param pt.size Point size for geom_violin
# @param cols Colors to use for plotting
# @param group.by Group (color) cells in different ways (for example, orig.ident)
# @param split.by A variable to split the plot by
# @param log plot Y axis on log scale
# @param slot Slot to pull expression data from (e.g. "counts" or "data")
# @param stack Horizontally stack plots for multiple feature
# @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
# ggplot object. If \code{FALSE}, return a list of ggplot objects
# @param fill.by Color violins/ridges based on either 'feature' or 'ident'
# @param flip flip plot orientation (identities on x-axis)
# @param add.noise determine if adding a small noise for plotting
# @param raster Convert points to raster format, default is \code{NULL} which
# automatically rasterizes if plotting more than 100,000 cells
#
# @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
# \code{combine = TRUE}; otherwise, a list of ggplot objects
#
#' @importFrom scales hue_pal
#' @importFrom ggplot2 xlab ylab
#' @importFrom patchwork wrap_plots
#
ExIPlot <- function(
    object,
    features,
    type = 'violin',
    idents = NULL,
    ncol = NULL,
    sort = FALSE,
    assay = NULL,
    y.max = NULL,
    same.y.lims = FALSE,
    adjust = 1,
    cols = NULL,
    pt.size = 0,
    group.by = NULL,
    split.by = NULL,
    log = FALSE,
    slot = 'data',
    stack = FALSE,
    combine = TRUE,
    fill.by = NULL,
    flip = FALSE,
    add.noise = TRUE,
    raster = NULL
) {
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
    if (!is.null(x = y.max)) {
      warning(
        "'y.max' is ignored when 'stack' is TRUE",
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
  pt.size <- pt.size %||% AutoPointSize(data = object)
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
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    # set `split` to the split by metadata column, subsetted for the cells in
    # the object, and convert to a factor if it is not already
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    # Set up colors for a split plot
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    } else if (length(x = cols) == 1 && cols == 'interaction') {
      # Splits each split.by group into secondary violins by ident class and
      # colors each ident-split.by combination separately
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
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  if (isTRUE(x = stack)) {
    return(MultiExIPlot(
      type = type,
      data = data,
      idents = idents,
      split = split,
      sort = sort,
      same.y.lims = same.y.lims,
      adjust = adjust,
      cols = cols,
      pt.size = pt.size,
      log = log,
      fill.by = fill.by,
      add.noise = add.noise,
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
        y.max = y.max,
        adjust = adjust,
        cols = cols,
        pt.size = pt.size,
        log = log,
        add.noise = add.noise,
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

#' Plot a single expression by identity on a plot
#'
#' @param data Data to plot
#' @param idents Idents to use
#' @param split Use a split violin plot
#' @param type Make either a \dQuote{ridge} or \dQuote{violin} plot
#' @param sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum Y value to plot
#' @param adjust Adjust parameter for geom_violin
#' @param pt.size Size of points for violin plots
#' @param cols Colors to use for plotting
#' @param seed.use Random seed to use. If NULL, don't set a seed
#' @param log plot Y axis on log10 scale
#' @param add.noise determine if adding small noise for plotting
#' @param raster Convert points to raster format. Requires 'ggrastr' to be installed.
#' default is \code{NULL} which automatically rasterizes if ggrastr is installed and
#' number of points exceed 100,000.
#'
#' @return A ggplot-based Expression-by-Identity plot
#'
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes_string theme labs geom_violin geom_jitter
#' ylim position_jitterdodge scale_fill_manual scale_y_log10 scale_x_log10
#' scale_y_discrete scale_x_continuous waiver
#' @importFrom cowplot theme_cowplot
#'
#' @keywords internal
#' @export
#'
SingleExIPlot <- function(
    data,
    idents,
    split = NULL,
    type = 'violin',
    sort = FALSE,
    y.max = NULL,
    adjust = 1,
    pt.size = 0,
    cols = NULL,
    seed.use = 42,
    log = FALSE,
    add.noise = TRUE,
    raster = NULL
) {
  if (!is.null(x = raster) && isTRUE(x = raster)){
    if (!PackageCheck('ggrastr', error = FALSE)) {
      stop("Please install ggrastr from CRAN to enable rasterization.")
    }
  }
  if (PackageCheck('ggrastr', error = FALSE)) {
    # Set rasterization to true if ggrastr is installed and
    # number of points exceeds 100,000
    if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
      message("Rasterizing points since number of points exceeds 100,000.",
              "\nTo disable this behavior set `raster=FALSE`")
      # change raster to TRUE
      raster <- TRUE
    }
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
  }
  feature <- colnames(x = data)
  data$ident <- idents
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(
        x = tapply(
          X = data[, feature],
          INDEX = data$ident,
          FUN = mean
        ),
        decreasing = grepl(pattern = paste0('^', tolower(x = sort)), x = 'decreasing')
      )))
    )
  }
  if (log) {
    noise <- rnorm(n = length(x = data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- rnorm(n = length(x = data[, feature])) / 100000
  }
  if (!add.noise) {
    noise <-  noise * 0
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else{
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- 'Expression Level'
  y.max <- y.max %||% max(data[, feature][is.finite(x = data[, feature])])
  if (type == 'violin' && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_violin
    fill <- 'split'
  } else if (type == 'splitViolin' && !is.null(x = split )) {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- 'split'
    type <- 'violin'
  } else {
    vln.geom <- geom_violin
    fill <- 'ident'
  }
  switch(
    EXPR = type,
    'violin' = {
      x <- 'ident'
      y <- paste0("`", feature, "`")
      xlab <- 'Identity'
      ylab <- axis.label
      geom <- list(
        vln.geom(scale = 'width', adjust = adjust, trim = TRUE),
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      if (is.null(x = split)) {
        if (isTRUE(x = raster)) {
          jitter <- ggrastr::rasterize(geom_jitter(height = 0, size = pt.size, show.legend = FALSE))
        } else {
          jitter <- geom_jitter(height = 0, size = pt.size, show.legend = FALSE)
        }
      } else {
        if (isTRUE(x = raster)) {
          jitter <- ggrastr::rasterize(geom_jitter(
            position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
            size = pt.size,
            show.legend = FALSE
          ))
        } else {
          jitter <- geom_jitter(
            position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
            size = pt.size,
            show.legend = FALSE
          )
        }
      }
      log.scale <- scale_y_log10()
      axis.scale <- ylim
    },
    'ridge' = {
      x <- paste0("`", feature, "`")
      y <- 'ident'
      xlab <- axis.label
      ylab <- 'Identity'
      geom <- list(
        geom_density_ridges(scale = 4),
        theme_ridges(),
        scale_y_discrete(expand = c(0.01, 0)),
        scale_x_continuous(expand = c(0, 0))
      )
      jitter <- geom_jitter(width = 0, size = pt.size, show.legend = FALSE)
      log.scale <- scale_x_log10()
      axis.scale <- function(...) {
        invisible(x = NULL)
      }
    },
    stop("Unknown plot type: ", type)
  )
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]
  ) +
    labs(x = xlab, y = ylab, title = feature, fill = NULL) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))
  plot <- do.call(what = '+', args = list(plot, geom))
  plot <- plot + if (log) {
    log.scale
  } else {
    axis.scale(min(data[, feature]), y.max)
  }
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      } else {
        unlist(x = lapply(
          X = idents,
          FUN = function(pattern, x) {
            x.mod <- gsub(
              pattern = paste0(pattern, '.'),
              replacement = paste0(pattern, ': '),
              x = x,
              fixed = TRUE
            )
            x.keep <- grep(pattern = ': ', x = x.mod, fixed = TRUE)
            x.return <- x.mod[x.keep]
            names(x = x.return) <- x[x.keep]
            return(x.return)
          },
          x = unique(x = as.vector(x = data$split))
        ))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    } else {
      labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  return(plot)
}

#' Plot expression of multiple features by identity on a plot
#'
#' @param data Data to plot
#' @param idents Idents to use
#' @param type Make either a 'ridge' or 'violin' plot
#' @param sort Sort identity classes and features based on hierarchical clustering
#' @param same.y.lims Indicates whether to use the same ylim for each feature
#' @param adjust Adjust parameter for geom_violin
#' @param cols Colors to use for plotting
#' @param log plot Y axis on log scale
#' @param fill.by Color violins/ridges based on either 'feature' or 'ident'
#' @param seed.use Random seed to use. If NULL, don't set a seed
#' @param flip flip plot orientation (identities on x-axis)
#'
#' @return A ggplot-based Expression-by-Identity plot
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom utils globalVariables
#' @importFrom stats rnorm dist hclust
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes_string facet_grid theme labs geom_rect
#' geom_violin geom_jitter ylim position_jitterdodge scale_fill_manual
#' scale_y_log10 scale_x_log10 scale_y_discrete scale_x_continuous
#' scale_y_continuous waiver
#'
MultiExIPlot <- function(
    data,
    idents,
    split = NULL,
    type = 'violin',
    sort = FALSE,
    same.y.lims = same.y.lims,
    adjust = 1,
    pt.size = 0,
    cols = NULL,
    seed.use = 42,
    log = FALSE,
    fill.by = NULL,
    add.noise = TRUE,
    flip = NULL
) {
  if (!(fill.by %in% c("feature", "ident"))) {
    stop("`fill.by` must be either `feature` or `ident`")
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) < 2) {
    stop("MultiExIPlot requires a data frame with >1 column")
  }
  data <- Melt(x = data)
  data <- data.frame(
    feature = data$cols,
    expression = data$vals,
    ident = rep_len(x = idents, length.out = nrow(x = data))
  )
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$feature <- as.vector(x = data$feature)
    data$ident <- as.vector(x = data$ident)
    # build matrix of average expression (#-features by #-idents), lexical ordering
    avgs.matrix <- sapply(
      X = split(x = data, f = data$ident),
      FUN = function(df) {
        return(tapply(
          X = df$expression,
          INDEX = df$feature,
          FUN = mean
        ))
      }
    )
    idents.order <- hclust(d = dist(x = t(x = L2Norm(mat = avgs.matrix, MARGIN = 2))))$order
    avgs.matrix <- avgs.matrix[,idents.order]
    avgs.matrix <- L2Norm(mat = avgs.matrix, MARGIN = 1)
    # order feature clusters by position of their "rank-1 idents"
    position <- apply(X = avgs.matrix, MARGIN = 1, FUN = which.max)
    mat <- hclust(d = dist(x = avgs.matrix))$merge
    orderings <- list()
    for (i in 1:nrow(mat)) {
      x <- if (mat[i,1] < 0) -mat[i,1] else orderings[[mat[i,1]]]
      y <- if (mat[i,2] < 0) -mat[i,2] else orderings[[mat[i,2]]]
      x.pos <- min(x = position[x])
      y.pos <- min(x = position[y])
      orderings[[i]] <- if (x.pos < y.pos) {
        c(x, y)
      } else {
        c(y, x)
      }
    }
    features.order <- orderings[[length(x = orderings)]]
    data$feature <- factor(
      x = data$feature,
      levels = unique(x = sort(x = data$feature))[features.order]
    )
    data$ident <- factor(
      x = data$ident,
      levels = unique(x = sort(x = data$ident))[rev(x = idents.order)]
    )
  } else {
    data$feature <- factor(x = data$feature, levels = unique(x = data$feature))
  }
  if (log) {
    noise <- rnorm(n = nrow(x = data)) / 200
    data$expression <- data$expression + 1
  } else {
    noise <- rnorm(n = nrow(x = data)) / 100000
  }
  if (!add.noise) {
    noise <- noise*0
  }
  for (f in unique(x = data$feature)) {
    if (all(data$expression[(data$feature == f)] == data$expression[(data$feature == f)][1])) {
      warning(
        "All cells have the same value of ",
        f,
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      data$expression[(data$feature == f)] <- data$expression[(data$feature == f)] + noise[(data$feature == f)]
    }
  }
  if (type == 'violin' && !is.null(x = split)) {
    data$split <- rep_len(x = split, length.out = nrow(data))
    vln.geom <- geom_violin
    fill.by <- 'split'
  } else if (type == 'splitViolin' && !is.null(x = split)) {
    data$split <- rep_len(x = split, length.out = nrow(data))
    vln.geom <- geom_split_violin
    fill.by <- 'split'
    type <- 'violin'
  } else {
    vln.geom <- geom_violin
  }
  switch(
    EXPR = type,
    'violin' = {
      geom <- list(vln.geom(scale = 'width', adjust = adjust, trim = TRUE))
    },
    'ridge' = {
      geom <- list(
        geom_density_ridges(scale = 4),
        theme_ridges(),
        scale_y_discrete(expand = c(0.01, 0))
      )
    },
    stop("Unknown plot type: ", type)
  )
  if (flip) {
    x <- 'ident'
    x.label <- 'Identity'
    y <- 'expression'
    y.label <- 'Expression Level'
  } else {
    y <- 'ident'
    y.label <- 'Identity'
    x <- 'expression'
    x.label <- 'Expression Level'
  }
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = x, y = y, fill = fill.by)[c(2, 3, 1)]
  ) +
    labs(x = x.label, y = y.label, fill = NULL) +
    theme_cowplot()
  plot <- do.call(what = '+', args = list(plot, geom))
  if (flip) {
    plot <- plot +
      scale_y_continuous(
        expand = c(0, 0),
        labels = function(x) c(rep(x = '', times = length(x)-2), x[length(x) - 1], '')) +
      facet_grid(feature ~ ., scales = (if (same.y.lims) 'fixed' else 'free')) +
      FacetTheme(
        panel.spacing = unit(0, 'lines'),
        panel.background = element_rect(fill = NA, color = "black"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y.right = element_text(angle = 0))
  } else {
    plot <- plot +
      scale_x_continuous(
        expand = c(0, 0),
        labels = function(x) c(rep(x = '', times = length(x)-2), x[length(x) - 1], '')) +
      facet_grid(. ~ feature, scales = (if (same.y.lims) 'fixed' else 'free')) +
      FacetTheme(
        panel.spacing = unit(0, 'lines'),
        panel.background = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(size = 7),
        strip.text.x = element_text(angle = -90))
  }
  if (log) {
    plot <- plot + scale_x_log10()
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      } else {
        unlist(x = lapply(
          X = idents,
          FUN = function(pattern, x) {
            x.mod <- gsub(
              pattern = paste0(pattern, '.'),
              replacement = paste0(pattern, ': '),
              x = x,
              fixed = TRUE
            )
            x.keep <- grep(pattern = ': ', x = x.mod, fixed = TRUE)
            x.return <- x.mod[x.keep]
            names(x = x.return) <- x[x.keep]
            return(x.return)
          },
          x = unique(x = as.vector(x = data$split))
        ))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    } else {
      labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  return(plot)
}
