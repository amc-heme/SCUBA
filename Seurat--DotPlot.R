#' Dot plot visualization
#'
#' Intuitive way of visualizing how feature expression changes across different
#' identity classes (clusters). The size of the dot encodes the percentage of
#' cells within a class, while the color encodes the AverageExpression level
#' across all cells within a class (blue is high).
#'
#' @param object Seurat object
#' @param assay Name of assay to use, defaults to the active assay
#' @param features Input vector of features, or named list of feature vectors
#' if feature-grouped panels are desired (replicates the functionality of the
#' old SplitDotPlotGG)
#' @param cols Colors to plot: the name of a palette from
#' \code{RColorBrewer::brewer.pal.info}, a pair of colors defining a gradient,
#' or 3+ colors defining multiple gradients (if split.by is set)
#' @param col.min Minimum scaled average expression threshold (everything
#' smaller will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger
#' will be set to this)
#' @param dot.min The fraction of cells at which to draw the smallest dot
#' (default is 0). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param dot.scale Scale the size of the points, similar to cex
#' @param idents Identity classes to include in plot (default is all)
#' @param group.by Factor to group the cells by
#' @param split.by Factor to split the groups by (replicates the functionality
#' of the old SplitDotPlotGG);
#' see \code{\link{FetchData}} for more details
#' @param cluster.idents Whether to order identities by hierarchical clusters
#' based on given features, default is FALSE
#' @param scale Determine whether the data is scaled, TRUE for default
#' @param scale.by Scale the size of the points by 'size' or by 'radius'
#' @param scale.min Set lower limit for scaling, use NA for default
#' @param scale.max Set upper limit for scaling, use NA for default
#'
#' @return A ggplot object
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point scale_size scale_radius
#' theme element_blank labs scale_color_identity scale_color_distiller
#' scale_color_gradient guides guide_legend guide_colorbar
#' facet_grid unit
#' @importFrom scattermore geom_scattermore
#' @importFrom stats dist hclust
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @export
#' @concept visualization
#'
#' @aliases SplitDotPlotGG
#' @seealso \code{RColorBrewer::brewer.pal.info}
#'
#' @examples
#' data("pbmc_small")
#' cd_genes <- c("CD247", "CD3E", "CD9")
#' DotPlot(object = pbmc_small, features = cd_genes)
#' pbmc_small[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = pbmc_small), replace = TRUE)
#' DotPlot(object = pbmc_small, features = cd_genes, split.by = 'groups')
#'
DotPlotTest <- function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA
) {
  # Determine the assay to use for gene expression data (if not provided)
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  # Determine whether to use
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))

  # Determine which ggplot2 scale_* function to use for scaling the dots
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )

  # Split features into groups if a named list is passed to `features`
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }

  # CellsByIdentities is a SeuratObject method
  # Prints cell names by each level in the ident class
  # i.e. if B cells and BM monocytes are passed to `idents`, cell names that
  # are listed as B cells and BM monocytes (according to the metadata column
  # set as the ident class) are returned
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

  # Fetch data for each feature
  data.features <- FetchData(object = object, vars = features, cells = cells)

  # Add metadata to table
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  # Convert metadata column to a factor if it is not already
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  # Store the factor levels in `id.levels` and convert back to a vector
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)

  # Handle split.by metadata if defined
  if (!is.null(x = split.by)) {
    # pull split.by metadata
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    # Assign colors to each unique split
    # Only works when an RColorBrewer palette is passed
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    # Combine
    data.features$id <- paste(data.features$id, splits, sep = '_')
    print("unique(data.features$id)")
    print(unique(data.features$id))
    unique.splits <- unique(x = splits)
    # Create labels for each group.by and split.by level with an underscore
    # (each combo reads group.by[i]_split.by[j] for all combinations of i
    # group.by and j split.by levels), *for which there are cells*
    id.levels <-
      paste0(
        rep(x = id.levels, each = length(x = unique.splits)),
        "_",
        rep(x = unique(x = splits), times = length(x = id.levels))
        )
    print("id.levels")
    print(id.levels)
  }

  # Compute statistics for each group #
  # Store results in a list
  # Groups may be levels of the group.by variable, or group.by-split.by
  # combinations if split.by is enabled
  # If there are multiple features, each element will be a vector with the
  # expression of each feature for the given group
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      # Pull data for computing statistics
      # Exclude last column of data.features ("id" column)
      data.use <-
        data.features[
          data.features$id == ident,
          1:(ncol(x = data.features) - 1),
          drop = FALSE
          ]

      # Average expression
      avg.exp <-
        apply(
          X = data.use,
          # Apply function by column
          MARGIN = 2,
          FUN = function(x) {
            return(mean(x = expm1(x = x)))
            }
          )

      # Percent of cells in group expressing feature
      pct.exp <-
        apply(
          X = data.use,
          MARGIN = 2,
          # Uses Seurat::PercentAbove()
          FUN = PercentAbove,
          # "threshold" is an argument of PercentAbove
          threshold = 0
          )

      return(
        list(
          avg.exp = avg.exp,
          pct.exp = pct.exp
          )
        )
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)

  # Cluster groups if `cluster.idents` is TRUE
  # Clustering is performed based on expression statistics for the currently
  # plotted set of features
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }

  # Construct dataframe from list of expression statistics
  # Return is a list of dataframes, one dataframe per group
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      # Add feature name as a column
      data.use$features.plot <- rownames(x = data.use)
      # Add column for the group ID
      data.use$id <- x
      return(data.use)
    }
  )
  # Combine list of dataframes produced above
  data.plot <- do.call(what = 'rbind', args = data.plot)

  # Convert the new "id" column in the dataframe to a factor
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }

  # Compute number of groups, and warn if the number of groups is 1, or less than 5
  ngroup <- length(x = levels(x = data.plot$id))

  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  # Scale expression across features
  avg.exp.scaled <-
    sapply(
      X = unique(x = data.plot$features.plot),
      FUN = function(x) {
        # Pull average expression column and subset to feature x
        data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
        if (scale) {
          print("Data before scaling")
          print(data.use)
          # Scale data (z-score)
          data.use <- scale(x = data.use)
          print("Data after scaling")
          print(data.use)
          # Any values with a z-score outside of col.min and col.max will be
          # normalized to those values
          data.use <- MinMax(data = data.use, min = col.min, max = col.max)
          print("Data after min/max transformation")
          print(data.use)
        } else {
          # If `scale` is FALSE, perform a log(x+1) transformation
          data.use <- log1p(x = data.use)
          }
        return(data.use)
        }
      )

  print("avg.exp.sacaled before conversion to vector")
  print(avg.exp.scaled)

  # Transfrom average expression to a vector
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))

  print("avg.exp.sacaled after conversion to vector")
  print(avg.exp.scaled)

  # "Bin" the average expression into 20 bins for split plots
  if (split.colors) {
    # Max value will be 20, min value will be 1, and all other values will be
    # binned based on their value
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }

  # Add scaled data to the data.plot dataframe and convert to factor
  data.plot$avg.exp.scaled <- avg.exp.scaled
  print("Plot data with scaled expression column")
  print(data.plot)
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  # Add NAs for percent expression values beneath dot.min
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  # Convert percent expression from fraction to percentage
  data.plot$pct.exp <- data.plot$pct.exp * 100

  print("Value of data.plot$id")
  print(data.plot$id)

  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      # Return template: character vector of the same length as X
      FUN.VALUE = character(length = 1L),
      # Parameters of gsub()
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    print("Value of splits.use")
    print(splits.use)

    # Color for each dot: use colorRampPalette to determine the color (using
    # a continuous scale between grey and the color entered, with 20 steps)
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        print("Value")
        print(value)
        print("Color")
        print(colorRampPalette(colors = c('grey', color))(20)[value])
        return(
          colorRampPalette(colors = c('grey', color))(20)[value]
          )
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )

    print("Split data after operations")
    print(head(data.plot))
  }

  # Color scale: if there are splits, use 20-level binned scale based on average
  # expression (feature-wide).
  # If not (default behavior), apply colors using ggplot2 scales
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')

  # Apply min/max cutoffs for *percent expression*
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }

  # Add data for feature groups to the table (if defined by the user)
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }

  plot <-
    ggplot(
      data = data.plot,
      # aes_string is deprecated in the current ggplot2 version
      mapping = aes_string(x = 'features.plot', y = 'id')
      ) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot()

  # Create facets for feature groups, if present
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank()
    )
  }
  if (split.colors) {
    # scale_color_identity applies colors that were defined based on the
    # 20-bin scale
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  return(plot)
}
