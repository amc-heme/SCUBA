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
#' \code{RColorBrewer::brewer_pal_info}, a pair of colors defining a gradient,
#' or 3+ colors defining multiple gradients (if split_by is set)
#' @param col_min Minimum scaled average expression threshold (everything
#' smaller will be set to this)
#' @param col_max Maximum scaled average expression threshold (everything larger
#' will be set to this)
#' @param dot_min The fraction of cells at which to draw the smallest dot
#' (default is 0). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param dot_scale Scale the size of the points, similar to cex
#' @param idents Identity classes to include in plot (default is all)
#' @param group_by Factor to group the cells by
#' @param split_by Factor to split the groups by (replicates the functionality
#' of the old SplitDotPlotGG);
#' see \code{\link{FetchData}} for more details
#' @param cluster_idents Whether to order identities by hierarchical clusters
#' based on given features, default is FALSE
#' @param scale Determine whether the data is scaled, TRUE for default
#' @param scale_by Scale the size of the points by 'size' or by 'radius'
#' @param scale_min Set lower limit for scaling, use NA for default
#' @param scale_max Set upper limit for scaling, use NA for default
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
#' @importFrom RColorBrewer brewer_pal_info
#'
#' @export
#' @concept visualization
#'
#' @seealso \code{RColorBrewer::brewer_pal_info}
#'
#' @examples
#' data("pbmc_small")
#' cd_genes <- c("CD247", "CD3E", "CD9")
#' DotPlot(object = pbmc_small, features = cd_genes)
#' pbmc_small[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = pbmc_small), replace = TRUE)
#' DotPlot(object = pbmc_small, features = cd_genes, split_by = 'groups')
#'
plot_dot <-
  function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "blue"),
    col_min = -2.5,
    col_max = 2.5,
    dot_min = 0,
    dot_scale = 6,
    idents = NULL,
    group_by = NULL,
    split_by = NULL,
    cluster_idents = FALSE,
    scale = TRUE,
    scale_by = 'radius',
    scale_min = NA,
    scale_max = NA
    ){
    # Determine the assay to use for gene expression data (if not provided)
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    # Determine whether to use
    split_colors <-
      !is.null(x = split_by) && !any(cols %in% rownames(x = brewer_pal_info))

    # Determine which ggplot2 scale_* function to use for scaling the dots
    scale_func <- switch(
      EXPR = scale_by,
      'size' = scale_size,
      'radius' = scale_radius,
      stop("'scale_by' must be either 'size' or 'radius'")
    )

    # Split features into groups if a named list is passed to `features`
    feature_groups <- NULL
    if (is.list(features) | any(!is.na(names(features)))) {
      feature_groups <- unlist(x = sapply(
        X = 1:length(features),
        FUN = function(x) {
          return(rep(x = names(x = features)[x], each = length(features[[x]])))
        }
      ))
      if (any(is.na(x = feature_groups))) {
        warning(
          "Some feature groups are unnamed.",
          call. = FALSE,
          immediate. = TRUE
        )
      }
      features <- unlist(x = features)
      names(x = feature_groups) <- features
    }

    # CellsByIdentities is a SeuratObject method
    # Prints cell names by each level in the ident class
    # i.e. if B cells and BM monocytes are passed to `idents`, cell names that
    # are listed as B cells and BM monocytes (according to the metadata column
    # set as the ident class) are returned
    cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

    # Fetch data for each feature
    expr_data <-
      FetchData(
        object = object,
        vars = features,
        cells = cells
        )

    # Add metadata to table
    expr_data$id <- if (is.null(x = group_by)) {
      Idents(object = object)[cells, drop = TRUE]
    } else {
      object[[group_by, drop = TRUE]][cells, drop = TRUE]
    }
    # Convert metadata column to a factor if it is not already
    if (!is.factor(x = expr_data$id)) {
      expr_data$id <- factor(x = expr_data$id)
    }
    # Store the factor levels in `id_levels` and convert back to a vector
    id_levels <- levels(x = expr_data$id)
    expr_data$id <- as.vector(x = expr_data$id)

    # Handle split_by metadata if defined
    if (!is.null(x = split_by)) {
      # pull split_by metadata
      splits <- object[[split_by, drop = TRUE]][cells, drop = TRUE]
      # Assign colors to each unique split
      # Only works when an RColorBrewer palette is passed
      if (split_colors) {
        if (length(x = unique(x = splits)) > length(x = cols)) {
          stop("Not enough colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
      }
      # Combine
      expr_data$id <- paste(expr_data$id, splits, sep = '_')
      print("unique(expr_data$id)")
      print(unique(expr_data$id))
      unique_splits <- unique(x = splits)
      # Create labels for each group_by and split_by level with an underscore
      # (each combo reads group_by[i]_split_by[j] for all combinations of i
      # group_by and j split_by levels), *for which there are cells*
      id_levels <-
        paste0(
          rep(x = id_levels, each = length(x = unique_splits)),
          "_",
          rep(x = unique(x = splits), times = length(x = id_levels))
        )
      print("id_levels")
      print(id_levels)
    }

    # Compute statistics for each group #
    # Store results in a list
    # Groups may be levels of the group_by variable, or group_by-split_by
    # combinations if split_by is enabled
    # If there are multiple features, each element will be a vector with the
    # expression of each feature for the given group
    plot_data <- lapply(
      X = unique(x = expr_data$id),
      FUN = function(ident) {
        # Pull data for computing statistics
        # Exclude last column of expr_data ("id" column)
        data_use <-
          expr_data[
            expr_data$id == ident,
            1:(ncol(x = expr_data) - 1),
            drop = FALSE
          ]

        # Average expression
        avg_exp <-
          apply(
            X = data_use,
            # Apply function by column
            MARGIN = 2,
            FUN = function(x) {
              return(mean(x = expm1(x = x)))
            }
          )

        # Percent of cells in group expressing feature
        pct_exp <-
          apply(
            X = data_use,
            MARGIN = 2,
            # Uses Seurat::PercentAbove()
            FUN = PercentAbove,
            # "threshold" is an argument of PercentAbove
            threshold = 0
          )

        return(
          list(
            avg_exp = avg_exp,
            pct_exp = pct_exp
          )
        )
      }
    )
    names(x = plot_data) <- unique(x = expr_data$id)

    # Cluster groups if `cluster_idents` is TRUE
    # Clustering is performed based on expression statistics for the currently
    # plotted set of features
    if (cluster_idents) {
      mat <- do.call(
        what = rbind,
        args = lapply(X = plot_data, FUN = unlist)
      )
      mat <- scale(x = mat)
      id_levels <- id_levels[hclust(d = dist(x = mat))$order]
    }

    # Construct dataframe from list of expression statistics
    # Return is a list of dataframes, one dataframe per group
    plot_data <- lapply(
      X = names(x = plot_data),
      FUN = function(x) {
        data_use <- as.data.frame(x = plot_data[[x]])
        # Add feature name as a column
        data_use$features_plot <- rownames(x = data_use)
        # Add column for the group ID
        data_use$id <- x
        return(data_use)
      }
    )
    # Combine list of dataframes produced above
    plot_data <- do.call(what = 'rbind', args = plot_data)

    # Convert the new "id" column in the dataframe to a factor
    if (!is.null(x = id_levels)) {
      plot_data$id <- factor(x = plot_data$id, levels = id_levels)
    }

    # Compute number of groups, and warn if the number of groups is 1, or less than 5
    ngroup <- length(x = levels(x = plot_data$id))

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
    avg_exp_scaled <-
      sapply(
        X = unique(x = plot_data$features_plot),
        FUN = function(x) {
          # Pull average expression column and subset to feature x
          data_use <- plot_data[plot_data$features_plot == x, 'avg_exp']
          if (scale) {
            print("Data before scaling")
            print(data_use)
            # Scale data (z-score)
            data_use <- scale(x = data_use)
            print("Data after scaling")
            print(data_use)
            # Any values with a z-score outside of col_min and col_max will be
            # normalized to those values
            data_use <- MinMax(data = data_use, min = col_min, max = col_max)
            print("Data after min/max transformation")
            print(data_use)
          } else {
            # If `scale` is FALSE, perform a log(x+1) transformation
            data_use <- log1p(x = data_use)
          }
          return(data_use)
        }
      )

    print("avg_exp.sacaled before conversion to vector")
    print(avg_exp_scaled)

    # Transfrom average expression to a vector
    avg_exp_scaled <- as.vector(x = t(x = avg_exp_scaled))

    print("avg_exp.sacaled after conversion to vector")
    print(avg_exp_scaled)

    # "Bin" the average expression into 20 bins for split plots
    if (split_colors) {
      # Max value will be 20, min value will be 1, and all other values will be
      # binned based on their value
      avg_exp_scaled <- as.numeric(x = cut(x = avg_exp_scaled, breaks = 20))
    }

    # Add scaled data to the plot_data dataframe and convert to factor
    plot_data$avg_exp_scaled <- avg_exp_scaled
    print("Plot data with scaled expression column")
    print(plot_data)
    plot_data$features_plot <- factor(
      x = plot_data$features_plot,
      levels = features
    )
    # Add NAs for percent expression values beneath dot_min
    plot_data$pct_exp[plot_data$pct_exp < dot_min] <- NA
    # Convert percent expression from fraction to percentage
    plot_data$pct_exp <- plot_data$pct_exp * 100

    print("Value of plot_data$id")
    print(plot_data$id)

    if (split_colors) {
      splits_use <- vapply(
        X = as.character(x = plot_data$id),
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
      print("Value of splits_use")
      print(splits_use)

      # Color for each dot: use colorRampPalette to determine the color (using
      # a continuous scale between grey and the color entered, with 20 steps)
      plot_data$colors <- mapply(
        FUN = function(color, value) {
          print("Value")
          print(value)
          print("Color")
          print(colorRampPalette(colors = c('grey', color))(20)[value])
          return(
            colorRampPalette(colors = c('grey', color))(20)[value]
          )
        },
        color = cols[splits_use],
        value = avg_exp_scaled
      )

      print("Split data after operations")
      print(head(plot_data))
    }

    # Color scale: if there are splits, use 20-level binned scale based on average
    # expression (feature-wide).
    # If not (default behavior), apply colors using ggplot2 scales
    color_by <- ifelse(test = split_colors, yes = 'colors', no = 'avg_exp_scaled')

    # Apply min/max cutoffs for *percent expression*
    if (!is.na(x = scale_min)) {
      plot_data[plot_data$pct_exp < scale_min, 'pct_exp'] <- scale_min
    }
    if (!is.na(x = scale_max)) {
      plot_data[plot_data$pct_exp > scale_max, 'pct_exp'] <- scale_max
    }

    # Add data for feature groups to the table (if defined by the user)
    if (!is.null(x = feature_groups)) {
      plot_data$feature_groups <- factor(
        x = feature_groups[plot_data$features_plot],
        levels = unique(x = feature_groups)
      )
    }

    plot <-
      ggplot(
        data = plot_data,
        # aes_string is deprecated in the current ggplot2 version
        mapping = aes_string(x = 'features_plot', y = 'id')
      ) +
      geom_point(mapping = aes_string(size = 'pct_exp', color = color_by)) +
      scale_func(range = c(0, dot_scale), limits = c(scale_min, scale_max)) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      guides(size = guide_legend(title = 'Percent Expressed')) +
      labs(
        x = 'Features',
        y = ifelse(test = is.null(x = split_by), yes = 'Identity', no = 'Split Identity')
      ) +
      theme_cowplot()

    # Create facets for feature groups, if present
    if (!is.null(x = feature_groups)) {
      plot <- plot + facet_grid(
        facets = ~feature_groups,
        scales = "free_x",
        space = "free_x",
        switch = "y"
      ) + theme(
        panel.spacing = unit(x = 1, units = "lines"),
        strip.background = element_blank()
      )
    }
    if (split_colors) {
      # scale_color_identity applies colors that were defined based on the
      # 20-bin scale
      plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
      plot <- plot + scale_color_distiller(palette = cols)
    } else {
      plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (!split_colors) {
      plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
    }
    return(plot)
  }
