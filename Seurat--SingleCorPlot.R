#' A single correlation plot
#'
#' @param data A data frame with two columns to be plotted
#' @param col.by A vector or factor of values to color the plot by
#' @param cols An optional vector of colors to use
#' @param pt.size Point size for the plot
#' @param smooth Make a smoothed scatter plot
#' @param rows.highight A vector of rows to highlight (like cells.highlight in
#' \code{\link{SingleDimPlot}})
#' @param legend.title Optional legend title
#' @param raster Convert points to raster format, default is \code{NULL}
#' which will automatically use raster if the number of points plotted is
#' greater than 100,000
#' @param raster.dpi the pixel resolution for rastered plots, passed to geom_scattermore().
#' Default is c(512, 512)
#' @param plot.cor ...
#' @param jitter Jitter for easier visualization of crowded points
#'
#' @return A ggplot2 object
#'
#' @importFrom stats cor
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 ggplot aes_string geom_point labs scale_color_brewer
#' scale_color_manual guides stat_density2d aes scale_fill_continuous
#' @importFrom scattermore geom_scattermore
#'
#' @keywords internal
#'
#' @export
#'
SingleCorPlot <- function(
    data,
    col.by = NULL,
    cols = NULL,
    pt.size = NULL,
    smooth = FALSE,
    rows.highlight = NULL,
    legend.title = NULL,
    na.value = 'grey50',
    span = NULL,
    raster = NULL,
    raster.dpi = NULL,
    plot.cor = TRUE,
    jitter = TRUE
) {
  pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
  if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e5)
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
      stop("'raster.dpi' must be a two-length numeric vector")
  }
  orig.names <- colnames(x = data)
  names.plot <- colnames(x = data) <- gsub(
    pattern = '-',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  names.plot <- colnames(x = data) <- gsub(
    pattern = ':',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  if (ncol(x = data) < 2) {
    msg <- "Too few variables passed"
    if (ncol(x = data) == 1) {
      msg <- paste0(msg, ', only have ', colnames(x = data)[1])
    }
    stop(msg, call. = FALSE)
  }
  plot.cor <- if (isTRUE(x = plot.cor)) {
    round(x = cor(x = data[, 1], y = data[, 2]), digits = 2)
  }
  else(
    ""
  )
  if (!is.null(x = rows.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = rows.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size,
      cols.highlight = 'red',
      col.base = 'black',
      pt.size = pt.size
    )
    cols <- highlight.info$color
    col.by <- factor(
      x = highlight.info$highlight,
      levels = rev(x = highlight.info$plot.order)
    )
    plot.order <- order(col.by)
    data <- data[plot.order, ]
    col.by <- col.by[plot.order]
  }
  if (!is.null(x = col.by)) {
    data$colors <- col.by
  }
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = names.plot[1], y = names.plot[2])
  ) +
    labs(
      x = orig.names[1],
      y = orig.names[2],
      title = plot.cor,
      color = legend.title
    )
  if (smooth) {
    # density <- kde2d(x = data[, names.plot[1]], y = data[, names.plot[2]], h = Bandwidth(data = data[, names.plot]), n = 200)
    # density <- data.frame(
    #   expand.grid(
    #     x = density$x,
    #     y = density$y
    #   ),
    #   density = as.vector(x = density$z)
    # )
    plot <- plot + stat_density2d(
      mapping = aes(fill = ..density.. ^ 0.25),
      geom = 'tile',
      contour = FALSE,
      n = 200,
      h = Bandwidth(data = data[, names.plot])
    ) +
      # geom_tile(
      #   mapping = aes_string(
      #     x = 'x',
      #     y = 'y',
      #     fill = 'density'
      #   ),
      #   data = density
      # ) +
      scale_fill_continuous(low = 'white', high = 'dodgerblue4') +
      guides(fill = FALSE)
  }
  position <- NULL
  if (jitter) {
    position <- 'jitter'
  } else {
    position <- 'identity'
  }
  if (!is.null(x = col.by)) {
    if (raster) {
      plot <- plot + geom_scattermore(
        mapping = aes_string(color = 'colors'),
        position = position,
        pointsize = pt.size,
        pixels = raster.dpi
      )
    } else {
      plot <- plot + geom_point(
        mapping = aes_string(color = 'colors'),
        position = position,
        size = pt.size
      )
    }
  } else {
    if (raster) {
      plot <- plot + geom_scattermore(position = position, pointsize = pt.size, pixels = raster.dpi)
    } else {
      plot <- plot + geom_point(position = position, size = pt.size)
    }
  }
  if (!is.null(x = cols)) {
    cols.scale <- if (length(x = cols) == 1 && cols %in% rownames(x = brewer.pal.info)) {
      scale_color_brewer(palette = cols)
    } else {
      scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + cols.scale
    if (!is.null(x = rows.highlight)) {
      plot <- plot + guides(color = FALSE)
    }
  }
  plot <- plot + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(x = span)) {
    plot <- plot + geom_smooth(
      mapping = aes_string(x = names.plot[1], y = names.plot[2]),
      method = 'loess',
      span = span
    )
  }
  return(plot)
}
