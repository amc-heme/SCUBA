#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @inheritParams plot_ridge
#' @param pt_size Point size for geom_violin
#' @param split_by A variable to split the violin plots by,
#' @param split_plot  plot each group of the split violin plots by multiple or
#' single violin shapes.
#' @param adjust Adjust parameter for geom_violin
#' @param flip flip plot orientation (identities on x-axis)
#' @param add_noise determine if adding a small noise for plotting
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
#' VlnPlot(object = pbmc_small, features = 'LYZ', split_by = 'groups')
#'
plot_violin <- function(
    object,
    features,
    cols = NULL,
    pt_size = NULL,
    idents = NULL,
    sort = FALSE,
    assay = NULL,
    group_by = NULL,
    split_by = NULL,
    adjust = 1,
    y_max = NULL,
    same_y_lims = FALSE,
    log = FALSE,
    ncol = NULL,
    slot = 'data',
    split_plot = FALSE,
    stack = FALSE,
    combine = TRUE,
    fill_by = 'feature',
    flip = FALSE,
    add_noise = TRUE,
    raster = NULL
) {
  # if (
  #   !is.null(x = split_by) &
  #   getOption(x = 'Seurat.warn.vlnplot.split', default = TRUE)
  # ) {
  #   message(
  #     "The default behaviour of split_by has changed.\n",
  #     "Separate violin plots are now plotted side-by-side.\n",
  #     "To restore the old behaviour of a single split violin,\n",
  #     "set split_plot = TRUE.
  #     \nThis message will be shown once per session."
  #   )
  #   options(Seurat.warn.vlnplot.split = FALSE)
  # }
  return(
    expr_plot(
      object = object,
      type = ifelse(test = split_plot, yes = 'splitViolin', no = 'violin'),
      features = features,
      idents = idents,
      ncol = ncol,
      sort = sort,
      assay = assay,
      y_max = y_max,
      same_y_lims = same_y_lims,
      adjust = adjust,
      pt_size = pt_size,
      cols = cols,
      group_by = group_by,
      split_by = split_by,
      log = log,
      slot = slot,
      stack = stack,
      combine = combine,
      fill_by = fill_by,
      flip = flip,
      add_noise = add_noise,
      raster = raster
      )
    )
}
