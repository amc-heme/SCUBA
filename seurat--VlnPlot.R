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
