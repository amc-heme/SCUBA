#' Single cell ridge plot
#'
#' Draws a ridge plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData)
#' @param cols Colors to use for plotting
#' @param idents Which classes to include in the plot (default is all)
#' @param sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted, can also pass 'increasing' or 'decreasing' to change sort direction
#' @param assay Name of assay to use, defaults to the active assay
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param log plot the feature axis on log scale
#' @param ncol Number of columns if multiple plots are displayed
#' @param slot Slot to pull expression data from (e.g. "counts" or "data")
#' @param stack Horizontally stack plots for each feature
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot
#' @param fill.by Color violins/ridges based on either 'feature' or 'ident'
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @export
#' @concept visualization
#'
#' @examples
#' data("pbmc_small")
#' RidgePlot(object = pbmc_small, features = 'PC_1')
#'
RidgePlot <- function(
    object,
    features,
    cols = NULL,
    idents = NULL,
    sort = FALSE,
    assay = NULL,
    group.by = NULL,
    y.max = NULL,
    same.y.lims = FALSE,
    log = FALSE,
    ncol = NULL,
    slot = 'data',
    stack = FALSE,
    combine = TRUE,
    fill.by = 'feature'
) {
  return(ExIPlot(
    object = object,
    type = 'ridge',
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    cols = cols,
    group.by = group.by,
    log = log,
    slot = slot,
    stack = stack,
    combine = combine,
    fill.by = fill.by
  ))
}

