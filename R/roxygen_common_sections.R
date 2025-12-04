#' Object Parameter
#'
#' The object parameter is used in nearly every instance of documentation in the SCUBA package. The parameter is described here to avoid excessive copy/pasting of documentation. When support is added for new object classes, the parameter description here should be updated.
#'
#' @keywords internal
#'
#' @param object A single cell object. Currently, Seurat, SingleCellExperiment,
#' and anndata objects are supported.
#'
#' @name object_param
NULL

#' Features Parameter
#'
#' Common parameter for functions that accept feature names.
#'
#' @keywords internal
#'
#' @param features A character vector of feature names to retrieve from the
#'   object.
#'
#' @name features_param
NULL

#' Layer Parameter
#'
#' Common parameter for functions that accept a layer name.
#'
#' @keywords internal
#'
#' @param layer The layer to pull data from. Layers are referred to as "slots"
#'   in Seurat objects v4 and earlier, and "assays" in SingleCellExperiment
#'   objects. If `NULL`, the default layer will be used.
#'
#' @name layer_param
NULL

#' Cells Parameter
#'
#' Common parameter for functions that accept cell names.
#'
#' @keywords internal
#'
#' @param cells A character vector of cell names to include, as they are named
#'   in the object (i.e. according to `colnames(object)`). If `NULL`, data will
#'   be returned for all cells in the object.
#'
#' @name cells_param
NULL
