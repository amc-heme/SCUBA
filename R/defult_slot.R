#' Get default reduction from object
#'
#' Deprecated: this function was deprecated 
#' in version 0.10.0. Please use [`default_layer()`] instead. 
#' 
#' `r lifecycle::badge("deprecated")`
#'
#' @param object a single cell object supported by SCUBA. Currently, Seurat, SingleCellExperiment, and anndata objects are supported.
#'
#' @keywords internal
#'
#' @export
default_slot <-
  function(
    object
  ){
    lifecycle::deprecate_stop(
      when = "0.10.0",
      what = "default_slot()",
      details = "Please use `default_layer()` instead."
      )
  }