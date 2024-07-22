#' Get default reduction from object
#'
#' Returns the default slot for the object passed.
#' 
#' `r lifecycle::badge("deprecated")`
#' 
#' This function was deprecated in version 0.10.0 and renamed to 
#' `default_layer()`. The function will be removed in 1.0.0.
#'
#' @param object a single cell object supported by SCUBA. Currently, Seurat, SingleCellExperiment, and anndata objects are supported.
#'
#' @export
default_slot <-
  function(
    object
  ){
    lifecycle::deprecate_warn(
      when = "0.10.0",
      what = "default_slot()",
      details = 
        paste0(
          "Please use `default_layer()` instead. `default_slot()` will be ",
          "removed in 1.0.0."
          )
      )
    
    default_layer(object)
  }