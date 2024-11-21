#' Get default reduction from object
#'
#' Returns the default layer for the object passed. The default layer is chosen
#' based on the conventions used in each object to name the normalized counts
#' layer.
#'
#' @param object A single cell object. Currently, Seurat, SingleCellExpleriment, 
#' and anndata objects are supported.
#' @param ... Currently unused.
#'
#' @rdname default_layer
#'
#' @export
default_layer <-
  function(
    object,
    ...
  ){
    UseMethod("default_layer")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
default_layer.default <-
  function(
    object
  ){
    warning(
      paste0(
        "default_layer does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn default_layer Seurat objects
#' 
#' For Seurat objects, the default layer is "data".
#' 
#' @export
default_layer.Seurat <-
  function(
    object
  ){
    # Seurat objects: "data"
    "data"
  }

#' @describeIn default_layer SingleCellExperiment objects 
#' 
#' For SingleCellExperiment objects, the default layer is "logcounts".
#' 
#' @export
default_layer.SingleCellExperiment <-
  function(
    object
  ){
    # SingleCellExperiment: default layer is "logcounts"
    "logcounts"
  }

#' @describeIn default_layer Anndata objects 
#' 
#' For Anndata objects, the default layer is \code{NULL}, which will direct 
#' FetchData to pull feature epxression data from the X matrix.
#' 
#' @export
default_layer.AnnDataR6 <-
  function(
    object
  ){
    NULL
  }

#' @describeIn default_layer MuData objects
#' 
#'  For MuData objects, the default layer is \code{NULL}, which will direct 
#' FetchData to pull feature epxression data from the X matrix.
#' 
#' @export
default_layer.MuData <-
  function(
    object
  ){
    NULL
  }
