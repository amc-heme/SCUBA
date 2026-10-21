#' Return the default layer
#'
#' Returns the default layer for the object passed. The default layer is chosen
#' based on the conventions used in each object to name the normalized counts
#' layer.
#' 
#' This is a utility function that is most useful for defining defaults in plotting functions, to reduce the number of required parameters end users need to supply. [see our website]() for an example of usage in a function.
#'
#' @inheritParams object_param
#'
#' @rdname default_layer
#'
#' @export
#' 
#' @examples
#' # Seurat objects
#' default_layer(AML_Seurat)
#' 
#' # SingleCellExperiment objects
#' default_layer(AML_SCE())
#' 
#' # anndata objects
#' default_layer(AML_h5ad())
default_layer <-
  function(
    object
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
    object,
    ...
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
