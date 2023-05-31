#' Get names of reductions in object
#'
#' Returns the names of all reductions in a single-cell object.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param ... Currently unused.
#'
#' @rdname reduction_names
#'
#' @export
reduction_names <-
  function(
    object,
    ...
  ){
    UseMethod("reduction_names")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
reduction_names.default <-
  function(
    object
  ){
    warning(
      paste0(
        "reduction_names does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn reduction_names Seurat objects
#' @export
reduction_names.Seurat <-
  function(
    object
  ){
    # Seurat objects: access reductions directly
    names(object@reductions)
  }

#' @describeIn reduction_names SingleCellExperiment objects
#' @export
reduction_names.SingleCellExperiment <-
  function(
    object
  ){
    # SingleCellExperiment objects: use reducedDimNames
    reducedDimNames(object)
  }
