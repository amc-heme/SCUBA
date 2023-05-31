#' Get names of assays/experiments in object
#'
#' Returns the names of all assays in a single-cell object. For Seurat objects,
#' returns assay names, and for SingleCellExperiment objects "experiments", the
#' equivalent of assays, are returned.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param ... Currently unused.
#'
#' @rdname assay_names
#'
#' @export
assay_names <-
  function(
    object,
    ...
  ){
    UseMethod("assay_names")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
assay_names.default <-
  function(
    object
  ){
    warning(
      paste0(
        "assay_names does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn assay_names Seurat objects
#' @export
assay_names.Seurat <-
  function(
    object
  ){
    # Seurat objects: access assays directly
    names(object@assays)
  }

#' @describeIn assay_names SingleCellExperiment objects
#' @export
assay_names.SingleCellExperiment <-
  function(
    object
  ){
    # SingleCellExperiment objects: return names of main, alternate experiments
    c(mainExpName(object), altExpNames(object))
  }
