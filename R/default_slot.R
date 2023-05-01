#' Get default reduction from object
#'
#' Returns the default slot for the object passed.
#'
#' @param object a single cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param ... Currently unused.
#'
#' @rdname default_slot
#'
#' @export
default_slot <-
  function(
    object,
    ...
  ){
    UseMethod("default_slot")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
default_slot.default <-
  function(
    object
  ){
    warning(
      paste0(
        "default_slot does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn default_slot Seurat objects (default slot is "data")
#' @export
default_slot.Seurat <-
  function(
    object
  ){
    # Seurat objects: "data"
    "data"
  }

#' @describeIn default_slot SingleCellExperiment objects (default slot is "logcounts")
#' @export
default_slot.SingleCellExperiment <-
  function(
    object
  ){
    # SingleCellExperiment: default slot is "logcounts"
    "logcounts"
  }
