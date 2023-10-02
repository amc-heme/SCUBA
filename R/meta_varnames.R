#' Object metadata summary
#'
#' Returns the names of all metadata variables in an object.
#'
#' @param object a single-cell object. Currently, Seurat, Anndata, and
#' SingleCellExperiment objects are supported.
#'
#' @rdname meta_varnames
#'
#' @export
meta_varnames <-
  function(
    object
  ){
    UseMethod("meta_varnames")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
meta_varnames.default <-
  function(
    object
  ){
    warning(
      paste0(
        "meta_varnames does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn meta_varnames Seurat objects
#' @export
meta_varnames.Seurat <-
  function(
    object
  ){
    colnames(object@meta.data)
  }

#' @describeIn meta_varnames SingleCellExperiment objects
#' @export
meta_varnames.SingleCellExperiment <-
  function(
    object
  ){
    colnames(colData(object))
  }

#' @describeIn meta_varnames Anndata objects
#' @export
meta_varnames.AnnDataR6 <-
  function(
    object
  ){
    object$obs_keys()
  }
