#' Summarize object metadata
#'
#' Returns the names of all metadata variables in an object.
#'
#' @inheritParams object_param
#' 
#' @rdname meta_varnames
#'
#' @export
#' 
#' @examples
#' # Seurat objects
#' meta_varnames(AML_Seurat)
#' 
#' # SingleCellExperiment objects
#' meta_varnames(AML_SCE())
#' 
#' # anndata objects
#' meta_varnames(AML_h5ad())
#' 
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
    .scuba_unsupported_class(object, "meta_varnames")
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
    .anndata_obs_columns(object)
  }
