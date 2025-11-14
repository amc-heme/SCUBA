#' Summarize object metadata
#'
#' Returns the names of all metadata variables in an object, or a single
#' modality of a mudata object.
#'
#' @inheritParams object_param
#' 
#' @param ... Additional arguments passed to S3 methods
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
    object,
    ...
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
    object,
    ...
  ){
    warning(
      paste0(
        "meta_varnames does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat, SingleCellExperiment and Anndata (AnnDataR6)."
      )
    )
  }

#' @describeIn meta_varnames Seurat objects
#' @export
meta_varnames.Seurat <-
  function(
    object,
    ...
  ){
    colnames(object@meta.data)
  }

#' @describeIn meta_varnames SingleCellExperiment objects
#' @export
meta_varnames.SingleCellExperiment <-
  function(
    object,
    ...
  ){
    colnames(colData(object))
  }

#' @describeIn meta_varnames Anndata objects
#' @export
meta_varnames.AnnDataR6 <-
  function(
    object,
    ...
  ){
    object$obs_keys()
  }

#' @describeIn meta_varnames Mudata objects
#' @export
meta_varnames.md._core.mudata.MuData <-
  function(
    object,
    # mod = NULL,
    ...
  ){
    # MuData: simply report column names of full table from fetch_metadata
    # Columns returned have an underscore between the modality and the 
    # metadata variable name instead of a ":", following MuData conventions.
    fetch_metadata(
      object,
      full_table = TRUE
      ) |> 
      colnames()
  }
  
#' @export
meta_varnames.mudata._core.mudata.MuData <-
  function(
    object,
    ...
  ){
    # mudata._core.mudata.MuData: possible class when loading 
    # Redirect to md._core.mudata.MuData method
    meta_varnames.md._core.mudata.MuData(
      object = object
    )
  }