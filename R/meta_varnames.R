#' Object metadata summary
#'
#' Returns the names of all metadata variables in an object, or a single
#' modality of a mudata object.
#'
#' @param object a single-cell object. Currently, Seurat, Anndata,
#' SingleCellExperiment, and Mudata objects are supported.
#' @param modality a mudata modality (e.g. "RNA", "ADT") to specify a specific
#' set of metadata variables
#'
#' @rdname meta_varnames
#'
#' @export
meta_varnames <-
  function(
    object,
    modality
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
        ". Currently supported classes: Seurat, SingleCellExperiment and Anndata (AnnDataR6)."
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

#' @describeIn meta_varnames Mudata objects
#' @export
meta_varnames.md._core.mudata.MuData <-
  function(
    object,
    modality = NULL
  ){
    library(reticulate)
    
    # Source fetch_reduction python script for mudata
    python_path =
      system.file(
        "extdata",
        "Python",
        "fetch_meta_varnames_mudata.py",
        package = "SCUBA"
      )
    
    reticulate::source_python(python_path)
  
    py$fetch_meta_varnames_mudata(
      obj = object, 
      modal = modality
    )
  }
  
#previous implementation for mudata (delete when additional logic is complete):
  # function(
  #   object
  # ){
  #   object$obs_keys()
  # }
