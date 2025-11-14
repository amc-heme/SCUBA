#' Metadata for multiple single cell object classes
#'
#' Updates the metadata of either a Seurat, SingleCellExperiment or Anndata 
#' object. This method is intended to be used with fetch_metadata to pull a  
#' table from an object, modify it, and then store the updated table. 
#'
#' @param object a single-cell object. Currently, Seurat, SingleCellExperiment
#' and AnnData objects are supported.
#' @param table a modified metadata table.
#' @param mod For MuData objects, the obs table of a specific modality can be
#' updated by supplying the mod parameter. If NULL (the default), the main obs 
#' table will be updated.
#'
#' @return the object passed to \code{object} with the modified metadata table.
#' @noRd
update_object_metadata <-
  function(
    object,
    table,
    ...
  ){
    UseMethod("update_object_metadata")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @export
#' @noRd
update_object_metadata.default <-
  function(
    object,
    table,
    ...
  ){
    warning(
      paste0(
        "update_object_metadata does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat, SingleCellExperiment and AnnData."
      )
    )
  }

#' @describeIn update_object_metadata Seurat objects
#' @export
#' @noRd
update_object_metadata.Seurat <-
  function(
    object,
    table,
    ...
  ){
    object@meta.data <- table
    
    object
  }

#' @describeIn update_object_metadata SingleCellExperiment objects
#' @importFrom SummarizedExperiment colData<-
#' @export
#' @noRd
update_object_metadata.SingleCellExperiment <-
  function(
    object,
    table,
    ...
  ){
    colData(object) <- table
    
    object
  }

#' @describeIn update_object_metadata Anndata objects
#' @export
#' @noRd
update_object_metadata.AnnDataR6 <-
  function(
    object,
    table,
    ...
  ){
    object$obs <- table
    
    object
  }


#' @describeIn update_object_metadata Anndata objects
#' @export
#' @noRd
update_object_metadata.md._core.mudata.MuData <-
  function(
    object,
    table,
    mod = NULL,
    ...
  ){
    if (is.null(mod)){
      object$obs <- table
    } else {
      object[[mod]]$obs <- table
    }
    
    object
  }

#' @export
#' @noRd
update_object_metadata.mudata._core.mudata.MuData <-
  function(
    object,
    table,
    mod = NULL,
    ...
    ){
    update_object_metadata.md._core.mudata.MuData(
      object = object,
      table = table,
      mod = mod
      )
    }