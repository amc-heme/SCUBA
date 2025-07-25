#' Get names of all features in an assay/experiment/modality
#'
#' Returns the names of all features in a modality. This utility function can be used in several applications:
#' - In Shiny apps, return available features for passage to a selection menu
#' - Before using [`fetch_data()`], generate a list of features in the object for searching, or test if a feature is present before requesting data
#' 
#' @inheritParams object_param
#' @param assay the name of an assay/modality for which to view features.
#'
#' @rdname features_in_assay
#' 
#' @export
#' 
#' @examples
#' features_in_assay(AML_Seurat, assay = "RNA") |> str()
#' 
#' # Check if a feature is present in an assay
#' "MEIS1" %in% features_in_assay(AML_Seurat, assay = "RNA")
features_in_assay <-
  function(
    object,
    assay
  ){
    if (length(assay) > 1){
      stop("Only one assay/experiment name may be passed to `assay`.")
    }
    
    UseMethod("features_in_assay")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
features_in_assay.default <-
  function(
    object,
    assay
  ){
    warning(
      paste0(
        "features_in_assay does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn features_in_assay Seurat objects
#' @export
features_in_assay.Seurat <-
  function(
    object,
    assay
  ){
    if (!assay %in% names(object@assays)){
      stop("Assay ", assay, " not found in the current object.")
    }
    
    # Return rownames of assay entered
    rownames(object@assays[[assay]])
  }

#' @describeIn features_in_assay SingleCellExperiment objects
#' @export
features_in_assay.SingleCellExperiment <-
  function(
    object,
    assay
  ){
    # Check if assay ("experiment" in SingleCellExperiment) is a valid  
    if (!assay %in% c(mainExpName(object), altExpNames(object))){
      stop("Assay ", assay, " not found in the current object.")
    }
    
    # If the assay is the main experiment, use rownames(object)
    if (assay == mainExpName(object)){
      rownames(object)
    } else {
      # For other experiments, use rownames on the alternate experiment 
      # corresponding to the assay name
      rownames(altExps(object)[[assay]])
    }
  }

#' @describeIn features_in_assay Anndata objects
#' @export
features_in_assay.AnnDataR6 <-
  function(
    object,
    assay
  ){
    # Check that assay is valid (either X, or a feature from obsm)
    if (!assay %in% c("X", object$obsm_keys())){
      stop("Assay ", assay, " not found in the obsm slot of the current object.")
    }
    
    if (assay == "X"){
      # For the X matrix, feature names are in var_names
      object$var_names
    } else {
      # Otherwise, use the column names of the matrix stored in obsm
      colnames(object$obsm[[assay]])
    }
  }