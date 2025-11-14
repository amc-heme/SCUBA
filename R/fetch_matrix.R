#' Fetch matrix or metadata from single-cell containers
#'
#' Generic function to retrieve a matrix (e.g. expression or embedding)
#' or metadata (cells or features) from various single-cell containers:
#' Seurat, SingleCellExperiment, Anndata, and MuData objects.
#'
#' @param object A single-cell object. Supported classes include:
#'   \itemize{
#'     \item Seurat
#'     \item SingleCellExperiment
#'     \item Anndata
#'     \item MuData
#'   }
#' @param matrix_location A string specifying what to fetch. For MuData, valid keys are:
#'   \itemize{
#'     \item \code{"obs"}: cell metadata of the main MuData object
#'     \item \code{"var"}: feature metadata of the main MuData object
#'     \item \code{"<modality>"}: the primary expression matrix \code{X} of a modality
#'     \item \code{"<modality>_<layer>"}: a named layer in a modality
#'     \item \code{"<modality>_obs"}: cell metadata of a modality
#'     \item \code{"<modality>_var"}: feature metadata of a modality
#'     \item \code{"<modality>_<obsm_key>"}: embedding or other \code{obsm} matrix for a modality
#'   }
#' @param layer Optional. If \code{NULL}, the matrix corresponding to the default layer (via `default_layer()`) is used. This is the matrix corresponding to normalized counts, or X in anndata and MuData objects.
#'
#' @return Depending on the request:
#'   \itemize{
#'     \item For expression or embedding matrices: a dense or sparse matrix.
#'     \item For metadata queries: a \code{data.frame} or DataFrame of cell/feature metadata.
#'   }
#'
#' @examples
#' # Seurat objects:
#' # expr <- fetch_matrix(seurat_obj, "counts")
#'
#' # MuData objects:
#' # umap <- fetch_matrix(mudata_obj, "rna_umap")
#' # genes <- fetch_matrix(mudata_obj, "rna_var")
#'
#' @keywords internal
#'
#' @export
#' @rdname fetch_matrix
fetch_matrix <-
  function(
    object,
    matrix_location,
    layer
  ){
    UseMethod("fetch_matrix")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
fetch_matrix.default <-
  function(
    object,
    matrix_location,
    layer
  ){
    warning(
      paste0(
        "fetch_matrix does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat, SingleCellExperiment, ",
        "and anndata."
      )
    )
  }

# # describeIn fetch_matrix Seurat objects
# # export
# fetch_matrix.Seurat <-
#   function(
#     object,
#     matrix_location
#   ){
#
#   }
#
# # describeIn fetch_matrix SingleCellExperiment objects
# # export
# fetch_matrix.SingleCellExperiment <-
#   function(
#     object,
#     matrix_location
#     ){
#
#     }



#' @describeIn fetch_matrix AnnDataR6 objects
#' @export
fetch_matrix.AnnDataR6 <-
  function(
    object,
    matrix_location,
    layer
  ){
    stop("fetch_matrix not yet implemented for anndata objects.")
    }

#' @describeIn fetch_matrix MuData objects
#' @export
fetch_matrix.md._core.mudata.MuData <-
  function(
    object,
    matrix_location,
    layer
  ){
    # Determine the keys of valid matrices 
    # Valid matrix "keys"
    # 1. "obs"
    # 2. "var"
    # 3. Set of modality names
    ## 3.a. Modality name, plus an obsm name (mod_obsm keys)
    ## 3.b. Modality name, plus "obs"
    ## 3.c. Modality name, plus "var"
    
    # All_keys returns "obs", modality, and mod_obsm keys  
    matrix_keys <- 
      all_keys(object = object, obsm_key_format = "mod_obsm") |> 
      unname()
    
    # Add obs, var keys for each modality, and "var" from main MuData object
    mod_plus_obs <- 
      sapply(
        object$mod_names,
        function(mod){
          paste(mod, "obs", sep = "_")
        }) |> unname()
    
    mod_plus_var <-
      sapply(
        object$mod_names,
        function(mod){
          paste(mod, "var", sep = "_")
        }) |> unname()
    
    matrix_keys <-
      c(matrix_keys,
        mod_plus_obs,
        mod_plus_var,
        "var"
        )
    
    # Determine the modality of the matrix associated with the key, if any
    mod_key_match <- match_key(matrix_location, keys = object$mod_names)
    
    mod <- mod_key_match$key
    
    if (!is.null(mod)){
      # If the key matches a modality, test if there is a suffix specifying
      # a sub-location
      # (i.e. mod+"obs", mod+"var", or mod_obsm keys)
      sub_location <- mod_key_match$suffix
      
      if (!is.null(sub_location)){
        # Suffix specifying sub-location is present
        # Return the matrix at the specified sub-location
        mod_obsm_keys <- object[[mod]]$obsm_keys()
        
        if (sub_location == "obs"){
          # Return cell metadata from main MuData object
          return(object[[mod]]$obs)
        } else if (sub_location == "var"){
          # Return feature metadata from main MuData object
          return(object[[mod]]$var)
        } else if (sub_location %in% mod_obsm_keys){
          return(object[[mod]]$obsm[[sub_location]])
        } else {
          # Throw an error for un-supported keys
          stop(
            paste0(
              "The requested matrix ", 
              sub_location, 
              " is not present in the modality ", 
              mod, 
              "."
              )
            )
        }
      } else {
        # No suffix: matrix location is the expression matrix 
        if (is.null(layer)){
          # If layer is not defined, return the X matrix from the modality
          return(object[[mod]]$X)
        } else {
          # If defined, return the matrix for the specified layer
          return(object[[mod]]$layers[[layer]])
        }
      }
    } else {
      # Matrix does not match a modality
      if (matrix_location == "obs"){
        # Cell metadata from main MuData object
        return(object$obs)
      } else if (matrix_location == "var"){
        # Feature metadata from main MuData object
        return(object$var)
      } else {
        stop(
          paste0(
            "Unrecognized input for matrix_location: ", 
            matrix_location, 
            "."
            )
          )
      }
    }
  }

#' @export
fetch_matrix.mudata._core.mudata.MuData <-
  function(
    object,
    matrix_location,
    layer
  ){
    # mudata._core.mudata.MuData: possible class when loading 
    # Redirect to md._core.mudata.MuData method
    fetch_matrix.md._core.mudata.MuData(
      object = object,
      matrix_location = matrix_location,
      layer = layer
    )
  }