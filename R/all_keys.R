#' Get keys for all assays/reductions in an object
#'
#' Returns the "keys" of all reductions and modalities/assays/experiments in an object, which are used to fetch data via the `vars` parameter of `fetch_data`. To fetch features from an object, use the key representing the modality the feature was recorded in, plus an underscore and the feature name. To fetch reduction coordinates, use the key of the reduction, plus an underscore, and a number representing the dimension for which to retrieve coordinates. 
#'
#' @param object a single-cell object. Currently, Seurat, SingleCellExperiment, and anndata objects are supported.
#' 
#' @param obsm_key_format For MuData objects, this determines the format in which keys are returned. If "obsm" (the default), this will return a set of the obsm keys from each modality. If "mod_obsm", the obsm keys will be formatted as modality-obsm pairs (i.e. a obsm matrix in "RNA" named "UMAP" will be returned as "RNA_UMAP"). Using keys in the "mod_obsm" format avoids issues with ambiguous obsm matrices, where an obsm matrix with the same name exists in multiple modalities. Ambiguous obsm entries will case an error in `fetch_data` and `fetch_reduction`
#'
#' @rdname all_keys
#'
#' @returns a named character vector. The names of the vector are names of the modalities and reductions in the object, and the values are the corresponding keys to be passed to fetch_data. For Seurat objects, a key for metadata will also be displayed.
#'
#' @export
#' 
#' @examples
#' ## View keys ##
#' # Seurat objects
#' all_keys(AML_Seurat)
#' 
#' # SingleCellExperiment objects 
#' all_keys(AML_SCE())
#' 
#' # anndata objects 
#' all_keys(AML_h5ad())
#' 
#' ## Use of keys to construct fetch_data query
#' # Fetch a feature from the "protein" 
#' # modality using its key from above
#' fetch_data(
#'   AML_h5ad(), 
#'   vars = "protein_CD9-AB"
#'   ) |> str()
#' 
#' # Fetch reduction coordinates using 
#' # the key for the UMAP reduction
#' fetch_data(
#'   AML_h5ad(), 
#'   vars = c("X_umap_1", "X_umap_2")
#'   ) |> str()
all_keys <-
  function(
    object,
    ...
  ){
    UseMethod("all_keys")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
all_keys.default <-
  function(
    object,
    ...
  ){
    warning(
      paste0(
        "all_keys does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat, SingleCellExperiment, and anndata (AnnDataR6)."
        )
      )
    }

#' @describeIn all_keys Seurat objects
#' 
#' @export
all_keys.Seurat <-
  function(
    object
  ){
    # Seurat objects: run SeuratObject.Key method
    Key(object)
  }

#' @describeIn all_keys SingleCellExperiment objects
#'  
#' @export
all_keys.SingleCellExperiment <-
  function(
    object
  ){
    # SingleCellExperiment objects: return names of experiments and reductions
    keys <-
      c(mainExpName(object),
      altExpNames(object),
      reducedDimNames(object)
      )

    # Add names to vector returned for consistency with SeuratObject.Key() method
    names(keys) <- keys

    keys
  }

#' @describeIn all_keys SingleCellExperiment objects
#' @export
all_keys.AnnDataR6 <-
  function(
    object
  ){
    # Since anndata objects don't have defined locations for modalities and
    # reductions, downstream functions calling this method may not work
    # in the same way they do for Seurat and SingleCellExperiment objects.

    # Anndata objects: return names of experiments and reductions
    keys <-
      c("X",
        object$obsm_keys()
        )

    # Add names to vector returned for consistency with
    # SeuratObject.Key() method
    names(keys) <- keys

    keys
  }

#' @describeIn all_keys MuData objects
#' @export
all_keys.md._core.mudata.MuData <-
  function(
    object,
    obsm_key_format = "obsm"
  ){
    # Form keys of different types
    # Modality keys
    mod_keys <- object$mod_names 
    
    names(mod_keys) <- mod_keys 
    
    # Obsm keys
    obsm_keys <- c()
    
    for (mod in object$mod_names){
      mod_reductions <- object[[mod]]$obsm_keys()
      
      if (length(mod_reductions) > 0){
        # Construct keys from each obsm key in the modality, 
        # if there are keys for that modality
        if (obsm_key_format == "obsm"){
          # "obsm": return the obsm keys for each modality, as they
          # appear in each modality
          mod_obsm_keys <- mod_reductions
          
          obsm_keys <- c(obsm_keys, mod_obsm_keys)
        } else if (obsm_key_format == "mod_obsm"){
          # "mod_obsm": prepend the modality key to obsm keys. This 
          # shows how to format the key to pull obsm data for the 
          # specific modality when a obsm matrix with the same name
          # appears in multiple modalities
          mod_obsm_keys <- paste(mod, mod_reductions, sep = "_")
          # For mod_obsm keys, show more descriptive names for obsm keys,
          # indicating the modality and reduction specified by the key
          names(mod_obsm_keys) <- 
            paste0(mod_reductions, " matrix, from ", mod, " modality")
          
          obsm_keys <- c(obsm_keys, mod_obsm_keys)
        }
      }
    }
    
    # When returning obsm key names without the modality prepended, 
    # show the obsm keys as a set of unique values 
    if (obsm_key_format == "obsm"){
      obsm_keys <- unique(obsm_keys)
      # Add names to unique keys (same as values)
      # Names are added for consistency with other methods
      names(obsm_keys) <- obsm_keys
    }
    
    # Add "obs" key for metadata
    metadata_key <- c("Metadata" = "obs")
    
    keys <- 
      c(metadata_key,
        mod_keys,
        obsm_keys
        )
    
    keys
  }

#' @export
all_keys.mudata._core.mudata.MuData <-
  function(
    object,
    obsm_key_format = "obsm"
  ){
    # mudata._core.mudata.MuData: possible class when loading 
    # Redirect to md._core.mudata.MuData method
    all_keys.md._core.mudata.MuData(
      object = object,
      obsm_key_format = obsm_key_format
    )
  }