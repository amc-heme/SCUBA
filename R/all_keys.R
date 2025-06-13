#' Get keys for all assays/reductions in an object
#'
#' Returns the "keys" of all reductions and modalities/assays/experiments in an object, which are used to fetch data via the `vars` parameter of `fetch_data`. To fetch features from an object, use the key representing the modality the feature was recorded in, plus an underscore and the feature name. To fetch reduction coordinates, use the key of the reduction, plus an underscore, and a number representing the dimension for which to retrieve coordinates. 
#'
#' @param object a single-cell object. Currently, Seurat, SingleCellExperiment, and anndata objects are supported.
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
    object
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
    object
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
    object
  ) {
    
     keys <- 
       c("X",
         object$obsm_keys()
         )
     
     names(keys) <- keys
     
     keys
  }
#all_keys(SCUBA:::AML_h5mu())
