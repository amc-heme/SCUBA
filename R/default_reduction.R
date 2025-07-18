#' Get default reduction from object
#'
#' Returns the default reduction for the single-cell object passed. `default_reduction` will check for reductions in the order below. If a reduction exists, it will be returned. If not, the function will check for the next reduction on the list. If none of the reductions in the list exist, the function will return an error.
#' 1. UMAP
#' 2. t-SNE
#' 3. PCA
#' 
#' @inherit default_layer details 
#' @inheritParams object_param
#'
#' @rdname default_reduction
#'
#' @export
#' 
#' @examples
#' # Seurat objects
#' default_reduction(AML_Seurat)
#' 
#' # SingleCellExperiment objects
#' default_reduction(AML_SCE())
#' 
#' # anndata objects
#' default_reduction(AML_h5ad())
default_reduction <-
  function(
    object
  ){
    UseMethod("default_reduction")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
default_reduction.default <-
  function(
    object
  ){
    warning(
      paste0(
        "default_reduction does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: AnnDataR6, Seurat, SingleCellExperiment and Mudata."
      )
    )
  }

#' @describeIn default_reduction Seurat objects
#' 
#' In Seurat objects, `default_reduction` is a wrapper for [SeuratObject::DefaultDimReduc()].
#' 
#' @export
default_reduction.Seurat <-
  function(
    object
  ){
    # Seurat objects: use existing SeuratObject method
    DefaultDimReduc(object)
  }

#' @describeIn default_reduction SingleCellExperiment objects
#' 
#' The search order for SingleCellExperiment objects is below. `fetch_reduction` 
#' will search in the order described for a reduction named exactly as described.
#' 1. UMAP
#' 2. TSNE
#' 3. PCA
#' 
#' @export
default_reduction.SingleCellExperiment <-
  function(
    object
  ){
    reductions <- reducedDimNames(object)

    # SingleCellExperiment objects
    # Conditional structure reflects priorities of reductions
    # 1. UMAP
    # 2. TSNE
    # 3. PCA
    if ("UMAP" %in% reductions){
      "UMAP"
    } else if ("TSNE" %in% reductions){
      "TSNE"
    } else if ("PCA" %in% reductions){
      "PCA"
    } else {
      stop('Unable to find a reduction matching "UMAP", "TSNE", or "PCA". Please specify the reduction to use via the `reduction` parameter.')
    }
  }

#' @describeIn default_reduction AnnDataR6 objects
#' 
#' The search order for anndata objects is below. `fetch_reduction` will search 
#' in the order described for a reduction named exactly as described.
#' 1. X_umap
#' 2. X_tsne
#' 3. X_pca
#' 
#' @export
default_reduction.AnnDataR6 <-
  function(
    object
  ){
    # All reductions should be stored in obsm. Obsm matrix names are fetched
    # using obsm_keys().
    reductions <- object$obsm_keys()

    # AnnData objects
    # Conditional structure reflects priorities of reductions
    # 1. UMAP
    # 2. TSNE
    # 3. PCA
    if ("X_umap" %in% reductions){
      "X_umap"
    } else if ("X_tsne" %in% reductions){
      "X_tsne"
    } else if ("X_pca" %in% reductions){
      "X_pca"
    } else {
      stop('Unable to find a reduction matching "X_umap", "X_tsne", or "X_pca". Please specify the reduction to use via the `reduction` parameter.')
    }
  }


#' @describeIn default_reduction Mudata objects
#' @export
default_reduction.md._core.mudata.MuData <-
  function(
    object
  ){
    # Priorities of reductions
    # 1. UMAP
    # 2. TSNE
    # 3. PCA
    # Construct list storing modalities for which the reduction is found, if any
    mod_reduction_matches <- list()
    
    # For each of the "default" embeddings, record the modalities 
    # in which it was found
    for (embedding in c("X_umap", "X_tsne", "X_pca")){
      # Initialize an empty vector before looping through modalities
      mod_reduction_matches[[embedding]] <- c()
      for (modality in object$mod_names){
        if (embedding %in% object$mod[modality]$obsm_keys()){
          # If the embedding is found, record the name of the modality
          # where it was found
          mod_reduction_matches[[embedding]] <- 
            c(mod_reduction_matches[[embedding]], modality)
        }
      }
    }
    
    # Process the mapping of reductions to modalities
    # Look for each default reduction in order of priority
    for (embedding in c("X_umap", "X_tsne", "X_pca")){
      if (embedding %in% names(mod_reduction_matches)){
        # When the preferred reduction is found, check if it is present in
        # one modality or multiple modalities
        if (length(mod_reduction_matches[[embedding]]) == 1){
          # Stop iteration through matches when the first match is identified,
          # and return the matching reduction
          return(embedding)
        } else {
          # If the reduction is present in multiple modalities, throw an 
          # error (user will need to say which modality to retrieve reduction
          # from)
          stop(
            paste0(
              "The reduction ", embedding, " was identified as the ",
              "preferred reduction, but it is present in multiple ",
              "modalities. It is not possible to determine which ",
              "modality-reduction combo should be considered the default. ",
              "Please specify the modality-reduction combo to use via the ",
              "`reduction` parameter of the function from which ",
              "default_reduction was called."
              )
          )
        }
      }
    }
    
    # Return an error if iteration completes and no matches are found
    stop(
      'Unable to find a reduction matching "X_umap", "X_tsne", or "X_pca" ',
      'in any modality of the MuData object. Please specify the reduction ',
      'to use via the `reduction` parameter of the function from which ',
      'default_reduction was called.'
    )
  }

#' @export
default_reduction.mudata._core.mudata.MuData <-
  function(
    object
    ){
    default_reduction.md._core.mudata.MuData(
      object = object
    )
  }