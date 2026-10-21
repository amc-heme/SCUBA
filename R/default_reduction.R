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
        ". Currently supported classes: AnnDataR6, Seurat and SingleCellExperiment."
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
