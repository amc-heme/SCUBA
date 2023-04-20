#' Get default reduction from object
#'
#' Returns a reduction from the object to use as the default.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param ... Currently unused.
#'
#' @rdname default_reduction
#'
#' @export
default_reduction <-
  function(
    object,
    ...
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
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn default_reduction Seurat objects
#' @export
default_reduction.Seurat <-
  function(
    object
  ){
    # Seurat objects: use existing SeuratObject method
    DefaultDimReduc(object)
  }

#' @describeIn default_reduction SingleCellExperiment objects
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
      stop('Unable to find a reduction matching "UMAP", "TSNE", or "PCA". Please specify the reudction to use via the `reduction` parameter.')
    }
  }
