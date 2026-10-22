#' Fetch an expression matrix
#'
#' Returns the expression values of a set of features (all features of the assay by default) for a
#' set of cells as a features x cells matrix, sparse (`dgCMatrix`) by default. On the Julia engine
#' the block is read natively (BPCells directories, h5ad files); otherwise the layer is sliced in R.
#'
#' For AnnData (reticulate) objects the legacy signature `fetch_matrix(object, matrix_location)`
#' still returns the raw slot (`"X"`, `c("obsm", "X_umap")`, ...).
#'
#' @param object a single-cell object (Seurat, SingleCellExperiment, `scuba_source`, AnnData).
#' @param features feature names (bare, or keyed with the assay key); `NULL` for all features.
#' @param ... Passed to methods (`matrix_location` for AnnData objects).
#' @param assay assay / modality (default: the default assay).
#' @param layer layer to read (default: [default_layer()]).
#' @param cells cells to include, in the returned order (default: all cells).
#' @param sparse return a `dgCMatrix` (default) or a dense matrix.
#'
#' @returns A features x cells matrix with dimnames.
#' @rdname fetch_matrix
#' @export
fetch_matrix <-
  function(
    object,
    features = NULL,
    ...,
    assay = NULL,
    layer = NULL,
    cells = NULL,
    sparse = TRUE
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
    features = NULL,
    ...
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
    features = NULL,
    ...,
    matrix_location = features
  ){
    if (length(matrix_location) == 1){
      return(object[[matrix_location]])
      } else if (length(matrix_location) == 2){
      return(object[[matrix_location[1]]][[matrix_location[2]]])
      } else {
      stop('Unexpected format for matrix_location. This should be a ',
           'character vector of length one or two, giving the location of ',
           'the desired matrix in the object. For example, to pull from ',
           'the "obs" slot, use "obs". To pull from the "X_umap" matrix ',
           'in "obsm", use c("obsm", "X_umap").')
      }
    }

#' @rdname fetch_matrix
#' @export
fetch_matrix.Seurat <-
  function(
    object,
    features = NULL,
    ...,
    assay = NULL,
    layer = NULL,
    cells = NULL,
    sparse = TRUE
  ){
    if (.scuba_julia_object(object)) {
      return(.jl_fetch_matrix(object, features, assay = assay, layer = layer, cells = cells, sparse = sparse))
    }
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    layer <- layer %||% default_layer(object)
    m <- SeuratObject::LayerData(object[[assay]], layer = layer)
    .r_matrix_slice(m, features, cells, sparse)
  }

#' @rdname fetch_matrix
#' @export
fetch_matrix.SingleCellExperiment <-
  function(
    object,
    features = NULL,
    ...,
    assay = NULL,
    layer = NULL,
    cells = NULL,
    sparse = TRUE
  ){
    if (.scuba_julia_object(object)) {
      return(.jl_fetch_matrix(object, features, assay = assay, layer = layer, cells = cells, sparse = sparse))
    }
    layer <- layer %||% default_layer(object)
    exp <- if (is.null(assay) || identical(assay, SingleCellExperiment::mainExpName(object))) object else SingleCellExperiment::altExp(object, assay)
    m <- SummarizedExperiment::assay(exp, layer)
    .r_matrix_slice(m, features, cells, sparse)
  }

#' @rdname fetch_matrix
#' @export
fetch_matrix.scuba_source <-
  function(
    object,
    features = NULL,
    ...,
    assay = NULL,
    layer = NULL,
    cells = NULL,
    sparse = TRUE
  ){
    .scuba_source_engine()
    .jl_fetch_matrix(object$path, features, assay = assay, layer = layer, cells = .src_cells(object, cells), sparse = sparse)
  }

# Slice a features x cells layer in R (any matrix class that supports `[`), as dgCMatrix or dense.
.r_matrix_slice <- function(m, features, cells, sparse) {
  f <- features %||% rownames(m)
  missing_f <- setdiff(f, rownames(m))
  if (length(missing_f)) stop("Features not found: ", paste(utils::head(missing_f, 10), collapse = ", "), call. = FALSE)
  c <- cells %||% colnames(m)
  missing_c <- setdiff(c, colnames(m))
  if (length(missing_c)) stop("Cells not found: ", paste(utils::head(missing_c, 10), collapse = ", "), call. = FALSE)
  out <- m[f, c, drop = FALSE]
  if (sparse) {
    if (!inherits(out, "dgCMatrix")) out <- methods::as(methods::as(as.matrix(out), "CsparseMatrix"), "generalMatrix")
  } else {
    out <- as.matrix(out)
  }
  dimnames(out) <- list(f, c)
  out
}
