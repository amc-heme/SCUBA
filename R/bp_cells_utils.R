## BPCells utility helpers (internal)
## Detection & fast feature access for BPCells-backed Seurat assays
##
## These functions enable SCUBA plotting functions to bypass SCUBA::fetch_data()
## for large feature queries when the underlying assay layers are BPCells-backed
## disk matrices. They are intentionally conservative: failures or missing
## dependencies fall back to standard FetchData behavior.

#' Detect BPCells-backed assay layer in a Seurat object
#'
#' Internal helper: returns TRUE if the indicated assay/layer appears to be
#' backed by a BPCells on-disk matrix (heuristic: matrix for the layer is
#' an IterableMatrix).
#'
#' @param object A Seurat object to check.
#' @param assay The assay to check. If `NULL`, the default assay will be used.
#' @param layer The layer to check within the assay. If `NULL`, the "data" layer
#'   will be checked.
#'
#' @return A logical scalar: `TRUE` if the specified assay/layer is backed by
#'   a BPCells IterableMatrix, `FALSE` otherwise.
#'
#' @keywords internal
#' @export
is_bp_cells_seurat <- function(
  object,
  assay = NULL,
  layer = NULL
) {
  if (!inherits(object, "Seurat")) {
    return(FALSE)
  }
  assay <- assay %||% DefaultAssay(object)
  layer <- layer %||% "data"

  # Attempt to access the matrix and check if it is an IterableMatrix
  # (BPCells Matrix class). Falls back to FALSE if layer access fails
  tryCatch(
    {
      matrix <- object@assays[[assay]]@layers[[layer]]@matrix

      inherits(matrix, "IterableMatrix")
    },
    error = function(e) FALSE
  )
}

#' Fast BPCells feature retrieval
#'
#' Internal helper: subset BPCells-backed matrix directly by accessing the
#' underlying IterableMatrix. Falls back to `fetch_data()` on any error to
#' preserve robustness. This function bypasses SCUBA's standard data retrieval
#' pipeline for large feature queries when the assay layer is BPCells-backed,
#' improving performance by leveraging direct matrix subsetting.
#'
#' @param object A Seurat object with a BPCells-backed assay layer.
#' @param features A character vector of feature names to retrieve.
#' @param assay The assay containing the BPCells-backed layer. If `NULL`, the
#'   default assay will be used.
#' @param layer The layer to retrieve from within the assay. If `NULL`, the
#'   "data" layer will be used.
#' @param cells A character vector of cell names to include. If `NULL`, all
#'   cells in the object will be included.
#'
#' @return A data.frame with cells as rows and features as columns, with feature
#'   values extracted from the BPCells-backed matrix. Missing features (those
#'   not present in the matrix) are included as columns of `NA` values.
#'
#' @keywords internal
#' @export
bp_cells_fetch_features <- function(
  object,
  features,
  assay = NULL,
  layer = NULL,
  cells = NULL
) {
  assay <- assay %||% DefaultAssay(object)
  layer <- layer %||% "data"
  cells <- cells %||% colnames(object)

  # Try direct matrix access, fall back to fetch_data on error
  expr_df <-
    tryCatch(
      {
        matrix <- object@assays[[assay]]@layers[[layer]]@matrix
        # Filter feature list to those present
        features_present <- intersect(features, rownames(matrix))
        if (length(features_present) == 0) {
          stop("No requested features present in BPCells-backed matrix.")
        }
        sub_mat <- matrix[features_present, cells, drop = FALSE]
        # Ensure dense for downstream operations
        df <- as.data.frame(t(as.matrix(sub_mat)))
        # Reintroduce any missing (non-present) features as NA
        # columns to match request
        missing <- setdiff(features, features_present)
        for (m in missing) {
          df[[m]] <- NA_real_
        }
        # Order columns to requested order
        df[, features, drop = FALSE]
      },
      error = function(e) {
        message(
          "BPCells direct access failed: ",
          conditionMessage(e),
          "\n",
          "Falling back to fetch_data."
        )
        SCUBA::fetch_data(
          object = object,
          vars = features,
          layer = layer,
          cells = cells
        )
      }
    )
  expr_df
}
