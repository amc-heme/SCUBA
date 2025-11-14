## BPCells utility helpers (internal)
## Detection & fast feature access for BPCells-backed Seurat assays
##
## These functions enable SCUBA plotting functions to bypass SeuratObject::FetchData
## for large feature queries when the underlying assay layers are BPCells-backed
## disk matrices. They are intentionally conservative: failures or missing
## dependencies fall back to standard FetchData behavior.

#' Detect BPCells-backed assay layer in a Seurat object
#'
#' Internal helper: returns TRUE if the indicated assay/layer appears to be
#' backed by a BPCells on-disk matrix (heuristic: presence of @matrix@dir slot).
#'
#' @keywords internal
is_bp_cells_seurat <- function(object, assay = NULL, layer = NULL){
  if (!inherits(object, "Seurat")) return(FALSE)
  assay <- assay %||% DefaultAssay(object)
  layer <- layer %||% "data"
  out <- try({
    lyr <- object@assays[[assay]]@layers[[layer]]
    dir_val <- try(lyr@matrix@dir, silent = TRUE)
    !inherits(dir_val, "try-error") && !is.null(dir_val)
  }, silent = TRUE)
  if (inherits(out, "try-error")) return(FALSE) else return(isTRUE(out))
}

#' Fast BPCells feature retrieval
#'
#' Internal helper: subset BPCells-backed matrix directly. Falls back to
#' FetchData on any error to preserve robustness.
#'
#' @keywords internal
bp_cells_fetch_features <- function(object, features, assay = NULL, layer = NULL, cells = NULL){
  assay <- assay %||% DefaultAssay(object)
  layer <- layer %||% "data"
  cells <- cells %||% colnames(object)
  # Try direct matrix access
  expr_df <- try({
    mat_obj <- object@assays[[assay]]@layers[[layer]]@matrix
    # Filter feature list to those present
    features_present <- intersect(features, rownames(mat_obj))
    if (length(features_present) == 0){
      stop("No requested features present in BPCells-backed matrix.")
    }
    sub_mat <- mat_obj[features_present, cells, drop = FALSE]
    # Ensure dense for downstream operations
    df <- as.data.frame(t(as.matrix(sub_mat)))
    # Reintroduce any missing (non-present) features as NA columns to match request
    missing <- setdiff(features, features_present)
    for (m in missing){ df[[m]] <- NA_real_ }
    # Order columns to requested order
    df <- df[, features, drop = FALSE]
    df
  }, silent = TRUE)
  if (inherits(expr_df, "try-error")){
    # Fallback
    expr_df <- FetchData(object = object, vars = features, layer = layer, cells = cells)
  }
  expr_df
}
