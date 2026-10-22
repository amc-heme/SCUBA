# File-backed sources served by the Julia engine --------------------------------------------
#
# `open_h5ad()` returns a lightweight `scuba_source` object: the file path, the full metadata
# table (read once from the file) and the schema (modalities, features, layers, reductions).
# Expression and reduction values are fetched through SCUBA.jl's native reader, keyed by path, so
# modified copies of the object (metadata updates, cell subsets) never invalidate the Julia-side
# source. Metadata lives in R, where cellDIVER-style `fetch_metadata(full_table = TRUE)` ->
# modify -> `update_object_metadata()` round trips cost nothing. R's copy-on-modify semantics
# make metadata updates and subsets per-object, unlike reticulate AnnData objects, which are
# mutated in place.

#' Open an .h5ad or .h5mu file as a SCUBA object (Julia engine)
#'
#' Reads the schema and the cell metadata of an AnnData `.h5ad` file (or a MuData `.h5mu` file
#' with `open_h5mu()`) and returns a `scuba_source` object that every SCUBA accessor accepts.
#' Expression values and reductions are read from the file on demand by the Julia engine
#' (SCUBA.jl), so no Python environment is needed. Requires `scuba_engine()` to be `"julia"`.
#'
#' Keys follow SCUBA.jl's rules. For `.h5ad`: the main matrix is modality `"X"` (features
#' addressed bare or as `X_<gene>`, layers `"X"` and `layers/*`), each `obsm` data frame is a
#' modality named after it (e.g. `protein_CD117-AB`), each numeric `obsm` array is a reduction
#' (`X_umap_1`), and metadata columns may be given bare or as `obs_<column>`. For `.h5mu`:
#' modalities are named after the MuData modalities (`rna_GAPDH`), reductions as
#' `<modality>_<key>` (`rna_X_umap_1`, bare `X_umap_1` when unique), metadata is the global
#' `obs` plus each modality's `obs` as `<modality>_<column>`; cells absent from a modality are
#' `NaN` / `NA`.
#'
#' @param path Path to an `.h5ad` (`open_h5ad()`) or `.h5mu` (`open_h5mu()`) file.
#' @returns An object of class `scuba_source`.
#' @seealso [update_object_metadata()], [scuba_source_modalities()]
#' @export
open_h5ad <- function(path) .open_scuba_source(path, "h5ad", "open_h5ad()")

#' @rdname open_h5ad
#' @export
open_h5mu <- function(path) .open_scuba_source(path, "h5mu", "open_h5mu()")

.open_scuba_source <- function(path, format, what) {
  path <- path.expand(path)
  if (!file.exists(path)) stop(what, ": file not found: ", path, call. = FALSE)
  path <- normalizePath(path)
  .scuba_source_engine(what)
  if (!isTRUE(JuliaCall::julia_eval("isdefined(ScubaR, :source_open)"))) {
    stop(what, " needs a SCUBA.jl with file-backed glue (ScubaR.source_open); update the SCUBA.jl checkout",
         call. = FALSE)
  }
  info <- .jl_call("source_open", path)
  cells <- as.character(info$cells)
  if (anyDuplicated(cells)) {
    stop(what, ": duplicated cell names (obs_names) in ", path, "; SCUBA objects need unique cell names",
         call. = FALSE)
  }
  meta <- .jl_table(info$obs)
  rownames(meta) <- cells
  mods <- as.character(info$modalities)
  as_list <- function(x, nm) {
    # arrives as a named list (Julia Dict); order it like `modalities`
    x <- lapply(nm, function(m) as.character(x[[m]]))
    names(x) <- nm
    x
  }
  reds <- as.character(info$reductions)
  structure(list(
    path = path,
    format = format,
    cells = NULL,                                   # NULL = every cell; a character vector after `[`
    meta = meta,
    info = list(
      modalities = mods,
      features = as_list(info$features, mods),
      layers = as_list(info$layers, mods),
      reductions = stats::setNames(as.integer(info$red_ndims), reds),
      default_reduction = if (nzchar(info$default_reduction)) info$default_reduction else NA_character_,
      default_layer = as.character(info$default_layer),
      uns_assays = if (length(info$uns_assays)) as.character(info$uns_assays) else NULL,
      uns_reductions = if (length(info$uns_reductions)) as.character(info$uns_reductions) else NULL
    )
  ), class = "scuba_source")
}

.scuba_source_engine <- function(what = "scuba_source objects") {
  if (!identical(scuba_engine(), "julia")) {
    stop(what, " need", if (what == "scuba_source objects") "" else "s", " the Julia engine (options(SCUBA.engine) is \"",
         getOption("SCUBA.engine", "auto"), "\" and Julia is ", if (isTRUE(scuba_julia_available(quiet = TRUE))) "available" else "not available",
         "). Install Julia and JuliaCall, or load the file with anndata::read_h5ad() for the reticulate path.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Cells of the object (all, or the subset recorded by `[`), optionally validated against a request.
.src_cells <- function(object, cells = NULL) {
  all <- object$cells %||% rownames(object$meta)
  if (is.null(cells)) return(all)
  cells <- as.character(cells)
  bad <- setdiff(cells, all)
  if (length(bad)) {
    stop("The following cells were not found in the object: ", paste(utils::head(bad, 10), collapse = ", "),
         if (length(bad) > 10) paste0(" (and ", length(bad) - 10, " more)"), call. = FALSE)
  }
  cells
}

#' @export
`[.scuba_source` <- function(x, i, j, ..., drop = FALSE) {
  if (!missing(j)) warning("Column subsetting is ignored for scuba_source objects", call. = FALSE)
  if (missing(i)) return(x)
  all <- get_all_cells(x)
  cells <- if (is.character(i)) .src_cells(x, i) else all[i]
  if (anyNA(cells)) stop("subscript out of bounds", call. = FALSE)
  x$cells <- cells
  x$meta <- x$meta[cells, , drop = FALSE]
  x
}

#' @export
dim.scuba_source <- function(x) c(length(x$info$features[[1]]), length(get_all_cells(x)))

#' @export
dimnames.scuba_source <- function(x) list(x$info$features[[1]], get_all_cells(x))

#' @export
print.scuba_source <- function(x, ...) {
  cat("SCUBA ", x$format, " source: ", x$path, "\n", sep = "")
  cat("  ", nrow(x), " features (", x$info$modalities[1], ") x ", ncol(x), " cells",
      if (!is.null(x$cells)) " (subset)" else "", "\n", sep = "")
  cat("  modalities: ", paste(x$info$modalities, collapse = ", "), "\n", sep = "")
  cat("  reductions: ", if (length(x$info$reductions)) paste(names(x$info$reductions), collapse = ", ") else "none", "\n", sep = "")
  cat("  metadata: ", ncol(x$meta), " columns\n", sep = "")
  invisible(x)
}

#' Modalities and reductions of a file-backed SCUBA object
#'
#' @param object A `scuba_source` object from [open_h5ad()].
#' @returns `scuba_source_modalities()`: a character vector of modality keys (main first).
#'   `scuba_source_reductions()`: a named integer vector of reduction keys and their dimension counts.
#' @export
scuba_source_modalities <- function(object) object$info$modalities

#' @rdname scuba_source_modalities
#' @export
scuba_source_reductions <- function(object) object$info$reductions

# --- accessors --------------------------------------------------------------------------------

#' @export
fetch_data.scuba_source <- function(object, vars, layer = NULL, cells = NULL, slot = NULL, ...) {
  .scuba_source_engine()
  vars <- as.character(vars)
  if (!length(vars)) stop("`vars` must contain at least one variable.", call. = FALSE)
  if (anyDuplicated(vars)) {
    warning("Duplicate entries in `vars` were dropped: ", paste(unique(vars[duplicated(vars)]), collapse = ", "), call. = FALSE)
    vars <- unique(vars)
  }
  cells <- .src_cells(object, cells)
  meta_names <- names(object$meta)
  bare <- sub("^obs_", "", vars)
  is_meta <- vars %in% meta_names | (startsWith(vars, "obs_") & bare %in% meta_names)
  out_meta <- NULL
  if (any(is_meta)) {
    out_meta <- object$meta[cells, bare[is_meta], drop = FALSE]
    names(out_meta) <- bare[is_meta]
    shadow <- intersect(bare[is_meta], object$info$features[[1]])
    if (length(shadow)) {
      warning("The following variables are both metadata columns and features of the main modality and were ",
              "returned as metadata: ", paste(shadow, collapse = ", "), ". Use the modality key for the feature (e.g. \"",
              object$info$modalities[1], "_", shadow[1], "\").", call. = FALSE)
    }
  }
  out_expr <- NULL; inputs <- character()
  if (any(!is_meta)) {
    res <- .jl_call("fetch_data", object$path, vars[!is_meta], slot %||% layer, cells, .scuba_julia_eltype(), any(is_meta))
    out_expr <- .jl_table(res)
    inputs <- as.character(res$inputs)
  }
  out <- if (is.null(out_expr)) out_meta else if (is.null(out_meta)) out_expr else cbind(out_expr, out_meta)
  requested <- c(inputs, vars[is_meta])
  out <- out[, order(match(requested, vars)), drop = FALSE]
  rownames(out) <- cells
  out
}

#' @export
fetch_feature.scuba_source <- function(object, features, assay = NULL, layer = NULL, cells = NULL) {
  .scuba_source_engine()
  cells <- .src_cells(object, cells)
  res <- .jl_call("fetch_feature", object$path, as.character(features), layer, cells, assay, .scuba_julia_eltype())
  out <- .jl_table(res)
  rownames(out) <- cells
  out
}

#' @export
fetch_metadata.scuba_source <- function(object, vars = NULL, cells = NULL, full_table = FALSE,
                                        return_class = "dataframe") {
  if (!return_class %in% c("dataframe", "vector")) {
    stop('return_class must be either "dataframe" or "vector".')
  }
  if (full_table == FALSE && is.null(vars)) {
    stop("`vars` must be defined, unless `full_table` is TRUE.")
  }
  if (full_table == TRUE && !is.null(vars)) {
    warning("`full_table` is TRUE, ignoring entries in `vars`.")
  }
  if (full_table == TRUE) {
    return(object$meta[.src_cells(object, cells), , drop = FALSE])
  }
  vars_present <- names(object$meta)
  not_found <- !vars %in% vars_present
  if (any(not_found)) {
    if (all(not_found)) {
      stop(paste0("\nNone of the variables entered in `vars` were found\nin the object. Variables entered: \"",
                  paste(vars, collapse = '", "'), "\". To view available entries for your object,\nrun SCUBA::meta_varnames()."))
    } else {
      stop(paste0("\nThe following variables entered in `vars` were not found in \nthe object: ",
                  paste(vars[not_found]), ". \n\nTo view available entries for your object, run SCUBA::meta_varnames()."))
    }
  }
  cells <- .src_cells(object, cells)
  data <- object$meta[cells, vars]
  if (return_class == "dataframe") {
    data <- as.data.frame(data)
    rownames(data) <- cells
    colnames(data) <- vars
  } else {
    names(data) <- cells
  }
  data
}

#' @export
fetch_reduction.scuba_source <- function(object, reduction, cells = NULL, dims = c(1, 2)) {
  .scuba_source_engine()
  if (length(dims) != 2) stop("'dims' must be a two-length vector")
  reds <- names(object$info$reductions)
  if (!reduction %in% reds) {
    stop(paste0('\nThe reduction "', reduction, '" was not found in the object passsed. \nReductions present in object: ',
                paste(reds, collapse = ", "), "."))
  }
  cells <- .src_cells(object, cells)
  out <- .jl_fetch_reduction(object$path, reduction, cells, dims)
  rownames(out) <- cells
  out
}

#' @export
get_all_cells.scuba_source <- function(object) object$cells %||% rownames(object$meta)

#' @export
fetch_cells.scuba_source <- function(object, meta_var, meta_levels) {
  col <- object$meta[[meta_var]]
  if (is.null(col)) stop("fetch_cells: metadata variable ", meta_var, " not found in object.", call. = FALSE)
  if (!inherits(col, c("character", "factor", "logical"))) {
    stop(meta_var, " is not a categorical variable. Only categorical variables are supported by fetch_cells at this time.")
  }
  rownames(object$meta)[col %in% meta_levels]
}

#' @export
meta_varnames.scuba_source <- function(object) names(object$meta)

#' @export
all_keys.scuba_source <- function(object) {
  k <- c(object$info$modalities, names(object$info$reductions))
  stats::setNames(k, k)
}

#' @export
features_in_assay.scuba_source <- function(object, assay) {
  if (!assay %in% object$info$modalities) stop("Assay ", assay, " not found in the current object.")
  object$info$features[[assay]]
}

#' @export
default_layer.scuba_source <- function(object) object$info$default_layer

#' @export
default_reduction.scuba_source <- function(object) {
  if (is.na(object$info$default_reduction)) {
    stop('Unable to find a reduction matching "umap", "tsne", or "pca". Please specify the reduction to use via the `reduction` parameter.')
  }
  object$info$default_reduction
}

#' @export
reduction_dimnames.scuba_source <- function(object, reduction, dims) {
  reds <- object$info$reductions
  if (!reduction %in% names(reds)) stop("Reduction ", reduction, " not found in the current object.", call. = FALSE)
  if (any(dims > reds[[reduction]])) stop("Reduction ", reduction, " has ", reds[[reduction]], " dimensions.", call. = FALSE)
  paste0(reduction, "_", dims)
}

#' @export
relabel_axis.scuba_source <- function(object, feature, ...) {
  reds <- names(object$info$reductions)
  for (r in reds) {
    if (grepl(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", r), "_[0-9]+$"), feature)) return("Embeddings Value")
  }
  key <- sub("_.*$", "", feature)
  if (key %in% object$info$modalities) return(NULL)
  if (feature %in% object$info$features[[1]]) return(NULL)
  ""
}

#' @export
update_object_metadata.scuba_source <- function(object, table) {
  table <- as.data.frame(table)
  all <- get_all_cells(object)
  if (nrow(table) != length(all)) {
    stop("`table` has ", nrow(table), " rows; the object has ", length(all), " cells.", call. = FALSE)
  }
  rn <- rownames(table)
  if (is.null(rn) || !setequal(rn, all)) {
    stop("The row names of `table` must be the object's cell names (see get_all_cells()).", call. = FALSE)
  }
  object$meta <- table[all, , drop = FALSE]
  object
}
