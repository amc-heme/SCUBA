# Pseudobulk aggregation ---------------------------------------------------------------------
#
# `pseudobulk()` sums (or averages) a counts layer over profiles = sample x group, for every
# object class SCUBA supports, without ever materialising the full features x cells matrix:
#
#   - Julia engine (Seurat, in-memory SingleCellExperiment, scuba_source): SCUBA.jl streams the
#     layer into a group accumulator (`ScubaR.pseudobulk`); only the small features x profiles
#     matrix crosses back into R.
#   - R engine: dgCMatrix / dense layers use one sparse indicator-matrix product, BPCells layers
#     use `BPCells::pseudobulk_matrix()`, DelayedArray / HDF5-backed layers use
#     `DelayedArray::colsum()` (block-wise), anything else a bounded block loop.
#   - AnnData (reticulate): `inst/extdata/Python/aggregate_anndata.py` aggregates in cell chunks,
#     in memory or backed.
#
# The grouping is always resolved in R (`.pseudobulk_codes()`), from metadata columns or from
# cell-named label vectors, and handed to the backend as integer codes. That is what lets
# R-only metadata columns (e.g. cellDIVER's `metacluster` on a `scuba_source`) drive the
# aggregation on every engine, and it fixes the profile order independently of the backend.

#' Pseudobulk profiles: aggregate a counts layer over groups of cells
#'
#' Sums (or averages) every feature of a layer over *profiles*, the combinations of a sample
#' identifier and, optionally, a group label (cell type, condition, cluster). The result is the
#' small features x profiles matrix that pseudobulk differential expression methods (edgeR,
#' DESeq2) consume, together with a table describing each profile.
#'
#' The full features x cells matrix is never built: on the Julia engine SCUBA.jl streams the
#' layer into a group accumulator, and on the R engine each matrix class is aggregated with its
#' own block-wise primitive (sparse indicator product for `dgCMatrix` and dense layers,
#' `BPCells::pseudobulk_matrix()` for BPCells layers, `DelayedArray::colsum()` for
#' HDF5-backed layers). AnnData objects (reticulate) are aggregated in cell chunks by a Python
#' helper, in memory or in backed mode.
#'
#' Cells whose sample or group value is `NA` or blank are left out (counted in
#' `n_cells_dropped`). Profiles are ordered sample-major, then by group; levels follow the
#' factor levels of the metadata column when it is a factor, and sorted order (C locale)
#' otherwise, so every engine and backend returns the same columns.
#'
#' @inheritParams object_param
#' @inheritParams cells_param
#' @param sample_by The sample (biological replicate) of each cell: the name of a metadata
#'   column, or a character/factor vector named by cell (one label per cell of the object, or
#'   at least per cell in `cells`). A vector lets a caller aggregate by labels that are not
#'   stored in the object.
#' @param group_by Optional group of each cell (same forms as `sample_by`). With `NULL` the
#'   profiles are the samples themselves.
#' @param ... Passed to methods.
#' @param assay Assay / modality to aggregate (default: the default assay of the object; the
#'   main modality `"X"` for AnnData and `.h5ad` objects).
#' @param layer Layer to aggregate. Defaults to `"counts"`, *not* to [default_layer()]: the
#'   default layer holds normalised values, which must not be summed. AnnData objects whose raw
#'   counts live in `X` need `layer = "X"` (or `"raw"` for `adata.raw`).
#' @param features Features to include (default: every feature of the assay).
#' @param fun `"sum"` (default) or `"mean"`; the mean divides each profile's sum by the number of
#'   cells that contributed to it.
#'
#' @returns An object of class `scuba_pseudobulk`: a list with
#'   - `counts`: features x profiles numeric matrix (sums, or means with `fun = "mean"`),
#'   - `profiles`: data.frame with one row per profile (`profile_id`, `sample`, `group`,
#'     `n_cells`), in the column order of `counts`; `group` is `NA` when `group_by` is `NULL`,
#'   - `assay`, `layer`, `fun`, `sample_by`, `group_by` (the column names, or `NA` for vectors),
#'   - `n_cells_dropped`: cells left out for a missing sample or group label,
#'   - `engine`: `"julia"`, `"r"` or `"python"`.
#'
#' @seealso [pseudobulk_design()] for the metadata-only preview of the profiles,
#'   [layer_names()] to find a counts layer, [fetch_matrix()].
#' @export
#'
#' @examples
#' pb <- pseudobulk(AML_Seurat, sample_by = "Batch", group_by = "condensed_cell_type")
#' pb
#' head(pb$profiles)
#' pb$counts[1:3, 1:3]
#'
#' # Samples only, as means
#' pseudobulk(AML_SCE(), sample_by = "Batch", fun = "mean")$profiles
pseudobulk <-
  function(
    object,
    sample_by,
    group_by = NULL,
    ...,
    assay = NULL,
    layer = "counts",
    cells = NULL,
    features = NULL,
    fun = c("sum", "mean")
  ) {
    UseMethod("pseudobulk")
  }

#' @export
pseudobulk.default <-
  function(
    object,
    sample_by,
    group_by = NULL,
    ...
  ) {
    stop(
      "pseudobulk() does not know how to handle an object of class ",
      paste(class(object), collapse = ", "),
      ". Supported classes: Seurat, SingleCellExperiment, AnnDataR6 and scuba_source ",
      "(from open_h5ad() / open_h5mu()).",
      call. = FALSE
    )
  }

#' @rdname pseudobulk
#' @export
pseudobulk.Seurat <-
  function(
    object,
    sample_by,
    group_by = NULL,
    ...,
    assay = NULL,
    layer = "counts",
    cells = NULL,
    features = NULL,
    fun = c("sum", "mean")
  ) {
    fun <- match.arg(fun)
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    if (!assay %in% SeuratObject::Assays(object)) {
      stop("Assay '", assay, "' not found; assays: ", paste(SeuratObject::Assays(object), collapse = ", "), call. = FALSE)
    }
    assay_object <- object[[assay]]
    layers <- .seurat_counts_layers(assay_object, layer)
    features <- .pseudobulk_features(features, rownames(assay_object))
    design <- .pseudobulk_codes(object, sample_by, group_by, cells)
    n_profiles <- nrow(design$profiles)
    total <- matrix(0, length(features), n_profiles)
    n <- integer(n_profiles)
    engine <- "r"
    use_julia <- .scuba_julia_object(object)
    for (one_layer in layers) {
      if (use_julia) {
        # Split layers (counts.1, counts.2, ...) hold disjoint cells: aggregate each over the
        # cells it holds and add the accumulators. A single layer is sent whole (NULL cells =
        # object order) unless the caller subset the cells.
        if (length(layers) == 1) {
          layer_cells <- design$cells
          jl_cells <- if (is.null(cells)) NULL else design$cells
        } else {
          layer_cells <- intersect(design$cells, SeuratObject::Cells(assay_object, layer = one_layer))
          jl_cells <- layer_cells
        }
        if (!length(layer_cells)) next
        codes <- design$codes[match(layer_cells, design$cells)]
        res <- .jl_pseudobulk(object, codes, n_profiles, assay = assay, layer = one_layer,
                              cells = jl_cells, features = features)
        res$sum <- res$sum[match(features, res$features), , drop = FALSE]
        engine <- "julia"
      } else {
        res <- .r_pseudobulk_matrix(SeuratObject::LayerData(assay_object, layer = one_layer), design, features,
                                    partial = length(layers) > 1)
      }
      total <- total + res$sum
      n <- n + res$n
    }
    .pseudobulk_result(total, n, design, features, assay, layer, fun, sample_by, group_by, engine)
  }

#' @rdname pseudobulk
#' @export
pseudobulk.SingleCellExperiment <-
  function(
    object,
    sample_by,
    group_by = NULL,
    ...,
    assay = NULL,
    layer = "counts",
    cells = NULL,
    features = NULL,
    fun = c("sum", "mean")
  ) {
    fun <- match.arg(fun)
    main <- SingleCellExperiment::mainExpName(object)
    experiments <- c(main %||% "main", SingleCellExperiment::altExpNames(object))
    is_main <- is.null(assay) || identical(assay, main) || (is.null(main) && identical(assay, "main"))
    if (!is_main && !assay %in% SingleCellExperiment::altExpNames(object)) {
      stop("Assay '", assay, "' not found; experiments: ", paste(experiments, collapse = ", "), call. = FALSE)
    }
    experiment <- if (is_main) object else SingleCellExperiment::altExp(object, assay)
    assay <- if (is_main) main %||% "main" else assay
    available <- SummarizedExperiment::assayNames(experiment)
    if (!layer %in% available) {
      stop("Layer '", layer, "' not found in experiment '", assay, "'; available: ",
           paste(available, collapse = ", "), ". Pass the layer holding raw counts.", call. = FALSE)
    }
    m <- SummarizedExperiment::assay(experiment, layer)
    features <- .pseudobulk_features(features, rownames(experiment))
    design <- .pseudobulk_codes(object, sample_by, group_by, cells)
    # HDF5-backed assays stay on the R engine even when Julia is selected: Julia's R-source path
    # densifies every slice in R, whereas DelayedArray::colsum() streams the HDF5 chunks.
    if (.scuba_julia_object(object) && !methods::is(m, "DelayedArray")) {
      res <- .jl_pseudobulk(object, design$codes, nrow(design$profiles), assay = assay, layer = layer,
                            cells = if (is.null(cells)) NULL else design$cells, features = features)
      res$sum <- res$sum[match(features, res$features), , drop = FALSE]
      engine <- "julia"
    } else {
      res <- .r_pseudobulk_matrix(m, design, features)
      engine <- "r"
    }
    .pseudobulk_result(res$sum, res$n, design, features, assay, layer, fun, sample_by, group_by, engine)
  }

#' @rdname pseudobulk
#' @export
pseudobulk.scuba_source <-
  function(
    object,
    sample_by,
    group_by = NULL,
    ...,
    assay = NULL,
    layer = "counts",
    cells = NULL,
    features = NULL,
    fun = c("sum", "mean")
  ) {
    fun <- match.arg(fun)
    .scuba_source_engine("pseudobulk()")
    assay <- assay %||% object$info$modalities[1]
    if (!assay %in% object$info$modalities) {
      stop("Modality '", assay, "' not found; modalities: ", paste(object$info$modalities, collapse = ", "), call. = FALSE)
    }
    available <- object$info$layers[[assay]]
    if (!layer %in% available) {
      stop("Layer '", layer, "' not found in modality '", assay, "'; available: ",
           paste(available, collapse = ", "), ". Pass the layer holding raw counts (e.g. layer = \"X\" ",
           "when X holds counts).", call. = FALSE)
    }
    features <- .pseudobulk_features(features, object$info$features[[assay]])
    design <- .pseudobulk_codes(object, sample_by, group_by, cells)
    # A whole object is aggregated in file order (NULL cells); a subset object or a `cells`
    # request sends the cell names, which Julia resolves in that order (matching `codes`).
    jl_cells <- if (is.null(cells)) object$cells else design$cells
    res <- .jl_pseudobulk(object$path, design$codes, nrow(design$profiles), assay = assay, layer = layer,
                          cells = jl_cells, features = features)
    res$sum <- res$sum[match(features, res$features), , drop = FALSE]
    .pseudobulk_result(res$sum, res$n, design, features, assay, layer, fun, sample_by, group_by, "julia")
  }

#' @rdname pseudobulk
#' @export
pseudobulk.AnnDataR6 <-
  function(
    object,
    sample_by,
    group_by = NULL,
    ...,
    assay = NULL,
    layer = "counts",
    cells = NULL,
    features = NULL,
    fun = c("sum", "mean")
  ) {
    fun <- match.arg(fun)
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop('Package "reticulate" must be installed to use pseudobulk() with anndata objects.', call. = FALSE)
    }
    design <- .pseudobulk_codes(object, sample_by, group_by, cells)
    script <- system.file("extdata", "Python", "aggregate_anndata.py", package = "SCUBA")
    py <- reticulate::py_run_file(script)
    res <- py$aggregate_anndata(
      obj = object,
      codes = as.integer(design$codes),
      n_groups = as.integer(nrow(design$profiles)),
      modality = assay,
      layer = layer,
      cells = if (is.null(cells)) NULL else as.list(design$cells),
      features = if (is.null(features)) NULL else as.list(as.character(features)),
      chunk_cells = as.integer(getOption("SCUBA.pseudobulk_chunk", 5000L))
    )
    total <- t(as.matrix(res$sum))
    if (!is.matrix(total) || nrow(total) != length(res$features)) {
      total <- matrix(as.numeric(res$sum), nrow = length(res$features), ncol = nrow(design$profiles), byrow = TRUE)
    }
    features <- as.character(res$features)
    .pseudobulk_result(total, as.integer(res$n), design, features, assay %||% "X", layer, fun, sample_by, group_by, "python")
  }

# anndata >= 0.8 (R package) returns plain reticulate objects of class "anndata.AnnData" instead
# of the AnnDataR6 wrappers; the same implementation serves both.
#' @rdname pseudobulk
#' @export
pseudobulk.anndata.AnnData <- pseudobulk.AnnDataR6

.is_anndata <- function(object) inherits(object, c("AnnDataR6", "anndata.AnnData"))

# Helpers that work on both anndata R classes: AnnDataR6 (anndata <= 0.7, R6 wrappers that hand
# back R values) and anndata.AnnData (anndata >= 0.8, reticulate objects created without
# conversion, so every attribute is a Python object until converted explicitly).
.py_chr <- function(x) {
  if (is.character(x)) return(x)
  as.character(unlist(reticulate::py_to_r(x)))
}
.anndata_cells <- function(object) {
  ids <- object$obs_names
  if (!is.character(ids)) ids <- ids$to_list()
  .py_chr(ids)
}
.anndata_keys <- function(mapping) {
  keys <- tryCatch(mapping$keys(), error = function(e) NULL)
  if (is.character(keys)) return(keys[!is.na(keys)])
  builtins <- reticulate::import_builtins(convert = TRUE)
  # a Mapping iterates over its keys; `keys()` is not callable on every reticulate wrapper
  keys <- if (is.null(keys)) builtins$list(mapping) else builtins$list(keys)
  keys <- .py_chr(keys)
  keys[!is.na(keys) & nzchar(keys)]
}
# One obs column as a cell-named R vector (a factor with the pandas categories when the column
# is categorical), for the requested cells.
.anndata_metadata_vector <- function(object, var, cells) {
  obs <- object$obs
  if (is.data.frame(obs)) {
    if (!var %in% names(obs)) stop("Metadata column '", var, "' not found; columns: ", paste(utils::head(names(obs), 20), collapse = ", "), call. = FALSE)
    v <- obs[[var]]
    names(v) <- rownames(obs)
    return(v[cells])
  }
  columns <- .py_chr(obs$columns$to_list())
  if (!var %in% columns) stop("Metadata column '", var, "' not found; columns: ", paste(utils::head(columns, 20), collapse = ", "), call. = FALSE)
  col <- obs[[var]]
  is_na <- as.logical(.py_chr(col$isna()$to_list()) == "TRUE")
  values <- .py_chr(col$astype("str")$to_list())
  values[is_na] <- NA
  categories <- tryCatch(.py_chr(col$cat$categories$to_list()), error = function(e) NULL)
  v <- if (!is.null(categories)) factor(values, levels = categories) else values
  names(v) <- .anndata_cells(object)
  v[cells]
}

# --- design preview ------------------------------------------------------------------------

#' Preview the pseudobulk profiles of an object
#'
#' Metadata-only companion of [pseudobulk()]: resolves the same sample x group profiles (same
#' rows, same order) and reports how many cells each one has, without reading any expression
#' data. Use it to check replication before aggregating, e.g. whether every group has enough
#' samples, or whether samples contribute cells to several groups (a paired design).
#'
#' @inheritParams pseudobulk
#' @returns An object of class `scuba_pseudobulk_design`: a list with `profiles` (as in
#'   [pseudobulk()]), `table` (samples x groups matrix of cell counts; one column when
#'   `group_by` is `NULL`), `paired` (`TRUE` when any sample has cells in more than one group),
#'   `samples_per_group` (samples with at least one cell per group), `samples_in_all_groups`,
#'   `n_cells_dropped`, `n_cells` (cells assigned to a profile), `sample_by` and `group_by`.
#' @export
#'
#' @examples
#' d <- pseudobulk_design(AML_Seurat, "Batch", "condensed_cell_type")
#' d$table
#' d$paired
pseudobulk_design <- function(object, sample_by, group_by = NULL, cells = NULL) {
  design <- .pseudobulk_codes(object, sample_by, group_by, cells)
  profiles <- design$profiles
  samples <- unique(profiles$sample)
  if (is.null(group_by)) {
    tab <- matrix(profiles$n_cells, ncol = 1, dimnames = list(profiles$sample, "all"))
  } else {
    groups <- unique(profiles$group)
    tab <- matrix(0L, length(samples), length(groups), dimnames = list(samples, groups))
    tab[cbind(match(profiles$sample, samples), match(profiles$group, groups))] <- profiles$n_cells
  }
  present <- tab > 0
  structure(list(
    profiles = profiles,
    table = tab,
    paired = any(rowSums(present) > 1),
    samples_per_group = colSums(present),
    samples_in_all_groups = sum(rowSums(present) == ncol(tab)),
    n_cells_dropped = design$n_cells_dropped,
    n_cells = sum(profiles$n_cells),
    sample_by = .pseudobulk_label_name(sample_by),
    group_by = if (is.null(group_by)) NA_character_ else .pseudobulk_label_name(group_by)
  ), class = "scuba_pseudobulk_design")
}

#' @export
print.scuba_pseudobulk_design <- function(x, ...) {
  grouped <- !all(is.na(x$profiles$group))
  cat("SCUBA pseudobulk design: ", nrow(x$profiles), " profiles (", nrow(x$table), " samples",
      if (grouped) paste0(" x ", ncol(x$table), " groups") else "", ")\n", sep = "")
  cat("  sample_by: ", x$sample_by, if (grouped) paste0(", group_by: ", x$group_by) else "", "\n", sep = "")
  cat("  cells assigned: ", x$n_cells, "; dropped: ", sum(x$n_cells_dropped), "\n", sep = "")
  if (grouped) {
    cat("  design: ", if (x$paired) "paired (samples contribute cells to several groups)" else "independent samples", "\n", sep = "")
    cat("  samples per group: ", paste(names(x$samples_per_group), x$samples_per_group, sep = " = ", collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

# --- S3 methods for the result -------------------------------------------------------------

#' @export
print.scuba_pseudobulk <- function(x, ...) {
  cat("SCUBA pseudobulk: ", nrow(x$counts), " features x ", ncol(x$counts), " profiles (",
      if (x$fun == "sum") "sum" else "mean", " of layer '", x$layer, "', assay '", x$assay, "', ", x$engine, " engine)\n", sep = "")
  cat("  samples: ", length(unique(x$profiles$sample)),
      if (!all(is.na(x$profiles$group))) paste0("; groups: ", length(unique(x$profiles$group))) else "",
      "; cells: ", sum(x$profiles$n_cells), " (dropped ", sum(x$n_cells_dropped), ")\n", sep = "")
  invisible(x)
}

#' @export
dim.scuba_pseudobulk <- function(x) dim(x$counts)

#' @export
dimnames.scuba_pseudobulk <- function(x) dimnames(x$counts)

#' @export
as.matrix.scuba_pseudobulk <- function(x, ...) x$counts

# --- shared internals ----------------------------------------------------------------------

# Is `x` the name of a metadata column (a single unnamed string) rather than a label vector?
.pseudobulk_is_column <- function(x) is.character(x) && length(x) == 1 && is.null(names(x))
.pseudobulk_label_name <- function(x) if (.pseudobulk_is_column(x)) x else NA_character_

# Per-cell labels for `sample_by` / `group_by`: a metadata column fetched for `cells`, or a
# cell-named vector indexed by `cells`.
.pseudobulk_labels <- function(object, x, cells, arg) {
  if (.pseudobulk_is_column(x)) {
    if (.is_anndata(object)) return(.anndata_metadata_vector(object, x, cells))
    v <- fetch_metadata(object, vars = x, cells = cells, return_class = "vector")
    if (length(v) != length(cells)) stop("fetch_metadata() returned ", length(v), " values for ", length(cells), " cells", call. = FALSE)
    return(v)
  }
  if (is.null(names(x))) {
    stop("`", arg, "` must be the name of a metadata column or a vector of labels named by cell.", call. = FALSE)
  }
  missing_cells <- setdiff(cells, names(x))
  if (length(missing_cells)) {
    stop("`", arg, "` has no label for ", length(missing_cells), " cell(s): ",
         paste(utils::head(missing_cells, 5), collapse = ", "), call. = FALSE)
  }
  x[cells]
}

# Observed levels in a deterministic order: factor levels (unused dropped) or sorted (C locale).
.pseudobulk_levels <- function(x, observed) {
  if (is.factor(x)) return(levels(x)[levels(x) %in% observed])
  sort(unique(observed), method = "radix")
}

.pseudobulk_features <- function(features, available) {
  if (is.null(features)) return(as.character(available))
  features <- as.character(features)
  missing_f <- setdiff(features, available)
  if (length(missing_f)) stop("Features not found: ", paste(utils::head(missing_f, 10), collapse = ", "), call. = FALSE)
  features
}

# Resolve the cells, their sample / group labels and the profile code of every cell.
# Returns list(cells, codes (integer, 0 = dropped), profiles, n_cells_dropped).
.pseudobulk_codes <- function(object, sample_by, group_by, cells) {
  all_cells <- if (.is_anndata(object)) .anndata_cells(object) else as.character(get_all_cells(object))
  if (is.null(cells)) {
    cells <- all_cells
  } else {
    cells <- as.character(cells)
    missing_cells <- setdiff(cells, all_cells)
    if (length(missing_cells)) {
      stop("Cells not found: ", paste(utils::head(missing_cells, 10), collapse = ", "), call. = FALSE)
    }
  }
  sample <- .pseudobulk_labels(object, sample_by, cells, "sample_by")
  s <- as.character(sample)
  drop_s <- is.na(s) | !nzchar(trimws(s))
  if (is.null(group_by)) {
    drop_g <- rep(FALSE, length(cells))
    g <- NULL
  } else {
    group <- .pseudobulk_labels(object, group_by, cells, "group_by")
    g <- as.character(group)
    drop_g <- is.na(g) | !nzchar(trimws(g))
  }
  keep <- !(drop_s | drop_g)
  sample_levels <- .pseudobulk_levels(sample, s[keep])
  n_groups <- 1L
  group_levels <- NULL
  if (!is.null(g)) {
    group_levels <- .pseudobulk_levels(group, g[keep])
    n_groups <- max(1L, length(group_levels))
  }
  si <- match(s, sample_levels)
  gi <- if (is.null(g)) rep(1L, length(cells)) else match(g, group_levels)
  combo <- (si - 1L) * n_groups + gi
  combo[!keep] <- NA
  observed <- sort(unique(combo[keep]))          # sample-major, then group
  codes <- match(combo, observed)
  codes[is.na(codes)] <- 0L
  codes <- as.integer(codes)
  p_sample <- sample_levels[(observed - 1L) %/% n_groups + 1L]
  p_group <- if (is.null(g)) rep(NA_character_, length(observed)) else group_levels[(observed - 1L) %% n_groups + 1L]
  profile_id <- if (is.null(g)) make.unique(p_sample) else make.unique(paste(p_sample, p_group, sep = "__"))
  profiles <- data.frame(
    profile_id = profile_id,
    sample = p_sample,
    group = p_group,
    n_cells = tabulate(codes, length(observed)),
    row.names = profile_id,
    stringsAsFactors = FALSE
  )
  list(
    cells = cells,
    codes = codes,
    profiles = profiles,
    n_cells_dropped = c(missing_sample = sum(drop_s), missing_group = sum(drop_g & !drop_s))
  )
}

.pseudobulk_result <- function(total, n, design, features, assay, layer, fun, sample_by, group_by, engine) {
  profiles <- design$profiles
  n <- as.integer(n)
  if (length(n) != nrow(profiles)) stop("the backend returned ", length(n), " cell counts for ", nrow(profiles), " profiles", call. = FALSE)
  profiles$n_cells <- n
  total <- as.matrix(total)
  storage.mode(total) <- "double"
  if (fun == "mean") total <- sweep(total, 2, pmax(n, 1L), "/")
  dimnames(total) <- list(features, profiles$profile_id)
  structure(list(
    counts = total,
    profiles = profiles,
    assay = assay,
    layer = layer,
    fun = fun,
    sample_by = .pseudobulk_label_name(sample_by),
    group_by = if (is.null(group_by)) NA_character_ else .pseudobulk_label_name(group_by),
    n_cells_dropped = design$n_cells_dropped,
    engine = engine
  ), class = "scuba_pseudobulk")
}

# The layer(s) of a Seurat assay holding `layer`: the layer itself, or the split layers
# `layer.*` of a v5 assay (LayerData(layer = "counts") would silently read only the first).
.seurat_counts_layers <- function(assay_object, layer) {
  if (methods::is(assay_object, "Assay5")) {
    layers <- SeuratObject::Layers(assay_object)
    if (layer %in% layers) return(layer)
    split <- SeuratObject::Layers(assay_object, search = layer)
    split <- split[startsWith(split, paste0(layer, "."))]
    if (length(split)) return(split)
  } else {
    layers <- c("counts", "data", "scale.data")
    present <- vapply(layers, function(l) {
      m <- tryCatch(SeuratObject::GetAssayData(assay_object, layer = l), error = function(e) NULL)
      !is.null(m) && length(m) > 0 && nrow(m) > 0
    }, logical(1))
    layers <- layers[present]
    if (layer %in% layers) return(layer)
  }
  stop("Layer '", layer, "' not found in the assay; available: ", paste(layers, collapse = ", "),
       ". Pass the layer holding raw counts.", call. = FALSE)
}

# --- R engine --------------------------------------------------------------------------------

# Aggregate a features x cells layer `m` by the profile codes of `design`. The layer is never
# subset by column (no copies, and BPCells rejects column-subset matrices): cells outside the
# request get code 0. `partial`: the layer holds only some of the design's cells (split layers).
.r_pseudobulk_matrix <- function(m, design, features, partial = FALSE) {
  cn <- colnames(m); rn <- rownames(m)
  if (is.null(cn) || is.null(rn)) stop("the layer has no row or column names", call. = FALSE)
  n_profiles <- nrow(design$profiles)
  pos <- match(design$cells, cn)
  if (anyNA(pos)) {
    if (!partial) stop("Cells not found in the layer: ", paste(utils::head(design$cells[is.na(pos)], 10), collapse = ", "), call. = FALSE)
    keep <- !is.na(pos); pos <- pos[keep]; design_codes <- design$codes[keep]
  } else {
    design_codes <- design$codes
  }
  codes <- integer(ncol(m))
  codes[pos] <- design_codes
  fidx <- match(features, rn)
  if (anyNA(fidx)) stop("Features not found: ", paste(utils::head(features[is.na(fidx)], 10), collapse = ", "), call. = FALSE)
  total <- if (inherits(m, "IterableMatrix")) {
    .r_pseudobulk_bpcells(m, codes, n_profiles)
  } else if (methods::is(m, "DelayedArray") || methods::is(m, "SparseArray")) {
    .r_pseudobulk_colsum(m, codes, n_profiles)
  } else if (inherits(m, "dgCMatrix") || is.matrix(m)) {
    .r_pseudobulk_indicator(m, codes, n_profiles)
  } else if (inherits(m, "sparseMatrix")) {
    .r_pseudobulk_indicator(methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix"), codes, n_profiles)
  } else {
    .r_pseudobulk_blocks(m, codes, n_profiles)
  }
  total <- as.matrix(total)[fidx, , drop = FALSE]
  list(sum = total, n = tabulate(codes, n_profiles))
}

# cells x profiles sparse indicator; `m %*% ind` sums the columns of each profile in one pass.
.r_pseudobulk_indicator_matrix <- function(codes, n_profiles) {
  keep <- which(codes > 0L)
  Matrix::sparseMatrix(i = keep, j = codes[keep], x = 1, dims = c(length(codes), n_profiles))
}
.r_pseudobulk_indicator <- function(m, codes, n_profiles) {
  as.matrix(m %*% .r_pseudobulk_indicator_matrix(codes, n_profiles))
}

# BPCells: its own streaming pseudobulk. The factor needs every level used (unused levels error)
# and no NA (NA folds into the first level), so dropped cells get an explicit sentinel level.
.r_pseudobulk_bpcells <- function(m, codes, n_profiles) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("The layer is a BPCells matrix; install BPCells (https://bnprks.github.io/BPCells/) to aggregate it.", call. = FALSE)
  }
  level_names <- as.character(seq_len(n_profiles))
  dropped <- any(codes == 0L)
  groups <- factor(ifelse(codes > 0L, as.character(codes), ".scuba_dropped"),
                   levels = c(level_names, if (dropped) ".scuba_dropped"))
  res <- tryCatch(
    BPCells::pseudobulk_matrix(m, groups, method = "sum", threads = as.integer(getOption("SCUBA.threads", 1L))),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    # transformed / queued BPCells matrices: fall back to the sparse product (also streamed)
    first <- conditionMessage(res)
    return(tryCatch(
      .r_pseudobulk_indicator(m, codes, n_profiles),
      error = function(e) stop("BPCells::pseudobulk_matrix() failed (", first, ") and so did the sparse product fallback (",
                               conditionMessage(e), ")", call. = FALSE)
    ))
  }
  res <- as.matrix(res)
  if (!is.null(colnames(res)) && all(level_names %in% colnames(res))) res[, level_names, drop = FALSE]
  else res[, seq_len(n_profiles), drop = FALSE]
}

# DelayedArray (HDF5Array) and SparseArray layers: block-wise column sums by group.
.r_pseudobulk_colsum <- function(m, codes, n_profiles) {
  if (!requireNamespace("DelayedArray", quietly = TRUE)) {
    stop("The layer is a DelayedArray; install DelayedArray to aggregate it.", call. = FALSE)
  }
  res <- as.matrix(DelayedArray::colsum(m, group = codes, reorder = TRUE))
  wanted <- as.character(seq_len(n_profiles))
  if (!is.null(colnames(res)) && all(wanted %in% colnames(res))) return(res[, wanted, drop = FALSE])
  # no dimnames: columns follow sort(unique(codes)); drop the leading 0 column when present
  present <- sort(unique(codes))
  out <- matrix(0, nrow(res), n_profiles)
  hit <- present > 0L
  out[, present[hit]] <- res[, which(hit), drop = FALSE]
  out
}

# Any other matrix class: bounded blocks of cells, densified one block at a time.
.r_pseudobulk_blocks <- function(m, codes, n_profiles, block = 5000L) {
  total <- matrix(0, nrow(m), n_profiles)
  for (idx in split(seq_len(ncol(m)), ceiling(seq_len(ncol(m)) / block))) {
    k <- idx[codes[idx] > 0L]
    if (!length(k)) next
    ind <- Matrix::sparseMatrix(i = seq_along(k), j = codes[k], x = 1, dims = c(length(k), n_profiles))
    total <- total + as.matrix(as.matrix(m[, k, drop = FALSE]) %*% ind)
  }
  total
}
