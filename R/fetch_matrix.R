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
#' @param matrix_location (AnnData objects) the legacy slot to return, for
#'   example `"X"` or `c("obsm", "X_umap")`; defaults to `features`.
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
    layer <- layer %||% default_layer(object)
    # Split layers (counts.1, counts.2, ... after split() or merge()) hold disjoint cells and
    # possibly different features: they are assembled on the R path on either engine (SCUBA.jl
    # indexes cells against the whole object, not the layer). Other requests on the Julia engine
    # go to SCUBA.jl, which resolves keys and modality names itself.
    if (.scuba_julia_object(object) && !.seurat_split_request(object, features, assay, layer)) {
      return(.jl_fetch_matrix(object, features, assay = assay, layer = layer, cells = cells, sparse = sparse))
    }
    request <- .r_matrix_request(object, features, assay)
    assay_object <- object[[request$assay]]
    pieces <- .seurat_layer_pieces(assay_object, layer)
    if (length(pieces) > 1) {
      return(.r_split_layer_slice(assay_object, pieces, request, cells, sparse))
    }
    m <- SeuratObject::LayerData(assay_object, layer = pieces)
    .r_matrix_slice(m, request$features, cells, sparse, row_names = request$requested)
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
    request <- .r_matrix_request(object, features, assay)
    layer <- layer %||% default_layer(object)
    exp <- if (identical(request$assay, SingleCellExperiment::mainExpName(object)) ||
               is.null(SingleCellExperiment::mainExpName(object)) && identical(request$assay, "main")) {
      object
    } else {
      SingleCellExperiment::altExp(object, request$assay)
    }
    m <- SummarizedExperiment::assay(exp, layer)
    .r_matrix_slice(m, request$features, cells, sparse, row_names = request$requested)
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
    .jl_fetch_matrix(object$path, features, assay = assay, layer = layer, cells = .src_cells_for_julia(object, cells),
                     sparse = sparse, cellnames = .src_cells(object, cells))
  }

# The layer itself, or its split pieces ("counts" -> counts.1, counts.2, ...) in a Seurat v5
# assay. Seurat v3 assays have no split layers.
.seurat_layer_pieces <- function(assay_object, layer) {
  if (!methods::is(assay_object, "Assay5")) return(layer)
  layers <- SeuratObject::Layers(assay_object)
  if (layer %in% layers) return(layer)
  split <- layers[startsWith(layers, paste0(layer, "."))]
  if (length(split)) return(split)
  stop("Layer '", layer, "' not found in the assay; available: ", paste(layers, collapse = ", "), call. = FALSE)
}

# Whether a Seurat fetch_matrix() request reads a split layer. Lenient on purpose: an assay name
# R does not recognise is left to SCUBA.jl (which also accepts modality keys) and reported there.
.seurat_split_request <- function(object, features, assay, layer) {
  assays <- SeuratObject::Assays(object)
  target <- if (!is.null(assay)) {
    if (!assay %in% assays) return(FALSE)
    assay
  } else {
    keys <- vapply(assays, function(a) SeuratObject::Key(object[[a]]), character(1))
    keyed <- assays[vapply(keys, function(k) any(startsWith(as.character(features), k)), logical(1))]
    if (length(keyed)) keyed[1] else SeuratObject::DefaultAssay(object)
  }
  assay_object <- object[[target]]
  methods::is(assay_object, "Assay5") && !layer %in% SeuratObject::Layers(assay_object) &&
    any(startsWith(SeuratObject::Layers(assay_object), paste0(layer, ".")))
}

# Which assay (Seurat) or experiment (SingleCellExperiment) a fetch_matrix() request reads, and
# the bare feature names to slice with. Keyed names ("rna_GAPDH", "AB_CD117-AB") select their
# assay and lose the key; bare names use `assay`, or the default assay. Features must come from
# one assay, as on the Julia engine. `requested` keeps the names as given, for the row names.
.r_matrix_request <- function(object, features, assay = NULL) {
  if (inherits(object, "Seurat")) {
    assays <- SeuratObject::Assays(object)
    keys <- vapply(assays, function(a) SeuratObject::Key(object[[a]]), character(1))
    default <- SeuratObject::DefaultAssay(object)
  } else {
    main <- SingleCellExperiment::mainExpName(object) %||% "main"
    assays <- c(main, SingleCellExperiment::altExpNames(object))
    keys <- paste0(assays, "_")
    default <- main
  }
  names(keys) <- assays
  if (!is.null(assay) && !assay %in% assays) {
    stop("Assay '", assay, "' not found; available: ", paste(assays, collapse = ", "), call. = FALSE)
  }
  if (is.null(features)) return(list(assay = assay %||% default, features = NULL, requested = NULL))
  features <- as.character(features)
  feature_assay <- character(length(features))
  bare <- features
  for (i in seq_along(features)) {
    candidates <- if (is.null(assay)) assays else assay
    hit <- candidates[startsWith(features[i], keys[candidates])]
    if (length(hit)) {
      feature_assay[i] <- hit[1]
      bare[i] <- substring(features[i], nchar(keys[hit[1]]) + 1L)
    } else {
      feature_assay[i] <- assay %||% default
    }
  }
  used <- unique(feature_assay)
  if (length(used) > 1) {
    stop("features span several assays (", paste(used, collapse = ", "),
         "); fetch_matrix takes features from one assay", call. = FALSE)
  }
  list(assay = used, features = bare, requested = features)
}

# Slice a features x cells layer in R (any matrix class that supports `[`), as dgCMatrix or dense.
# `row_names`: names for the returned rows (the features as requested, possibly keyed).
.r_matrix_slice <- function(m, features, cells, sparse, row_names = features) {
  f <- features %||% rownames(m)
  missing_f <- setdiff(f, rownames(m))
  if (length(missing_f)) stop("Features not found: ", paste(utils::head(missing_f, 10), collapse = ", "), call. = FALSE)
  c <- cells %||% colnames(m)
  missing_c <- setdiff(c, colnames(m))
  if (length(missing_c)) stop("Cells not found: ", paste(utils::head(missing_c, 10), collapse = ", "), call. = FALSE)
  out <- m[f, c, drop = FALSE]
  out <- if (sparse) .as_dgcmatrix(out) else as.matrix(out)
  dimnames(out) <- list(row_names %||% f, c)
  out
}

# Coerce a (sliced) matrix to dgCMatrix without going through a dense copy where the class allows:
# other Matrix sparse classes convert directly, and DelayedArray (HDF5-backed SCE) and BPCells
# matrices have their own block-wise coercions. Only plain or unknown matrices become dense first.
.as_dgcmatrix <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  if (methods::is(x, "sparseMatrix")) {
    return(methods::as(methods::as(methods::as(x, "CsparseMatrix"), "generalMatrix"), "dMatrix"))
  }
  if (inherits(x, "IterableMatrix") || methods::is(x, "DelayedArray")) {
    return(methods::as(x, "dgCMatrix"))
  }
  methods::as(methods::as(as.matrix(x), "CsparseMatrix"), "generalMatrix")
}

# Place the rows of a sparse matrix into a zero matrix with rows `all_rows` (by name), for layers
# that hold only some of an assay's features.
.zero_fill_rows <- function(m, all_rows) {
  if (identical(rownames(m), all_rows)) return(m)
  triplets <- methods::as(m, "TsparseMatrix")
  Matrix::sparseMatrix(
    i = match(rownames(m), all_rows)[triplets@i + 1L], j = triplets@j + 1L, x = triplets@x,
    dims = c(length(all_rows), ncol(m)), dimnames = list(all_rows, colnames(m))
  )
}

# fetch_matrix() over split Seurat v5 layers: read each piece for the requested features and the
# cells it holds, zero-fill features a piece lacks, bind the cells, and return them in request
# order.
.r_split_layer_slice <- function(assay_object, pieces, request, cells, sparse) {
  f <- request$features %||% rownames(assay_object)
  known <- unique(unlist(lapply(pieces, function(p) rownames(SeuratObject::LayerData(assay_object, layer = p)))))
  missing_f <- setdiff(f, known)
  if (length(missing_f)) stop("Features not found: ", paste(utils::head(missing_f, 10), collapse = ", "), call. = FALSE)
  blocks <- lapply(pieces, function(p) {
    m <- SeuratObject::LayerData(assay_object, layer = p)
    keep_cells <- if (is.null(cells)) colnames(m) else intersect(cells, colnames(m))
    if (!length(keep_cells)) return(NULL)
    present <- intersect(f, rownames(m))
    .zero_fill_rows(.as_dgcmatrix(m[present, keep_cells, drop = FALSE]), f)
  })
  out <- do.call(cbind, Filter(Negate(is.null), blocks))
  c <- cells %||% colnames(assay_object)
  missing_c <- setdiff(c, colnames(out))
  if (length(missing_c)) stop("Cells not found: ", paste(utils::head(missing_c, 10), collapse = ", "), call. = FALSE)
  out <- out[, c, drop = FALSE]
  if (!sparse) out <- as.matrix(out)
  dimnames(out) <- list(request$requested %||% f, c)
  out
}
