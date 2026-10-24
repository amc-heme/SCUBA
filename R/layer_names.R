#' Layers of an assay
#'
#' Names of the layers (Seurat), assays (SingleCellExperiment experiment), matrices (AnnData:
#' `"X"`, every entry of `layers`, and `"raw"` when `adata.raw` is set) or layers (`.h5ad` /
#' `.h5mu` modality) available for an assay. Use it to find a raw-count layer before calling
#' [pseudobulk()].
#'
#' @inheritParams object_param
#' @param assay Assay / modality to inspect (default: the default assay).
#' @returns A character vector of layer names (empty when the assay has none).
#' @export
#'
#' @examples
#' layer_names(AML_Seurat)
#' layer_names(AML_SCE(), assay = "AB")
layer_names <- function(object, assay = NULL) UseMethod("layer_names")

#' @export
layer_names.default <- function(object, assay = NULL) {
  warning("layer_names does not know how to handle object of class ", paste(class(object), collapse = ", "),
          ". Currently supported classes: Seurat, SingleCellExperiment, AnnDataR6 and scuba_source.", call. = FALSE)
  character()
}

#' @rdname layer_names
#' @export
layer_names.Seurat <- function(object, assay = NULL) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) {
    stop("Assay '", assay, "' not found; assays: ", paste(SeuratObject::Assays(object), collapse = ", "), call. = FALSE)
  }
  assay_object <- object[[assay]]
  if (methods::is(assay_object, "Assay5")) return(SeuratObject::Layers(assay_object))
  layers <- c("counts", "data", "scale.data")
  layers[vapply(layers, function(l) {
    m <- tryCatch(SeuratObject::GetAssayData(assay_object, layer = l), error = function(e) NULL)
    !is.null(m) && length(m) > 0 && nrow(m) > 0
  }, logical(1))]
}

#' @rdname layer_names
#' @export
layer_names.SingleCellExperiment <- function(object, assay = NULL) {
  main <- SingleCellExperiment::mainExpName(object)
  if (is.null(assay) || identical(assay, main) || (is.null(main) && identical(assay, "main"))) {
    return(SummarizedExperiment::assayNames(object))
  }
  if (!assay %in% SingleCellExperiment::altExpNames(object)) {
    stop("Assay '", assay, "' not found; experiments: ",
         paste(c(main %||% "main", SingleCellExperiment::altExpNames(object)), collapse = ", "), call. = FALSE)
  }
  SummarizedExperiment::assayNames(SingleCellExperiment::altExp(object, assay))
}

#' @rdname layer_names
#' @export
layer_names.AnnDataR6 <- function(object, assay = NULL) {
  raw <- object$raw
  has_raw <- !is.null(raw) && !inherits(raw, "python.builtin.NoneType")
  if (is.null(assay) || identical(assay, "X")) {
    layers <- c("X", .anndata_keys(object$layers))
    if (has_raw) layers <- c(layers, "raw")
    return(layers)
  }
  if (identical(assay, "raw")) {
    if (!has_raw) stop("The object has no raw matrix", call. = FALSE)
    return("X")
  }
  obsm <- .anndata_keys(object$obsm)
  # only obsm data frames (e.g. protein panels) are modalities; numeric arrays are reductions
  obsm <- obsm[vapply(obsm, function(k) {
    x <- object$obsm[[k]]
    is.data.frame(x) || any(grepl("DataFrame", class(x), fixed = TRUE))
  }, logical(1))]
  if (assay %in% obsm) return("X")
  stop("Modality '", assay, "' not found; modalities: ", paste(c("X", obsm), collapse = ", "), call. = FALSE)
}

#' @rdname layer_names
#' @exportS3Method layer_names anndata.AnnData
layer_names.anndata.AnnData <- layer_names.AnnDataR6
#' @rdname layer_names
#' @exportS3Method layer_names anndata._core.anndata.AnnData
layer_names.anndata._core.anndata.AnnData <- layer_names.AnnDataR6

#' @rdname layer_names
#' @export
layer_names.scuba_source <- function(object, assay = NULL) {
  assay <- assay %||% object$info$modalities[1]
  if (!assay %in% object$info$modalities) {
    stop("Modality '", assay, "' not found; modalities: ", paste(object$info$modalities, collapse = ", "), call. = FALSE)
  }
  as.character(object$info$layers[[assay]])
}
