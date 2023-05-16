#' Subset cells based on metadata
#'
#' Returns a character vector with the cell names matching the levels/classes
#' of a specified metadata variable.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param meta_var a metadata variable used as the basis for subsetting cells.
#' @param meta_levels the levels of the specified metadata variable in
#' \code{meta_var} to include in the
#' @param ... Currently unused.
#'
#' @rdname fetch_cells
#'
#' @export
fetch_cells <-
  function(
    object,
    meta_var,
    meta_levels,
    ...
  ){
    UseMethod("fetch_cells")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
fetch_cells.default <-
  function(
    object,
    meta_var,
    meta_levels
  ){
    warning(
      paste0(
        "fetch_cells does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn fetch_cells Seurat objects
#' @export
fetch_cells.Seurat <-
  function(
    object,
    meta_var,
    meta_levels
  ){
    # Seurat objects: subset meta.data object, extract rownames
    object@meta.data[object@meta.data[[meta_var]] %in% meta_levels, ] |>
      rownames()
  }

#' @describeIn fetch_cells SingleCellExperiment objects
#' @export
fetch_cells.SingleCellExperiment <-
  function(
    object,
    meta_var,
    meta_levels
  ){
    # SingleCellExperiment objects: subset colData, extract rownames
    colData(object)[colData(object)[[meta_var]] %in% meta_levels, ] |>
      rownames()
  }
