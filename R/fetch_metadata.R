#' Fetch metadata
#'
#' Returns object metadata for a specified set of cells.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param vars metadata variables to pull from object
#' @param cells cells to include in the returned metadata
#' @param ... Currently unused.
#'
#' @rdname fetch_metadata
#'
#' @export
fetch_metadata <-
  function(
    object,
    vars,
    cells,
    ...
  ){
    UseMethod("fetch_metadata")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
fetch_metadata.default <-
  function(
    object,
    vars,
    cells
  ){
    warning(
      paste0(
        "fetch_metadata does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn fetch_metadata Seurat objects
#' @export
fetch_metadata.Seurat <-
  function(
    object,
    vars,
    cells
  ){
    # Seurat objects: pull metadata with `[[`
    data <-
      as.data.frame(
        object[[vars]][cells, ]
        )

    colnames(data) <- vars

    data
  }

#' @describeIn fetch_metadata SingleCellExperiment objects
#' @export
fetch_metadata.SingleCellExperiment <-
  function(
    object,
    vars,
    cells
  ){
    # SingleCellExperiment objects: use colData
    data <-
      as.data.frame(
        colData(object)[cells, vars]
        )

    colnames(data) <- vars

    data
  }
