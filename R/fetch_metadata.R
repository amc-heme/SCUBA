#' Fetch metadata
#'
#' Returns object metadata for a specified set of cells.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param vars metadata variables to pull from object.
#' @param cells cells to include in the returned metadata.
#' @param return_class class of data returned. Set to "dataframe" by default to
#' return a data.frame, and may also be set to "vector" to yield a vector of
#' values.
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
    return_class = "dataframe",
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
    cells,
    return_class = "dataframe"
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
    cells,
    return_class = "dataframe"
  ){
    # Seurat objects: pull metadata with `[[`
    data <- object[[vars]][cells, ]

    # Use as.data.frame() if return_class is equal to "dataframe"
    # (if "vector", no additional operations neccessary)
    if (return_class == "dataframe"){
      data <-
        as.data.frame(
          data
        )

      colnames(data) <- vars
    }

    data
  }

#' @describeIn fetch_metadata SingleCellExperiment objects
#' @export
fetch_metadata.SingleCellExperiment <-
  function(
    object,
    vars,
    cells,
    return_class = "dataframe"
  ){
    # SingleCellExperiment objects: use colData
    data <-
      colData(object)[cells, vars]

    # Use as.data.frame() if return_class is equal to "dataframe"
    # (if "vector", no additional operations neccessary)
    if (return_class == "dataframe"){
      data <-
        as.data.frame(
          data
        )

      colnames(data) <- vars
    }

    data
  }
