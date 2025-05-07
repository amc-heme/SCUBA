#' Return cell IDs for a subset based on metadata
#'
#' Returns a character vector with the cell names matching the levels/classes
#' of a specified metadata variable.
#'
#' @inheritParams object_param
#' @param meta_var a metadata variable used as the basis for subsetting cells.
#' @param meta_levels the levels of the specified metadata variable in
#' `meta_var` to include in the
#'
#' @rdname fetch_cells
#'
#' @export
fetch_cells <-
  function(
    object,
    meta_var,
    meta_levels
  ){
    # Return an error if meta_levels is equal to NULL
    if (is.null(meta_levels)){
      stop("Meta_levels must not be NULL.")
    }

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
        ". Currently supported classes: Seurat, SingleCellExperiment, ",
        "and Anndata."
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
    # Test if meta_var is a categorical variable. Numeric variables are not
    # supported
    if (
      !inherits(
        object@meta.data[[meta_var]],
        c("character", "factor", "logical")
        )
      ){
      stop(
        meta_var, " is not a categorical variable. Only categorical variables ",
        "are supported by fetch_cells at this time."
        )
    }

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
    # Test if meta_var is a categorical variable. Numeric variables are not
    # supported
    if (
      !inherits(
        colData(object)[[meta_var]],
        c("character", "factor", "logical")
      )
    ){
      stop(
        meta_var, " is not a categorical variable. Only categorical variables ",
        "are supported by fetch_cells at this time."
      )
    }

    # SingleCellExperiment objects: subset colData, extract rownames
    colData(object)[colData(object)[[meta_var]] %in% meta_levels, ] |>
      rownames()
  }

#' @describeIn fetch_cells Anndata objects
#' @export
fetch_cells.AnnDataR6 <-
  function(
    object,
    meta_var,
    meta_levels
  ){
    # Test if meta_var is a categorical variable. Numeric variables are not
    # supported
    if (
      !inherits(
        object$obs[[meta_var]],
        c("character", "factor", "logical")
      )
    ){
      stop(
        meta_var, " is not a categorical variable. Only categorical variables ",
        "are supported by fetch_cells at this time."
      )
    }

    # Anndata objects: use obs
    object$obs[object$obs[[meta_var]] %in% meta_levels, ] |>
      rownames()
  }
