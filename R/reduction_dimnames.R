#' Get names of reduction keys
#'
#' Given the name of a reduction and a set of dims, this function will return
#' the names of each dim as it appears in the reduction matrix. The output of
#' this function is passed to FetchData to pull information for reductions, and
#' it is also used to label reductions on DimPlots and feature plots.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param reduction the reduction from which names should be formed
#' @param dims a two-element integer vector with the dimensions for which
#' names should be returned
#' @param ... Currently unused.
#'
#' @rdname reduction_dimnames
#'
#' @export
reduction_dimnames <-
  function(
    object,
    reduction,
    dims,
    ...
  ){
    UseMethod("reduction_dimnames")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
reduction_dimnames.default <-
  function(
    object,
    reduction,
    dims
  ){
    warning(
      paste0(
        "reduction_dimnames does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn reduction_dimnames Seurat objects
#' @export
reduction_dimnames.Seurat <-
  function(
    object,
    reduction,
    dims
  ){
    if (length(x = dims) != 2) {
      stop("'dims' must be a two-length vector")
    }

    # Seurat objects: fetch the key for each reduction, and append the dim
    paste0(Key(object[[reduction]]), dims)
  }

#' @describeIn reduction_dimnames SingleCellExperiment objects
#' @export
reduction_dimnames.SingleCellExperiment <-
  function(
    object,
    reduction,
    dims
  ){
    if (length(x = dims) != 2) {
      stop("'dims' must be a two-length vector")
    }

    # SingleCellExperiment objects: reduction passed does not to be converted.
    # Simply append each dim after an underscore (<reduction>_<dim>)
    paste0(reduction, "_", dims)
  }


#' @describeIn reduction_dimnames AnnDataR6 objects
#' @export
reduction_dimnames.AnnDataR6 <-
  function(
    object,
    reduction,
    dims
  ){
    if (length(x = dims) != 2) {
      stop("'dims' must be a two-length vector")
    }
    
    #AnnData
    paste0(reduction, "_", dims)
  }
