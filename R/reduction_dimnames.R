#' From names of reduction keys for fetch_data
#'
#' Given the name of a reduction and a set of dimensions, this function will 
#' return the names of each dimension as it appears in the reduction matrix.
#' The output of this function can be passed to the `vars` parameter of 
#' [`fetch_data()`] to facilitate the specification of reduction coordinates to 
#' return from this function.
#'
#' @inheritParams object_param
#' @param reduction the name of the reduction.
#' @param dims a numeric vector with the dimensions for which
#' names should be returned. 
#'
#' @rdname reduction_dimnames
#'
#' @export
#' 
#' @examples
#' # From names for first and second UMAP dimensions and 
#' # pass to fetch_data
#' dimnames <- reduction_dimnames(
#'   AML_Seurat,
#'   reduction = "umap",
#'   dims = c(1, 2)
#'   )
#'   
#' dimnames
#' 
#' fetch_data(
#'   AML_Seurat,
#'   vars = dimnames
#'   ) |> str()
#'   
#'   
#' # Form names for first five PCA dimensions and 
#' # pass to fetch_data
#' dimnames <- reduction_dimnames(
#'   AML_Seurat,
#'   reduction = "pca",
#'   dims = c(1:5)
#'   )
#'   
#' dimnames
#' 
#' fetch_data(
#'   AML_Seurat,
#'   vars = dimnames
#'   ) |> str()
reduction_dimnames <-
  function(
    object,
    reduction,
    dims
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
    #AnnData: same as for SingleCellExperiment objects
    paste0(reduction, "_", dims)
  }

#' @describeIn reduction_dimnames MuData objects
#' @export
reduction_dimnames.md._core.mudata.MuData <-
  function(
    object,
    reduction,
    dims
  ){
    # Simply the reduction and dimensions with an underscore
    paste0(reduction, "_", dims)
  }


#' @describeIn reduction_dimnames MuData objects
#' @export
reduction_dimnames.mudata._core.mudata.MuData <-
  function(
    object,
    reduction,
    dims
  ){
    # mudata._core.mudata.MuData: possible class when loading 
    # Redirect to md._core.mudata.MuData method
    reduction_dimnames.md._core.mudata.MuData(
      object = object,
      reduction = reduction,
      dims = dims
    )
  }
