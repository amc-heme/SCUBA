#' Fetch reduction coordinates from single-cell objects
#'
#' Returns the reduction coordinates matching the name and dimensions supplied.
#'
#' See our GitHub.io website for additional information and examples.
#'
#' @inheritParams object_param
#' @param reduction the name of the reduction to pull coordinates from.
#' @param cells cell IDs for which to pull reduction data. If `NULL`, 
#' coordinates will be returned from all cells in the object. Cell IDs can be 
#' generated with [`fetch_cells()`].
#' @param dims a numeric vector indicating the dimensions to pull. Currently, 
#' only two dimensions are supported, but [`fetch_data()`] supports more than 
#' two dimensions. For instructions on pulling more than two dimensions at 
#' once, see the examples of [`fetch_data()`].
#'
#' @rdname fetch_reduction
#'
#' @export
#' 
#' @examples
#' # Return the first and second dimensions from the UMAP reduction
#' fetch_reduction(
#'   AML_Seurat,
#'   reduction = "umap",
#'   dims = c(1, 2)
#'   ) |> str()
fetch_reduction <-
  function(
    object,
    reduction,
    cells = NULL,
    dims = c(1, 2)
  ){
    UseMethod("fetch_reduction")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
fetch_reduction.default <-
  function(
    object,
    reduction,
    cells = NULL,
    dims = c(1, 2)
  ){
    warning(
      paste0(
        "fetch_reduction does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn fetch_reduction Seurat objects
#' @export
fetch_reduction.Seurat <-
  function(
    object,
    reduction,
    cells = NULL,
    dims = c(1, 2)
  ){
    if (length(x = dims) != 2) {
      stop("'dims' must be a two-length vector")
    }
    
    # Throw an error if the reduction entered does not exist
    if (!reduction %in% names(object@reductions)){
      # Collapse available reductions into a single-length character 
      # vector for display
      reductions_str <- paste(names(object@reductions), collapse = ", ")
      
      stop(
        paste0(
          '\nThe reduction "', reduction, 
          '" was not found in the object passsed. \n',
          'Reductions present in object: ', 
          reductions_str, 
          '.'
          )
        )
      }

    # Cells: if NULL, use all cells in the object
    cells <- cells %||% get_all_cells(object)

    # Fetch reduction coordinates using SeuratObject Embeddings method
    as.data.frame(
      Embeddings(object = object[[reduction]])[cells, dims]
      )
  }

#' @describeIn fetch_reduction SingleCellExperiment objects
#' @export
fetch_reduction.SingleCellExperiment <-
  function(
    object,
    reduction,
    cells = NULL,
    dims = c(1, 2)
  ){
    if (length(x = dims) != 2) {
      stop("'dims' must be a two-length vector")
    }
    
    # Throw an error if the reduction entered does not exist
    if (!reduction %in% reducedDimNames(object)){
      # Collapse available reductions into a single-length character 
      # vector for display
      reductions_str <- paste(reducedDimNames(object), collapse = ", ")
      
      stop(
        paste0(
          '\nThe reduction "', reduction, 
          '" was not found in the object passsed. \n',
          'Reductions present in object: ', 
          reductions_str, '.'
        )
      )
    }

    # Cells: if NULL, use all cells in the object
    cells <- cells %||% get_all_cells(object)

    # SingleCellExperiment objects: pull reduction matrix, then subset for
    # cells and dims
    data <- as.data.frame(
      reducedDims(object)[[reduction]][cells, dims]
    )
    
    # Special case: reduction does not have column names
    # Data will be returned, but it will have no column names. 
    # This will not cause an error with fetch_reduction, but behavior will 
    # be inconsistent with fetch_data
    if (is.null(colnames(reducedDims(object)[[reduction]]))){
      # Append the dims found to the reduction name with an underscore
      # to match input (i.e. if the reduction `key` is UMAP, this will
      # be UMAP_1, UMAP_2, etc.)
      colnames(data) <- paste(reduction, dims, sep = "_")
    }
    
    data
  }


#' @describeIn fetch_reduction AnnDataR6 objects
#' @export
fetch_reduction.AnnDataR6 <-
  function(
    object,
    reduction,
    cells = NULL,
    dims = c(1, 2)
  ){
    if (length(x = dims) != 2) {
      stop("'dims' must be a two-length vector")
    }

    # Throw an error if the reduction entered does not exist
    if (!reduction %in% object$obsm_keys()){
      # Collapse available reductions into a single-length character 
      # vector for display
      reductions_str <- paste(object$obsm_keys(), collapse = ", ")
      
      stop(
        paste0(
          '\nThe reduction "', reduction, 
          '" was not found in the obsm slot of object passsed. \n',
          'Available obsm matrices: ', 
          reductions_str, '.',
          '\n\n(This list includes all reductions, but is not limited to \n',
          'them since anndata objects do not have a location specific to \n',
          'reductions only. obsm matrices that are not reductions can be \n',
          'entered, but unexpected results may occur.)'
        )
      )
    }
    
    # Cells: if NULL, use all cells in the object
    cells <- cells %||% get_all_cells(object)

    fetch_data(
        object,
        vars = reduction_dimnames(object, reduction, dims),
        cells = cells
        )
  }

#' @describeIn fetch_reduction MuData objects
#' @export
fetch_reduction.md._core.mudata.MuData <-
  function(
    object,
    reduction,
    cells = NULL,
    dims = c(1, 2)
    ){
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop(
        paste0(
          'Package "reticulate" must be installed to use this ',
          'function with anndata objects.'
        ),
        call. = FALSE
      )
    }
    
    library(reticulate)
    
    # Establish Python package dependencies
    # Reticulate will automatically manage a Python environment with these 
    # packages, installing each if they are not already present
    try({
      # py_require statements called more than once per session 
      # may result in an error
      py_require("anndata>=0.11.4")
      py_require("pandas>=2.0.0")
      py_require("numpy")
      py_require("scipy>=1.14.0")
      py_require("mudata>=0.3.1")
    })
    
    # Source fetch_reduction python script for mudata
    python_path =
      system.file(
        "extdata",
        "Python",
        "fetch_data.py",
        package = "SCUBA"
        )
    
    py_objs <- reticulate::py_run_file(python_path)
    
    # Cells: if NULL, use all cells in the object
    cells <- cells %||% get_all_cells(object)
    
    # Convert dims to a character vector if it isn't already
    dims <- as.character(dims)
    
    # Pull reduction info via python function
    data <- py_objs$fetch_reduction_mudata(
      obj = object, 
      reduction = reduction, 
      cells = cells, 
      dims = dims
      )
    
    # Remove Pandas index artifact from return data.frame
    if (!is.null(idx <- attr(data, "pandas.index"))){
      # idx is the Pandas index, if it exists
      rownames(data) <- 
        idx$to_list() |> 
        reticulate::py_to_r()
      # drop Pandas index after setting rownames
      attr(data, "pandas.index") <- NULL
    }
    
    data
  }

#' @export
fetch_reduction.mudata._core.mudata.MuData <-
  function(
    object,
    reduction,
    cells = NULL,
    dims = c(1, 2)
    ){
    # mudata._core.mudata.MuData: possible class when loading 
    # Redirect to fetch_data.md._core.mudata.MuData method
    fetch_reduction.md._core.mudata.MuData(
      object = object,
      reduction = reduction,
      cells = cells,
      dims = dims
      )
    }