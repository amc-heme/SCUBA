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
    as.data.frame(
      reducedDims(object)[[reduction]][cells, dims]
    )
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

