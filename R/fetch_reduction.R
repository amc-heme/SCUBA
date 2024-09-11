#' Get names of reduction keys
#'
#' Given the name of a reduction and a set of dims, this function will return
#' the corresponding reduction data.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param reduction the reduction to pull coordinates from
#' @param cells cells for which to pull reduction data
#' @param dims a two-element integer vector with the dimensions for which
#' data should be returned
#' @param ... Currently unused.
#'
#' @rdname fetch_reduction
#'
#' @export
fetch_reduction <-
  function(
    object,
    reduction,
    cells = NULL,
    dims = c(1, 2),
    ...
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

    FetchData(
        object,
        vars = reduction_dimnames(object, reduction, dims),
        cells = cells
        )
  }

