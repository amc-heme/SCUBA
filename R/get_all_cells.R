#' Plotting utility function to return all cell IDs
#'
#' Returns a character vector with all cell names in the object. This is a utility function used to set defaults for plotting functions created by SCUBA
#' 
#' For additional information, see our GitHub.io website ("User Guide" article, "Data Visualization" section)
#'
#' @inheritParams object_param
#'
#' @rdname get_all_cells
#'
#' @export
#' 
#' @examples
#' get_all_cells(AML_Seurat) |> str()
get_all_cells <-
  function(
    object
  ){
    UseMethod("get_all_cells")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
get_all_cells.default <-
  function(
    object
  ){
    warning(
      paste0(
        "get_all_cells does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
        )
      )
    }

#' @describeIn get_all_cells Seurat objects
#' @export
get_all_cells.Seurat <-
  function(
    object
  ){
    # Seurat objects: cell names are column names of object
    colnames(object)
  }

#' @describeIn get_all_cells SingleCellExperiment objects
#' @export
get_all_cells.SingleCellExperiment <-
  function(
    object
  ){
    # SingleCellExperiment objects: cell names are column names of object
    colnames(object)
  }

#' @describeIn get_all_cells SingleCellExperiment objects
#' @export
get_all_cells.AnnDataR6 <-
  function(
    object
  ){
    # Use obs_names method from anndata
    object$obs_names
  }

#' @describeIn get_all_cells BPCells objects
#' @export
get_all_cells.BPCells <-
  function(
    object
  ){
    if (!is.null(object$cells)){
      return(object$cells)
    }
    # Attempt to infer from first reduction matrix
    if (!is.null(object$reductions) && length(object$reductions) >= 1){
      first_red <- object$reductions[[1]]
      if (!is.null(rownames(first_red))){
        return(rownames(first_red))
      }
    }
    stop("Unable to determine cell IDs for BPCells object (no 'cells' vector and reduction matrices lack rownames).")
  }

