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

