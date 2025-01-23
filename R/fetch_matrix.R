#' Fetch matrix
#'
#' Returns object metadata for a specified set of cells.
#'
#' @param object a single-cell object (Seurat, SingleCellExperiment, Anndata).
#' @param matrix_location the
#' @param ... Currently unused.
#'
#' @keywords internal
#'
#' @rdname fetch_matrix
#'
#' @export
fetch_matrix <-
  function(
    object,
    matrix_location,
    ...
  ){
    UseMethod("fetch_matrix")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
fetch_matrix.default <-
  function(
    object,
    matrix_location
  ){
    warning(
      paste0(
        "fetch_matrix does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat, SingleCellExperiment, ",
        "and anndata."
      )
    )
  }

# # describeIn fetch_matrix Seurat objects
# # export
# fetch_matrix.Seurat <-
#   function(
#     object,
#     matrix_location
#   ){
#
#   }
#
# # describeIn fetch_matrix SingleCellExperiment objects
# # export
# fetch_matrix.SingleCellExperiment <-
#   function(
#     object,
#     matrix_location
#     ){
#
#     }



#' @describeIn fetch_matrix AnnDataR6 objects
#' @export
fetch_matrix.AnnDataR6 <-
  function(
    object,
    matrix_location,
    densify = FALSE
  ){
    if (length(matrix_location) == 1){
      return(object[[matrix_location]])
      } else if (length(matrix_location) == 2){
      return(object[[matrix_location[1]]][[matrix_location[2]]])
      } else {
      stop('Unexpected format for matrix_location. This should be a ',
           'character vector of length one or two, giving the location of ',
           'the desired matrix in the object. For example, to pull from ',
           'the "obs" slot, use "obs". To pull from the "X_umap" matrix ',
           'in "obsm", use c("obsm", "X_umap").')
      }
    }
