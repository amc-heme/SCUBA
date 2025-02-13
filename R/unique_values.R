#' Summarize values in a metadata variable
#'
#' Returns the unique values for the metadata variable provided. `unique_values`
#' is a utility function for summarizing metadata in an object.
#'
#' @inheritParams object_param
#' @param var the metadata variable for which to return unique values.
#' 
#' @returns a character vector giving the unique values in the specified
#' metadata variable.
#'
#' @export
#' 
#' @examples
#' unique_values(AML_Seurat, var = "Batch")
#' 
#' unique_values(AML_Seurat, var = "condensed_cell_type")
unique_values <-
  function(
    object,
    var
  ){
    if (length(var) > 1){
      stop('Only one variable may be passed to "var"')
    }

    if (!class(var) %in% c("character", "factor")){
      stop('"var" is not a character vector. Please pass the name of one metadata variable to "var".')
    }

    # Handle errors caused by passing a var not found in the object
    # smetadata
    if (!var %in% SCUBA::meta_varnames(object)){
      stop("unique_values: var ", var, " not found in object.")
    }

    fetch_metadata(
      object = object,
      vars = var,
      cells = get_all_cells(object),
      return_class = "vector"
      ) |>
      unique()
  }
