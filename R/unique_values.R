#' Unique metadata values
#'
#' Returns the unique values for a given metadata variable in a single cell object.
#' The function supports all object types that can be passed \code{\link[fetch_metadata]{fetch_metadata}}.
#'
#' @param object any object type supported by
#' \code{\link[fetch_metadata]{fetch_metadata}}.
#' @param var the metadata variable for which to return unique values.
#' @return a character vector giving the unique values in the specified
#' metadata variable.
#'
#' @export
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
