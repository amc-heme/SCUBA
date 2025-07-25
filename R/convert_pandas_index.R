#' Remove Pandas index attribute and set row names
#'
#' For MuData objects, queries such as object$obs_names will return a Pandas
#' Index, which will not be properly converted to a character vector when
#' accessed in R. `convert_pandas_index` converts a Pandas index to a Python
#' list so it can be properly converted to an R character vector by Reticulate.
#'
#' @param data a Python object returned via Reticulate operations that may be 
#' a Pandas index. 
#'
#' @return If `data` is a Pandas Index, it is converted to a list in Python 
#' and returned as a character vector in R. Otherwise, it is returned as-is.
#'
#' @keywords internal
#'
#' @export
convert_pandas_index <-
  function(
    data
    ){
    # Ensure data is returned as an R character vector 
    # instead of a Python index
    if (inherits(data, "pandas.core.indexes.base.Index")){
      # If a Pandas Index, convert to a list in Python
      data$to_list()
    } else {
      data
    }
  }