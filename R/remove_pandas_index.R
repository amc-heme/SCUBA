#' Remove Pandas index attribute and set row names
#'
#' When Pandas Dataframes are converted to an R data.frame, they may still 
#' contain a "pandas.index" attribute. In SCUBA, this happens for data extracted
#' via the fetch_* functions, which construct and return a Pandas DataFrame. 
#' This function is used to fully convert data aggregated from Python object 
#' formats to R syntax.
#'
#' @param table A `data.frame` object possibly containing a `pandas.index` attribute. 
#'
#' @return Returns the modified `data.frame` with row names set
#' to the original Pandas index and the `pandas.index` attribute removed.
#'
#' @keywords internal
#'
#' @export
remove_pandas_index <-
  function(
    table
    ){
    if (!is.null(idx <- attr(table, "pandas.index"))){
      # idx is the Pandas index, if it exists
      rownames(table) <- 
        idx$to_list() |> 
        reticulate::py_to_r()
      # drop Pandas index after setting rownames
      attr(table, "pandas.index") <- NULL
    }
    
    table
  }