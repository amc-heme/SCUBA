#' Object metadata: number of unique values
#'
#' Returns the number of unique values in a single cell object, for the specified
#' metadata variable. Any object type supported by
#' \code{\link[unique_values]{unique_values}} may be used.
#'
#' @param object A single cell object.
#' @param meta_var Name of a metadata variable.
#'
#' @return Number of unique classes (integer).
#'
#' @keywords internal
#'
#' @noRd
n_unique <-
  function(
    object,
    meta_var
  ){
    SCUBA::unique_values(
      object = object,
      var = meta_var
      ) |>
      length()
  }
