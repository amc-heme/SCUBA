#' Match a string against a set of keys and extract the suffix
#'
#' This function mirrors the behavior of the Python `_match_key` helper in 
#' inst/extdata/Python/fetch_data.py. It checks whether
#' a given string `var` matches exactly one of the provided `keys`, or begins with a key
#' followed by an underscore, returning the matched key and the trailing content.
#'
#' @param var A single character string to test.
#' @param keys A character vector of possible keys.
#' @param sep The separator between each key and the following content. This is 
#' usually "_", but for MuData objects, this may be ":", as the colon is the 
#' convention with MuData objects.
#' 
#' @return A list with two elements:
#'   \describe{
#'     \item{key}{The matched key, or NULL if no match was found.}
#'     \item{suffix}{The substring after the first underscore following the key, or NULL if none.}
#'   }
#' @examples
#' match_key("MOFA_UMAP_feature", c("MOFA_UMAP", "MOFA", "XYZ"))
#' #> list(key = "MOFA_UMAP", suffix = "feature")
#'
#' @keywords internal
#' @export
match_key <- function(var, keys, sep = "_"){
  # Ensure keys is a character vector (single string remains length 1)
  if (!is.character(keys)) {
    stop("`keys` must be a character vector")
  }
  
  # Sort keys by decreasing length to avoid premature matches
  keys_sorted <- keys[order(nchar(keys), decreasing = TRUE)]
  
  for (k in keys_sorted) {
    # Exact match
    if (identical(var, k)) {
      return(list(key = k, suffix = NULL))
    }
    # Define prefix to match var against
    # Prefix is constructed with a separator (usually an underscore)
    prefix <- paste0(k, sep)
    if (startsWith(var, prefix)) {
      # Extract content after key and underscore
      suffix <- substring(var, nchar(prefix) + 1)
      return(list(key = k, suffix = suffix))
    }
  }
  
  # No match found
  return(list(key = NULL, suffix = NULL))
}