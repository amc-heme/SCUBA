# reticulate converters for Python anndata >= 0.13 and pandas >= 3 ####
#
# Stopgap until the anndataR migration (amc-heme/SCUBA#173): delete this file
# and the call to .scuba_register_py_converters() in .onLoad() (R/zzz.R).
#
# reticulate turns a Python object into an R object by dispatching py_to_r()
# on the Python class name ("<module>.<class>"). The superseded anndata R
# package registers its converters under the names older Python releases
# used. anndata 0.13 renamed `anndata._core.anndata.AnnData` to
# `anndata.AnnData`, and pandas 3 renamed `pandas.core.indexes.base.Index` to
# `pandas.Index`, so without these two methods anndata::read_h5ad() returns
# an unwrapped Python object (not AnnDataR6) whose obs_names stay a Python
# Index.

#' Wrap a Python AnnData object as AnnDataR6
#'
#' @param x A Python `anndata.AnnData` object (reticulate reference).
#'
#' @return An `AnnDataR6` object from the anndata R package.
#'
#' @noRd
.py_to_r_anndata <- function(x) {
  anndata::AnnDataR6$new(x)
}

#' Convert a pandas Index to an R vector
#'
#' Mirrors the anndata R package's method for the pre-pandas-3 class name:
#' the values as an R vector, with the Index name as the `name` attribute.
#'
#' @param x A Python `pandas.Index` object (reticulate reference).
#'
#' @return An R vector of the index values.
#'
#' @noRd
.py_to_r_pandas_index <- function(x) {
  out <- reticulate::import_builtins()$list(x)
  attr(out, "name") <- reticulate::py_to_r(x$name)
  out
}

#' Register SCUBA's py_to_r converters where none exists yet
#'
#' Each converter is registered only when reticulate has no method for that
#' class yet. So an anndata R release that fixes the names itself, or
#' anndataR (which registers its own AnnData converter), keeps its method when
#' it was loaded first; a package loaded after SCUBA replaces SCUBA's method,
#' as S3 registration always does. SCUBA never overrides another package.
#'
#' @param converters Named list of functions; names are Python class names.
#'
#' @return The class names registered, invisibly.
#'
#' @noRd
.scuba_register_py_converters <- function(
    converters = list(
      anndata.AnnData = .py_to_r_anndata,
      pandas.Index = .py_to_r_pandas_index
    )
) {
  reticulate_namespace <- asNamespace("reticulate")
  registered <- character()
  for (class_name in names(converters)) {
    existing <- utils::getS3method(
      "py_to_r", class_name,
      optional = TRUE, envir = reticulate_namespace
    )
    if (is.null(existing)) {
      registerS3method(
        "py_to_r", class_name, converters[[class_name]],
        envir = reticulate_namespace
      )
      registered <- c(registered, class_name)
    }
  }
  invisible(registered)
}
