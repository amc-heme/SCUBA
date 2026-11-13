# Errors for objects no SCUBA method handles ####

#' Python packages SCUBA declares for AnnData support
#'
#' The requirements `.onLoad()` passes to `reticulate::py_require()`. They live
#' here so the unsupported-class hint can quote exactly the same versions.
#'
#' - anndata >= 0.11.2: 0.10.x and 0.11.0-0.11.1 read backed sparse matrices
#'   through scipy's private `_validate_indices()`, whose return value changed
#'   in scipy 1.15 (backed reads then fail); 0.11.2 handles both forms.
#'   < 0.14: AnnData support rests on the superseded anndata R package plus
#'   the converters in R/anndata_shim.R, tested through the 0.13 series.
#' - pandas >= 2.1 (the oldest series CI tests, and the first with Python
#'   3.12 wheels); < 4: the next major release is untested. pandas 3 also
#'   needs reticulate >= 1.45.0, which first converts its renamed
#'   `pandas.DataFrame`; `.onLoad()` adds `pandas<3` for older reticulate.
#' - scipy >= 1.14: the oldest series CI tests; anndata bounds it otherwise.
#' - Python >= 3.10: required by anndata 0.11.
#'
#' Lift the anndata and pandas caps only after CI passes on the new series.
#'
#' @noRd
.scuba_python_requirements <- list(
  packages = c(
    "anndata>=0.11.2,<0.14", "pandas>=2.1.0,<4", "numpy", "scipy>=1.14.0"
  ),
  python_version = ">=3.10"
)

#' Classes a SCUBA generic has methods for
#'
#' Read from SCUBA's own S3 methods table, so the list in error messages always
#' matches what NAMESPACE registers. The `default` method is left out, and so
#' are the raw Python AnnData class names some generics register as aliases of
#' their AnnDataR6 method (users know those objects as AnnData).
#'
#' @param generic Name of the generic, e.g. `"fetch_data"`.
#'
#' @return A sorted character vector of class names.
#'
#' @noRd
.scuba_supported_classes <- function(generic) {
  methods_table <- asNamespace("SCUBA")[[".__S3MethodsTable__."]]
  prefix <- paste0(generic, ".")
  method_names <- ls(methods_table, all.names = TRUE)
  classes <- substring(
    method_names[startsWith(method_names, prefix)],
    nchar(prefix) + 1L
  )
  classes <- setdiff(
    classes,
    c("default", "anndata.AnnData", "anndata._core.anndata.AnnData")
  )
  sort(classes, method = "radix")
}

#' Extra explanation for Python objects reached through reticulate
#'
#' @param object The object a SCUBA generic could not handle.
#'
#' @return `""`, or a paragraph (starting with a blank line) explaining why a
#'   Python object was not dispatched and what to do about it.
#'
#' @noRd
.scuba_unsupported_hint <- function(object) {
  anndata_classes <- c("anndata.AnnData", "anndata._core.anndata.AnnData")
  if (inherits(object, anndata_classes)) {
    return(paste0(
      "\n\nThis is a Python AnnData object that the R anndata package did not ",
      "wrap as AnnDataR6. Either another package (such as anndataR) has ",
      "registered its own reticulate converter for this class, or the Python ",
      "environment reticulate is using is outside the stack SCUBA supports: ",
      paste(.scuba_python_requirements$packages, collapse = ", "),
      " (Python ", .scuba_python_requirements$python_version, "). ",
      "SCUBA declares these with reticulate::py_require(), but they are ",
      "ignored when RETICULATE_PYTHON, use_python() or use_condaenv() selects ",
      "another Python. With the Julia engine, open_h5ad() reads .h5ad files ",
      "without Python."
    ))
  }
  if (inherits(object, "python.builtin.object")) {
    return(paste0(
      "\n\nThis is a Python object reached through reticulate. Open .h5ad / ",
      ".h5mu files with open_h5ad() / open_h5mu(), or convert the object to ",
      "one of the supported classes."
    ))
  }
  ""
}

#' Stop with an informative error for an unsupported object
#'
#' Every `.default` method of a SCUBA generic calls this, so all of them fail
#' the same way: an error (not a warning that returns the message text) naming
#' the generic, the object's class and the classes SCUBA does support.
#'
#' @param object The object passed to the generic.
#' @param generic Name of the generic, e.g. `"fetch_data"`.
#'
#' @return Does not return; always signals an error.
#'
#' @noRd
.scuba_unsupported_class <- function(object, generic) {
  classes <- .scuba_supported_classes(generic)
  classes[classes == "scuba_source"] <-
    "scuba_source (from open_h5ad() / open_h5mu())"
  stop(
    generic, "() does not know how to handle an object of class '",
    paste(class(object), collapse = "', '"), "'. Supported classes: ",
    paste(classes, collapse = ", "), ".",
    .scuba_unsupported_hint(object),
    call. = FALSE
  )
}
