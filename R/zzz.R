.onLoad <- function(libname, pkgname) {
  # Establish Python package dependencies only when reticulate is available.
  # reticulate is in Suggests, so it may not be installed. Python-backed
  # functionality (AnnData objects) will fail with an informative error at
  # call time when reticulate or anndata are absent.
  if (requireNamespace("reticulate", quietly = TRUE)) {
    reticulate::py_require("anndata>=0.11.4")
    reticulate::py_require("pandas>=2.0.0")
    reticulate::py_require("numpy")
    reticulate::py_require("scipy>=1.14.0")
  }
}
