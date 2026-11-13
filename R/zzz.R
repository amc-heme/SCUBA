.onLoad <- function(libname, pkgname) {
  # Establish Python package dependencies only when reticulate is available.
  # reticulate is in Suggests, so it may not be installed. Python-backed
  # functionality (AnnData objects) will fail with an informative error at
  # call time when reticulate or anndata are absent.
  if (requireNamespace("reticulate", quietly = TRUE)) {
    tryCatch(
      {
        # The anndata R package (superseded) wraps Python AnnData objects as
        # AnnDataR6 only for Python anndata < 0.11; from 0.11 the class is
        # exposed as `anndata.AnnData` and SCUBA's AnnData methods are not
        # dispatched. That series needs pandas < 3 and scipy < 1.15 (backed
        # matrices call scipy's removed `_validate_indices`), and scipy 1.14
        # has no wheels for Python 3.14. Freeze that stack for
        # reticulate-managed environments until the anndataR migration
        # (https://github.com/amc-heme/SCUBA/issues/173); environments the
        # user manages (RETICULATE_PYTHON) are not affected. The pins live in
        # .scuba_python_requirements (R/unsupported_class.R) so error
        # messages can quote them.
        reticulate::py_require(
          .scuba_python_requirements$packages,
          python_version = .scuba_python_requirements$python_version
        )
      },
      error = function(e) {
        warning(
          paste(
            "Optional Python dependencies for AnnData support could not be",
            "configured during package load. AnnData functionality may be",
            "unavailable until the Python environment is set up.",
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
  }
}
