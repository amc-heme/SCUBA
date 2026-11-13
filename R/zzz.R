.onLoad <- function(libname, pkgname) {
  # Establish Python package dependencies only when reticulate is available.
  # reticulate is in Suggests, so it may not be installed. Python-backed
  # functionality (AnnData objects) will fail with an informative error at
  # call time when reticulate or anndata are absent.
  if (requireNamespace("reticulate", quietly = TRUE)) {
    tryCatch(
      {
        # Python stack for reticulate-managed environments. The bounds and
        # the reason for each live in .scuba_python_requirements
        # (R/unsupported_class.R), so error messages can quote them.
        # Environments the user manages (RETICULATE_PYTHON, use_python(),
        # use_condaenv()) ignore py_require(); the converters registered
        # below cover those too.
        packages <- .scuba_python_requirements$packages
        # reticulate converts pandas 3 data frames only from 1.45.0; keep
        # an older reticulate on pandas 2 rather than hand it objects it
        # cannot convert.
        if (utils::packageVersion("reticulate") < "1.45.0") {
          packages <- c(packages, "pandas<3")
        }
        reticulate::py_require(
          packages,
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
    # Converters for the class names Python anndata >= 0.13 and pandas >= 3
    # use (R/anndata_shim.R). Only when the anndata R package is installed:
    # the converters wrap objects in its AnnDataR6 class. system.file() checks
    # for it without loading it.
    if (nzchar(system.file(package = "anndata"))) {
      tryCatch(
        .scuba_register_py_converters(),
        error = function(e) {
          warning(
            paste(
              "SCUBA could not register its reticulate converters for",
              "Python anndata >= 0.13 / pandas >= 3; AnnData objects read",
              "with those versions may not be wrapped as AnnDataR6.",
              conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )
    }
  }
}
