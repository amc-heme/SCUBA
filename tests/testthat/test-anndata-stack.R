# AnnData through reticulate on every supported Python stack: objects are
# wrapped as AnnDataR6 (R/anndata_shim.R covers anndata >= 0.13 / pandas >= 3),
# every X storage format reads the same in memory and backed, a realistically
# sized backed file aggregates correctly, and SCUBA's Python code raises no
# deprecation warnings.

skip_if_no_python_anndata <- function() {
  skip_if_not_installed("reticulate")
  skip_if_not_installed("anndata")
  skip_if_not(
    tryCatch(reticulate::py_module_available("anndata"), error = function(e) FALSE),
    "python anndata not available"
  )
}

# Write CSR, CSC and dense copies of the bundled fixture (X and every layer)
# to `directory`, with Python anndata. Returns the file paths, by format.
write_format_copies <- function(directory) {
  reticulate::py_run_string(sprintf(r"(
import anndata as _ad, numpy as _np, scipy.sparse as _sp
_source = _ad.read_h5ad(r"%s")
_convert = {
    "csr": _sp.csr_matrix,
    "csc": _sp.csc_matrix,
    "dense": lambda m: m.toarray() if _sp.issparse(m) else _np.asarray(m),
}
for _format, _function in _convert.items():
    _copy = _source.copy()
    _copy.X = _function(_copy.X)
    for _layer in [k for k in _copy.layers.keys() if isinstance(k, str)]:
        _copy.layers[_layer] = _function(_copy.layers[_layer])
    _copy.write_h5ad(r"%s/AML_" + _format + ".h5ad")
)", system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA"), directory))
  stats::setNames(
    file.path(directory, paste0("AML_", c("csr", "csc", "dense"), ".h5ad")),
    c("csr", "csc", "dense")
  )
}

test_that("read_h5ad() returns AnnDataR6 with character obs and var names", {
  skip_if_no_python_anndata()
  for (object in list(AML_h5ad(), AML_h5ad_backed())) {
    expect_s3_class(object, "AnnDataR6")
    expect_type(object$obs_names, "character")
    expect_type(object$var_names, "character")
    expect_identical(get_all_cells(object), object$obs_names)
  }
})

test_that("converters are registered only for classes without a method", {
  skip_if_not_installed("reticulate")
  reticulate_namespace <- asNamespace("reticulate")
  methods_table <- reticulate_namespace[[".__S3MethodsTable__."]]
  test_classes <- c("scuba.TestTaken", "scuba.TestFree")
  withr::defer(rm(list = intersect(paste0("py_to_r.", test_classes), ls(methods_table)),
                  envir = methods_table))
  existing <- function(x) "existing"
  registerS3method("py_to_r", "scuba.TestTaken", existing, envir = reticulate_namespace)
  ours <- function(x) "ours"
  registered <- SCUBA:::.scuba_register_py_converters(
    list(scuba.TestTaken = ours, scuba.TestFree = ours)
  )
  expect_identical(registered, "scuba.TestFree")
  expect_identical(
    utils::getS3method("py_to_r", "scuba.TestTaken", envir = reticulate_namespace),
    existing
  )
  expect_identical(
    utils::getS3method("py_to_r", "scuba.TestFree", envir = reticulate_namespace),
    ours
  )
})

test_that("CSR, CSC and dense X read the same in memory and backed", {
  skip_if_no_python_anndata()
  paths <- write_format_copies(withr::local_tempdir())
  features <- c("GAPDH", "FLT3", "MEIS1", "CD34", "MPO")
  reference <- AML_h5ad()
  expected_data <- fetch_data(reference, c(features, "condensed_cell_type"))
  expected_counts <- fetch_data(reference, paste0("X_", features), layer = "counts")
  expected_pseudobulk <- pseudobulk(reference, "Batch", "condensed_cell_type")
  cells <- rev(get_all_cells(reference))[1:60]
  for (format in names(paths)) {
    for (backed in c(FALSE, TRUE)) {
      label <- paste(format, if (backed) "backed" else "in memory")
      object <- if (backed) {
        anndata::read_h5ad(paths[[format]], backed = "r")
      } else {
        anndata::read_h5ad(paths[[format]])
      }
      expect_equal(fetch_data(object, c(features, "condensed_cell_type")),
                   expected_data, ignore_attr = TRUE, label = label)
      expect_equal(fetch_data(object, paste0("X_", features), layer = "counts"),
                   expected_counts, ignore_attr = TRUE, label = label)
      expect_equal(fetch_data(object, features, cells = cells),
                   expected_data[cells, features], ignore_attr = TRUE, label = label)
      result <- pseudobulk(object, "Batch", "condensed_cell_type")
      expect_equal(result$counts, expected_pseudobulk$counts, label = label)
      expect_identical(result$profiles$n_cells, expected_pseudobulk$profiles$n_cells)
    }
  }
})

test_that("a realistically sized backed file aggregates and reads correctly", {
  skip_on_cran()
  skip_if_no_python_anndata()
  path <- file.path(withr::local_tempdir(), "synthetic.h5ad")
  # 20,000 cells x 2,000 genes at 5% density (about 2 million non-zeros) in
  # 8 samples; the Python reference sums are an indicator-matrix product.
  reticulate::py_run_string(sprintf(r"(
import anndata as _ad, numpy as _np, pandas as _pd, scipy.sparse as _sp
_rng = _np.random.default_rng(173)
_counts = _sp.random(20000, 2000, density=0.05, format="csr", random_state=_rng,
                     data_rvs=lambda n: _rng.integers(1, 20, n)).astype(_np.float32)
_samples = _np.array(["s%%d" %% i for i in range(8)])[_rng.integers(0, 8, 20000)]
_obs = _pd.DataFrame({"sample": _pd.Categorical(_samples)},
                     index=["cell%%d" %% i for i in range(20000)])
_var = _pd.DataFrame(index=["gene%%d" %% i for i in range(2000)])
_ad.AnnData(X=_counts, obs=_obs, var=_var, layers={"counts": _counts}).write_h5ad(r"%s")
_indicator = _sp.csr_matrix((_np.ones(20000),
                             (_obs["sample"].cat.codes.to_numpy(), _np.arange(20000))),
                            shape=(8, 20000))
scuba_reference_sums = _np.asarray((_indicator @ _counts).todense()).T
scuba_reference_samples = list(_obs["sample"].cat.categories)
scuba_reference_genes = _counts[:, [0, 999, 1999]].toarray()
)", path))
  python_main <- reticulate::py
  object <- anndata::read_h5ad(path, backed = "r")
  result <- pseudobulk(object, "sample")
  expected <- python_main$scuba_reference_sums
  dimnames(expected) <- list(paste0("gene", 0:1999), python_main$scuba_reference_samples)
  expect_equal(as.matrix(result$counts)[, colnames(expected)], expected, ignore_attr = TRUE)
  expect_identical(sum(result$profiles$n_cells), 20000L)
  genes <- fetch_data(object, c("gene0", "gene999", "gene1999"), layer = "counts")
  expect_equal(unname(as.matrix(genes)), python_main$scuba_reference_genes,
               ignore_attr = TRUE)
})

test_that("SCUBA's AnnData code raises no Python deprecation warnings", {
  skip_if_no_python_anndata()
  # Record every Python warning raised while SCUBA's AnnData functions run.
  # A deprecation attributed to SCUBA's scripts, or to "<sys>" (a call made
  # from R through reticulate), is SCUBA's to fix.
  reticulate::py_run_string(r"(
import warnings as _warnings
_scuba_recorder = _warnings.catch_warnings(record=True)
scuba_recorded_warnings = _scuba_recorder.__enter__()
_warnings.simplefilter("always")
)")
  withr::defer(reticulate::py_run_string("_scuba_recorder.__exit__(None, None, None)"))
  for (object in list(AML_h5ad(), AML_h5ad_backed())) {
    fetch_data(object, c("X_GAPDH", "MEIS1", "condensed_cell_type", "X_umap_1", "protein_CD34-AB"))
    meta_varnames(object)
    all_keys(object)
    default_reduction(object)
    fetch_reduction(object, "X_umap")
    features_in_assay(object, "protein")
    fetch_matrix(object, c("obsm", "X_umap"))
    pseudobulk(object, "Batch", "condensed_cell_type")
    layer_names(object)
  }
  script_directory <- normalizePath(system.file("extdata", "Python", package = "SCUBA"))
  reticulate::py_run_string(sprintf(r"(
scuba_deprecations = [
    "%%s:%%s %%s: %%s" %% (w.filename, w.lineno, w.category.__name__, w.message)
    for w in scuba_recorded_warnings
    if issubclass(w.category, (DeprecationWarning, FutureWarning))
    and (w.filename.startswith(r"%s") or w.filename == "<sys>")
]
)", script_directory))
  expect_identical(reticulate::py$scuba_deprecations, list())
})
