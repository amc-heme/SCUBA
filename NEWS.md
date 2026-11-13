# SCUBA 2.0.0 (development)

## Unsupported classes are an error

* Calling a SCUBA generic on an object it has no method for now raises an error instead of
  warning and returning the message text (so `get_all_cells()` no longer returns that
  sentence where cell names were expected). The message names the right function
  (`fetch_data()` used to report itself as `fetch_metadata`), the object's class and every
  supported class (AnnDataR6, Seurat, SingleCellExperiment, scuba_source), read from the
  registered methods so it cannot drift. A Python AnnData object that the R anndata package
  did not wrap as AnnDataR6 (Python anndata >= 0.11 in a user-managed environment) gets a
  hint naming the Python stack SCUBA supports; other reticulate objects point to
  `open_h5ad()` / `open_h5mu()`.

## Fixes from the review of the 2.0 branch

* `fetch_matrix()` on the R engine accepts keyed names (`"rna_GAPDH"`, `"AB_CD117-AB"`),
  which select their assay as on the Julia engine; features from several assays are an error.
* `fetch_matrix()` on Seurat v5 objects with split layers (`counts.1`, `counts.2`, ... after
  `split()` or `merge()`) returns every cell, zero-filling features a layer does not hold,
  instead of the first layer only.
* `fetch_matrix(sparse = TRUE)` converts other sparse classes, DelayedArray (HDF5-backed) and
  BPCells matrices to `dgCMatrix` without a dense copy.
* `pseudobulk()` on split layers always runs on the R path (SCUBA.jl cannot yet index cells
  within a layer), and aggregates layers holding different feature sets.
* `pseudobulk()` on an AnnData `obsm` modality requires `layer = "X"` instead of silently
  aggregating it under the default `layer = "counts"`.
* `fetch_data()` rejects conflicting `layer` and `slot` values on every engine and class, and
  keeps the very-wide-query warning on the Julia engine. A metadata column literally named
  `obs_<name>` on file-backed objects is no longer read as `<name>`.
* SeuratObject 4 stays supported: `fetch_data()` passes `slot =` to its `FetchData()`, and
  Seurat objects use the R engine under SeuratObject 4.
* The Julia engine no longer falls back to an embedded copy of its glue code for old SCUBA.jl
  versions (it did not match the R side); such versions are rejected, so `"auto"` uses the R
  engine. Paths and URLs reach Julia as values, so Windows paths in `SCUBA_JL_PATH` work.

## Julia engine

* SCUBA.jl is installed from its GitHub repository when the registry lookup fails (it is not
  registered yet); previously the engine was unavailable without a `SCUBA_JL_PATH` checkout.
  `options(SCUBA.jl_url)` / `SCUBA_JL_URL` override the source.

## Python environment

* reticulate-managed Python environments are pinned to `anndata>=0.10.9,<0.11`,
  `pandas<3`, `scipy<1.15` and Python < 3.14 (`py_require()` in `.onLoad`). The
  anndata R package, now superseded, only wraps Python AnnData objects as
  `AnnDataR6` for Python anndata < 0.11; later versions expose them as plain
  `anndata.AnnData` objects that SCUBA's AnnData methods do not handle. This is
  an interim measure until AnnData support moves to anndataR (#173).
  Environments set through `RETICULATE_PYTHON` are not affected.
* `pseudobulk()` and `layer_names()` accept unwrapped Python AnnData objects under
  both class names (`anndata.AnnData`, `anndata._core.anndata.AnnData`).
* `fetch_data(slot = )` on Seurat objects is an alias of `layer`; SeuratObject >= 5.4
  no longer accepts `slot` in `FetchData()`.

## Pseudobulk aggregation

* New `pseudobulk(object, sample_by, group_by)` sums (or averages) a counts layer over
  sample x group profiles for Seurat (including BPCells layers and split v5 layers),
  SingleCellExperiment (in memory or HDF5-backed), AnnData (in memory or backed, through
  reticulate) and file-backed `.h5ad` / `.h5mu` objects (`open_h5ad()` / `open_h5mu()`),
  returning the small features x profiles matrix and a table of the profiles. The full
  features x cells matrix is never materialised: the Julia engine streams every backend into a
  group accumulator, and the R engine uses the class's own block-wise primitive
  (`BPCells::pseudobulk_matrix()`, `DelayedArray::colsum()`, a sparse indicator product).
  `sample_by` / `group_by` may be metadata column names or cell-named label vectors.
* New `pseudobulk_design()` previews the same profiles from metadata alone (cells per
  sample x group, replicates per group, paired designs) without reading expression data.
* New `layer_names()` lists the layers of an assay for every supported class.

## Julia engine

* SCUBA's accessors (`fetch_data()`, `fetch_metadata()`, `fetch_reduction()`, `fetch_cells()`,
  `get_all_cells()`, `meta_varnames()`, `all_keys()`, `features_in_assay()`,
  `reduction_dimnames()`, `default_layer()`, `default_reduction()`, `unique_values()`,
  `n_unique()`, `fetch_feature()`) can now run on the Julia package
  [SCUBA.jl](https://github.com/amc-heme/SCUBA.jl), which implements the variable-resolution
  rules once for R, Python and Julia. The engine is chosen with `options(SCUBA.engine =
  "auto" | "julia" | "r")`; `"auto"` (the default) uses Julia when the `JuliaCall` package and
  a Julia installation are available and falls back to the legacy R implementation otherwise.
  See `scuba_engine()`, `scuba_julia_available()` and `with_engine()`.
* Seurat and SingleCellExperiment objects are read in place through RCall; BPCells layers are
  read natively from disk by SCUBA.jl. AnnData objects continue to use the reticulate path.
* The Julia engine normalises behaviour that differed between object classes in 1.x: keys match
  only as `key_`, duplicates always collapse with a warning (including `GAPDH` + `rna_GAPDH`),
  partial misses produce one aggregated warning on every class, a single reduction dimension
  works on SingleCellExperiment objects, and unknown cells raise an error.

## Breaking changes

* `FetchData()` methods for SingleCellExperiment and AnnData objects (deprecated in 1.1.2) were
  removed; use `fetch_data()`. `plot_dot()`, `plot_scatter()`, `plot_violin()` and
  `plot_ridge()` now use `fetch_data()` internally.
* `default_slot()` (deprecated in 0.10.0) now errors; use `default_layer()`.
