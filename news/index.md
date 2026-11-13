# Changelog

## SCUBA 2.0.0 (development)

### Split Seurat v5 layers on the Julia engine

- With a SCUBA.jl that reads split layers (`counts.<group>` after
  [`split()`](https://rdrr.io/r/base/split.html), `counts.1` /
  `counts.2` after [`merge()`](https://rdrr.io/r/base/merge.html)),
  every accessor on such objects runs on the Julia engine: the logical
  layer names (`counts`, `data`) join the pieces,
  [`pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.md)
  aggregates them in one pass, and BPCells-backed pieces are read
  natively. An explicit piece (`layer = "counts.1"`) returns every
  requested cell, `NaN` outside the piece, where the R engine drops
  those cells.
- With an older SCUBA.jl (no split-layer support), Seurat objects with
  split layers use the R engine for every accessor. Previously every
  Julia-engine accessor failed on them, including
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md)
  with the default layer.

### Fixes from the review of the 2.0 branch

- [`fetch_matrix()`](https://amc-heme.github.io/SCUBA/reference/fetch_matrix.md)
  on the R engine accepts keyed names (`"rna_GAPDH"`, `"AB_CD117-AB"`),
  which select their assay as on the Julia engine; features from several
  assays are an error.
- [`fetch_matrix()`](https://amc-heme.github.io/SCUBA/reference/fetch_matrix.md)
  on Seurat v5 objects with split layers (`counts.1`, `counts.2`, …
  after [`split()`](https://rdrr.io/r/base/split.html) or
  [`merge()`](https://rdrr.io/r/base/merge.html)) returns every cell,
  zero-filling features a layer does not hold, instead of the first
  layer only.
- `fetch_matrix(sparse = TRUE)` converts other sparse classes,
  DelayedArray (HDF5-backed) and BPCells matrices to `dgCMatrix` without
  a dense copy.
- [`pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.md)
  on split layers always runs on the R path (SCUBA.jl cannot yet index
  cells within a layer), and aggregates layers holding different feature
  sets.
- [`pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.md)
  on an AnnData `obsm` modality requires `layer = "X"` instead of
  silently aggregating it under the default `layer = "counts"`.
- [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md)
  rejects conflicting `layer` and `slot` values on every engine and
  class, and keeps the very-wide-query warning on the Julia engine. A
  metadata column literally named `obs_<name>` on file-backed objects is
  no longer read as `<name>`.
- SeuratObject 4 stays supported:
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md)
  passes `slot =` to its `FetchData()`, and Seurat objects use the R
  engine under SeuratObject 4.
- The Julia engine no longer falls back to an embedded copy of its glue
  code for old SCUBA.jl versions (it did not match the R side); such
  versions are rejected, so `"auto"` uses the R engine. Paths and URLs
  reach Julia as values, so Windows paths in `SCUBA_JL_PATH` work.

### Julia engine

- SCUBA.jl is installed from its GitHub repository when the registry
  lookup fails (it is not registered yet); previously the engine was
  unavailable without a `SCUBA_JL_PATH` checkout.
  `options(SCUBA.jl_url)` / `SCUBA_JL_URL` override the source.

### Python environment

- reticulate-managed Python environments are pinned to
  `anndata>=0.10.9,<0.11`, `pandas<3`, `scipy<1.15` and Python \< 3.14
  (`py_require()` in `.onLoad`). The anndata R package, now superseded,
  only wraps Python AnnData objects as `AnnDataR6` for Python anndata \<
  0.11; later versions expose them as plain `anndata.AnnData` objects
  that SCUBA’s AnnData methods do not handle. This is an interim measure
  until AnnData support moves to anndataR (#173). Environments set
  through `RETICULATE_PYTHON` are not affected.
- [`pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.md)
  and
  [`layer_names()`](https://amc-heme.github.io/SCUBA/reference/layer_names.md)
  accept unwrapped Python AnnData objects under both class names
  (`anndata.AnnData`, `anndata._core.anndata.AnnData`).
- `fetch_data(slot = )` on Seurat objects is an alias of `layer`;
  SeuratObject \>= 5.4 no longer accepts `slot` in `FetchData()`.

### Pseudobulk aggregation

- New `pseudobulk(object, sample_by, group_by)` sums (or averages) a
  counts layer over sample x group profiles for Seurat (including
  BPCells layers and split v5 layers), SingleCellExperiment (in memory
  or HDF5-backed), AnnData (in memory or backed, through reticulate) and
  file-backed `.h5ad` / `.h5mu` objects
  ([`open_h5ad()`](https://amc-heme.github.io/SCUBA/reference/open_h5ad.md)
  /
  [`open_h5mu()`](https://amc-heme.github.io/SCUBA/reference/open_h5ad.md)),
  returning the small features x profiles matrix and a table of the
  profiles. The full features x cells matrix is never materialised: the
  Julia engine streams every backend into a group accumulator, and the R
  engine uses the class’s own block-wise primitive
  ([`BPCells::pseudobulk_matrix()`](https://bnprks.github.io/BPCells/reference/pseudobulk_matrix.html),
  [`DelayedArray::colsum()`](https://rdrr.io/pkg/S4Arrays/man/rowsum.html),
  a sparse indicator product). `sample_by` / `group_by` may be metadata
  column names or cell-named label vectors.
- New
  [`pseudobulk_design()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk_design.md)
  previews the same profiles from metadata alone (cells per sample x
  group, replicates per group, paired designs) without reading
  expression data.
- New
  [`layer_names()`](https://amc-heme.github.io/SCUBA/reference/layer_names.md)
  lists the layers of an assay for every supported class.

### Julia engine

- SCUBA’s accessors
  ([`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md),
  [`fetch_metadata()`](https://amc-heme.github.io/SCUBA/reference/fetch_metadata.md),
  [`fetch_reduction()`](https://amc-heme.github.io/SCUBA/reference/fetch_reduction.md),
  [`fetch_cells()`](https://amc-heme.github.io/SCUBA/reference/fetch_cells.md),
  [`get_all_cells()`](https://amc-heme.github.io/SCUBA/reference/get_all_cells.md),
  [`meta_varnames()`](https://amc-heme.github.io/SCUBA/reference/meta_varnames.md),
  [`all_keys()`](https://amc-heme.github.io/SCUBA/reference/all_keys.md),
  [`features_in_assay()`](https://amc-heme.github.io/SCUBA/reference/features_in_assay.md),
  [`reduction_dimnames()`](https://amc-heme.github.io/SCUBA/reference/reduction_dimnames.md),
  [`default_layer()`](https://amc-heme.github.io/SCUBA/reference/default_layer.md),
  [`default_reduction()`](https://amc-heme.github.io/SCUBA/reference/default_reduction.md),
  [`unique_values()`](https://amc-heme.github.io/SCUBA/reference/unique_values.md),
  [`n_unique()`](https://amc-heme.github.io/SCUBA/reference/n_unique.md),
  [`fetch_feature()`](https://amc-heme.github.io/SCUBA/reference/fetch_feature.md))
  can now run on the Julia package
  [SCUBA.jl](https://github.com/amc-heme/SCUBA.jl), which implements the
  variable-resolution rules once for R, Python and Julia. The engine is
  chosen with `options(SCUBA.engine = "auto" | "julia" | "r")`; `"auto"`
  (the default) uses Julia when the `JuliaCall` package and a Julia
  installation are available and falls back to the legacy R
  implementation otherwise. See
  [`scuba_engine()`](https://amc-heme.github.io/SCUBA/reference/scuba_engine.md),
  [`scuba_julia_available()`](https://amc-heme.github.io/SCUBA/reference/scuba_julia_available.md)
  and
  [`with_engine()`](https://amc-heme.github.io/SCUBA/reference/with_engine.md).
- Seurat and SingleCellExperiment objects are read in place through
  RCall; BPCells layers are read natively from disk by SCUBA.jl. AnnData
  objects continue to use the reticulate path.
- The Julia engine normalises behaviour that differed between object
  classes in 1.x: keys match only as `key_`, duplicates always collapse
  with a warning (including `GAPDH` + `rna_GAPDH`), partial misses
  produce one aggregated warning on every class, a single reduction
  dimension works on SingleCellExperiment objects, and unknown cells
  raise an error.

### Breaking changes

- `FetchData()` methods for SingleCellExperiment and AnnData objects
  (deprecated in 1.1.2) were removed; use
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md).
  [`plot_dot()`](https://amc-heme.github.io/SCUBA/reference/plot_dot.md),
  [`plot_scatter()`](https://amc-heme.github.io/SCUBA/reference/plot_scatter.md),
  [`plot_violin()`](https://amc-heme.github.io/SCUBA/reference/plot_violin.md)
  and
  [`plot_ridge()`](https://amc-heme.github.io/SCUBA/reference/plot_ridge.md)
  now use
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md)
  internally.
- [`default_slot()`](https://amc-heme.github.io/SCUBA/reference/default_slot.md)
  (deprecated in 0.10.0) now errors; use
  [`default_layer()`](https://amc-heme.github.io/SCUBA/reference/default_layer.md).
