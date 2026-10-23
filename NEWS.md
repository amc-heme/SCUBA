# SCUBA 2.0.0 (development)

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
