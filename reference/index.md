# Package index

## Data exploration

- [`meta_varnames()`](https://amc-heme.github.io/SCUBA/reference/meta_varnames.md)
  : Summarize object metadata
- [`unique_values()`](https://amc-heme.github.io/SCUBA/reference/unique_values.md)
  : Summarize values in a metadata variable

## Data access

- [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md)
  : Fetch feature expression data, reduction coordinates, or metadata
  from single-cell objects
- [`fetch_metadata()`](https://amc-heme.github.io/SCUBA/reference/fetch_metadata.md)
  : Fetch metadata from single-cell objects
- [`fetch_reduction()`](https://amc-heme.github.io/SCUBA/reference/fetch_reduction.md)
  : Fetch reduction coordinates from single-cell objects
- [`fetch_matrix()`](https://amc-heme.github.io/SCUBA/reference/fetch_matrix.md)
  : Fetch an expression matrix

## Aggregation

- [`pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.md)
  : Pseudobulk profiles: aggregate a counts layer over groups of cells
- [`pseudobulk_design()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk_design.md)
  : Preview the pseudobulk profiles of an object
- [`layer_names()`](https://amc-heme.github.io/SCUBA/reference/layer_names.md)
  : Layers of an assay

## Utility functions

- [`reduction_dimnames()`](https://amc-heme.github.io/SCUBA/reference/reduction_dimnames.md)
  : From names of reduction keys for fetch_data
- [`fetch_cells()`](https://amc-heme.github.io/SCUBA/reference/fetch_cells.md)
  : Return cell IDs for a subset based on metadata
- [`all_keys()`](https://amc-heme.github.io/SCUBA/reference/all_keys.md)
  : Get keys for all assays/reductions in an object
- [`features_in_assay()`](https://amc-heme.github.io/SCUBA/reference/features_in_assay.md)
  : Get names of all features in an assay/experiment/modality
- [`default_layer()`](https://amc-heme.github.io/SCUBA/reference/default_layer.md)
  : Return the default layer
- [`default_reduction()`](https://amc-heme.github.io/SCUBA/reference/default_reduction.md)
  : Get default reduction from object
- [`get_all_cells()`](https://amc-heme.github.io/SCUBA/reference/get_all_cells.md)
  : Plotting utility function to return all cell IDs
- [`n_unique()`](https://amc-heme.github.io/SCUBA/reference/n_unique.md)
  : Number of unique values in a metadata variable
- [`update_object_metadata()`](https://amc-heme.github.io/SCUBA/reference/update_object_metadata.md)
  : Metadata for multiple single cell object classes

## Julia engine and file-backed objects

- [`scuba_engine()`](https://amc-heme.github.io/SCUBA/reference/scuba_engine.md)
  : Which engine SCUBA is using
- [`with_engine()`](https://amc-heme.github.io/SCUBA/reference/with_engine.md)
  : Run an expression with a given SCUBA engine
- [`scuba_julia_available()`](https://amc-heme.github.io/SCUBA/reference/scuba_julia_available.md)
  : Is the Julia engine available?
- [`scuba_julia_warmup()`](https://amc-heme.github.io/SCUBA/reference/scuba_julia_warmup.md)
  : Warm up the Julia engine
- [`scuba_build_sysimage()`](https://amc-heme.github.io/SCUBA/reference/scuba_build_sysimage.md)
  : Build a Julia system image with SCUBA.jl compiled in
- [`open_h5ad()`](https://amc-heme.github.io/SCUBA/reference/open_h5ad.md)
  [`open_h5mu()`](https://amc-heme.github.io/SCUBA/reference/open_h5ad.md)
  : Open an .h5ad or .h5mu file as a SCUBA object (Julia engine)
- [`scuba_source_modalities()`](https://amc-heme.github.io/SCUBA/reference/scuba_source_modalities.md)
  [`scuba_source_reductions()`](https://amc-heme.github.io/SCUBA/reference/scuba_source_modalities.md)
  : Modalities and reductions of a file-backed SCUBA object

## Datasets

- [`AML_Seurat`](https://amc-heme.github.io/SCUBA/reference/AML_Seurat.md)
  : Reference dataset used for testing and demonstration (Seurat)
- [`AML_SCE()`](https://amc-heme.github.io/SCUBA/reference/AML_SCE.md) :
  Reference dataset used for testing and demonstration
  (SingleCellExperiment)
- [`AML_h5ad()`](https://amc-heme.github.io/SCUBA/reference/anndata_example.md)
  [`AML_h5ad_backed()`](https://amc-heme.github.io/SCUBA/reference/anndata_example.md)
  : Reference dataset used for testing and demonstration (anndata)
- [`AML_BPCells`](https://amc-heme.github.io/SCUBA/reference/AML_BPCells.md)
  : Reference dataset used for testing and demonstration (Seurat with
  BPCells-backed matrices)
