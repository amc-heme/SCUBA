# Fetch feature expression data from single-cell objects

**This function is still in development.** Currently, it only has a full
implementation for Seurat objects. For all other object types, it acts
as a wrapper around
[`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md).

This function retrieves feature expression data from single-cell
objects. For Seurat objects with BPCells-backed assay layers, this
function uses optimized direct matrix access for improved performance.
For all other cases, it delegates to
[`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md).

Unlike
[`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md),
this function is specifically designed for feature retrieval and does
not support fetching metadata or reduction coordinates.

See our GitHub.io website for additional information and examples.

## Usage

``` r
fetch_feature(object, features, assay = NULL, layer = NULL, cells = NULL)

# S3 method for class 'Seurat'
fetch_feature(object, features, assay = NULL, layer = NULL, cells = NULL)

# S3 method for class 'SingleCellExperiment'
fetch_feature(object, features, assay = NULL, layer = NULL, cells = NULL)

# S3 method for class 'AnnDataR6'
fetch_feature(object, features, assay = NULL, layer = NULL, cells = NULL)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExperiment, and
  anndata objects are supported.

- features:

  A character vector of feature names to retrieve from the object.
  Features should be specified with the assay key prefix if pulling from
  a non-default assay (e.g., "rna_FLT3"). To determine the key that
  corresponds to the assay to pull data from, run
  [`all_keys()`](https://amc-heme.github.io/SCUBA/reference/all_keys.md).

- assay:

  For Seurat objects, the assay to pull data from. If `NULL`, the
  default assay will be used.

- layer:

  The layer to pull data from. Layers are referred to as "slots" in
  Seurat objects v4 and earlier, and "assays" in SingleCellExperiment
  objects. If `NULL`, the default layer will be used.

- cells:

  A character vector of cell names to include, as they are named in the
  object (i.e. according to `colnames(object)`). If `NULL`, data will be
  returned for all cells in the object.

## Value

A data.frame with the requested `features` as columns and the cells as
rows.

## Methods (by class)

- `fetch_feature(Seurat)`: Seurat objects. Uses optimized BPCells access
  when available, otherwise delegates to
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md).

- `fetch_feature(SingleCellExperiment)`: SingleCellExperiment objects.
  Wrapper for
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md).

- `fetch_feature(AnnDataR6)`: AnnDataR6 objects. Wrapper for
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md).

## Examples

``` r
# Fetch feature expression data
fetch_feature(AML_Seurat, features = c("FLT3", "NPM1")) |> head()
#>          FLT3     NPM1
#> 487013_1    0 3.056462
#> 39207_1     0 2.801581
#> 861619_1    0 1.548350
#> 561110_1    0 2.260502
#> 283967_1    0 4.935888
#> 422573_1    0 2.530927

# Fetch from a specific layer
fetch_feature(AML_Seurat, features = "FLT3", layer = "data") |> head()
#>          FLT3
#> 487013_1    0
#> 39207_1     0
#> 861619_1    0
#> 561110_1    0
#> 283967_1    0
#> 422573_1    0
```
