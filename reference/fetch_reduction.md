# Fetch reduction coordinates from single-cell objects

Returns the reduction coordinates matching the name and dimensions
supplied.

## Usage

``` r
fetch_reduction(object, reduction, cells = NULL, dims = c(1, 2))

# S3 method for class 'Seurat'
fetch_reduction(object, reduction, cells = NULL, dims = c(1, 2))

# S3 method for class 'SingleCellExperiment'
fetch_reduction(object, reduction, cells = NULL, dims = c(1, 2))

# S3 method for class 'AnnDataR6'
fetch_reduction(object, reduction, cells = NULL, dims = c(1, 2))
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

- reduction:

  the name of the reduction to pull coordinates from.

- cells:

  cell IDs for which to pull reduction data. If `NULL`, coordinates will
  be returned from all cells in the object. Cell IDs can be generated
  with
  [`fetch_cells()`](https://amc-heme.github.io/SCUBA/reference/fetch_cells.md).

- dims:

  a numeric vector indicating the dimensions to pull. Currently, only
  two dimensions are supported, but
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md)
  supports more than two dimensions. For instructions on pulling more
  than two dimensions at once, see the examples of
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md).

## Details

See our GitHub.io website for additional information and examples.

## Methods (by class)

- `fetch_reduction(Seurat)`: Seurat objects

- `fetch_reduction(SingleCellExperiment)`: SingleCellExperiment objects

- `fetch_reduction(AnnDataR6)`: AnnDataR6 objects

## Examples

``` r
# Return the first and second dimensions from the UMAP reduction
fetch_reduction(
  AML_Seurat,
  reduction = "umap",
  dims = c(1, 2)
  ) |> str()
#> 'data.frame':    250 obs. of  2 variables:
#>  $ UMAP_1: num  -1.64 -1.5 -1.45 -1.38 -1.41 ...
#>  $ UMAP_2: num  9.9 10.13 10.21 10.51 3.39 ...
```
