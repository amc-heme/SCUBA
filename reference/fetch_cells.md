# Return cell IDs for a subset based on metadata

Returns a character vector with the cell names matching the
levels/classes of a specified metadata variable.

## Usage

``` r
fetch_cells(object, meta_var, meta_levels)

# S3 method for class 'Seurat'
fetch_cells(object, meta_var, meta_levels)

# S3 method for class 'SingleCellExperiment'
fetch_cells(object, meta_var, meta_levels)

# S3 method for class 'AnnDataR6'
fetch_cells(object, meta_var, meta_levels)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

- meta_var:

  a metadata variable used as the basis for subsetting cells.

- meta_levels:

  the levels of the specified metadata variable in `meta_var` to include
  in the

## Methods (by class)

- `fetch_cells(Seurat)`: Seurat objects

- `fetch_cells(SingleCellExperiment)`: SingleCellExperiment objects

- `fetch_cells(AnnDataR6)`: Anndata objects
