# Plotting utility function to return all cell IDs

Returns a character vector with all cell names in the object. This is a
utility function used to set defaults for plotting functions created by
SCUBA

## Usage

``` r
get_all_cells(object)

# S3 method for class 'Seurat'
get_all_cells(object)

# S3 method for class 'SingleCellExperiment'
get_all_cells(object)

# S3 method for class 'AnnDataR6'
get_all_cells(object)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

## Details

For additional information, see our GitHub.io website ("User Guide"
article, "Data Visualization" section)

## Methods (by class)

- `get_all_cells(Seurat)`: Seurat objects

- `get_all_cells(SingleCellExperiment)`: SingleCellExperiment objects

- `get_all_cells(AnnDataR6)`: SingleCellExperiment objects

## Examples

``` r
get_all_cells(AML_Seurat) |> str()
#>  chr [1:250] "487013_1" "39207_1" "861619_1" "561110_1" "283967_1" ...
```
