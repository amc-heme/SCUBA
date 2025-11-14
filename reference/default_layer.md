# Return the default layer

Returns the default layer for the object passed. The default layer is
chosen based on the conventions used in each object to name the
normalized counts layer.

## Usage

``` r
default_layer(object)

# S3 method for class 'Seurat'
default_layer(object)

# S3 method for class 'SingleCellExperiment'
default_layer(object)

# S3 method for class 'AnnDataR6'
default_layer(object)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

## Details

This is a utility function that is most useful for defining defaults in
plotting functions, to reduce the number of required parameters end
users need to supply. [see our
website](https://amc-heme.github.io/SCUBA/reference/see%20our%20website)
for an example of usage in a function.

## Methods (by class)

- `default_layer(Seurat)`: Seurat objects

  For Seurat objects, the default layer is "data".

- `default_layer(SingleCellExperiment)`: SingleCellExperiment objects

  For SingleCellExperiment objects, the default layer is "logcounts".

- `default_layer(AnnDataR6)`: Anndata objects

  For Anndata objects, the default layer is `NULL`, which will direct
  FetchData to pull feature epxression data from the X matrix.

## Examples

``` r
# Seurat objects
default_layer(AML_Seurat)
#> [1] "data"

# SingleCellExperiment objects
default_layer(AML_SCE())
#> [1] "logcounts"

# anndata objects
default_layer(AML_h5ad())
#> NULL
```
