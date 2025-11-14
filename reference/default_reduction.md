# Get default reduction from object

Returns the default reduction for the single-cell object passed.
`default_reduction` will check for reductions in the order below. If a
reduction exists, it will be returned. If not, the function will check
for the next reduction on the list. If none of the reductions in the
list exist, the function will return an error.

1.  UMAP

2.  t-SNE

3.  PCA

## Usage

``` r
default_reduction(object)

# S3 method for class 'Seurat'
default_reduction(object)

# S3 method for class 'SingleCellExperiment'
default_reduction(object)

# S3 method for class 'AnnDataR6'
default_reduction(object)
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

- `default_reduction(Seurat)`: Seurat objects

  In Seurat objects, `default_reduction` is a wrapper for
  [`SeuratObject::DefaultDimReduc()`](https://satijalab.github.io/seurat-object/reference/DefaultDimReduc.html).

- `default_reduction(SingleCellExperiment)`: SingleCellExperiment
  objects

  The search order for SingleCellExperiment objects is below.
  `fetch_reduction` will search in the order described for a reduction
  named exactly as described.

  1.  UMAP

  2.  TSNE

  3.  PCA

- `default_reduction(AnnDataR6)`: AnnDataR6 objects

  The search order for anndata objects is below. `fetch_reduction` will
  search in the order described for a reduction named exactly as
  described.

  1.  X_umap

  2.  X_tsne

  3.  X_pca

## Examples

``` r
# Seurat objects
default_reduction(AML_Seurat)
#> [1] "umap"

# SingleCellExperiment objects
default_reduction(AML_SCE())
#> [1] "UMAP"

# anndata objects
default_reduction(AML_h5ad())
#> [1] "X_umap"
```
