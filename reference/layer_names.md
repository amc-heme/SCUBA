# Layers of an assay

Names of the layers (Seurat), assays (SingleCellExperiment experiment),
matrices (AnnData: `"X"`, every entry of `layers`, and `"raw"` when
`adata.raw` is set) or layers (`.h5ad` / `.h5mu` modality) available for
an assay. Use it to find a raw-count layer before calling
[`pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.md).

## Usage

``` r
layer_names(object, assay = NULL)

# S3 method for class 'Seurat'
layer_names(object, assay = NULL)

# S3 method for class 'SingleCellExperiment'
layer_names(object, assay = NULL)

# S3 method for class 'AnnDataR6'
layer_names(object, assay = NULL)

# S3 method for class 'anndata.AnnData'
layer_names(object, assay = NULL)

# S3 method for class 'anndata._core.anndata.AnnData'
layer_names(object, assay = NULL)

# S3 method for class 'scuba_source'
layer_names(object, assay = NULL)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExperiment, and
  anndata objects are supported.

- assay:

  Assay / modality to inspect (default: the default assay).

## Value

A character vector of layer names (empty when the assay has none).

## Examples

``` r
layer_names(AML_Seurat)
#> [1] "counts"     "data"       "scale.data"
layer_names(AML_SCE(), assay = "AB")
#> [1] "counts"    "logcounts"
```
