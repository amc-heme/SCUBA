# Fetch an expression matrix

Returns the expression values of a set of features (all features of the
assay by default) for a set of cells as a features x cells matrix,
sparse (`dgCMatrix`) by default. On the Julia engine the block is read
natively (BPCells directories, h5ad files); otherwise the layer is
sliced in R.

## Usage

``` r
fetch_matrix(
  object,
  features = NULL,
  ...,
  assay = NULL,
  layer = NULL,
  cells = NULL,
  sparse = TRUE
)

# S3 method for class 'AnnDataR6'
fetch_matrix(object, features = NULL, ..., matrix_location = features)

# S3 method for class 'Seurat'
fetch_matrix(
  object,
  features = NULL,
  ...,
  assay = NULL,
  layer = NULL,
  cells = NULL,
  sparse = TRUE
)

# S3 method for class 'SingleCellExperiment'
fetch_matrix(
  object,
  features = NULL,
  ...,
  assay = NULL,
  layer = NULL,
  cells = NULL,
  sparse = TRUE
)

# S3 method for class 'scuba_source'
fetch_matrix(
  object,
  features = NULL,
  ...,
  assay = NULL,
  layer = NULL,
  cells = NULL,
  sparse = TRUE
)
```

## Arguments

- object:

  a single-cell object (Seurat, SingleCellExperiment, `scuba_source`,
  AnnData).

- features:

  feature names (bare, or keyed with the assay key); `NULL` for all
  features.

- ...:

  Passed to methods (`matrix_location` for AnnData objects).

- assay:

  assay / modality (default: the default assay).

- layer:

  layer to read (default:
  [`default_layer()`](https://amc-heme.github.io/SCUBA/reference/default_layer.md)).

- cells:

  cells to include, in the returned order (default: all cells).

- sparse:

  return a `dgCMatrix` (default) or a dense matrix.

- matrix_location:

  (AnnData objects) the legacy slot to return, for example `"X"` or
  `c("obsm", "X_umap")`; defaults to `features`.

## Value

A features x cells matrix with dimnames.

## Details

For AnnData (reticulate) objects the legacy signature
`fetch_matrix(object, matrix_location)` still returns the raw slot
(`"X"`, `c("obsm", "X_umap")`, ...).

## Methods (by class)

- `fetch_matrix(AnnDataR6)`: AnnDataR6 objects
