# Fast BPCells feature retrieval

Internal helper: subset BPCells-backed matrix directly by accessing the
underlying IterableMatrix. Falls back to
[`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md)
on any error to preserve robustness. This function bypasses SCUBA's
standard data retrieval pipeline for large feature queries when the
assay layer is BPCells-backed, improving performance by leveraging
direct matrix subsetting.

## Usage

``` r
bp_cells_fetch_features(
  object,
  features,
  assay = NULL,
  layer = NULL,
  cells = NULL
)
```

## Arguments

- object:

  A Seurat object with a BPCells-backed assay layer.

- features:

  A character vector of feature names to retrieve.

- assay:

  The assay containing the BPCells-backed layer. If `NULL`, the default
  assay will be used.

- layer:

  The layer to retrieve from within the assay. If `NULL`, the "data"
  layer will be used.

- cells:

  A character vector of cell names to include. If `NULL`, all cells in
  the object will be included.

## Value

A data.frame with cells as rows and features as columns, with feature
values extracted from the BPCells-backed matrix. Only features that are
present in the matrix are included in the output.
