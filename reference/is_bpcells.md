# Detect BPCells-backed assay layer in a Seurat object

Internal helper: returns TRUE if the indicated assay/layer appears to be
backed by a BPCells on-disk matrix (heuristic: matrix for the layer is
an IterableMatrix).

## Usage

``` r
is_bpcells(object, assay = NULL, layer = NULL)
```

## Arguments

- object:

  A Seurat object to check.

- assay:

  The assay to check. If `NULL`, the default assay will be used.

- layer:

  The layer to check within the assay. If `NULL`, the "data" layer will
  be checked.

## Value

A logical scalar: `TRUE` if the specified assay/layer is backed by a
BPCells IterableMatrix, `FALSE` otherwise.
