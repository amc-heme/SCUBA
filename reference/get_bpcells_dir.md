# Get BPCells directory path

Internal helper: retrieves the directory path for a BPCells-backed
matrix in a Seurat object. This function validates that the specified
assay/layer is BPCells-backed before attempting to access the directory
path.

## Usage

``` r
get_bpcells_dir(object, assay = NULL, layer = NULL)
```

## Arguments

- object:

  A Seurat object containing BPCells-backed assay layers.

- assay:

  The name of the assay containing the BPCells-backed layer. If `NULL`,
  the default assay will be used.

- layer:

  The name of the layer within the assay to query. If `NULL`, the "data"
  layer will be used.

## Value

A character string containing the directory path of the BPCells matrix,
or an error if the layer is not BPCells-backed.
