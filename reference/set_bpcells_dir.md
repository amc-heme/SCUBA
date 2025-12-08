# Set BPCells directory path

Internal helper: sets the directory path for a BPCells-backed matrix in
a Seurat object. This function validates that the specified assay/layer
is BPCells-backed before attempting to modify the directory path.

## Usage

``` r
set_bpcells_dir(object, assay, layer, dirname)
```

## Arguments

- object:

  A Seurat object containing BPCells-backed assay layers.

- assay:

  The name of the assay containing the BPCells-backed layer.

- layer:

  The name of the layer within the assay to modify.

- dirname:

  The target directory path to set for the BPCells matrix. Can be
  specified as either an absolute or relative path; the path will be
  normalized to an absolute path before storage.

## Value

The modified Seurat object (invisibly). The object is returned invisibly
so it will not print when called interactively. To apply the changes,
you must assign the output back to the object:
`obj <- set_bpcells_dir(obj, "RNA", "data", "/path")`. Alternatively,
use the pipe operator:
`obj <- obj |> set_bpcells_dir("RNA", "data", "/path")`.
