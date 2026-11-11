# Metadata for multiple single cell object classes

Updates the metadata of either a Seurat, SingleCellExperiment or Anndata
object. This method is intended to be used with fetch_metadata to pull a
table from an object, modify it, and then store the updated table.

## Usage

``` r
update_object_metadata(object, table)
```

## Arguments

- object:

  a single-cell object. Currently, Seurat, SingleCellExperiment and
  AnnData objects are supported.

- table:

  a modified metadata table.

## Value

the object passed to `object` with the modified metadata table.

## Examples

``` r
table <- fetch_metadata(AML_Seurat, full_table = TRUE)
table$flag <- table$nCount_RNA > stats::median(table$nCount_RNA)
updated <- update_object_metadata(AML_Seurat, table)
"flag" %in% meta_varnames(updated)
#> [1] TRUE
```
