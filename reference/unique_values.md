# Summarize values in a metadata variable

Returns the unique values for the metadata variable provided.
`unique_values` is a utility function for summarizing metadata in an
object.

## Usage

``` r
unique_values(object, var)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

- var:

  the metadata variable for which to return unique values.

## Value

a character vector giving the unique values in the specified metadata
variable.

## Examples

``` r
unique_values(AML_Seurat, var = "Batch")
#> [1] "BM_200AB"   "PBMC_200AB"

unique_values(AML_Seurat, var = "condensed_cell_type")
#>  [1] "Plasma cells"                 "Primitive"                   
#>  [3] "Dendritic cells"              "Plasmacytoid dendritic cells"
#>  [5] "BM Monocytes"                 "NK Cells"                    
#>  [7] "CD8+ T Cells"                 "B Cells"                     
#>  [9] "CD4+ T Cells"                 "PBMC Monocytes"              
```
