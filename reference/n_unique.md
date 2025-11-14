# Number of unique values in a metadata variable

Returns the number of unique values in a single cell object, for the
specified metadata variable. This is a utility function useful for
directing the behavior of plotting functions created with SCUBA. For
example, you may have a plotting function that changes the number of
columns in a legend or the palette used when plotting a metadata
variable with a high number of values.

## Usage

``` r
n_unique(object, meta_var)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

- meta_var:

  Name of a metadata variable.

## Value

An integer giving the number of unique values in the specified metadata
variable.

## Examples

``` r
n_unique(AML_Seurat, meta_var = "ct")
#> [1] 27

metadata_variable <- "ct"

if (n_unique(AML_Seurat, metadata_variable) > 25){
  print("Execute alternate plotting code when a variable has more than 
  25 values.")
}
#> [1] "Execute alternate plotting code when a variable has more than \n  25 values."
```
