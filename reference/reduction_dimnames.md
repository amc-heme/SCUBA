# From names of reduction keys for fetch_data

Given the name of a reduction and a set of dimensions, this function
will return the names of each dimension as it appears in the reduction
matrix. The output of this function can be passed to the `vars`
parameter of
[`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md)
to facilitate the specification of reduction coordinates to return from
this function.

## Usage

``` r
reduction_dimnames(object, reduction, dims)

# S3 method for class 'Seurat'
reduction_dimnames(object, reduction, dims)

# S3 method for class 'SingleCellExperiment'
reduction_dimnames(object, reduction, dims)

# S3 method for class 'AnnDataR6'
reduction_dimnames(object, reduction, dims)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

- reduction:

  the name of the reduction.

- dims:

  a numeric vector with the dimensions for which names should be
  returned.

## Methods (by class)

- `reduction_dimnames(Seurat)`: Seurat objects

- `reduction_dimnames(SingleCellExperiment)`: SingleCellExperiment
  objects

- `reduction_dimnames(AnnDataR6)`: AnnDataR6 objects

## Examples

``` r
# From names for first and second UMAP dimensions and 
# pass to fetch_data
dimnames <- reduction_dimnames(
  AML_Seurat,
  reduction = "umap",
  dims = c(1, 2)
  )
  
dimnames
#> [1] "UMAP_1" "UMAP_2"

fetch_data(
  AML_Seurat,
  vars = dimnames
  ) |> str()
#> 'data.frame':    250 obs. of  2 variables:
#>  $ UMAP_1: num  -1.64 -1.5 -1.45 -1.38 -1.41 ...
#>  $ UMAP_2: num  9.9 10.13 10.21 10.51 3.39 ...
  
  
# Form names for first five PCA dimensions and 
# pass to fetch_data
dimnames <- reduction_dimnames(
  AML_Seurat,
  reduction = "pca",
  dims = c(1:5)
  )
  
dimnames
#> [1] "PC_1" "PC_2" "PC_3" "PC_4" "PC_5"

fetch_data(
  AML_Seurat,
  vars = dimnames
  ) |> str()
#> 'data.frame':    250 obs. of  5 variables:
#>  $ PC_1: num  3.91 3.67 2.97 2.88 42.94 ...
#>  $ PC_2: num  1.676 -0.157 0.716 1.103 -11.258 ...
#>  $ PC_3: num  1.86 3.03 2.6 3.68 -4 ...
#>  $ PC_4: num  1.308 0.649 -0.415 0.734 25.255 ...
#>  $ PC_5: num  -1.09 -2.1 -1.68 -1.16 8.28 ...
```
