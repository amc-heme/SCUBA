# Get keys for all assays/reductions in an object

Returns the "keys" of all reductions and modalities/assays/experiments
in an object, which are used to fetch data via the `vars` parameter of
`fetch_data`. To fetch features from an object, use the key representing
the modality the feature was recorded in, plus an underscore and the
feature name. To fetch reduction coordinates, use the key of the
reduction, plus an underscore, and a number representing the dimension
for which to retrieve coordinates.

## Usage

``` r
all_keys(object)

# S3 method for class 'Seurat'
all_keys(object)

# S3 method for class 'SingleCellExperiment'
all_keys(object)

# S3 method for class 'AnnDataR6'
all_keys(object)
```

## Arguments

- object:

  a single-cell object. Currently, Seurat, SingleCellExperiment, and
  anndata objects are supported.

## Value

a named character vector. The names of the vector are names of the
modalities and reductions in the object, and the values are the
corresponding keys to be passed to fetch_data. For Seurat objects, a key
for metadata will also be displayed.

## Methods (by class)

- `all_keys(Seurat)`: Seurat objects

- `all_keys(SingleCellExperiment)`: SingleCellExperiment objects

- `all_keys(AnnDataR6)`: SingleCellExperiment objects

## Examples

``` r
## View keys ##
# Seurat objects
all_keys(AML_Seurat)
#> meta.data       RNA        AB       pca      umap 
#>     "md_"    "rna_"     "ab_"     "PC_"   "UMAP_" 

# SingleCellExperiment objects 
all_keys(AML_SCE())
#>    RNA     AB    PCA   UMAP 
#>  "RNA"   "AB"  "PCA" "UMAP" 

# anndata objects 
all_keys(AML_h5ad())
#>         X     X_pca    X_umap   protein 
#>       "X"   "X_pca"  "X_umap" "protein" 

## Use of keys to construct fetch_data query
# Fetch a feature from the "protein" 
# modality using its key from above
fetch_data(
  AML_h5ad(), 
  vars = "protein_CD9-AB"
  ) |> str()
#> 'data.frame':    250 obs. of  1 variable:
#>  $ protein_CD9-AB: num  1.247 0.982 2.661 0.565 0.674 ...
#>  - attr(*, "pandas.index")=Index(['487013_1', '39207_1', '861619_1', '561110_1', '283967_1', '422573_1',
#>        '453256_1', '531766_1', '796968_1', '624345_1',
#>        ...
#>        '883406_2', '662718_2', '691696_2', '431743_2', '158371_2', '679107_2',
#>        '844492_2', '729807_2', '545562_2', '849364_2'],
#>       dtype='object', length=250)

# Fetch reduction coordinates using 
# the key for the UMAP reduction
fetch_data(
  AML_h5ad(), 
  vars = c("X_umap_1", "X_umap_2")
  ) |> str()
#> 'data.frame':    250 obs. of  2 variables:
#>  $ X_umap_1: num  -1.64 -1.5 -1.45 -1.38 -1.41 ...
#>  $ X_umap_2: num  9.9 10.13 10.21 10.51 3.39 ...
#>  - attr(*, "pandas.index")=Index(['487013_1', '39207_1', '861619_1', '561110_1', '283967_1', '422573_1',
#>        '453256_1', '531766_1', '796968_1', '624345_1',
#>        ...
#>        '883406_2', '662718_2', '691696_2', '431743_2', '158371_2', '679107_2',
#>        '844492_2', '729807_2', '545562_2', '849364_2'],
#>       dtype='object', length=250)
```
