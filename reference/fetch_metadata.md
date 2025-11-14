# Fetch metadata from single-cell objects

Returns object metadata for a specified set of cells.

## Usage

``` r
fetch_metadata(
  object,
  vars = NULL,
  cells = NULL,
  full_table = FALSE,
  return_class = "dataframe"
)

# S3 method for class 'Seurat'
fetch_metadata(
  object,
  vars = NULL,
  cells = NULL,
  full_table = FALSE,
  return_class = "dataframe"
)

# S3 method for class 'SingleCellExperiment'
fetch_metadata(
  object,
  vars = NULL,
  cells = NULL,
  full_table = FALSE,
  return_class = "dataframe"
)

# S3 method for class 'AnnDataR6'
fetch_metadata(
  object,
  vars = NULL,
  cells = NULL,
  full_table = FALSE,
  return_class = "dataframe"
)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

- vars:

  metadata variables to pull from object. This must be defined, unless
  "full_table" is set to `TRUE`.

- cells:

  cell IDs for which to pull metadata. If `NULL`, coordinates will be
  returned from all cells in the object. Cell IDs can be generated with
  [`fetch_cells()`](https://amc-heme.github.io/SCUBA/reference/fetch_cells.md).

- full_table:

  if `TRUE`, return the entire metadata table. This is `FALSE` by
  default.

- return_class:

  class of data returned. Set to "dataframe" by default to return a
  data.frame, and may also be set to "vector" to yield a vector of
  values. This is ignored if "full_table" is set to `TRUE`.

## Details

See our GitHub.io website for additional information and examples.

## Methods (by class)

- `fetch_metadata(Seurat)`: Seurat objects

- `fetch_metadata(SingleCellExperiment)`: SingleCellExperiment objects

- `fetch_metadata(AnnDataR6)`: AnnDataR6 objects

## Examples

``` r
# Return several metadata variables as a data.frame
fetch_metadata(
  AML_Seurat, 
  vars = c("condensed_cell_type", "Batch", "nCount_RNA")
  ) |> str()
#> 'data.frame':    250 obs. of  3 variables:
#>  $ condensed_cell_type: chr  "Plasma cells" "Plasma cells" "Plasma cells" "Plasma cells" ...
#>  $ Batch              : chr  "BM_200AB" "BM_200AB" "BM_200AB" "BM_200AB" ...
#>  $ nCount_RNA         : num  10863 8403 8100 8151 8828 ...
  
# Return data for a single metadata variable as a vector
fetch_metadata(
  AML_Seurat, 
  vars = "condensed_cell_type",
  return_class = "vector"
  ) |> str()
#>  Named chr [1:250] "Plasma cells" "Plasma cells" "Plasma cells" ...
#>  - attr(*, "names")= chr [1:250] "487013_1" "39207_1" "861619_1" "561110_1" ...

# Return all metadata 
fetch_metadata(
  AML_Seurat,
  full_table = TRUE
  ) |> str()
#> 'data.frame':    250 obs. of  29 variables:
#>  $ orig.ident         : chr  "SeuratProject" "SeuratProject" "SeuratProject" "SeuratProject" ...
#>  $ nCount_RNA         : num  10863 8403 8100 8151 8828 ...
#>  $ nFeature_RNA       : int  228 210 196 179 242 147 264 232 229 246 ...
#>  $ nCount_AB          : num  25709 31367 28166 14440 8203 ...
#>  $ nFeature_AB        : int  195 195 195 194 191 192 194 192 193 194 ...
#>  $ nCount_BOTH        : num  36572 39770 36266 22591 17031 ...
#>  $ nFeature_BOTH      : int  423 405 391 373 433 339 458 424 422 440 ...
#>  $ BOTH_snn_res.0.9   : chr  "17" "17" "17" "17" ...
#>  $ seurat_clusters    : Factor w/ 15 levels "0","1","2","3",..: 3 3 3 3 6 3 6 6 6 6 ...
#>  $ Prediction_Ind     : chr  "Plasma Cells" "Plasma Cells" "Plasma Cells" "Plasma Cells" ...
#>  $ BOTH_snn_res.1     : chr  "12" "26" "26" "26" ...
#>  $ ClusterID          : chr  "12" "26" "26" "26" ...
#>  $ Batch              : chr  "BM_200AB" "BM_200AB" "BM_200AB" "BM_200AB" ...
#>  $ x                  : num  -9.56 -9.53 -9.56 -9.53 -4.81 ...
#>  $ y                  : num  1.49 1.42 1.56 1.47 -2.44 ...
#>  $ x_mean             : num  -9.56 -9.51 -9.53 -9.48 -4.74 ...
#>  $ y_mean             : num  1.49 1.43 1.5 1.43 -2.49 ...
#>  $ cor                : num  0.852 0.856 0.878 0.855 0.923 ...
#>  $ ct                 : Factor w/ 38 levels "Plasma cells",..: 1 1 1 1 2 1 7 6 6 15 ...
#>  $ prop               : num  1 1 1 1 1 1 1 0.8 0.8 1 ...
#>  $ meandist           : num  0.045 0.0835 0.1062 0.1264 0.4902 ...
#>  $ cDC                : num  NaN NaN NaN NaN NaN ...
#>  $ B.cells            : num  NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
#>  $ Myelocytes         : num  NaN NaN NaN NaN NaN ...
#>  $ Erythroid          : num  NaN NaN NaN NaN 7.98 ...
#>  $ Megakaryocte       : num  NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
#>  $ Ident              : Factor w/ 28 levels "Plasma cells",..: 1 1 1 1 2 1 4 4 4 4 ...
#>  $ RNA_snn_res.0.4    : Factor w/ 15 levels "0","1","2","3",..: 3 3 3 3 6 3 6 6 6 6 ...
#>  $ condensed_cell_type: chr  "Plasma cells" "Plasma cells" "Plasma cells" "Plasma cells" ...
```
