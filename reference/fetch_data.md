# Fetch feature expression data, reduction coordinates, or metadata from single-cell objects

This function extends the behavior of SeuratObject's FetchData to other
single-cell objects, allowing for expression data, metdadata, or
reduction coordinates to be pulled using consistent syntax.

## Usage

``` r
fetch_data(object, vars = NULL, layer = NULL, cells = NULL, ...)

# S3 method for class 'Seurat'
fetch_data(object, vars = NULL, layer = NULL, cells = NULL, slot = NULL, ...)

# S3 method for class 'SingleCellExperiment'
fetch_data(object, vars, layer = NULL, cells = NULL, ...)

# S3 method for class 'AnnDataR6'
fetch_data(object, vars, layer = NULL, cells = NULL, ...)
```

## Arguments

- object:

  A single-cell object. Currently, Seurat, SingleCellExperiment, and
  anndata objects are supported. If a Seruat object is passed to this
  generic, the FetchData method from the SeruatObject method will be
  ran. For all other objects, methods that replicate the behavior of
  FetchData in that object will be ran.

- vars:

  A character vector with desired features, metadata variables, or
  reduction dimensions to pull from the object. By default, features are
  returned from the default assay (or "experiment" in
  SingleCellExperiment objects, or "modality" in anndata objects). To
  pull feature expression data, the assay to pull data from should be
  defined using the "key" of the assay before the feature name. To
  determine the key that corresponds to the assay to pull data from, run
  `all_keys`. For more information, see [our user
  guide](https://amc-heme.github.io/SCUBA/reference/our%20user%20guide).

- layer:

  For feature expression data, the layer to pull data from. Layers are
  referred to as "slots" in Seurat objects v4 and earlier, and "assays"
  in SingleCellExperiment objects.

- cells:

  A character vector of cells to include, as they are named in the
  object (i.e. according to colNames(object)). If `NULL`, data will be
  returned for all cells in the object.

- ...:

  Additional parameters, beyond the ones listed above, to be passed to
  S3 methods. This includes the following, all of which are documented
  in the parameter entries below:

  - fetch_data.Seurat: `slot` parameter

- slot:

  This parameter is added for backwards compatability with Seruat v4 and
  earlier. This was deprecated in Seurat version 5.0.0. *If you have
  Seruat v5.0.0 or later, this should not be used. This parameter should
  also not be used for SingleCellExperiment objects or anndata objects
  and will not work at all for these object classes. Use* `layer`
  *instead.*

## Value

A data.frame with the requested `vars` as columns and the cells as rows.

## Details

See our GitHub.io website for additional information and examples.

## Methods (by class)

- `fetch_data(Seurat)`: Seurat objects. This will run FetchData from the
  SeuratObject package.

- `fetch_data(SingleCellExperiment)`: SingleCellExperiment objects

- `fetch_data(AnnDataR6)`: anndata Objects

## Examples

``` r
# Feature expression data
fetch_data(AML_Seurat, vars = "rna_FLT3") |> str()
#> 'data.frame':    250 obs. of  1 variable:
#>  $ rna_FLT3: num  0 0 0 0 0 ...

# Reduction coordinates
fetch_data(AML_Seurat, vars = c("UMAP_1", "UMAP_2")) |> str()
#> 'data.frame':    250 obs. of  2 variables:
#>  $ UMAP_1: num  -1.64 -1.5 -1.45 -1.38 -1.41 ...
#>  $ UMAP_2: num  9.9 10.13 10.21 10.51 3.39 ...
fetch_data(AML_Seurat, vars = c("PC_1", "PC_2", "PC_3")) |> str()
#> 'data.frame':    250 obs. of  3 variables:
#>  $ PC_1: num  3.91 3.67 2.97 2.88 42.94 ...
#>  $ PC_2: num  1.676 -0.157 0.716 1.103 -11.258 ...
#>  $ PC_3: num  1.86 3.03 2.6 3.68 -4 ...

# Metadata
fetch_data(AML_Seurat, vars = c("condensed_cell_type", "Batch", "nCount_RNA")) |> str()
#> 'data.frame':    250 obs. of  3 variables:
#>  $ condensed_cell_type: chr  "Plasma cells" "Plasma cells" "Plasma cells" "Plasma cells" ...
#>  $ Batch              : chr  "BM_200AB" "BM_200AB" "BM_200AB" "BM_200AB" ...
#>  $ nCount_RNA         : num  10863 8403 8100 8151 8828 ...
```
