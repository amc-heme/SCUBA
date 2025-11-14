# Get names of all features in an assay/experiment/modality

Returns the names of all features in a modality. This utility function
can be used in several applications:

- In Shiny apps, return available features for passage to a selection
  menu

- Before using
  [`fetch_data()`](https://amc-heme.github.io/SCUBA/reference/fetch_data.md),
  generate a list of features in the object for searching, or test if a
  feature is present before requesting data

## Usage

``` r
features_in_assay(object, assay)

# S3 method for class 'Seurat'
features_in_assay(object, assay)

# S3 method for class 'SingleCellExperiment'
features_in_assay(object, assay)

# S3 method for class 'AnnDataR6'
features_in_assay(object, assay)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

- assay:

  the name of an assay/modality for which to view features.

## Methods (by class)

- `features_in_assay(Seurat)`: Seurat objects

- `features_in_assay(SingleCellExperiment)`: SingleCellExperiment
  objects

- `features_in_assay(AnnDataR6)`: Anndata objects

## Examples

``` r
features_in_assay(AML_Seurat, assay = "RNA") |> str()
#>  chr [1:462] "ACTG1" "ADGRG1" "AHSP" "AIF1" "ANK1" "ANKRD28" "ANLN" ...

# Check if a feature is present in an assay
"MEIS1" %in% features_in_assay(AML_Seurat, assay = "RNA")
#> [1] TRUE
```
