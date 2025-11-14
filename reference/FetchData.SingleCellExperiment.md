# FetchData Equivalent for SingleCellExperiment Objects

The SingleCellExperiment equivalent of the SeuratObject FetchData
method.

## Usage

``` r
# S3 method for class 'SingleCellExperiment'
FetchData(object, vars, layer = NULL, cells = NULL, ...)
```

## Arguments

- object:

  A SingleCellExperiment object.

- vars:

  A character vector with desired features or metadata variables to pull
  from the object. To include features from an experiment other than the
  main experiment, use the name of the experiment as a prefix (i.e.
  AB_CD4 for a feature in the experiment "AB" named "CD4".)

- layer:

  The assay (equivalent of layer in Seruat objects (slot in v4 and
  earlier)) to pull data from. To view the list of available assays in
  your object, use `assayNames({your object})`.

- cells:

  A character vector of cells to include, as they are named in the
  object (i.e. according to colNames(object)). If `NULL`, data will be
  returned for all cells in the object.

- ...:

  parameter provided for consistency with S3 generic/methods

## Value

A data.frame object containing the requested expression data or
metadata.
