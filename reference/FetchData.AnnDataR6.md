# FetchData Equivalent for Anndata Objects

The anndata equivalent of the SeuratObject FetchData method.

## Usage

``` r
# S3 method for class 'AnnDataR6'
FetchData(object, vars, layer = NULL, cells = NULL, ...)
```

## Arguments

- object:

  An anndata object.

- vars:

  A character vector with desired features or metadata variables to pull
  from the object. Any combination of entries in the genes matrix (X),
  metadata (obs), or obsm matrices can be provided here. To include a
  feature from layers, use the `layers` parameter. It is greatly
  preferred to specify the matrix a variable is in with an underscore.
  for example, to pull the FIS1 gene from the genes matrix (X), specify
  "X_FIS1" instead of "FIS1". To pull metadata, use "obs\_", and to pull
  data from a matrix in obsm, use the name of that matrix, not obsm. For
  example, if a matrix in obsm is named "protein", use "protein\_" to
  pull data from that matrix. Variables that are entered without a key
  can still be found, as long as there is only one matrix in the object
  with that variable name. Variables that do not have a valid key (X\_,
  obs\_, and a key from obj.obsm_names()) will be ignored, as will
  duplicate variables.

- layer:

  The layer to pull data from. If unspecified, the feature will be
  pulled from the X matrix. To view a list of available layers, run
  `object$layers`. Layers will not work with alternate modalities stored
  in `obsm`, only the main modality (the modality in the X matrix of the
  object).

- cells:

  A character vector of cells to include, as they are named in the
  object (i.e. according to colNames(object)). If `NULL`, data will be
  returned for all cells in the object.

- ...:

  parameter provided for consistency with S3 generic/methods

## Value

A data.frame object containing the requested expression data or
metadata.
