# Modalities and reductions of a file-backed SCUBA object

Modalities and reductions of a file-backed SCUBA object

## Usage

``` r
scuba_source_modalities(object)

scuba_source_reductions(object)
```

## Arguments

- object:

  A `scuba_source` object from
  [`open_h5ad()`](https://amc-heme.github.io/SCUBA/reference/open_h5ad.md).

## Value

`scuba_source_modalities()`: a character vector of modality keys (main
first). `scuba_source_reductions()`: a named integer vector of reduction
keys and their dimension counts.
