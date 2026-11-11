# Pseudobulk profiles: aggregate a counts layer over groups of cells

Sums (or averages) every feature of a layer over *profiles*, the
combinations of a sample identifier and, optionally, a group label (cell
type, condition, cluster). The result is the small features x profiles
matrix that pseudobulk differential expression methods (edgeR, DESeq2)
consume, together with a table describing each profile.

## Usage

``` r
pseudobulk(
  object,
  sample_by,
  group_by = NULL,
  ...,
  assay = NULL,
  layer = "counts",
  cells = NULL,
  features = NULL,
  fun = c("sum", "mean")
)

# S3 method for class 'Seurat'
pseudobulk(
  object,
  sample_by,
  group_by = NULL,
  ...,
  assay = NULL,
  layer = "counts",
  cells = NULL,
  features = NULL,
  fun = c("sum", "mean")
)

# S3 method for class 'SingleCellExperiment'
pseudobulk(
  object,
  sample_by,
  group_by = NULL,
  ...,
  assay = NULL,
  layer = "counts",
  cells = NULL,
  features = NULL,
  fun = c("sum", "mean")
)

# S3 method for class 'scuba_source'
pseudobulk(
  object,
  sample_by,
  group_by = NULL,
  ...,
  assay = NULL,
  layer = "counts",
  cells = NULL,
  features = NULL,
  fun = c("sum", "mean")
)

# S3 method for class 'AnnDataR6'
pseudobulk(
  object,
  sample_by,
  group_by = NULL,
  ...,
  assay = NULL,
  layer = "counts",
  cells = NULL,
  features = NULL,
  fun = c("sum", "mean")
)

# S3 method for class 'anndata.AnnData'
pseudobulk(
  object,
  sample_by,
  group_by = NULL,
  ...,
  assay = NULL,
  layer = "counts",
  cells = NULL,
  features = NULL,
  fun = c("sum", "mean")
)

# S3 method for class 'anndata._core.anndata.AnnData'
pseudobulk(
  object,
  sample_by,
  group_by = NULL,
  ...,
  assay = NULL,
  layer = "counts",
  cells = NULL,
  features = NULL,
  fun = c("sum", "mean")
)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExperiment, and
  anndata objects are supported.

- sample_by:

  The sample (biological replicate) of each cell: the name of a metadata
  column, or a character/factor vector named by cell (one label per cell
  of the object, or at least per cell in `cells`). A vector lets a
  caller aggregate by labels that are not stored in the object.

- group_by:

  Optional group of each cell (same forms as `sample_by`). With `NULL`
  the profiles are the samples themselves.

- ...:

  Passed to methods.

- assay:

  Assay / modality to aggregate (default: the default assay of the
  object; the main modality `"X"` for AnnData and `.h5ad` objects).

- layer:

  Layer to aggregate. Defaults to `"counts"`, *not* to
  [`default_layer()`](https://amc-heme.github.io/SCUBA/reference/default_layer.md):
  the default layer holds normalised values, which must not be summed.
  AnnData objects whose raw counts live in `X` need `layer = "X"` (or
  `"raw"` for `adata.raw`). An AnnData `obsm` modality (for example a
  protein panel) is a single matrix with no layers: pass `layer = "X"`;
  any other layer is an error rather than a silent substitution.

- cells:

  A character vector of cell names to include, as they are named in the
  object (i.e. according to `colnames(object)`). If `NULL`, data will be
  returned for all cells in the object.

- features:

  Features to include (default: every feature of the assay).

- fun:

  `"sum"` (default) or `"mean"`; the mean divides each profile's sum by
  the number of cells that contributed to it.

## Value

An object of class `scuba_pseudobulk`: a list with

- `counts`: features x profiles numeric matrix (sums, or means with
  `fun = "mean"`),

- `profiles`: data.frame with one row per profile (`profile_id`,
  `sample`, `group`, `n_cells`), in the column order of `counts`;
  `group` is `NA` when `group_by` is `NULL`,

- `assay`, `layer`, `fun`, `sample_by`, `group_by` (the column names, or
  `NA` for vectors),

- `n_cells_dropped`: cells left out for a missing sample or group label,

- `engine`: `"julia"`, `"r"` or `"python"`.

## Details

The full features x cells matrix is never built: on the Julia engine
SCUBA.jl streams the layer into a group accumulator, and on the R engine
each matrix class is aggregated with its own block-wise primitive
(sparse indicator product for `dgCMatrix` and dense layers,
[`BPCells::pseudobulk_matrix()`](https://bnprks.github.io/BPCells/reference/pseudobulk_matrix.html)
for BPCells layers,
[`DelayedArray::colsum()`](https://rdrr.io/pkg/S4Arrays/man/rowsum.html)
for HDF5-backed layers). AnnData objects (reticulate) are aggregated in
cell chunks by a Python helper, in memory or in backed mode.

Cells whose sample or group value is `NA` or blank are left out (counted
in `n_cells_dropped`). Profiles are ordered sample-major, then by group;
levels follow the factor levels of the metadata column when it is a
factor, and sorted order (C locale) otherwise, so every engine and
backend returns the same columns.

## See also

[`pseudobulk_design()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk_design.md)
for the metadata-only preview of the profiles,
[`layer_names()`](https://amc-heme.github.io/SCUBA/reference/layer_names.md)
to find a counts layer,
[`fetch_matrix()`](https://amc-heme.github.io/SCUBA/reference/fetch_matrix.md).

## Examples

``` r
pb <- pseudobulk(AML_Seurat, sample_by = "Batch", group_by = "condensed_cell_type")
pb
#> SCUBA pseudobulk: 462 features x 18 profiles (sum of layer 'counts', assay 'RNA', r engine)
#>   samples: 2; groups: 10; cells: 250 (dropped 0)
head(pb$profiles)
#>                                          profile_id   sample           group
#> BM_200AB__B Cells                 BM_200AB__B Cells BM_200AB         B Cells
#> BM_200AB__BM Monocytes       BM_200AB__BM Monocytes BM_200AB    BM Monocytes
#> BM_200AB__CD4+ T Cells       BM_200AB__CD4+ T Cells BM_200AB    CD4+ T Cells
#> BM_200AB__CD8+ T Cells       BM_200AB__CD8+ T Cells BM_200AB    CD8+ T Cells
#> BM_200AB__Dendritic cells BM_200AB__Dendritic cells BM_200AB Dendritic cells
#> BM_200AB__NK Cells               BM_200AB__NK Cells BM_200AB        NK Cells
#>                           n_cells
#> BM_200AB__B Cells              14
#> BM_200AB__BM Monocytes         25
#> BM_200AB__CD4+ T Cells          8
#> BM_200AB__CD8+ T Cells         17
#> BM_200AB__Dendritic cells      15
#> BM_200AB__NK Cells             18
pb$counts[1:3, 1:3]
#>        BM_200AB__B Cells BM_200AB__BM Monocytes BM_200AB__CD4+ T Cells
#> ACTG1                 70                    548                     57
#> ADGRG1                 1                      2                      0
#> AHSP                   0                      0                      0

# Samples only, as means
pseudobulk(AML_SCE(), sample_by = "Batch", fun = "mean")$profiles
#>            profile_id     sample group n_cells
#> BM_200AB     BM_200AB   BM_200AB  <NA>     152
#> PBMC_200AB PBMC_200AB PBMC_200AB  <NA>      98
```
