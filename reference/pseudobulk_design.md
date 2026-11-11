# Preview the pseudobulk profiles of an object

Metadata-only companion of
[`pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.md):
resolves the same sample x group profiles (same rows, same order) and
reports how many cells each one has, without reading any expression
data. Use it to check replication before aggregating, e.g. whether every
group has enough samples, or whether samples contribute cells to several
groups (a paired design).

## Usage

``` r
pseudobulk_design(object, sample_by, group_by = NULL, cells = NULL)
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

- cells:

  A character vector of cell names to include, as they are named in the
  object (i.e. according to `colnames(object)`). If `NULL`, data will be
  returned for all cells in the object.

## Value

An object of class `scuba_pseudobulk_design`: a list with `profiles` (as
in
[`pseudobulk()`](https://amc-heme.github.io/SCUBA/reference/pseudobulk.md)),
`table` (samples x groups matrix of cell counts; one column when
`group_by` is `NULL`), `paired` (`TRUE` when any sample has cells in
more than one group), `samples_per_group` (samples with at least one
cell per group), `samples_in_all_groups`, `n_cells_dropped`, `n_cells`
(cells assigned to a profile), `sample_by` and `group_by`.

## Examples

``` r
d <- pseudobulk_design(AML_Seurat, "Batch", "condensed_cell_type")
d$table
#>            B Cells BM Monocytes CD4+ T Cells CD8+ T Cells Dendritic cells
#> BM_200AB        14           25            8           17              15
#> PBMC_200AB      11            0           17            8              10
#>            NK Cells Plasma cells Plasmacytoid dendritic cells Primitive
#> BM_200AB         18           12                           19        24
#> PBMC_200AB        7           13                            6         1
#>            PBMC Monocytes
#> BM_200AB                0
#> PBMC_200AB             25
d$paired
#> [1] TRUE
```
