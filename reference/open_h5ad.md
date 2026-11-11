# Open an .h5ad or .h5mu file as a SCUBA object (Julia engine)

Reads the schema and the cell metadata of an AnnData `.h5ad` file (or a
MuData `.h5mu` file with `open_h5mu()`) and returns a `scuba_source`
object that every SCUBA accessor accepts. Expression values and
reductions are read from the file on demand by the Julia engine
(SCUBA.jl), so no Python environment is needed. Requires
[`scuba_engine()`](https://amc-heme.github.io/SCUBA/reference/scuba_engine.md)
to be `"julia"`.

## Usage

``` r
open_h5ad(path)

open_h5mu(path)
```

## Arguments

- path:

  Path to an `.h5ad` (`open_h5ad()`) or `.h5mu` (`open_h5mu()`) file.

## Value

An object of class `scuba_source`.

## Details

Keys follow SCUBA.jl's rules. For `.h5ad`: the main matrix is modality
`"X"` (features addressed bare or as `X_<gene>`, layers `"X"` and
`layers/*`), each `obsm` data frame is a modality named after it (e.g.
`protein_CD117-AB`), each numeric `obsm` array is a reduction
(`X_umap_1`), and metadata columns may be given bare or as
`obs_<column>`. For `.h5mu`: modalities are named after the MuData
modalities (`rna_GAPDH`), reductions as `<modality>_<key>`
(`rna_X_umap_1`, bare `X_umap_1` when unique), metadata is the global
`obs` plus each modality's `obs` as `<modality>_<column>`; cells absent
from a modality are `NaN` / `NA`.

## See also

[`update_object_metadata()`](https://amc-heme.github.io/SCUBA/reference/update_object_metadata.md),
[`scuba_source_modalities()`](https://amc-heme.github.io/SCUBA/reference/scuba_source_modalities.md)
