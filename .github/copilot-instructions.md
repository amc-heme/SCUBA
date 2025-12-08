# SCUBA Copilot Instructions

Purpose: Enable AI agents to work productively in the SCUBA R package (unified single-cell data access API for Seurat, SingleCellExperiment, and anndata objects).

Please view the conventions in `.github/agent_conventions.md` for coding best practices before proceeding.

## Core Architecture & Data Flow
- All user‑facing access goes through S3 generics named `fetch_*`, `default_*`, `plot_*`, and small utilities (`features_in_assay`, `all_keys`, `meta_varnames`, etc.).
- Pattern: each generic defined as `fn <- function(object, ...) UseMethod("fn")`; class methods: `.Seurat`, `.SingleCellExperiment`, `.AnnDataR6`; `.default` issues warning.
- `fetch_data()` is the central API: returns a dense `data.frame` with cells as rows and requested variables as columns independent of storage format.
- Feature/reduction disambiguation relies on "keys" prepended to variable names (e.g. `RNA_FLT3`, `UMAP_1`); SCE alternate experiments and reductions, and AnnData modalities are unified via these keys.
- Large queries (≥1000 features) intentionally warn; do not suppress unless explicitly asked.
- AnnData access (`fetch_data.AnnDataR6`) delegates to Python `inst/extdata/Python/fetch_anndata.py` via `reticulate::py_run_file` and `py_require()`; Python deps are auto‑managed (anndata, pandas, numpy, scipy).

## Object Class Support & Extensibility
- Adding a new object class: implement methods for core generics: `fetch_data`, `fetch_cells`, `fetch_metadata`, `fetch_reduction`, `default_layer`, `default_reduction`, plotting helpers; add a `.default` fallback only in the generic file.
- Update shared param docs in `R/roxygen_common_sections.R` (e.g. `object_param`) to reflect new class; then run `devtools::document()`.
- Ensure naming/keying semantics remain consistent: features must be uniquely addressable by optional key prefix; reductions use dimension suffix (`UMAP_1`).

## Conventions & Naming
- Retrieval: `fetch_*` (returns vectors/data.frames), metadata utilities (`meta_varnames`, `unique_values`, `all_keys`).
- Defaults: `default_layer()` (preferred) vs deprecated `default_slot()` (file `defult_slot.R` intentionally left for backward compatibility; do not rename typo file unless doing a deprecation cleanup).
- Plotting: `plot_*` functions accept unified data structures from `fetch_data` and default reduction/layer helpers.
- Tests mirror function names: `tests/testthat/test-fetch_data.R`, `test-fetch_cells.R`, etc.; add analogous tests when extending.
- Warnings purposefully informative (duplicate features, ambiguous metadata vs feature names, missing vars across multiple experiments); preserve logic and messaging.

## Feature Naming & Key Prefixes

**Overview**: The `vars` parameter of `fetch_data()` accepts feature names (for expression data) and reduction dimensions (for coordinates). Both use key prefixes to disambiguate their source, but the keying conventions differ.

### Feature Keys (for expression data)

Feature keys identify which assay/modality/experiment a feature belongs to. The syntax is `<assay_key>_<feature_name>`.

**By object class**:
- **Seurat**: Keys retrieved via `SCUBA::all_keys()` (wrapper for `SeuratObject::Key()`; e.g., `rna_`, `adt_`). Example: `rna_FLT3`, `adt_CD19`.
- **SingleCellExperiment**: Keys are experiment names—`mainExpName(object)` for the main experiment, plus `altExpNames(object)` for alternate experiments. Example: `RNA_FLT3`, `protein_CD9`.
- **anndata**: The main matrix uses key `X`. Other modalities stored in `obsm` use their `obsm` key name. Example: `X_GAPDH`, `protein_CD9` (if protein modality is keyed as `protein` in `obsm`).

**Default behavior** (feature requested without key):
- SCUBA searches the default layer of the main/primary experiment (or `X` for anndata).
- If the feature exists in metadata *and* as a feature, metadata takes precedence and a warning directs the user to add the key prefix.
- If a feature exists in multiple alternate experiments, the query fails with a warning requesting the user to specify the experiment key.

**Examples**:
```r
# Fetch feature with explicit assay key
fetch_data(obj, vars = "rna_FLT3")
fetch_data(obj, vars = c("adt_CD9", "adt_CD19"))

# Fetch from default assay (no key needed if unambiguous)
fetch_data(obj, vars = "GAPDH")

# anndata: fetch from main matrix X
fetch_data(adata, vars = "X_GAPDH")  
fetch_data(adata, vars = "GAPDH")    # also works for X

# anndata: fetch from obsm modality
fetch_data(adata, vars = "protein_CD9")  # 'protein' is the obsm key
```

### Reduction Keys (for coordinates)

Reduction keys identify which dimensionality reduction to pull coordinates from. The syntax is `<reduction_name>_<dimension_number>`. This pattern is consistent across all object classes.

**By object class**:
- **Seurat**: Reduction names from the object (e.g., `umap`, `pca`). Example: `UMAP_1`, `PC_1`.
- **SingleCellExperiment**: Reduction names from `reducedDimNames(object)`. Example: `UMAP_1`, `PCA_2`.
- **anndata**: Reduction keys from `obsm_keys()`, typically prefixed with `X_`. Example: `X_umap_1`, `X_pca_2`.

**Examples**:
```r
# Seurat / SCE reductions
fetch_data(obj, vars = c("UMAP_1", "UMAP_2"))
fetch_data(obj, vars = c("PC_1", "PC_2", "PC_3"))

# anndata reductions
fetch_data(adata, vars = c("X_umap_1", "X_umap_2"))
```

### Discovering Keys

Use `SCUBA::all_keys(object)` to retrieve all available keys (both assay and reduction keys together). This function is a SCUBA wrapper that works uniformly across Seurat, SingleCellExperiment, and anndata objects. Users are expected to know which keys correspond to assays vs. reductions based on their object structure. **Note**: In internal SCUBA code, always prefer `SCUBA::all_keys()` over `SeuratObject::Key()` for consistency across object classes.

## Common Workflows (Local Dev)
- Install (development): `devtools::load_all()`; regenerate docs: `devtools::document()`; run tests: `devtools::test()`; full check: `devtools::check()`.
- Example feature query (Seurat): `fetch_data(AML_Seurat, vars = c("rna_FLT3", "UMAP_1", "condensed_cell_type"))`.
- Feature listing for UI scaffolding: `features <- features_in_assay(AML_Seurat, assay = "RNA")`.
- Subsetting cells by metadata: `cells <- fetch_cells(AML_Seurat, meta_var = "Batch", meta_levels = c("B1","B2"))` then pass `cells` to other `fetch_*` functions.

## Python Integration Notes
- Do not manually manage Python env in code paths; rely on `py_require()` calls inside AnnData methods.
- Any new AnnData functionality should follow pattern: source script under `inst/extdata/Python/` and call via `reticulate::py_run_file`; avoid embedding long Python strings inline.

## Performance & Safety
- All matrix extraction coerces to dense `matrix` then transposes—be cautious adding very large feature batches (>1000) as memory may spike.
- Maintain feature name integrity after transposition (`colnames` set explicitly before conversion to list/data.frame).
- Avoid altering order semantics: after retrieval, columns are re-ordered to match user `vars` input using `pmatch`; preserve this behavior.

## Adding Tests
- Test files: `tests/testthat/test-<function>.R`; follow existing assertion style (presence/absence, warning expectation, dimensionality checks).
- For new warning conditions, use `expect_warning()` with partial message matching (keep messages stable & user-helpful).

## Documentation
- Always add roxygen documentation for new R functions and methods.
- Use `R/roxygen_common_sections.R` for shared parameter documentation (e.g., `object_param`, `features_param`, `layer_param`, `cells_param`).
- Use `@inheritParams` to reference shared parameter docs instead of duplicating descriptions.
- Run `devtools::document()` after adding or modifying roxygen comments.

## Deprecations & Backward Compatibility
- `default_slot()` deprecated (warns, delegates to `default_layer()`); keep until version threshold specified in file.
- Backward Seurat v4 support uses `slot` arg; do not remove without updating README & docs.

## When Implementing Changes
- Prefer minimal, targeted edits: add method files under `R/` matching existing naming style (snake_case).
- Update documentation via roxygen then run tests before proposing changes.
- If touching key resolution logic (e.g. in `fetch_data.SingleCellExperiment`), replicate existing warning flows for duplicates/ambiguities.

## Quick Validation Commands
```r
# Load dev version
devtools::load_all()
# Regenerate docs
devtools::document()
# Run tests
devtools::test()
# Example smoke test
fetch_data(AML_Seurat, vars = c("rna_FLT3", "UMAP_1")) |> head()
```

Feedback welcome: Please indicate any unclear sections (e.g. key naming, adding new class methods) for refinement.
