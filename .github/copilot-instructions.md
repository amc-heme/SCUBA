# SCUBA Copilot Instructions

Purpose: Enable AI agents to work productively in the SCUBA R package (unified single-cell data access API for Seurat, SingleCellExperiment, and anndata objects).

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
