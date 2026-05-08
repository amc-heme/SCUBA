test_that("AML_h5ad() errors informatively when anndata is not installed", {
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) FALSE,
    .package = "SCUBA"
  )
  expect_error(
    AML_h5ad(),
    regexp = 'Package "anndata" must be installed to use this function\\. Install it with: install\\.packages\\("anndata"\\)',
    fixed = FALSE
  )
})

test_that("AML_h5ad_backed() errors informatively when anndata is not installed", {
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) FALSE,
    .package = "SCUBA"
  )
  expect_error(
    AML_h5ad_backed(),
    regexp = 'Package "anndata" must be installed to use this function\\. Install it with: install\\.packages\\("anndata"\\)',
    fixed = FALSE
  )
})
