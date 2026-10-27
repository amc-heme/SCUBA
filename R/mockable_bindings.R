# Placeholder bindings so tests can mock base functions
#
# testthat::local_mocked_bindings() can only replace bindings that exist in
# the package namespace, and `requireNamespace()` lives in base. The tests
# for the "package not installed" errors (AML_h5ad(), fetch_data() on AnnData
# objects) mock it with `.package = "SCUBA"`, which needs this binding.
# R skips non-function bindings when looking up a function to call, so every
# `requireNamespace(...)` call in the package still reaches the base
# function. See ?testthat::local_mocked_bindings, "Base functions".
requireNamespace <- NULL
