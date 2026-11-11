# Classification of requested names for engine routing. Pure R: no Julia needed.

test_that("feature, metadata and reduction names are told apart", {
  kinds <- SCUBA:::.scuba_var_kinds(
    AML_Seurat, c("rna_GAPDH", "GAPDH", "PC_1", "UMAP_2", "nCount_RNA", "md_nCount_RNA", "PC_X")
  )
  expect_identical(kinds$is_meta, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE))
  # a reduction coordinate needs a trailing dimension number
  expect_identical(kinds$is_red, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))

  expect_true(SCUBA:::.scuba_feature_only(AML_Seurat, c("rna_GAPDH", "ab_CD117-AB")))
  expect_false(SCUBA:::.scuba_feature_only(AML_Seurat, "PC_1"))
  expect_false(SCUBA:::.scuba_feature_only(AML_Seurat, c("rna_GAPDH", "md_nCount_RNA")))
  expect_false(SCUBA:::.scuba_feature_only(AML_Seurat, "nCount_RNA"))

  # explicit metadata names and reduction keys (as used for file-backed objects)
  expect_false(SCUBA:::.scuba_feature_only(
    NULL, "X_umap_1", meta = "nCount_RNA", reductions = c("X_pca_", "X_umap_")
  ))
  expect_true(SCUBA:::.scuba_feature_only(
    NULL, "X_GAPDH", meta = "nCount_RNA", reductions = c("X_pca_", "X_umap_")
  ))
})
