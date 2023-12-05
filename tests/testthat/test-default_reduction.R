suppressMessages(library(SCUBA, quietly = TRUE))

# Load Seurat Object
sobj <- AML_Seurat

# Load SCE Object (with DelayedArray integration)
sce <- AML_SCE()

test_that("default_reduction works", {
  # Test if UMAP is returned (objects have UMAP and PCA reductions)
  expect_equal(
    object = default_reduction(sobj),
    expected = "umap"
  )

  expect_equal(
    object = default_reduction(sce),
    expected = "UMAP"
  )

  # Delete the UMAP assay from both objects and test if PCA is returned
  reducedDims(sce)[["UMAP"]] <- NULL
  sobj@reductions$umap <- NULL

  expect_equal(
    object = default_reduction(sobj),
    expected = "pca"
  )

  expect_equal(
    object = default_reduction(sce),
    expected = "PCA"
  )

  # Check if an error is returned if no reductions exist
  #reducedDims(sce)[["PCA"]] <- NULL
  #sobj@reductions$pca <- NULL

  expect_error(
    object = default_reduction(sobj)
  )

  expect_error(
    object = default_reduction(sce)
  )
})
