test_that("fetch_reduction returns the same results for Anndata R6, SingleCellExperiment and Seurat objects.
    ", {
  seurat_reduction <-
    fetch_reduction(
      AML_Seurat, 
      reduction = default_reduction(AML_Seurat),
      cells = get_all_cells(AML_Seurat)
    )
  sce_reduction <-
    fetch_reduction(
      AML_SCE(), 
      reduction = default_reduction(AML_SCE()),
      cells = get_all_cells(AML_SCE())
    )
    
  h5ad_reduction <-
    fetch_reduction(
      AML_h5ad(), 
      reduction = default_reduction(AML_h5ad()),
      cells = get_all_cells(AML_h5ad())
      
    )
  
  #Check rowSums and colSums of data
  expect_equal(
    object = h5ad_reduction,
    expected = seurat_reduction,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(
    object = seurat_reduction,
    expected = sce_reduction,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(
    object = h5ad_reduction,
    expected = sce_reduction,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})
