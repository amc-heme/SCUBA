test_that("meta_varnames returns the same results for Anndata R6, SingleCellExperiment and Seurat objects.
    ", {
      
  seurat_meta_vars <- meta_varnames(AML_Seurat)
  sce_meta_vars <- meta_varnames(AML_SCE())
  h5ad_meta_vars <- meta_varnames(AML_h5ad())
  
  #These delete orig.ident (AML_Seurat+AML_SCE()) and ident (AML_SCE()) 
  #Once the test object are uniform this should be removed.
  seurat_meta_vars <- seurat_meta_vars[-1]
  sce_meta_vars <- sce_meta_vars[-1]
  sce_meta_vars <- sce_meta_vars[-29]
  
  #Check rowSums and colSums of data
  expect_equal(
    object = seurat_meta_vars,
    expected = sce_meta_vars,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(
    object = seurat_meta_vars,
    expected = h5ad_meta_vars,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})
