test_that("fetch_cells returns the same results for Anndata R6, SingleCellExperiment and Seurat objects.
    ", {
  # Subset for cells matching a combination of levels within a variable
  seurat_cells <-
    fetch_cells(
      AML_Seurat, 
      meta_var = "seurat_clusters",
      meta_levels = unique_values(AML_Seurat, var="seurat_clusters")
    )
  sce_cells <-
    fetch_cells(
      AML_SCE(),  
      meta_var = "seurat_clusters",
      meta_levels = unique_values(AML_SCE(), var="seurat_clusters")
    )
    
  h5ad_cells <-
    fetch_cells(
      AML_h5ad(),  
      meta_var = "seurat_clusters",
      meta_levels = unique_values(AML_h5ad(), var="seurat_clusters")
    )
  
  h5mu_cells <-
    fetch_cells(
      AML_h5mu(),  
      meta_var = "seurat_clusters",
      meta_levels = unique_values(AML_h5ad(), var="seurat_clusters")
    )
  
  #Check rowSums and colSums of data
  expect_equal(
    object = h5ad_cells,
    expected = seurat_cells,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(
    object = seurat_cells,
    expected = sce_cells,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(
    object = h5mu_cells,
    expected = seurat_cells,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  
  #Null value for meta_levels should return error
  expect_error(
    object = 
      fetch_cells(
        AML_Seurat,
        meta_var = "seurat_clusters",
        meta_levels = NULL
      )
  )
  
  #Only categorical values permitted for meta_levels
  expect_error(
    object = 
      fetch_cells(
        AML_Seurat,
        meta_var = "nCount_RNA",
        meta_levels = 
          unique_values(
            AML_Seurat, 
            var = "nCount_RNA"
            )
        )
    )
})
