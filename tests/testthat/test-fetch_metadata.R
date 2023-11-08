test_that("fetch_metadata returns the same results for Anndata R6, SingleCellExperiment and Seurat objects.
    ", {
  seurat_metadata <-
    fetch_metadata(
      AML_Seurat, 
      full_table=T
#      vars = "condensed_cell_type",
#      cells = get_all_cells(AML_Seurat)
    )
  sce_metadata <-
    fetch_metadata(
      AML_SCE(), 
      full_table = TRUE
    )
    
  h5ad_metadata <-
    fetch_metadata(
      AML_h5ad(), 
      full_table = TRUE
    )
  
  seurat_metadata$orig.ident <- NULL
  sce_metadata$orig.ident <- NULL
  sce_metadata$ident <- NULL
  #Check rowSums and colSums of data
  expect_equal(
    object = as.matrix(h5ad_metadata),
    expected = as.matrix(seurat_metadata),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(
    object = as.data.frame(seurat_metadata),
    expected = as.data.frame(sce_metadata),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})
