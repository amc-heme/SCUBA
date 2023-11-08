test_that("unique_values returns the same results for Anndata R6, SingleCellExperiment and Seurat objects.
    ", {
      
  h5ad_meta_vars <- meta_varnames(AML_h5ad())
  
  #Run through and test every single variable in meta_varnames.
  #Test comparing all three object types to each other.
  #Cast to matrix required for h5ad objects required until test
  #object and functions have been updated.
  for(temp_meta_var in h5ad_meta_vars){
    expect_equal(
      object = unique_values(AML_Seurat, temp_meta_var),
      expected = unique_values(AML_SCE(), temp_meta_var),
      tolerance = 1e-6,
      ignore_attr = TRUE
    )
    expect_equal(
      object = as.matrix(unique_values(AML_Seurat, temp_meta_var)),
      expected = as.matrix(unique_values(AML_h5ad(), temp_meta_var)),
      tolerance = 1e-6,
      ignore_attr = TRUE
    )
    expect_equal(
      object = as.matrix(unique_values(AML_h5ad(), temp_meta_var)),
      expected = as.matrix(unique_values(AML_SCE(), temp_meta_var)),
      tolerance = 1e-6,
      ignore_attr = TRUE
    )
  }
})
