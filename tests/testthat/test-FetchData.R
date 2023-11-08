test_that("Fetchdata.SingleCellExperiment, Fetchdata.SeuratObject and FetchData AnnData all return the same data.", {
  seurat_data <-
    FetchData(
      AML_Seurat,
      slot = "data",
      vars =
        c("ab_CD117-AB",
          "ab_CD123-AB",
          "ab_CD11c-AB",
          "rna_GAPDH",
          "rna_MEIS1",
          # Reductions
          "UMAP_1",
          "UMAP_2",
          # Metadata
          "nCount_RNA",
          "nFeature_RNA",
          # "Ambiguous" feature not in RNA assay
          "CD11a-AB"
        )
    )

  sce_data <-
    FetchData(
      AML_SCE(),
      slot = "logcounts",
      vars =
        c("AB_CD117-AB",
          "AB_CD123-AB",
          "AB_CD11c-AB",
          "RNA_GAPDH",
          "RNA_MEIS1",
          # Reductions
          "UMAP_1",
          "UMAP_2",
          # Metadata
          "nCount_RNA",
          "nFeature_RNA",
          # "Ambiguous" feature not in RNA assay
          "CD11a-AB"
        )
    )
  
  h5ad_data <-  
    FetchData(
      AML_h5ad(),
      vars =
        c("protein_CD117-AB",
          "protein_CD123-AB",
          "protein_CD11c-AB",
          "X_GAPDH",
          "X_MEIS1",
          # Reductions
          "X_umap_1",
          "X_umap_2",
          # Metadata
          "nCount_RNA",
          "nFeature_RNA",
          # "Ambiguous" feature not in RNA assay
          "CD11a-AB"
      )
  )
  
  # Check rowSums and colSums of data
  expect_equal(
    object = rowSums(sce_data),
    expected = rowSums(seurat_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  expect_equal(
    object = colSums(sce_data),
    expected = colSums(seurat_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  
  expect_equal(
    object = rowSums(sce_data),
    expected = rowSums(h5ad_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  
  expect_equal(
    object = colSums(sce_data),
    expected = colSums(h5ad_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  
  
  #Test for error on missing feature on Seurat
  expect_error(
    FetchData(
      AML_Seurat,
      slot = "data",
      vars =
        # Nonexistent feature
        c("ab_CD900")
    )
  )
  #Test for error on missing feature on SingleCellExperiment
  expect_error(
    FetchData(
      AML_SCE(),
      slot = "logcounts",
      vars =
        # Nonexistent feature
        c("AB_CD900")
    )
  )
  #Test for error on missing feature on AnnData
  expect_error(
    FetchData(
      AML_h5ad(),
      vars =
        # Nonexistent feature
        c("ab_CD900")
    )
  )
  #Test for error 
  expect_error(
    FetchData(
      AML_Seurat,
      slot = "data",
      vars =
        # Nonexistent feature
        c("ab_CD900")
    )
  )
})
