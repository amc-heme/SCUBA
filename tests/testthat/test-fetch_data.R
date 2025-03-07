# test_that scripts
# A "normal" run with features/metadata/reductions in the object
test_that(
  'fetch_data: "standard" run is consistent in return structure',
  {
    # 1. "standard" run completes without errors
    # Results from each object format are saved for comparison in a later test
    expect_no_error({
      seurat_data <-
        SCUBA::fetch_data(
          SCUBA::AML_Seurat,
          layer = "data",
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
              "nFeature_RNA"
              )
          )
      
      sce_data <-
        SCUBA::fetch_data(
          SCUBA::AML_SCE(),
          layer = "logcounts",
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
        SCUBA::fetch_data(
          SCUBA::AML_h5ad(),
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
      })
    
    # 2. Data returned from all object classes matches the expected format
    # Load from csv in tests/testthat/test_data
    expected_output <-
      utils::read.csv(
        test_path("test_data/fetch_data_standard.csv"),
        # Row names were stored in the first column
        row.names = 1
        )
    
    # Check rowSums and colSums of data
    expect_equal(
      object = rowSums(seurat_data),
      expected = rowSums(expected_output),
      tolerance = 1e-6,
      ignore_attr = TRUE
      )
    
    expect_equal(
      object = rowSums(sce_data),
      expected = rowSums(expected_output),
      tolerance = 1e-6,
      ignore_attr = TRUE
      )
    
    expect_equal(
      object = rowSums(h5ad_data),
      expected = rowSums(expected_output),
      tolerance = 1e-6,
      ignore_attr = TRUE
      )
    })

test_that(
  "fetch_data: raw counts are returned properly",
  {
    # 1. No errors when pulling from counts
    expect_no_error({
      seurat_data <-
        SCUBA::fetch_data(
          SCUBA::AML_Seurat,
          # Raw counts
          layer = "counts",
          vars =
            # Features from default (genes) assay
            c("rna_MEIS1",
              "rna_GAPDH",
              # Features from AB modality
              "ab_CD117-AB",
              "ab_CD123-AB",
              "ab_CD11c-AB"
              )
        )
      
      sce_data <-
        SCUBA::fetch_data(
          SCUBA::AML_SCE(),
          layer = "counts",
          vars =
            # Features from default (genes) assay 
            c("RNA_MEIS1",
              "RNA_GAPDH",
              # Features from AB modality
              "AB_CD117-AB",
              "AB_CD123-AB",
              "AB_CD11c-AB"
              )
          )
      
      h5ad_data <-  
        SCUBA::fetch_data(
          SCUBA::AML_h5ad(),
          layer = "counts",
          vars =
            # Features from default (genes) assay
            c("X_MEIS1",
              "X_GAPDH",
              # Features from AB modality
              "protein_CD117-AB",
              "protein_CD123-AB",
              "protein_CD11c-AB"
              )
          )
    })
    
    # 2. Output is as expected
    expected_output <-
      utils::read.csv(
        test_path("test_data/fetch_data_raw_counts.csv"),
        # Row names were stored in the first column
        row.names = 1
        )
    
    for (df in c(seurat_data, sce_data, h5ad_data)){
      expect_equal(
        object = df,
        expected = rowSums(expected_output),
        tolerance = 1e-6,
        ignore_attr = TRUE
        )
      }
  })


# "Ambiguous" features
# Features from non-default assays without an assay key return a warning,
# But not an error
test_that(
  "fetch_data: ambiguous features are returned properly, with a warning",
  {
    expect_warning({
      ambiguous_seurat <- 
        SCUBA::fetch_data(
          SCUBA::AML_Seurat,
          # "Ambiguous" feature not in RNA assay
          vars = "CD11a-AB"
          )
        }) |> 
      expect_no_error()
    
    expect_warning({
      ambiguous_sce <- 
        SCUBA::fetch_data(
          SCUBA::AML_Seurat,
          # "Ambiguous" feature not in RNA assay
          vars = "CD11a-AB"
        )
      }) |> 
      expect_no_error()
    
    expect_warning({
      ambiguous_h5ad <- 
        SCUBA::fetch_data(
          SCUBA::AML_Seurat,
          # "Ambiguous" feature not in RNA assay
          vars = "CD11a-AB"
        )
      }) |> 
      expect_no_error()
  })

# Test that the class returned from the ambiguous feature query in the last is
# a data.frame, and that the inputs are as expected
test_that(
  "fetch_data: ambiguous features return a data.frame with the correct values",
  {
    expected_output <-
      utils::read.csv(
        test_path("test_data/fetch_data_ambiguous_feature.csv"),
        # Row names were stored in the first column
        row.names = 1
        )
    
    # Test each of the tables created in the previous test
    for (df in c(ambiguous_seurat, ambiguous_sce, ambiguous_h5ad)){
      expect_true(inherits(df, "data.frame"))
      
      expect_equal(
        object = rowSums(df),
        expected = rowSums(expected_output),
        tolerance = 1e-6,
        ignore_attr = TRUE
      )
    }
  })

# Nonsensical vars return an error
test_that(
  "fetch_data: nonsensical vars return an error", 
  {
    expect_error(
      SCUBA::fetch_data(
        SCUBA::AML_Seurat,
        layer = "data",
        vars =
          # Nonexistent feature
          c("ab_CD900")
        )
      )
    
    expect_error(
      SCUBA::fetch_data(
        SCUBA::AML_SCE(),
        layer = "logcounts",
        vars =
          # Nonexistent feature
          c("AB_CD900")
        )
      )
    
    expect_error(
      SCUBA::fetch_data(
        SCUBA::AML_h5ad(),
        vars =
          # Nonexistent feature
          c("protein_CD900")
        )
      )
    })
