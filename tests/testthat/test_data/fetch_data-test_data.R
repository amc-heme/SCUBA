# Script to generate CSV file used to test fetch_data outputs
# Run with standard inputs
standard_expected <-
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

# Feature expression data from Raw counts assay 
raw_counts <-
  SCUBA::fetch_data(
    SCUBA::AML_Seurat,
    layer = "counts",
    vars =
      # Features from deafult (genes) assay
      c("rna_MEIS1", 
        "rna_GAPDH",
        # Features from AB modality
        "ab_CD117-AB",
        "ab_CD123-AB",
        "ab_CD11c-AB"
        )
    )

# "Ambiguous" feature in an alternate assay that is not specified with an 
# assay key
ambiguous <-
  SCUBA::fetch_data(
    SCUBA::AML_Seurat,
    # "Ambiguous" feature not in RNA assay
    vars = "CD11a-AB"
  )


# Write to CSV
utils::write.csv(standard_expected, file = "./tests/testthat/test_data/fetch_data_standard.csv")

utils::write.csv(raw_counts, file = "./tests/testthat/test_data/fetch_data_raw_counts.csv")

utils::write.csv(ambiguous, file = "./tests/testthat/test_data/fetch_data_ambiguous_feature.csv")



