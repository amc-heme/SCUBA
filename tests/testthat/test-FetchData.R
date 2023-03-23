library(SeuratObject, quietly = TRUE, warn.conflicts = FALSE)
library(Seurat, quietly = TRUE, warn.conflicts = FALSE)

# Dependencies of SingleCellExperiment and HDF5Array
# (these are loaded quietly to avoid messages in the testing section)
library(MatrixGenerics, quietly = TRUE, warn.conflicts = FALSE)
library(BiocGenerics, quietly = TRUE, warn.conflicts = FALSE)
library(S4Vectors, quietly = TRUE, warn.conflicts = FALSE)
library(IRanges, quietly = TRUE, warn.conflicts = FALSE)
library(Biobase, quietly = TRUE, warn.conflicts = FALSE)
library(SummarizedExperiment, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
library(DelayedArray, quietly = TRUE, warn.conflicts = FALSE)
library(rhdf5, quietly = TRUE, warn.conflicts = FALSE)

library(SingleCellExperiment, quietly = TRUE, warn.conflicts = FALSE)
library(HDF5Array, quietly = TRUE, warn.conflicts = FALSE)

# Load Seurat Object
sobj <- AML_Seurat

# Load SCE Object (with DelayedArray integration)
sce <- AML_SCE()

test_that("Fetchdata.SingleCellExperiment equivalent to Seurat method", {
  seurat_data <-
    FetchData(
      sobj,
      slot = "data",
      vars =
        c("ab_CD117-AB",
          "ab_CD123-AB",
          "ab_CD11c-AB",
          "rna_GAPDH",
          "rna_MEIS1",
          # Nonexistent features
          "ab_CD900",
          # Metadata
          "nCount_RNA",
          "nFeature_RNA",
          # "Ambiguous" feature not in RNA assay
          "CD11a-AB"
        )
    )

  sce_data <-
    FetchData(
      sce,
      slot = "logcounts",
      vars =
        c("AB_CD117-AB",
          "AB_CD123-AB",
          "AB_CD11c-AB",
          "RNA_GAPDH",
          "RNA_MEIS1",
          # Nonexistent features
          "AB_CD900",
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
})
