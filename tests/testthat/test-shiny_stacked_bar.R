# Load SingleCellExperiment object from HDF5 storage
sce <- AML_SCE()

# shiny_stacked_bar test: test if plot runs without error
test_that("Shiny_stacked_bar runs without error", {
  # Plotting with SingleCellExperiment object
  expect_no_error(
    SCUBA:::shiny_stacked_bar(
      sce,
      group_by = "seurat_clusters",
      split_by = "Batch"
    )
  )
  # Seurat object
  expect_no_error(
    SCUBA:::shiny_stacked_bar(
      AML_Seurat,
      group_by = "seurat_clusters",
      split_by = "Batch"
    )
  )
})

test_that("Shiny_stacked_bar creates doppelganger.", {
  vdiffr::expect_doppelganger(
    title = "shiny_stacked_bar.SingleCellExperiment test",
    SCUBA:::shiny_stacked_bar(
      sce,
      group_by = "seurat_clusters",
      split_by = "Batch"
    )
  )
  #Seurat Object
  vdiffr::expect_doppelganger(
    title = "shiny_stacked_bar.Seurat test",
    SCUBA:::shiny_stacked_bar(
      AML_Seurat,
      group_by = "seurat_clusters",
      split_by = "Batch"
    )
  )
})
