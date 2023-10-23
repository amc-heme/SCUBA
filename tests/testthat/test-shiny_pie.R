# Load SingleCellExperiment object from HDF5 storage
sce <- AML_SCE()

# shiny_pie test: test if plot runs without error
test_that("Shiny_pie runs without error", {
  # Plotting with SingleCellExperiment object
  expect_no_error(
    SCUBA:::shiny_pie(
      sce,
      patient_colname = "Batch",
      group_by = "seurat_clusters"
    )
  )
  # Seurat object
  expect_no_error(
    SCUBA:::shiny_pie(
      AML_Seurat,
      patient_colname = "Batch",
      group_by = "seurat_clusters"
    )
  )
})

test_that("Shiny_pie creates doppelganger.", {
  vdiffr::expect_doppelganger(
    title = "shiny_pie.SingleCellExperiment test",
    SCUBA:::shiny_pie(
      sce,
      patient_colname = "Batch",
      group_by = "seurat_clusters"
    )
  )
  #Seurat Object
  vdiffr::expect_doppelganger(
    title = "shiny_pie.Seurat test",
    SCUBA:::shiny_pie(
      AML_Seurat,
      patient_colname = "Batch",
      group_by = "seurat_clusters"
    )
  )
})
