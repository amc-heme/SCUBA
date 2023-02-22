# Load SingleCellExperiment object from HDF5 storage
sce <- AML_SCE()

# Dimplot test: test if plot runs without error
test_that("DimPlot runs without error", {
  # Plotting with SingleCellExperiment object
  expect_no_error(
    SCEPlots::DimPlot(
      sce,
      group_by = "condensed_cell_type",
      split_by = "Batch"
    )
  )
  # Seurat object
  expect_no_error(
    SCEPlots::DimPlot(
      AML_Seurat,
      group_by = "condensed_cell_type",
      split_by = "Batch"
    )
  )
})
