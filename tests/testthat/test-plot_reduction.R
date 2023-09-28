# Load SingleCellExperiment object from HDF5 storage
sce <- AML_SCE()
AML_anndata <- AML_h5ad()

# Dimplot test: test if plot runs without error
test_that("plot_reduction runs without error", {
  # Plotting with SingleCellExperiment object
  expect_no_error(
    SCUBA::plot_reduction(
      sce,
      group_by = "condensed_cell_type",
      split_by = "Batch"
    )
  )
  # Seurat object
  expect_no_error(
    SCUBA::plot_reduction(
      AML_Seurat,
      group_by = "condensed_cell_type",
      split_by = "Batch"
    )
  )
  #AnndataR6 object
  expect_no_error(
    SCUBA::plot_reduction(
      AML_anndata,
      group_by="condensed_cell_type",
      split_by="Batch"
    )
  )
})



test_that("plot_reduction creates doppelganger.", {
  vdiffr::expect_doppelganger(
    title = "plot_reduction.SingleCellExperiment test",
    SCUBA::plot_reduction(
      sce,
      group_by = "condensed_cell_type",
      split_by = "Batch"
    )
  )
  #Seurat Object
  vdiffr::expect_doppelganger(
    title = "plot_reduction.Seurat test",
    SCUBA::plot_reduction(
      AML_Seurat,
      group_by = "condensed_cell_type",
      split_by = "Batch"
    )
  )
  #AnnDataR6 Object
  vdiffr::expect_doppelganger(
    title = "plot_reduction.AnnDataR6 test",
    SCUBA::plot_reduction(
      AML_anndata,
      group_by = "condensed_cell_type",
      split_by = "Batch"
    )
  )
})
