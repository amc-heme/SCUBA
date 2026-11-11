# layer_names(): the layers of an assay, for every supported class.

test_that("layer_names on Seurat objects", {
  expect_identical(layer_names(AML_Seurat), c("counts", "data", "scale.data"))
  expect_identical(layer_names(AML_Seurat, assay = "AB"), c("counts", "data"))
  expect_error(layer_names(AML_Seurat, assay = "nope"), "Assay 'nope' not found")
  a5 <- AML_Seurat
  a5[["RNA"]] <- suppressWarnings(methods::as(a5[["RNA"]], "Assay5"))
  expect_identical(layer_names(a5), SeuratObject::Layers(a5[["RNA"]]))
})

test_that("layer_names on SingleCellExperiment objects", {
  sce <- AML_SCE()
  expect_identical(layer_names(sce), SummarizedExperiment::assayNames(sce))
  expect_identical(layer_names(sce, assay = SingleCellExperiment::mainExpName(sce)), SummarizedExperiment::assayNames(sce))
  expect_identical(layer_names(sce, assay = "AB"), SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, "AB")))
  expect_error(layer_names(sce, assay = "nope"), "Assay 'nope' not found")
})

test_that("layer_names on unsupported objects warns and returns nothing", {
  expect_warning(res <- layer_names(1:3), "does not know how to handle")
  expect_identical(res, character())
})

test_that("layer_names on anndata objects", {
  skip_if_not_installed("reticulate")
  skip_if_not_installed("anndata")
  skip_if_not(tryCatch(reticulate::py_module_available("anndata"), error = function(e) FALSE), "python anndata not available")
  ad <- AML_h5ad()
  expect_identical(layer_names(ad), c("X", "counts"))
  expect_identical(layer_names(ad, assay = "protein"), "X")
  expect_error(layer_names(ad, assay = "nope"), "Modality 'nope' not found")
})

test_that("layer_names on file-backed objects", {
  skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")
  h <- with_engine("julia", open_h5ad(system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA")))
  expect_identical(layer_names(h), c("X", "counts"))
  expect_identical(layer_names(h, assay = "protein"), "X")
  expect_error(layer_names(h, assay = "nope"), "Modality 'nope' not found")
})
