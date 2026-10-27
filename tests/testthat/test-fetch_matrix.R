# fetch_matrix(): features x cells blocks, sparse by default, on both engines and the file-backed class.
ref <- SeuratObject::LayerData(AML_Seurat[["RNA"]], layer = "data")[c("GAPDH", "MEIS1"), ]

test_that("fetch_matrix on Seurat (R engine) returns a features x cells dgCMatrix", {
  m <- with_engine("r", fetch_matrix(AML_Seurat, c("GAPDH", "MEIS1")))
  expect_s4_class(m, "dgCMatrix")
  expect_equal(dim(m), c(2L, ncol(AML_Seurat)))
  expect_equal(as.matrix(m), as.matrix(ref))
  d <- with_engine("r", fetch_matrix(AML_Seurat, c("GAPDH", "MEIS1"), sparse = FALSE, cells = colnames(AML_Seurat)[3:1]))
  expect_true(is.matrix(d))
  expect_equal(d, as.matrix(ref)[, colnames(AML_Seurat)[3:1]])
  expect_equal(dim(with_engine("r", fetch_matrix(AML_Seurat, layer = "counts"))), dim(AML_Seurat[["RNA"]]))
  expect_error(with_engine("r", fetch_matrix(AML_Seurat, "nope")), "not found")
})

test_that("fetch_matrix on SingleCellExperiment (R engine)", {
  sce <- AML_SCE()
  m <- with_engine("r", fetch_matrix(sce, c("GAPDH", "MEIS1")))
  expect_s4_class(m, "dgCMatrix")
  expect_equal(as.matrix(m), as.matrix(ref)[, colnames(sce)], ignore_attr = TRUE)
  expect_equal(nrow(with_engine("r", fetch_matrix(sce, assay = "AB"))), nrow(SingleCellExperiment::altExp(sce, "AB")))
})

test_that("fetch_matrix on the Julia engine matches the R engine", {
  skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")
  for (obj in list(AML_Seurat, AML_SCE())) {
    j <- with_engine("julia", fetch_matrix(obj, c("GAPDH", "MEIS1")))
    r <- with_engine("r", fetch_matrix(obj, c("GAPDH", "MEIS1")))
    expect_s4_class(j, "dgCMatrix")
    expect_equal(as.matrix(j), as.matrix(r), tolerance = 1e-6)
    jd <- with_engine("julia", fetch_matrix(obj, c("GAPDH", "MEIS1"), sparse = FALSE, cells = colnames(obj)[5:1]))
    expect_equal(jd, as.matrix(r)[, colnames(obj)[5:1]], tolerance = 1e-6)
  }
  skip_if_not(exists("open_h5ad", asNamespace("SCUBA")))
  h <- with_engine("julia", open_h5ad(system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA")))
  mh <- fetch_matrix(h, c("GAPDH", "MEIS1"))
  expect_s4_class(mh, "dgCMatrix")
  expect_equal(as.matrix(mh), as.matrix(ref)[, colnames(mh)], tolerance = 1e-6)
  expect_equal(dim(fetch_matrix(h, assay = "protein")), c(197L, 250L))
})
