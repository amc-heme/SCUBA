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


test_that("keyed feature names select their assay on the R engine (review #10)", {
  data_ab <- SeuratObject::LayerData(AML_Seurat[["AB"]], layer = "data")
  keyed <- with_engine("r", fetch_matrix(AML_Seurat, "rna_GAPDH"))
  expect_identical(rownames(keyed), "rna_GAPDH")
  expect_equal(unname(as.matrix(keyed)), unname(as.matrix(ref["GAPDH", , drop = FALSE])))
  protein <- with_engine("r", fetch_matrix(AML_Seurat, "ab_CD117-AB"))
  expect_equal(unname(as.matrix(protein)), unname(as.matrix(data_ab["CD117-AB", , drop = FALSE])))
  expect_identical(rownames(with_engine("r", fetch_matrix(AML_Seurat, c("GAPDH", "rna_MEIS1")))), c("GAPDH", "rna_MEIS1"))
  expect_error(with_engine("r", fetch_matrix(AML_Seurat, c("rna_GAPDH", "ab_CD117-AB"))), "several assays")
  sce <- AML_SCE()
  expect_identical(rownames(with_engine("r", fetch_matrix(sce, "RNA_GAPDH"))), "RNA_GAPDH")
  expect_equal(dim(with_engine("r", fetch_matrix(sce, "AB_CD117-AB"))), c(1L, ncol(sce)))
  skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")
  for (request in list("rna_GAPDH", "ab_CD117-AB", c("GAPDH", "rna_MEIS1"))) {
    expect_equal(with_engine("r", fetch_matrix(AML_Seurat, request)),
                 with_engine("julia", fetch_matrix(AML_Seurat, request)))
  }
})

test_that("split v5 layers are assembled for every cell (review #3)", {
  counts <- SeuratObject::LayerData(AML_Seurat[["RNA"]], layer = "counts")
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  object$Batch <- AML_Seurat$Batch
  object[["RNA"]] <- split(object[["RNA"]], f = object$Batch)
  engines <- if (scuba_julia_available(quiet = TRUE)) c("r", "julia") else "r"
  for (engine in engines) {
    full <- with_engine(engine, fetch_matrix(object, layer = "counts"))
    expect_equal(dim(full), dim(counts))
    expect_equal(unname(as.matrix(full)), unname(as.matrix(counts[, colnames(full)])))
    some <- colnames(counts)[c(250, 1, 130)]
    sub <- with_engine(engine, fetch_matrix(object, c("GAPDH", "MEIS1"), layer = "counts", cells = some))
    expect_identical(colnames(sub), some)
    expect_equal(unname(as.matrix(sub)), unname(as.matrix(counts[c("GAPDH", "MEIS1"), some])))
  }
})

test_that("sparse results do not go through a dense copy (review #4)", {
  m <- SeuratObject::LayerData(AML_Seurat[["RNA"]], layer = "counts")[1:20, 1:30]
  dense_calls <- 0L
  original_as_matrix <- base::as.matrix
  local_mocked_bindings(as.matrix = function(x, ...) { dense_calls <<- dense_calls + 1L; original_as_matrix(x, ...) },
                        .package = "base")
  inputs <- list(methods::as(m, "TsparseMatrix"), methods::as(m, "RsparseMatrix"))
  if (requireNamespace("DelayedArray", quietly = TRUE)) inputs <- c(inputs, list(DelayedArray::DelayedArray(m)))
  for (x in inputs) {
    out <- SCUBA:::.as_dgcmatrix(x)
    expect_s4_class(out, "dgCMatrix")
    expect_equal(unname(as(out, "matrix")), unname(as(m, "matrix")))
  }
  expect_identical(dense_calls, 0L)
})
