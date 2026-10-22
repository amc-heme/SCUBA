# The Julia engine must give the same results as the legacy R engine for Seurat and SCE
# objects, apart from the documented normalisations (see docs: Compatibility with SCUBA (R)).
skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")

objects <- list(seurat = AML_Seurat, sce = AML_SCE())
vars <- list(
  seurat = c("rna_GAPDH", "ab_CD117-AB", "UMAP_1", "UMAP_2", "nCount_RNA", "condensed_cell_type", "Batch"),
  sce = c("RNA_GAPDH", "AB_CD117-AB", "UMAP_1", "UMAP_2", "nCount_RNA", "condensed_cell_type", "Batch")
)

expect_same <- function(a, b, ...) {
  expect_equal(dim(a), dim(b))
  expect_equal(colnames(a), colnames(b))
  expect_equal(rownames(a), rownames(b))
  for (col in colnames(a)) {
    if (is.factor(a[[col]]) || is.factor(b[[col]])) {
      expect_equal(as.character(a[[col]]), as.character(b[[col]]), info = col)
    } else {
      expect_equal(a[[col]], b[[col]], tolerance = 1e-6, ignore_attr = TRUE, info = col)
    }
  }
}

test_that("routing: objects without reductions and non-character vars", {
  skip_if_not_installed("SeuratObject")
  bare <- SeuratObject::CreateSeuratObject(counts = Matrix::rsparsematrix(5, 4, 0.5, rand.x = function(n) rpois(n, 2) + 1,
                                                                            dimnames = list(paste0("g", 1:5), paste0("c", 1:4))))
  expect_true(with_engine("julia", SCUBA:::.scuba_route_fetch_data(bare, "g1")))
  expect_true(with_engine("julia", SCUBA:::.scuba_route_fetch_data(bare, factor("g1"))))
  expect_false(with_engine("julia", SCUBA:::.scuba_route_fetch_data(bare, "nCount_RNA")))
})

test_that("scuba_julia_warmup(objects) wraps the given objects and leaves results unchanged", {
  r <- scuba_julia_warmup(objects = list(a = AML_Seurat, b = AML_SCE()))
  expect_identical(r, c(a = TRUE, b = TRUE))
  expect_identical(unname(scuba_julia_warmup(objects = AML_Seurat)), TRUE)
  skip_if_not_installed("BPCells")
  # an object whose BPCells directory has vanished: the engine cannot read it
  m <- Matrix::rsparsematrix(6, 5, 0.5, rand.x = function(n) rpois(n, 2) + 1, dimnames = list(paste0("g", 1:6), paste0("c", 1:5)))
  d <- tempfile("gone"); bp <- BPCells::write_matrix_dir(methods::as(m, "IterableMatrix"), d)
  bad <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = bp)); bad[["RNA"]]$data <- bp
  unlink(d, recursive = TRUE)
  expect_warning(r2 <- scuba_julia_warmup(objects = list(bad = bad, good = AML_Seurat)), "warm-up failed for bad")
  expect_identical(r2, c(bad = FALSE, good = TRUE))
  expect_gte(JuliaCall::julia_eval("ScubaR.cache_size()"), 1)
  expect_same(with_engine("julia", fetch_data(AML_Seurat, "rna_GAPDH")), with_engine("r", fetch_data(AML_Seurat, "rna_GAPDH")))
})

test_that("BPCells layers stored transposed on disk (t(write_matrix_dir(t(m)))) match legacy", {
  skip_if_not_installed("BPCells")
  counts <- SeuratObject::LayerData(AML_Seurat[["RNA"]], layer = "counts")
  d <- tempfile("cxg")
  # a cells x genes matrix written as such (col storage = gene-major), then t() for the Seurat layer:
  # BPCells keeps the files and flips a `transpose` flag, so the on-disk shape is swapped
  cxg <- BPCells::write_matrix_dir(methods::as(Matrix::t(counts), "IterableMatrix"), d)
  expect_identical(readBin(file.path(d, "shape"), "integer", n = 4, size = 4)[-(1:2)], c(ncol(counts), nrow(counts)))
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = BPCells::t(cxg)))
  obj[["RNA"]]$data <- BPCells::t(cxg)
  expect_true(obj[["RNA"]]$data@transpose)   # BPCells sets the flag on the outer wrapper too
  expect_silent(with_engine("julia", fetch_data(obj, "GAPDH")))   # native path, no shape-mismatch warning
  # transpose_storage_order() output (row-major, transpose flag TRUE, logical shape on disk) is native too
  rowdir <- BPCells::transpose_storage_order(BPCells::write_matrix_dir(methods::as(counts, "IterableMatrix"), tempfile("gxc")), outdir = tempfile("row"), tmpdir = tempdir())
  obj2 <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = rowdir)); obj2[["RNA"]]$data <- rowdir
  expect_silent(with_engine("julia", fetch_data(obj2, "GAPDH")))
  expect_same(with_engine("julia", fetch_data(obj2, c("GAPDH", "FLT3"))), with_engine("r", fetch_data(obj2, c("GAPDH", "FLT3"))))
  for (g in c(rownames(obj)[1], "GAPDH", rownames(obj)[nrow(obj)])) {
    expect_same(with_engine("julia", fetch_data(obj, g)), with_engine("r", fetch_data(obj, g)))
  }
  expect_same(with_engine("julia", fetch_data(obj, c("GAPDH", "FLT3"), cells = colnames(obj)[c(5, 1, 3)])),
              with_engine("r", fetch_data(obj, c("GAPDH", "FLT3"), cells = colnames(obj)[c(5, 1, 3)])))
  unlink(d, recursive = TRUE)
})

test_that("fetch_feature: keyed names from other assays and metadata collisions match legacy", {
  obj <- AML_Seurat
  expect_same(with_engine("julia", fetch_feature(obj, c("ab_CD117-AB", "GAPDH"))), with_engine("r", fetch_feature(obj, c("ab_CD117-AB", "GAPDH"))))
  expect_same(with_engine("julia", fetch_feature(obj, "ab_CD117-AB", assay = "AB")), with_engine("r", fetch_feature(obj, "ab_CD117-AB", assay = "AB")))
  obj@meta.data[["GAPDH"]] <- seq_len(ncol(obj))   # bare name colliding with a metadata column
  expect_true(SCUBA:::.scuba_feature_collides(obj, "GAPDH"))
  expect_false(SCUBA:::.scuba_feature_collides(obj, "rna_GAPDH"))
  expect_same(suppressWarnings(with_engine("julia", fetch_feature(obj, "GAPDH"))), suppressWarnings(with_engine("r", fetch_feature(obj, "GAPDH"))))
  expect_same(with_engine("julia", fetch_feature(obj, "rna_GAPDH")), with_engine("r", fetch_feature(obj, "rna_GAPDH")))
})

for (form in names(objects)) {
  obj <- objects[[form]]
  v <- vars[[form]]

  test_that(paste("routing: metadata-only fetch_data stays in R, mixed requests go to Julia:", form), {
    expect_false(SCUBA:::.scuba_route_fetch_data(obj, c("nCount_RNA", "Batch")))
    expect_false(SCUBA:::.scuba_route_fetch_data(obj, c("UMAP_1", "UMAP_2")))
    expect_true(with_engine("julia", SCUBA:::.scuba_route_fetch_data(obj, c(v[1], "Batch"))))
    # two reduction dims: legacy SCE drops a single dim (documented divergence)
    expect_same(with_engine("julia", fetch_data(obj, c(v[1], "Batch", "UMAP_1", "UMAP_2"))), with_engine("r", fetch_data(obj, c(v[1], "Batch", "UMAP_1", "UMAP_2"))))
  })

  test_that(paste("fetch_data parity:", form), {
    expect_same(with_engine("julia", fetch_data(obj, v)), with_engine("r", fetch_data(obj, v)))
    cells <- rev(get_all_cells(obj)[1:20])
    expect_same(with_engine("julia", fetch_data(obj, v[c(1, 5, 7)], cells = cells)),
                with_engine("r", fetch_data(obj, v[c(1, 5, 7)], cells = cells)))
    expect_same(with_engine("julia", fetch_data(obj, v[1], layer = "counts")),
                with_engine("r", fetch_data(obj, v[1], layer = "counts")))
    expect_same(with_engine("julia", fetch_data(obj, c("GAPDH", "nCount_RNA"))),
                with_engine("r", fetch_data(obj, c("GAPDH", "nCount_RNA"))))
  })

  test_that(paste("warnings and errors parity:", form), {
    expect_warning(with_engine("julia", fetch_data(obj, "CD11a-AB")), "CD11a-AB")
    expect_warning(with_engine("julia", fetch_data(obj, c(v[1], "nonsense"))), "not found")
    expect_error(with_engine("julia", fetch_data(obj, paste0(sub("_.*", "", v[2]), "_CD900"))),
                 "None of the requested variables were found")
    expect_error(with_engine("julia", fetch_metadata(obj, vars = c("nCount_RNA", "nonsense"))))
    expect_error(with_engine("julia", fetch_cells(obj, "nCount_RNA", 1)), "categorical")
  })

  test_that(paste("metadata, cells, reductions, helpers parity:", form), {
    expect_same(with_engine("julia", fetch_metadata(obj, vars = c("nCount_RNA", "ct", "prop"))),
                with_engine("r", fetch_metadata(obj, vars = c("nCount_RNA", "ct", "prop"))))
    expect_equal(with_engine("julia", fetch_metadata(obj, vars = "nCount_RNA", return_class = "vector")),
                 with_engine("r", fetch_metadata(obj, vars = "nCount_RNA", return_class = "vector")))
    expect_equal(with_engine("julia", fetch_cells(obj, "Batch", "PBMC_200AB")),
                 with_engine("r", fetch_cells(obj, "Batch", "PBMC_200AB")))
    expect_equal(with_engine("julia", get_all_cells(obj)), with_engine("r", get_all_cells(obj)))
    expect_equal(with_engine("julia", meta_varnames(obj)), with_engine("r", meta_varnames(obj)))
    expect_equal(with_engine("julia", default_layer(obj)), with_engine("r", default_layer(obj)))
    expect_equal(with_engine("julia", default_reduction(obj)), with_engine("r", default_reduction(obj)))
    red <- default_reduction(obj)
    expect_same(as.data.frame(with_engine("julia", fetch_reduction(obj, red))),
                as.data.frame(with_engine("r", fetch_reduction(obj, red))))
    expect_equal(with_engine("julia", reduction_dimnames(obj, red, 1:2)),
                 with_engine("r", reduction_dimnames(obj, red, 1:2)))
    expect_equal(as.character(with_engine("julia", unique_values(obj, "condensed_cell_type"))),
                 as.character(with_engine("r", unique_values(obj, "condensed_cell_type"))))
    expect_equal(with_engine("julia", n_unique(obj, "Batch")), with_engine("r", n_unique(obj, "Batch")))
  })

  test_that(paste("features_in_assay and all_keys parity:", form), {
    main <- if (form == "seurat") "RNA" else SingleCellExperiment::mainExpName(obj)
    expect_equal(with_engine("julia", features_in_assay(obj, main)), with_engine("r", features_in_assay(obj, main)))
    # all_keys() is served by the R method under both engines (metadata-type calls are routed to R)
    expect_identical(with_engine("julia", all_keys(obj)), with_engine("r", all_keys(obj)))
  })
}
