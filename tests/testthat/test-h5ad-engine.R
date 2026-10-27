# File-backed .h5ad objects (class scuba_source) served by the Julia engine.
skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")
skip_if_not(isTRUE(JuliaCall::julia_eval("isdefined(ScubaR, :source_open)")), "SCUBA.jl without file-backed glue")

path <- system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA")
obj <- with_engine("julia", open_h5ad(path))
seu <- AML_Seurat
cells <- get_all_cells(obj)

same_values <- function(a, b, tolerance = 1e-6) {
  expect_equal(nrow(a), nrow(b))
  for (k in seq_len(ncol(a))) {
    x <- a[[k]]; y <- b[[k]]
    if (is.factor(x) || is.factor(y)) expect_equal(as.character(x), as.character(y))
    else expect_equal(unname(x), unname(y), tolerance = tolerance, ignore_attr = TRUE)
  }
}

test_that("open_h5ad reads the schema and metadata", {
  expect_s3_class(obj, "scuba_source")
  expect_equal(dim(obj), c(462L, 250L))
  expect_equal(rownames(obj)[1:2], features_in_assay(obj, "X")[1:2])
  expect_length(cells, 250)
  expect_setequal(meta_varnames(obj), setdiff(colnames(seu[[]]), "orig.ident"))
  expect_equal(unname(all_keys(obj)), c("X", "protein", "X_pca", "X_umap"))
  expect_equal(names(all_keys(obj)), unname(all_keys(obj)))
  expect_length(features_in_assay(obj, "protein"), 197)
  expect_error(features_in_assay(obj, "nope"), "not found")
  expect_equal(default_layer(obj), "X")
  expect_equal(default_reduction(obj), "X_umap")
  expect_equal(reduction_dimnames(obj, "X_umap", 1:2), c("X_umap_1", "X_umap_2"))
  expect_error(reduction_dimnames(obj, "X_umap", 1:3), "dimensions")
  expect_equal(unname(scuba_source_reductions(obj)), c(50L, 2L))
  expect_equal(scuba_source_modalities(obj), c("X", "protein"))
  expect_output(print(obj), "462 features")
  expect_error(with_engine("r", open_h5ad(path)), "Julia engine")
})

test_that("fetch_data matches the Seurat object and the R engine", {
  d <- fetch_data(obj, c("X_GAPDH", "protein_CD117-AB", "X_umap_1", "X_umap_2", "Batch", "nCount_RNA", "condensed_cell_type"))
  s <- with_engine("r", fetch_data(seu, c("rna_GAPDH", "ab_CD117-AB", "UMAP_1", "UMAP_2", "Batch", "nCount_RNA", "condensed_cell_type")))
  expect_equal(rownames(d), cells)
  same_values(d, s[cells, ])
  # bare names, obs_ prefix, mixed order restored
  d2 <- fetch_data(obj, c("Batch", "GAPDH", "obs_nCount_RNA", "CD117-AB"))
  expect_equal(names(d2), c("Batch", "GAPDH", "nCount_RNA", "protein_CD117-AB"))
  expect_equal(d2$GAPDH, d$X_GAPDH)
  # cell subset in the requested order; counts layer
  sub <- cells[c(20, 3, 7)]
  d3 <- fetch_data(obj, c("GAPDH", "Batch"), cells = sub)
  expect_equal(rownames(d3), sub)
  expect_equal(d3$GAPDH, d$X_GAPDH[match(sub, cells)])
  cnt <- fetch_data(obj, "GAPDH", layer = "counts")
  scnt <- with_engine("r", fetch_data(seu, "rna_GAPDH", layer = "counts"))
  expect_equal(cnt$GAPDH, unname(scnt[cells, "rna_GAPDH"]), tolerance = 1e-6)
  # metadata-only request never touches Julia; expression-only request has no metadata
  expect_equal(names(fetch_data(obj, c("Batch", "nCount_RNA"))), c("Batch", "nCount_RNA"))
  expect_equal(names(fetch_data(obj, "GAPDH")), "GAPDH")
})

test_that("fetch_data warnings and errors", {
  expect_warning(fetch_data(obj, c("GAPDH", "X_GAPDH")), "dropped|same data")
  expect_warning(d <- fetch_data(obj, c("GAPDH", "nope", "Batch")), "not found")
  expect_equal(names(d), c("GAPDH", "Batch"))
  expect_error(fetch_data(obj, "nope"), "None of the requested variables")
  expect_error(fetch_data(obj, "GAPDH", cells = c(cells[1], "not_a_cell")), "not_a_cell")
  expect_warning(fetch_data(obj, "CD117-AB"), "protein_CD117-AB")
})

test_that("metadata accessors", {
  full <- fetch_metadata(obj, full_table = TRUE)
  expect_equal(dim(full), c(250L, 28L))
  expect_equal(rownames(full), cells)
  expect_equal(as.character(full$Batch), as.character(seu[[]][cells, "Batch"]))
  v <- fetch_metadata(obj, vars = "Batch", cells = cells[3:1], return_class = "vector")
  expect_equal(names(v), cells[3:1])
  expect_equal(as.character(v), as.character(seu[[]][cells[3:1], "Batch"]))
  df <- fetch_metadata(obj, vars = c("Batch", "nCount_RNA"))
  expect_equal(names(df), c("Batch", "nCount_RNA"))
  expect_error(fetch_metadata(obj, vars = "nope"), "None of the variables")
  expect_setequal(unique_values(obj, "Batch"), unique(as.character(seu$Batch)))
  expect_equal(n_unique(obj, "condensed_cell_type"), length(unique(seu$condensed_cell_type)))
  expect_setequal(fetch_cells(obj, "Batch", "BM_200AB"), colnames(seu)[seu$Batch == "BM_200AB"])
  expect_error(fetch_cells(obj, "nCount_RNA", 1), "categorical")
})

test_that("fetch_feature resolves reduction coordinates and metadata through fetch_data", {
  request <- c("X_GAPDH", "X_umap_1", "nCount_RNA")
  expect_equal(fetch_feature(obj, request), fetch_data(obj, request))
  expect_equal(fetch_feature(obj, "X_umap_1"), fetch_data(obj, "X_umap_1"))
})

test_that("reductions", {
  r <- fetch_reduction(obj, "X_umap", cells = cells[1:3])
  expect_equal(names(r), c("X_umap_1", "X_umap_2"))
  expect_equal(rownames(r), cells[1:3])
  s <- with_engine("r", fetch_reduction(seu, "umap", cells = cells[1:3]))
  expect_equal(unname(as.matrix(r)), unname(as.matrix(s)), tolerance = 1e-6)
  expect_error(fetch_reduction(obj, "X_umap", dims = 1:3), "two-length")
  expect_error(fetch_reduction(obj, "nope"), "not found")
  expect_equal(relabel_axis(obj, "X_umap_1"), "Embeddings Value")
  expect_null(relabel_axis(obj, "X_GAPDH"))
  expect_null(relabel_axis(obj, "GAPDH"))
  expect_equal(relabel_axis(obj, "Batch"), "")
})

test_that("update_object_metadata is copy-on-modify and feeds fetch_data", {
  m <- fetch_metadata(obj, full_table = TRUE)
  m$new_col <- seq_len(nrow(m))
  m$Batch <- factor(m$Batch, levels = rev(levels(m$Batch)))
  obj2 <- update_object_metadata(obj, m)
  expect_false("new_col" %in% meta_varnames(obj))
  expect_equal(fetch_data(obj2, c("GAPDH", "new_col"))$new_col, seq_len(250))
  expect_equal(levels(fetch_metadata(obj2, "Batch")$Batch), rev(levels(fetch_metadata(obj, "Batch")$Batch)))
  expect_error(update_object_metadata(obj, m[1:10, ]), "rows")
  bad <- m; rownames(bad) <- NULL
  expect_error(update_object_metadata(obj, bad), "row names")
  shuffled <- m[rev(seq_len(nrow(m))), ]
  expect_equal(fetch_metadata(update_object_metadata(obj, shuffled), full_table = TRUE)$new_col, seq_len(250))
  expect_warning(fetch_data(obj2, "GAPDH")  , NA)
})

test_that("cell subsets", {
  keep <- cells[c(5, 1, 3)]
  sub <- obj[keep, ]
  expect_equal(get_all_cells(sub), keep)
  expect_equal(dim(sub), c(462L, 3L))
  expect_equal(rownames(fetch_data(sub, "GAPDH")), keep)
  expect_equal(rownames(fetch_metadata(sub, full_table = TRUE)), keep)
  expect_equal(rownames(fetch_reduction(sub, "X_umap")), keep)
  expect_error(fetch_data(sub, "GAPDH", cells = cells[2]), "not found")
  expect_equal(get_all_cells(obj[c(TRUE, rep(FALSE, 249))]), cells[1])
  expect_warning(obj[keep, 1], "ignored")
})

test_that("warm-up accepts file-backed objects and the Julia cache is keyed by path", {
  expect_identical(scuba_julia_warmup(objects = list(h = obj)), c(h = TRUE))
  n0 <- JuliaCall::julia_eval("ScubaR.path_cache_size()")
  m <- fetch_metadata(obj, full_table = TRUE); m$x <- 1
  fetch_data(update_object_metadata(obj, m), "GAPDH")
  fetch_data(obj[cells[1:5], ], "GAPDH")
  expect_equal(JuliaCall::julia_eval("ScubaR.path_cache_size()"), n0)
})

test_that("eltype option: single precision gather gives the same doubles in R", {
  d64 <- fetch_data(obj, c("GAPDH", "protein_CD117-AB"))
  d32 <- withr::with_options(list(SCUBA.julia_eltype = "single"), fetch_data(obj, c("GAPDH", "protein_CD117-AB")))
  expect_type(d32$GAPDH, "double")
  expect_equal(d32$GAPDH, d64$GAPDH, tolerance = 1e-6)
  s64 <- with_engine("julia", fetch_data(seu, "rna_GAPDH"))
  s32 <- withr::with_options(list(SCUBA.julia_eltype = "single"), with_engine("julia", fetch_data(seu, "rna_GAPDH")))
  expect_equal(s32$rna_GAPDH, s64$rna_GAPDH, tolerance = 1e-6)
})

test_that("open_h5mu serves MuData files through the same class", {
  h5mu <- file.path(Sys.getenv("SCUBA_JL_PATH", "../../../SCUBA.jl"), "fixtures", "data", "AML.h5mu")
  skip_if_not(file.exists(h5mu), "AML.h5mu fixture not available")
  mu <- with_engine("julia", open_h5mu(h5mu))
  expect_s3_class(mu, "scuba_source")
  expect_equal(mu$format, "h5mu")
  expect_equal(scuba_source_modalities(mu), c("rna", "prot"))
  expect_equal(default_reduction(mu), "rna_X_umap")
  d <- fetch_data(mu, c("rna_GAPDH", "prot_CD117-AB", "rna_X_umap_1", "Batch"))
  h <- fetch_data(obj, c("X_GAPDH", "protein_CD117-AB", "X_umap_1", "Batch"))
  expect_equal(unname(as.matrix(d[, 1:3])), unname(as.matrix(h[, 1:3])), tolerance = 1e-6)
  expect_equal(as.character(d$Batch), as.character(h$Batch))
  expect_true("rna_condensed_cell_type" %in% meta_varnames(mu))
  expect_equal(names(fetch_data(mu, "X_umap_1")), "rna_X_umap_1")
})
