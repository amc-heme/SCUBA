# pseudobulk(): sums of a counts layer per sample x group, identical on every class, backend and
# engine. The reference is computed in base R from the in-memory AML_Seurat counts.

counts_ref <- SeuratObject::LayerData(AML_Seurat[["RNA"]], layer = "counts")
meta_ref <- AML_Seurat[[]]
pb_reference <- function(m, sample, group = NULL, fun = "sum") {
  key <- if (is.null(group)) as.character(sample) else paste(sample, group, sep = "__")
  out <- t(rowsum(t(as.matrix(m)), key))
  if (fun == "mean") out <- sweep(out, 2, as.numeric(table(key)[colnames(out)]), "/")
  out
}
want_ref <- pb_reference(counts_ref, meta_ref$Batch, meta_ref$condensed_cell_type)
expect_pb_equal <- function(pb, want = want_ref, tolerance = 1e-8) {
  expect_s3_class(pb, "scuba_pseudobulk")
  expect_true(is.matrix(pb$counts) && is.double(pb$counts))
  expect_setequal(colnames(pb$counts), colnames(want))
  expect_equal(unname(pb$counts), unname(want[rownames(pb$counts), colnames(pb$counts), drop = FALSE]), tolerance = tolerance)
  expect_identical(colnames(pb$counts), pb$profiles$profile_id)
  expect_identical(rownames(pb$profiles), pb$profiles$profile_id)
}
# the shipped BPCells object points at the directories of another machine
aml_bpcells <- function() {
  obj <- AML_BPCells
  root <- system.file("extdata", "AML_BPCells", package = "SCUBA")
  for (a in c("RNA", "AB")) for (l in c("counts", "data")) obj <- set_bpcells_dir(obj, a, l, file.path(root, paste0(a, "_", l)))
  obj
}

test_that("pseudobulk on Seurat (R engine) matches the base-R reference", {
  pb <- with_engine("r", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type"))
  expect_pb_equal(pb)
  expect_identical(pb$engine, "r")
  expect_identical(pb$layer, "counts"); expect_identical(pb$assay, "RNA"); expect_identical(pb$fun, "sum")
  expect_identical(pb$sample_by, "Batch"); expect_identical(pb$group_by, "condensed_cell_type")
  expect_equal(pb$n_cells_dropped, c(missing_sample = 0L, missing_group = 0L))
  expect_equal(sum(pb$profiles$n_cells), ncol(AML_Seurat))
  # profiles are sample-major, then group, matching the counts columns
  expect_identical(pb$profiles$profile_id, paste(pb$profiles$sample, pb$profiles$group, sep = "__"))
  expect_false(is.unsorted(match(pb$profiles$sample, unique(pb$profiles$sample))))
  expect_equal(pb$profiles$n_cells, as.integer(table(paste(meta_ref$Batch, meta_ref$condensed_cell_type, sep = "__"))[pb$profiles$profile_id]))
  expect_equal(dim(pb), dim(pb$counts)); expect_identical(dimnames(pb), dimnames(pb$counts)); expect_identical(as.matrix(pb), pb$counts)
  expect_output(print(pb), "462 features x 18 profiles")
})

test_that("samples only, means, assays, feature and cell subsets (R engine)", {
  pb <- with_engine("r", pseudobulk(AML_Seurat, "Batch"))
  expect_pb_equal(pb, pb_reference(counts_ref, meta_ref$Batch))
  expect_true(all(is.na(pb$profiles$group)))
  expect_identical(pb$profiles$profile_id, c("BM_200AB", "PBMC_200AB"))
  expect_true(is.na(pb$group_by))
  pm <- with_engine("r", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type", fun = "mean"))
  expect_pb_equal(pm, pb_reference(counts_ref, meta_ref$Batch, meta_ref$condensed_cell_type, fun = "mean"))
  ab <- SeuratObject::LayerData(AML_Seurat[["AB"]], layer = "counts")
  expect_pb_equal(with_engine("r", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type", assay = "AB")),
                  pb_reference(ab, meta_ref$Batch, meta_ref$condensed_cell_type))
  pf <- with_engine("r", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type", features = c("MEIS1", "GAPDH")))
  expect_identical(rownames(pf$counts), c("MEIS1", "GAPDH"))
  expect_pb_equal(pf, want_ref[c("MEIS1", "GAPDH"), ])
  cells <- rev(colnames(AML_Seurat)[1:80])
  pc <- with_engine("r", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type", cells = cells))
  expect_pb_equal(pc, pb_reference(counts_ref[, cells], meta_ref[cells, "Batch"], meta_ref[cells, "condensed_cell_type"]))
  expect_equal(sum(pc$profiles$n_cells), 80L)
})

test_that("cell-named label vectors, NA / blank labels, factor levels and profile ids", {
  lymph <- stats::setNames(ifelse(meta_ref$condensed_cell_type %in% c("B Cells", "NK Cells"), "lymph", "other"), colnames(AML_Seurat))
  pv <- with_engine("r", pseudobulk(AML_Seurat, "Batch", lymph))
  expect_pb_equal(pv, pb_reference(counts_ref, meta_ref$Batch, lymph))
  expect_true(is.na(pv$group_by))
  # a vector may cover more cells than requested, in any order
  cells <- colnames(AML_Seurat)[10:1]
  expect_pb_equal(with_engine("r", pseudobulk(AML_Seurat, rev(lymph), cells = cells)),
                  pb_reference(counts_ref[, cells], lymph[cells]))
  expect_error(with_engine("r", pseudobulk(AML_Seurat, unname(lymph))), "named by cell")
  expect_error(with_engine("r", pseudobulk(AML_Seurat, lymph[1:5])), "no label for")
  # NA and blank labels are dropped and counted
  obj <- AML_Seurat
  obj$sample2 <- as.character(obj$Batch); obj$sample2[1:5] <- NA; obj$sample2[6:7] <- "  "
  obj$group2 <- as.character(obj$condensed_cell_type); obj$group2[8] <- ""
  pd <- with_engine("r", pseudobulk(obj, "sample2", "group2"))
  expect_equal(pd$n_cells_dropped, c(missing_sample = 7L, missing_group = 1L))
  expect_equal(sum(pd$profiles$n_cells), ncol(obj) - 8L)
  keep <- colnames(obj)[-(1:8)]
  expect_pb_equal(pd, pb_reference(counts_ref[, keep], meta_ref[keep, "Batch"], meta_ref[keep, "condensed_cell_type"]))
  # within a sample, profiles follow the factor levels, or the C-locale order for characters
  groups_in_order <- function(pb, levels) {
    all(vapply(split(pb$profiles$group, pb$profiles$sample), function(g) !is.unsorted(match(g, levels)) && !anyNA(match(g, levels)), logical(1)))
  }
  obj$grp_f <- factor(as.character(obj$condensed_cell_type), levels = rev(sort(unique(as.character(obj$condensed_cell_type)))))
  pf <- with_engine("r", pseudobulk(obj, "Batch", "grp_f"))
  expect_true(groups_in_order(pf, levels(obj$grp_f)))
  pc <- with_engine("r", pseudobulk(obj, "Batch", "condensed_cell_type"))
  expect_true(groups_in_order(pc, sort(unique(as.character(obj$condensed_cell_type)), method = "radix")))
  expect_false(groups_in_order(pc, levels(obj$grp_f)))
  # profile ids stay unique when labels contain the separator
  obj$odd <- ifelse(obj$Batch == "BM_200AB", "a__b", "a")
  obj$odd2 <- ifelse(obj$Batch == "BM_200AB", "c", "b__c")
  po <- with_engine("r", pseudobulk(obj, "odd", "odd2"))
  expect_false(anyDuplicated(po$profiles$profile_id) > 0)
})

test_that("matrix classes on the R engine: dense, SVT, HDF5-backed, BPCells, split v5 layers", {
  dense <- AML_Seurat
  dense[["RNA"]]$counts <- as.matrix(counts_ref)
  expect_pb_equal(with_engine("r", pseudobulk(dense, "Batch", "condensed_cell_type")))
  sce <- AML_SCE()
  expect_s4_class(SummarizedExperiment::assay(sce, "counts"), "DelayedMatrix")
  ph <- with_engine("r", pseudobulk(sce, "Batch", "condensed_cell_type"))
  expect_pb_equal(ph)
  expect_identical(ph$assay, SingleCellExperiment::mainExpName(sce))
  expect_pb_equal(with_engine("r", pseudobulk(sce, "Batch", "condensed_cell_type", assay = "AB")),
                  pb_reference(SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, "AB"), "counts"), meta_ref$Batch, meta_ref$condensed_cell_type))
  svt <- sce
  SummarizedExperiment::assay(svt, "counts") <- methods::as(as.matrix(SummarizedExperiment::assay(sce, "counts")), "SVT_SparseMatrix")
  expect_pb_equal(with_engine("r", pseudobulk(svt, "Batch", "condensed_cell_type")))
  skip_if_not_installed("BPCells")
  bp <- aml_bpcells()
  expect_true(inherits(SeuratObject::LayerData(bp[["RNA"]], layer = "counts"), "IterableMatrix"))
  expect_pb_equal(with_engine("r", pseudobulk(bp, "Batch", "condensed_cell_type")))
  # a cell subset exercises the sentinel level for dropped cells
  cells <- rev(colnames(bp)[1:100])
  expect_pb_equal(with_engine("r", pseudobulk(bp, "Batch", "condensed_cell_type", cells = cells)),
                  pb_reference(counts_ref[, cells], meta_ref[cells, "Batch"], meta_ref[cells, "condensed_cell_type"]))
  # row-major storage
  row <- BPCells::transpose_storage_order(SeuratObject::LayerData(bp[["RNA"]], layer = "counts"), outdir = tempfile("row"), tmpdir = tempdir())
  bpr <- bp; bpr[["RNA"]]$counts <- row
  expect_pb_equal(with_engine("r", pseudobulk(bpr, "Batch", "condensed_cell_type")))
  # split v5 layers (counts.BM_200AB, counts.PBMC_200AB) are summed layer by layer
  a5 <- AML_Seurat
  a5[["RNA"]] <- suppressWarnings(methods::as(a5[["RNA"]], "Assay5"))
  a5[["RNA"]] <- suppressMessages(split(a5[["RNA"]], f = a5$Batch))
  expect_true(all(c("counts.BM_200AB", "counts.PBMC_200AB") %in% SeuratObject::Layers(a5[["RNA"]])))
  p5 <- with_engine("r", pseudobulk(a5, "Batch", "condensed_cell_type"))
  expect_pb_equal(p5)
  expect_equal(sum(p5$profiles$n_cells), ncol(a5))
})

test_that("large counts do not overflow", {
  m <- Matrix::sparseMatrix(i = c(1, 1, 2), j = c(1, 2, 2), x = c(2e9, 2e9, 3e9), dims = c(2, 3),
                            dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = m, meta.data = data.frame(s = c("a", "a", "b"), row.names = colnames(m))))
  pb <- with_engine("r", pseudobulk(obj, "s"))
  expect_equal(unname(pb$counts["g1", ]), c(4e9, 0))
  expect_equal(unname(pb$counts["g2", ]), c(3e9, 0))
})

test_that("errors are informative", {
  expect_error(pseudobulk(1:3, "s"), "does not know how to handle")
  expect_error(with_engine("r", pseudobulk(AML_Seurat, "Batch", layer = "nope")), "Layer 'nope' not found")
  expect_error(with_engine("r", pseudobulk(AML_Seurat, "Batch", assay = "nope")), "Assay 'nope' not found")
  expect_error(with_engine("r", pseudobulk(AML_Seurat, "nope")), "nope")
  expect_error(with_engine("r", pseudobulk(AML_Seurat, "Batch", features = "ZZZ")), "Features not found")
  expect_error(with_engine("r", pseudobulk(AML_Seurat, "Batch", cells = "zz")), "Cells not found")
  expect_error(with_engine("r", pseudobulk(AML_SCE(), "Batch", layer = "logcounts2")), "Layer 'logcounts2' not found")
  expect_error(with_engine("r", pseudobulk(AML_Seurat, "Batch", fun = "max")))
})

test_that("pseudobulk_design previews the same profiles from metadata only", {
  d <- pseudobulk_design(AML_Seurat, "Batch", "condensed_cell_type")
  expect_s3_class(d, "scuba_pseudobulk_design")
  pb <- with_engine("r", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type"))
  expect_identical(d$profiles, pb$profiles)
  expect_equal(dim(d$table), c(2L, 10L))
  expect_equal(d$table["BM_200AB", "B Cells"], 14L)
  expect_equal(d$table["PBMC_200AB", "BM Monocytes"], 0L)
  expect_true(d$paired)
  expect_equal(unname(d$samples_per_group["B Cells"]), 2L)
  expect_equal(unname(d$samples_per_group["PBMC Monocytes"]), 1L)
  expect_equal(d$samples_in_all_groups, 0L)
  expect_equal(d$n_cells, 250L)
  expect_identical(d$sample_by, "Batch"); expect_identical(d$group_by, "condensed_cell_type")
  expect_output(print(d), "paired")
  # an independent design: samples nested in groups
  nested <- pseudobulk_design(AML_Seurat, "Batch", stats::setNames(as.character(meta_ref$Batch), colnames(AML_Seurat)))
  expect_false(nested$paired)
  expect_equal(unname(nested$samples_per_group), c(1L, 1L))
  expect_output(print(nested), "independent")
  only <- pseudobulk_design(AML_Seurat, "Batch")
  expect_equal(dim(only$table), c(2L, 1L)); expect_true(is.na(only$group_by))
})

test_that("anndata objects (reticulate) aggregate in memory and in backed mode", {
  skip_if_not_installed("reticulate")
  skip_if_not_installed("anndata")
  skip_if_not(tryCatch(reticulate::py_module_available("anndata"), error = function(e) FALSE), "python anndata not available")
  ad <- AML_h5ad()
  expect_true("counts" %in% layer_names(ad))
  pa <- pseudobulk(ad, "Batch", "condensed_cell_type")
  expect_identical(pa$engine, "python")
  expect_pb_equal(pa)
  expect_identical(pa$assay, "X")
  cells <- rev(colnames(AML_Seurat)[1:50])
  ps <- pseudobulk(ad, "Batch", "condensed_cell_type", cells = cells, features = c("GAPDH", "MEIS1"))
  expect_pb_equal(ps, pb_reference(counts_ref[c("GAPDH", "MEIS1"), cells], meta_ref[cells, "Batch"], meta_ref[cells, "condensed_cell_type"]))
  # X holds the normalised values; the protein obsm data frame is a modality of its own
  px <- pseudobulk(ad, "Batch", layer = "X")
  expect_false(all(px$counts == round(px$counts)))
  expect_equal(dim(pseudobulk(ad, "Batch", "condensed_cell_type", assay = "protein", layer = "X")), c(197L, 18L))
  expect_error(pseudobulk(ad, "Batch", layer = "nope"), "Layer 'nope' not found")
  adb <- AML_h5ad_backed()
  expect_pb_equal(pseudobulk(adb, "Batch", "condensed_cell_type"))
  expect_pb_equal(pseudobulk(adb, "Batch", "condensed_cell_type", cells = cells),
                  pb_reference(counts_ref[, cells], meta_ref[cells, "Batch"], meta_ref[cells, "condensed_cell_type"]))
})

test_that("the Julia engine gives the same profiles and sums as the R engine", {
  skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")
  skip_if_not(isTRUE(JuliaCall::julia_eval("isdefined(ScubaR, :pseudobulk)")), "SCUBA.jl without ScubaR.pseudobulk")
  pj <- with_engine("julia", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type"))
  pr <- with_engine("r", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type"))
  expect_identical(pj$engine, "julia")
  expect_pb_equal(pj)
  expect_identical(pj$profiles, pr$profiles)
  expect_equal(pj$counts, pr$counts)
  cells <- rev(colnames(AML_Seurat)[1:60])
  expect_equal(with_engine("julia", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type", cells = cells, features = c("MEIS1", "GAPDH")))$counts,
               with_engine("r", pseudobulk(AML_Seurat, "Batch", "condensed_cell_type", cells = cells, features = c("MEIS1", "GAPDH")))$counts)
  expect_equal(with_engine("julia", pseudobulk(AML_Seurat, "Batch", fun = "mean", assay = "AB"))$counts,
               with_engine("r", pseudobulk(AML_Seurat, "Batch", fun = "mean", assay = "AB"))$counts)
  # HDF5-backed SCE stays on the R primitive even on the Julia engine; an in-memory SCE goes to Julia
  sce <- AML_SCE()
  ph <- with_engine("julia", pseudobulk(sce, "Batch", "condensed_cell_type"))
  expect_identical(ph$engine, "r"); expect_pb_equal(ph)
  mem <- sce
  SummarizedExperiment::assay(mem, "counts") <- methods::as(as.matrix(SummarizedExperiment::assay(sce, "counts")), "dgCMatrix")
  pm <- with_engine("julia", pseudobulk(mem, "Batch", "condensed_cell_type"))
  expect_identical(pm$engine, "julia"); expect_pb_equal(pm)
  skip_if_not_installed("BPCells")
  bp <- aml_bpcells()
  pb <- with_engine("julia", pseudobulk(bp, "Batch", "condensed_cell_type"))
  expect_identical(pb$engine, "julia"); expect_pb_equal(pb)
  expect_pb_equal(with_engine("julia", pseudobulk(bp, "Batch", "condensed_cell_type", cells = cells)),
                  pb_reference(counts_ref[, cells], meta_ref[cells, "Batch"], meta_ref[cells, "condensed_cell_type"]))
})

test_that("file-backed h5ad / h5mu objects aggregate through the Julia engine", {
  skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")
  skip_if_not(isTRUE(JuliaCall::julia_eval("isdefined(ScubaR, :pseudobulk)")), "SCUBA.jl without ScubaR.pseudobulk")
  h <- with_engine("julia", open_h5ad(system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA")))
  expect_identical(layer_names(h), c("X", "counts"))
  ph <- pseudobulk(h, "Batch", "condensed_cell_type")
  expect_identical(ph$engine, "julia")
  expect_pb_equal(ph)
  expect_identical(ph$assay, "X")
  cells <- rev(colnames(AML_Seurat)[1:60])
  want_sub <- pb_reference(counts_ref[c("MEIS1", "GAPDH"), cells], meta_ref[cells, "Batch"], meta_ref[cells, "condensed_cell_type"])
  expect_pb_equal(pseudobulk(h, "Batch", "condensed_cell_type", cells = cells, features = c("MEIS1", "GAPDH")), want_sub)
  expect_pb_equal(pseudobulk(h[cells], "Batch", "condensed_cell_type", features = c("MEIS1", "GAPDH")), want_sub)
  # an R-side metadata column drives the grouping
  tab <- fetch_metadata(h, full_table = TRUE)
  tab$metacluster <- ifelse(tab$condensed_cell_type %in% c("B Cells", "NK Cells"), "lymph", "other")
  h2 <- update_object_metadata(h, tab)
  po <- pseudobulk(h2, "Batch", "metacluster")
  expect_pb_equal(po, pb_reference(counts_ref, meta_ref$Batch, tab[colnames(counts_ref), "metacluster"]))
  expect_identical(pseudobulk_design(h2, "Batch", "metacluster")$profiles, po$profiles)
  expect_equal(dim(pseudobulk(h, "Batch", "condensed_cell_type", assay = "protein", layer = "X", fun = "mean")), c(197L, 18L))
  expect_error(pseudobulk(h, "Batch", layer = "nope"), "Layer 'nope' not found")
  expect_error(with_engine("r", pseudobulk(h, "Batch")), "Julia engine")
  # MuData, through SCUBA.jl's fixtures when they are present
  h5mu <- tryCatch(JuliaCall::julia_eval('SCUBA.fixture_path("AML.h5mu")'), error = function(e) "")
  skip_if_not(nzchar(h5mu) && file.exists(h5mu), "AML.h5mu fixture not available")
  mu <- open_h5mu(h5mu)
  expect_identical(scuba_source_modalities(mu)[1], "rna")
  pmu <- pseudobulk(mu, "Batch", "condensed_cell_type")
  expect_pb_equal(pmu)
  expect_identical(pmu$assay, "rna")
  expect_equal(dim(pseudobulk(mu, "Batch", "condensed_cell_type", assay = "prot", layer = "X")), c(197L, 18L))
})
