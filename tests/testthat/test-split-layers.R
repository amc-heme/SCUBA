# Seurat v5 objects with split layers (counts.<group> after split(), counts.1 / counts.2 after
# merge()): the R and Julia engines agree, and a SCUBA.jl without split-layer support keeps such
# objects on the R engine.

split_fixture <- function() {
  counts <- SeuratObject::LayerData(AML_Seurat[["RNA"]], layer = "counts")
  whole <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  whole <- SeuratObject::SetAssayData(whole, layer = "data", new.data = log1p(counts))
  whole$batch <- as.character(AML_Seurat$Batch)
  split <- whole
  split[["RNA"]] <- split(split[["RNA"]], f = split$batch)
  list(whole = whole, split = split)
}

test_that("split layers are detected; a plain scale.data layer is not a split", {
  objects <- split_fixture()
  expect_true(SCUBA:::.seurat_has_split_layers(objects$split))
  expect_false(SCUBA:::.seurat_has_split_layers(objects$whole))
  expect_false(SCUBA:::.seurat_has_split_layers(AML_Seurat))
  scaled <- objects$whole
  SeuratObject::LayerData(scaled[["RNA"]], layer = "scale.data") <-
    as.matrix(SeuratObject::LayerData(scaled[["RNA"]], layer = "data"))
  expect_false(SCUBA:::.seurat_has_split_layers(scaled))
})

test_that("without SCUBA.jl split-layer support, split objects stay on the R engine", {
  objects <- split_fixture()
  local_mocked_bindings(scuba_engine = function() "julia", .scuba_jl_split_layers = function() FALSE)
  expect_false(SCUBA:::.scuba_julia_object(objects$split))
  expect_true(SCUBA:::.scuba_julia_object(objects$whole))
  # with support, both go to Julia
  local_mocked_bindings(.scuba_jl_split_layers = function() TRUE)
  expect_true(SCUBA:::.scuba_julia_object(objects$split))
})

test_that("split objects read the same on the Julia engine as on the R engine", {
  skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")
  skip_if_not(SCUBA:::.scuba_jl_split_layers(), "SCUBA.jl without split-layer support")
  objects <- split_fixture()
  split <- objects$split
  features <- c("GAPDH", "MEIS1", "FLT3")
  for (layer in list(NULL, "counts", "data")) {
    expect_equal(with_engine("julia", fetch_data(split, features, layer = layer)),
                 with_engine("r", fetch_data(split, features, layer = layer)), ignore_attr = TRUE)
  }
  cells <- rev(colnames(split))[1:40]
  expect_equal(with_engine("julia", fetch_data(split, features, cells = cells)),
               with_engine("r", fetch_data(split, features, cells = cells)), ignore_attr = TRUE)
  expect_equal(with_engine("julia", fetch_feature(split, features)),
               with_engine("r", fetch_feature(split, features)), ignore_attr = TRUE)
  expect_equal(with_engine("julia", fetch_matrix(split, layer = "counts")),
               with_engine("r", fetch_matrix(split, layer = "counts")), ignore_attr = TRUE)
  julia_pb <- with_engine("julia", pseudobulk(split, "batch"))
  expect_identical(julia_pb$engine, "julia")
  expect_equal(julia_pb$counts, with_engine("r", pseudobulk(split, "batch"))$counts)
  expect_identical(julia_pb$profiles$n_cells, with_engine("r", pseudobulk(split, "batch"))$profiles$n_cells)
  # an explicit piece: the R engine drops cells outside it; Julia returns them as NaN
  piece <- grep("^counts[.]", SeuratObject::Layers(split[["RNA"]]), value = TRUE)[1]
  r_piece <- suppressMessages(suppressWarnings(with_engine("r", fetch_data(split, "GAPDH", layer = piece))))
  julia_piece <- with_engine("julia", fetch_data(split, "GAPDH", layer = piece))
  expect_equal(nrow(julia_piece), ncol(split))
  expect_equal(julia_piece[rownames(r_piece), , drop = FALSE], r_piece, ignore_attr = TRUE)
  expect_true(all(is.nan(julia_piece[setdiff(colnames(split), rownames(r_piece)), 1])))
})
