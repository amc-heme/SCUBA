# Setup: Load BPCells test object and set paths to the correct context
# (BPCells matrix paths are absolute and must be set to the path to the test
# dataset matrices in whatever context the tests are being run)
BPCells_object <- AML_BPCells

# Set BPCells directory paths for RNA assay (counts and data layers)
SCUBA::set_bpcells_dir(
  object = BPCells_object,
  assay = "RNA",
  layer = "counts",
  dirname = system.file(
    "extdata",
    "AML_BPCells",
    "RNA_counts",
    package = "SCUBA"
  )
)

SCUBA::set_bpcells_dir(
  object = BPCells_object,
  assay = "RNA",
  layer = "data",
  dirname = system.file(
    "extdata",
    "AML_BPCells",
    "RNA_data",
    package = "SCUBA"
  )
)

# Set BPCells directory paths for AB assay (counts and data layers)
SCUBA::set_bpcells_dir(
  object = BPCells_object,
  assay = "AB",
  layer = "counts",
  dirname = system.file(
    "extdata",
    "AML_BPCells",
    "AB_counts",
    package = "SCUBA"
  )
)

SCUBA::set_bpcells_dir(
  object = BPCells_object,
  assay = "AB",
  layer = "data",
  dirname = system.file(
    "extdata",
    "AML_BPCells",
    "AB_data",
    package = "SCUBA"
  )
)

seurat_keyed_features <-
  c(
    "rna_GAPDH",
    "rna_FLT3",
    "rna_FIS1",
    "ab_CD19-AB",
    "ab_CD4-AB",
    "ab_CD8-AB"
  )

seurat_unkeyed_features <-
  c(
    "GAPDH",
    "FLT3",
    "FIS1"
  )

# Test 1: fetch_feature vs. fetch_data equivalence for standard Seurat object
test_that("fetch_feature returns same data as fetch_data for Seurat objects with keyed features", {
  data_from_fetch_data <-
    SCUBA::fetch_data(
      object = AML_Seurat,
      vars = seurat_keyed_features,
      layer = "data"
    )

  data_from_fetch_feature <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_keyed_features,
      layer = "data"
    )

  # Check that dimensions match
  expect_equal(
    object = dim(data_from_fetch_feature),
    expected = dim(data_from_fetch_data)
  )

  # Check that column names match
  expect_equal(
    object = colnames(data_from_fetch_feature),
    expected = colnames(data_from_fetch_data)
  )

  # Check that row names match
  expect_equal(
    object = rownames(data_from_fetch_feature),
    expected = rownames(data_from_fetch_data)
  )

  # Check that values are numerically equivalent
  expect_equal(
    object = data_from_fetch_feature,
    expected = data_from_fetch_data,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})

# Test 2: fetch_feature vs. fetch_data equivalence for standard Seurat with unkeyed features
test_that("fetch_feature returns same data as fetch_data for Seurat objects with unkeyed features", {
  data_from_fetch_data <-
    SCUBA::fetch_data(
      object = AML_Seurat,
      vars = seurat_unkeyed_features,
      layer = "data"
    )

  data_from_fetch_feature <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_unkeyed_features,
      layer = "data"
    )

  # Check that dimensions match
  expect_equal(
    object = dim(data_from_fetch_feature),
    expected = dim(data_from_fetch_data)
  )

  # Check that column names match
  expect_equal(
    object = colnames(data_from_fetch_feature),
    expected = colnames(data_from_fetch_data)
  )

  # Check that values are numerically equivalent
  expect_equal(
    object = data_from_fetch_feature,
    expected = data_from_fetch_data,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})

# Test 3: fetch_feature with cell subsetting
test_that("fetch_feature correctly subsets cells", {
  selected_cells <- colnames(AML_Seurat)[1:10]

  data_from_fetch_data <-
    SCUBA::fetch_data(
      object = AML_Seurat,
      vars = seurat_unkeyed_features,
      layer = "data",
      cells = selected_cells
    )

  data_from_fetch_feature <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_unkeyed_features,
      layer = "data",
      cells = selected_cells
    )

  # Check that only requested cells are returned
  expect_equal(
    object = nrow(data_from_fetch_feature),
    expected = 10
  )

  # Check that rownames match selected cells
  expect_equal(
    object = rownames(data_from_fetch_feature),
    expected = selected_cells
  )

  # Check equivalence with fetch_data
  expect_equal(
    object = data_from_fetch_feature,
    expected = data_from_fetch_data,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})

# Test 4: fetch_feature preserves column order
test_that("fetch_feature preserves requested feature order", {
  unordered_features <- c("ab_CD8-AB", "rna_GAPDH", "ab_CD4-AB", "rna_FLT3")

  data_from_fetch_feature <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = unordered_features,
      layer = "data"
    )

  # Check that column order matches requested order
  expect_equal(
    object = colnames(data_from_fetch_feature),
    expected = unordered_features
  )
})

# Test 5: fetch_feature with BPCells object equivalence
test_that("fetch_feature returns same data for BPCells and standard Seurat", {
  # Check that both objects have the same dimensions
  expect_equal(
    object = dim(AML_Seurat),
    expected = dim(BPCells_object)
  )

  data_from_standard <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_unkeyed_features,
      layer = "data"
    )

  data_from_bpcells <-
    SCUBA::fetch_feature(
      object = BPCells_object,
      features = seurat_unkeyed_features,
      layer = "data"
    )

  # Check dimensions match
  expect_equal(
    object = dim(data_from_bpcells),
    expected = dim(data_from_standard)
  )

  # Check that values are equivalent (BPCells may have minor floating point diff)
  expect_equal(
    object = data_from_bpcells,
    expected = data_from_standard,
    tolerance = 1e-4,
    ignore_attr = TRUE
  )
})

# Test 6: fetch_feature with default layer
test_that("fetch_feature uses default layer when layer is NULL", {
  data_with_null_layer <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_unkeyed_features,
      layer = NULL
    )

  data_with_default_layer <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_unkeyed_features,
      layer = SCUBA::default_layer(AML_Seurat)
    )

  expect_equal(
    object = data_with_null_layer,
    expected = data_with_default_layer,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})

# Test 7: fetch_feature with default cells (all cells)
test_that("fetch_feature returns all cells when cells is NULL", {
  data_with_null_cells <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_unkeyed_features,
      layer = "data",
      cells = NULL
    )

  data_with_all_cells <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_unkeyed_features,
      layer = "data",
      cells = colnames(AML_Seurat)
    )

  expect_equal(
    object = nrow(data_with_null_cells),
    expected = ncol(AML_Seurat)
  )

  expect_equal(
    object = data_with_null_cells,
    expected = data_with_all_cells,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})

# Test 8: fetch_feature returns data.frame
test_that("fetch_feature returns a data.frame", {
  result <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = seurat_unkeyed_features,
      layer = "data"
    )

  expect_s3_class(result, "data.frame")
})

# Test 9: fetch_feature handles mixed keyed and unkeyed features
test_that("fetch_feature handles mixed keyed and unkeyed features", {
  mixed_features <- c("rna_GAPDH", "FLT3", "ab_CD4-AB")

  data_from_fetch_data <-
    SCUBA::fetch_data(
      object = AML_Seurat,
      vars = mixed_features,
      layer = "data"
    )

  data_from_fetch_feature <-
    SCUBA::fetch_feature(
      object = AML_Seurat,
      features = mixed_features,
      layer = "data"
    )

  # Check dimensions match
  expect_equal(
    object = dim(data_from_fetch_feature),
    expected = dim(data_from_fetch_data)
  )

  # Check column order matches
  expect_equal(
    object = colnames(data_from_fetch_feature),
    expected = mixed_features
  )

  # Check values are equivalent
  expect_equal(
    object = data_from_fetch_feature,
    expected = data_from_fetch_data,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})

# Test 10: fetch_feature returns warning for unsupported object types
test_that("fetch_feature.default returns a warning for unsupported object types", {
  unsupported_object <- list(data = matrix(1:10, nrow = 2))

  expect_warning(
    SCUBA::fetch_feature(
      object = unsupported_object,
      features = c("feature1", "feature2")
    )
  )
})
