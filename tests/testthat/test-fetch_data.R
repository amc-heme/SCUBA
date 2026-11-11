# Fetch data for equivalence testing across objects
# The data below is used in two test_that statements, so is ran before
# both statements
seurat_data <-
  # Warnings are expected here due to ambiguous feature inputs
  # The warnings are not the subject of the tests based on this data, so they
  # are suppressed
  suppressWarnings(
    fetch_data(
      AML_Seurat,
      layer = "data",
      vars = c(
        "ab_CD117-AB",
        "ab_CD123-AB",
        "ab_CD11c-AB",
        "rna_GAPDH",
        "rna_MEIS1",
        # Reductions
        "UMAP_1",
        "UMAP_2",
        # Metadata
        "nCount_RNA",
        "nFeature_RNA",
        # "Ambiguous" feature not in RNA assay
        "CD11a-AB"
      )
    )
  )

sce_data <-
  # Warnings are expected here due to ambiguous feature inputs
  suppressWarnings(
    fetch_data(
      AML_SCE(),
      layer = "logcounts",
      vars = c(
        "AB_CD117-AB",
        "AB_CD123-AB",
        "AB_CD11c-AB",
        "RNA_GAPDH",
        "RNA_MEIS1",
        # Reductions
        "UMAP_1",
        "UMAP_2",
        # Metadata
        "nCount_RNA",
        "nFeature_RNA",
        # "Ambiguous" feature not in RNA assay
        "CD11a-AB"
      )
    )
  )

h5ad_data <-
  suppressWarnings(
    fetch_data(
      AML_h5ad(),
      vars = c(
        "protein_CD117-AB",
        "protein_CD123-AB",
        "protein_CD11c-AB",
        "X_GAPDH",
        "X_MEIS1",
        # Reductions
        "X_umap_1",
        "X_umap_2",
        # Metadata
        "nCount_RNA",
        "nFeature_RNA",
        # "Ambiguous" feature not in RNA assay
        "CD11a-AB"
      )
    )
  )

h5ad_backed_data <-
  suppressWarnings(
    fetch_data(
      AML_h5ad_backed(),
      vars = c(
        "protein_CD117-AB",
        "protein_CD123-AB",
        "protein_CD11c-AB",
        "X_GAPDH",
        "X_MEIS1",
        # Reductions
        "X_umap_1",
        "X_umap_2",
        # Metadata
        "nCount_RNA",
        "nFeature_RNA",
        # "Ambiguous" feature not in RNA assay
        "CD11a-AB"
      )
    )
  )

test_that("All fetch_data methods return the same data.", {
  # Check rowSums and colSums of data
  expect_equal(
    object = rowSums(sce_data),
    expected = rowSums(seurat_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  expect_equal(
    object = colSums(sce_data),
    expected = colSums(seurat_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  expect_equal(
    object = rowSums(sce_data),
    expected = rowSums(h5ad_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  expect_equal(
    object = colSums(sce_data),
    expected = colSums(h5ad_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})

test_that("Disk-backed outputs match in-memory outputs", {
  # Seurat vs. BPCells: incomplete

  # SingleCellExperiment: incomplete

  # anndata
  expect_equal(
    object = rowSums(h5ad_data),
    expected = rowSums(h5ad_backed_data),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})

test_that("fetch_data returns a warning, but no errors, for ambiguous feature inputs not in the main modality", {
  # Seurat
  expect_warning(
    fetch_data(
      AML_Seurat,
      vars = c(
        # Ambiguously entered feature in an alternate assay (AB)
        "CD11a-AB"
      )
    )
  )

  # SingleCellExperiment
  expect_warning(
    fetch_data(
      AML_SCE(),
      vars = c(
        # Ambiguously entered feature in an alternate experiment (AB)
        "CD11a-AB"
      )
    )
  )

  # anndata
  expect_warning(
    fetch_data(
      AML_h5ad(),
      vars = c(
        # Ambiguously entered feature not in an alternate modality (AB)
        "CD11a-AB"
      )
    )
  )

  # Disk-backed anndata
  expect_warning(
    fetch_data(
      AML_h5ad_backed(),
      vars = c(
        # Ambiguously entered feature not in an alternate modality (AB)
        "CD11a-AB"
      )
    )
  )
})

test_that(
  paste0(
    "fetch_data does not warn for ambiguously entered features that are ",
    "in the main modality or are metadata."
  ),
  {
    # No warning for an ambiguous feature in the genes assay
    expect_no_warning(
      fetch_data(AML_Seurat, vars = c("GAPDH", "nCount_RNA"))
    )

    expect_no_warning(
      fetch_data(AML_SCE(), vars = c("GAPDH", "nCount_RNA"))
    )

    expect_no_warning(
      fetch_data(AML_h5ad(), vars = c("GAPDH", "nCount_RNA"))
    )

    expect_no_warning(
      fetch_data(AML_h5ad_backed(), vars = c("GAPDH", "nCount_RNA"))
    )
  }
)

test_that("fetch_data returns a warning, but no errors, for mixed feature inputs", {
  # Inputs below are present features mixed with features not present
  # Seurat. The Julia engine aggregates the missing variables into one warning
  # (as the SCE and AnnData paths always did); SeuratObject::FetchData emits four.
  if (identical(scuba_engine(), "julia")) {
    expect_warning(
      fetch_data(
        AML_Seurat,
        vars = c("ab_CD110-AB", "rna_GAPDH", "rna_MEIS1", "UMAP_1",
                 "nonsense", "rna_MESS1", "ab_CD1110-AB", "UMAP_11")
      ),
      "not found"
    )
  } else expect_warning(
    fetch_data(
      AML_Seurat,
      vars = c(
        # Legitamate features
        "ab_CD110-AB",
        "rna_GAPDH",
        "rna_MEIS1",
        "UMAP_1",
        # Nonsensical input
        "nonsense",
        # Feature with genes key, but that doesn't exist
        "rna_MESS1",
        # Misspelled feature from AB modality
        "ab_CD1110-AB",
        # Reduction from nonsensical index (due to a typo)
        "UMAP_11"
      )
    )
  ) |>
    # Four warnings will be returned in this case. Each one needs to be caught
    # in an expect_warning statement for the test to succeed.
    expect_warning() |>
    expect_warning() |>
    expect_warning()

  # SingleCellExperiment
  expect_warning(
    fetch_data(
      AML_SCE(),
      vars = c(
        # Legitamate features
        "AB_CD110-AB",
        "RNA_GAPDH",
        "RNA_MEIS1",
        "UMAP_1",
        # Nonsensical input
        "nonsense",
        # Feature with genes key, but that doesn't exist
        "RNA_MESS1",
        # Misspelled feature from AB modality
        "AB_CD1110-AB",
        # Reduction from nonsensical index (due to a typo)
        "UMAP_11"
      )
    )
  )

  # anndata
  expect_warning(
    fetch_data(
      AML_h5ad(),
      vars = c(
        # Legitimate vars
        "protein_CD110-AB",
        "X_GAPDH",
        "X_MEIS1",
        "X_umap_1",
        # Nonsensical input
        "nonsense",
        # Feature with genes key, but that doesn't exist
        "X_MESS1",
        # Misspelled feature from AB modality
        "protein_CD1110-AB",
        # Reduction from nonsensical index (due to a typo)
        "X_umap_11"
      )
    )
  )

  # Disk-backed anndata
  expect_warning(
    fetch_data(
      AML_h5ad_backed(),
      vars = c(
        # Legitimate features
        "protein_CD110-AB",
        "X_GAPDH",
        "X_MEIS1",
        "X_umap_1",
        # Nonsensical input
        "nonsense",
        # Feature with genes key, but that doesn't exist
        "X_MESS1",
        # Misspelled feature from AB modality
        "protein_CD1110-AB",
        # Reduction from nonsensical index (due to a typo)
        "X_umap_11"
      )
    )
  )
})

test_that(
  paste0(
    "fetch_data returns a warning when duplicate features are entered, ",
    "but still returns one instance of the duplicate feature."
  ),
  {
    # 1. Seurat
    ## 1.a. Warning with duplicate features
    # Doesn't exist in the Seurat FetchData method. Nice to add but low priority
    # expect_warning(
    #   fetch_data(
    #     AML_Seurat,
    #     vars = c("rna_GAPDH", "rna_GAPDH")
    #     )
    # )

    ## 1.b. But normal feature return
    expect_equal(
      # Test that...
      fetch_data(
        AML_Seurat,
        vars = c("rna_GAPDH", "rna_GAPDH", "rna_MEIS1", "rna_MEIS1")
      ),
      # Is the same as...
      fetch_data(
        AML_Seurat,
        vars = c("rna_GAPDH", "rna_MEIS1")
      ),
      tolerance = 1e-6,
      ignore_attr = TRUE
    )

    ## 1.c. Duplicate features caused by entering a key with one but not the
    # other are caught
    # The Seurat method does not catch this.
    # expect_equal(
    #   # Test that...
    #   fetch_data(
    #     AML_Seurat,
    #     vars = c("rna_GAPDH", "GAPDH")
    #   ),
    #   # Is the same as...
    #   fetch_data(
    #     AML_Seurat,
    #     vars = c("rna_GAPDH")
    #   ),
    #   tolerance = 1e-6,
    #   ignore_attr = TRUE
    # )

    # 2. SingleCellExperiment
    ## 2.a. Proper return of duplicate features, but with a warning
    expect_equal(
      # Test that...
      fetch_data(
        AML_SCE(),
        vars = c("RNA_GAPDH", "RNA_GAPDH", "RNA_MEIS1", "RNA_MEIS1")
      ),
      # Is the same as...
      fetch_data(
        AML_SCE(),
        vars = c("RNA_GAPDH", "RNA_MEIS1")
      ),
      tolerance = 1e-6,
      ignore_attr = TRUE
    ) |>
      expect_warning()

    ## 2.b. Duplicate features caused by entering a key with one but not the
    # other are caught
    expect_equal(
      # Test that...
      fetch_data(
        AML_SCE(),
        vars = c("RNA_GAPDH", "GAPDH", "RNA_MEIS1", "MEIS1")
      ),
      # Is the same as...
      fetch_data(
        AML_SCE(),
        vars = c("RNA_GAPDH", "RNA_MEIS1")
      ),
      tolerance = 1e-6,
      ignore_attr = TRUE
    ) |>
      expect_warning()

    ## 2.c. Column names of return in 2.c. should be equal to the feature
    # specified with the key
    expect_equal(
      fetch_data(
        AML_SCE(),
        vars = c("RNA_GAPDH", "GAPDH", "RNA_MEIS1", "MEIS1")
      ) |>
        colnames(),
      # Is the same as...
      c("RNA_GAPDH", "RNA_MEIS1"),
      ignore_attr = TRUE
    ) |>
      expect_warning()

    # 3. anndata
    ## 3.a. Proper return of duplicate features, but with a warning
    expect_equal(
      # Test that...
      fetch_data(
        AML_h5ad(),
        vars = c("X_GAPDH", "X_GAPDH", "X_MEIS1", "X_MEIS1")
      ),
      # Is the same as...
      fetch_data(
        AML_h5ad(),
        vars = c("X_GAPDH", "X_MEIS1")
      ),
      tolerance = 1e-6,
      ignore_attr = TRUE
    ) |>
      expect_warning()

    ## 3.b. Duplicate features caused by entering a key with one but not the
    # other are caught
    expect_equal(
      # Test that...
      fetch_data(
        AML_h5ad(),
        vars = c("X_GAPDH", "GAPDH", "X_MEIS1", "MEIS1")
      ),
      # Is the same as...
      fetch_data(
        AML_h5ad(),
        vars = c("X_GAPDH", "X_MEIS1")
      ),
      tolerance = 1e-6,
      ignore_attr = TRUE
    ) |>
      expect_warning()

    ## 3.c. Column names of return in 3.c. should be equal to the feature
    # specified with the key
    expect_equal(
      fetch_data(
        AML_h5ad(),
        vars = c("X_GAPDH", "GAPDH", "X_MEIS1", "MEIS1")
      ) |>
        colnames(),
      # Is the same as...
      c("X_GAPDH", "X_MEIS1"),
      ignore_attr = TRUE
    ) |>
      expect_warning()

    # 4. anndata (disk backed)
    ## 4.a. Proper return of duplicate features, but with a warning
    expect_equal(
      # Test that...
      fetch_data(
        AML_h5ad_backed(),
        vars = c("X_GAPDH", "X_GAPDH", "X_MEIS1", "X_MEIS1")
      ),
      # Is the same as...
      fetch_data(
        AML_h5ad_backed(),
        vars = c("X_GAPDH", "X_MEIS1")
      ),
      tolerance = 1e-6,
      ignore_attr = TRUE
    ) |>
      expect_warning()

    ## 4.b. Duplicate features caused by entering a key with one but not the
    # other are caught
    expect_equal(
      # Test that...
      fetch_data(
        AML_h5ad_backed(),
        vars = c("X_GAPDH", "GAPDH", "X_MEIS1", "MEIS1")
      ),
      # Is the same as...
      fetch_data(
        AML_h5ad_backed(),
        vars = c("X_GAPDH", "X_MEIS1")
      ),
      tolerance = 1e-6,
      ignore_attr = TRUE
    ) |>
      expect_warning()

    ## 4.c. Column names of return in 4.c. should be equal to the feature
    # specified with the key
    expect_equal(
      fetch_data(
        AML_h5ad_backed(),
        vars = c("X_GAPDH", "GAPDH", "X_MEIS1", "MEIS1")
      ) |>
        colnames(),
      # Is the same as...
      c("X_GAPDH", "X_MEIS1"),
      ignore_attr = TRUE
    ) |>
      expect_warning()
  }
)

test_that("fetch_data returns an error when a nonsensical feature is entered", {
  #Test for error on missing feature on Seurat
  expect_error(
    fetch_data(
      AML_Seurat,
      layer = "data",
      # Nonexistent feature
      vars = c("ab_CD900")
    )
  )

  #Test for error on missing feature on SingleCellExperiment
  expect_error(
    fetch_data(
      AML_SCE(),
      layer = "logcounts",
      # Nonexistent feature
      vars = c("AB_CD900")
    )
  )

  # anndata
  expect_error(
    fetch_data(
      AML_h5ad(),
      # Nonexistent feature
      vars = c("ab_CD900")
    )
  )

  # Backed anndata
  expect_error(
    fetch_data(
      AML_h5ad_backed(),
      # Nonexistent feature
      vars = c("ab_CD900")
    )
  )
})

test_that("fetch_data.AnnDataR6 errors informatively when reticulate is not installed", {
  # Create a minimal mock AnnDataR6 object
  mock_anndata <- structure(list(), class = "AnnDataR6")
  
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) FALSE,
    .package = "SCUBA"
  )
  
  expect_error(
    fetch_data(mock_anndata, vars = "test_feature"),
    regexp = 'Package "reticulate" must be installed to use this function with anndata objects\\. Install it with: install\\.packages\\("reticulate"\\)',
    fixed = FALSE
  )
})


test_that("conflicting layer and slot fail on every engine; slot is an alias (review #8)", {
  engines <- if (scuba_julia_available(quiet = TRUE)) c("r", "julia") else "r"
  for (engine in engines) with_engine(engine, {
    expect_error(fetch_data(AML_Seurat, "rna_GAPDH", layer = "data", slot = "counts"), "Both `layer` and `slot`")
    expect_identical(fetch_data(AML_Seurat, "rna_GAPDH", slot = "counts"), fetch_data(AML_Seurat, "rna_GAPDH", layer = "counts"))
    sce <- AML_SCE()
    expect_identical(fetch_data(sce, "RNA_GAPDH", slot = "counts"), fetch_data(sce, "RNA_GAPDH", layer = "counts"))
  })
})

test_that("very wide queries warn once on every engine (review #9)", {
  vars <- paste0("rna_", rownames(AML_Seurat)[rep_len(seq_len(nrow(AML_Seurat)), 1000)])
  engines <- if (scuba_julia_available(quiet = TRUE)) c("r", "julia") else "r"
  for (engine in engines) {
    seen <- 0L
    withCallingHandlers(
      invisible(with_engine(engine, fetch_data(AML_Seurat, vars))),
      warning = function(w) {
        if (grepl("A very large number of features was requested", conditionMessage(w))) seen <<- seen + 1L
        invokeRestart("muffleWarning")
      }
    )
    expect_identical(seen, 1L, info = engine)
  }
})

test_that("SeuratObject 4 gets `slot =` and stays on the R engine (review #2)", {
  captured <- NULL
  local_mocked_bindings(.seuratobject_v5 = function() FALSE)
  local_mocked_bindings(FetchData = function(object, vars, ...) { captured <<- list(...); data.frame(x = 1) },
                        .package = "SeuratObject")
  with_engine("r", fetch_data(AML_Seurat, "rna_GAPDH", layer = "counts"))
  expect_identical(captured$slot, "counts")
  expect_null(captured$layer)
  with_engine("r", fetch_data(AML_Seurat, "rna_GAPDH"))
  expect_identical(captured$slot, "data")
  # with the Julia engine selected, Seurat objects still stay on R (no Julia needed to check)
  local_mocked_bindings(scuba_engine = function() "julia")
  expect_false(SCUBA:::.scuba_julia_object(AML_Seurat))
  expect_true(SCUBA:::.scuba_julia_object(AML_SCE()))
})
