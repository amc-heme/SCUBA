# List of test objects to iterate through
test_objects <- list(AML_Seurat, AML_SCE(), AML_h5ad(), AML_h5mu())

# Dummy sets of variables
# none_wrong: all variables are in the test objects
none_wrong <- c("nCount_RNA", "ct", "prop")
# some_wrong: one variable is not in the test object. This should return a 
# warning.
some_wrong <- c("nCount_RNA", "ct", "stuff")
# all_wrong: all variables are not in the test object.
all_wrong <- c("stuff", "things", "junk")

test_that(
  paste0(
    "fetch_metadata runs without errors when all variables entered are ",
    "present in object"
    ),
  {
    # Test on each supported object class
    for (object in test_objects){
      testthat::expect_no_error(
        fetch_metadata(
          object, 
          vars = none_wrong,
          cells = get_all_cells(object)
          )
        )
      }
    })

test_that(
  paste0(
    "fetch_metadata returns an error if some or all variables ",
    "entered are not found"
    ),
  {
    # Seurat objects
    
    # Test on each supported object class
    # Test for no error, but a warning 
    for (object in test_objects[1:3]){
      # MuData does not currently return errors
      testthat::expect_error(
        fetch_metadata(
          object, 
          vars = some_wrong,
          cells = get_all_cells(object)
        )
      )
      
      testthat::expect_error(
        fetch_metadata(
          object, 
          vars = all_wrong,
          cells = get_all_cells(object)
        )
      )
    }
  })


test_that(
  paste0(
    "fetch_metadata with full table = TRUE returns the same results for ",
    "Anndata R6, SingleCellExperiment and Seurat objects."
    ),
  {
    seurat_metadata <-
      fetch_metadata(
        AML_Seurat, 
        full_table = T
  #      vars = "condensed_cell_type",
  #      cells = get_all_cells(AML_Seurat)
  )
    
  sce_metadata <-
    fetch_metadata(
      AML_SCE(), 
      full_table = TRUE
    )
    
  h5ad_metadata <-
    fetch_metadata(
      AML_h5ad(), 
      full_table = TRUE
    )
  
  h5mu_metadata <-
    fetch_metadata(
      AML_h5mu(), 
      full_table = TRUE
    )
  
  # orig.ident and ident columns exist in some test objects but not others
  # This is not an issue with SCUBA
  # Once the test object are uniform this should be removed.
  seurat_metadata$orig.ident <- NULL
  sce_metadata$orig.ident <- NULL
  sce_metadata$ident <- NULL
  h5mu_metadata$RNA_orig.ident <- NULL
  
  #Check rowSums and colSums of data
  expect_equal(
    object = as.matrix(h5ad_metadata),
    expected = as.matrix(seurat_metadata),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(
    object = as.data.frame(seurat_metadata),
    expected = as.data.frame(sce_metadata),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  
  expect_equal(
    object = as.matrix(h5mu_metadata),
    expected = as.matrix(seurat_metadata),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
})
