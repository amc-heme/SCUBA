# Every SCUBA generic stops on an object it has no method for, naming itself,
# the object's class and every supported class. Python objects get a hint.

# One call per generic, with just enough arguments to reach dispatch.
unsupported_calls <- list(
  all_keys = function(x) all_keys(x),
  default_layer = function(x) default_layer(x),
  default_reduction = function(x) default_reduction(x),
  features_in_assay = function(x) features_in_assay(x, "RNA"),
  fetch_cells = function(x) fetch_cells(x, "group", "a"),
  fetch_data = function(x) fetch_data(x, "a"),
  fetch_feature = function(x) fetch_feature(x, "a"),
  fetch_matrix = function(x) fetch_matrix(x),
  fetch_metadata = function(x) fetch_metadata(x, "a"),
  fetch_reduction = function(x) fetch_reduction(x, "umap"),
  get_all_cells = function(x) get_all_cells(x),
  layer_names = function(x) layer_names(x),
  meta_varnames = function(x) meta_varnames(x),
  pseudobulk = function(x) pseudobulk(x, "sample"),
  reduction_dimnames = function(x) reduction_dimnames(x, "umap", 1:2),
  relabel_axis = function(x) relabel_axis(x, "a"),
  update_object_metadata = function(x) update_object_metadata(x, data.frame())
)

supported_classes <- c("AnnDataR6", "Seurat", "SingleCellExperiment", "scuba_source")

unsupported_message <- function(call, object) {
  condition <- tryCatch(call(object), error = identity, warning = identity)
  expect_s3_class(condition, "error")
  conditionMessage(condition)
}

test_that("every generic with a default method is covered here", {
  namespace_names <- ls(asNamespace("SCUBA"), all.names = TRUE)
  defaults <- sub("\\.default$", "", grep("^[a-z_]+\\.default$", namespace_names, value = TRUE))
  expect_setequal(names(unsupported_calls), defaults)
})

test_that("each generic errors with its own name and every supported class", {
  object <- structure(list(), class = "scuba_not_a_class")
  for (generic in names(unsupported_calls)) {
    message <- unsupported_message(unsupported_calls[[generic]], object)
    expect_true(
      startsWith(message, paste0(generic, "() does not know how to handle")),
      label = paste(generic, "names itself:", message)
    )
    expect_match(message, "'scuba_not_a_class'", fixed = TRUE)
    for (class_name in supported_classes) {
      expect_match(message, class_name, fixed = TRUE, label = paste(generic, "lists", class_name))
    }
    expect_no_match(message, "anndata._core", fixed = TRUE)
    expect_no_match(message, "did not wrap", fixed = TRUE)
  }
})

test_that("supported classes are read from the registered methods", {
  expect_identical(SCUBA:::.scuba_supported_classes("fetch_data"), supported_classes)
  # layer_names also registers the raw Python AnnData classes; they are aliases
  expect_identical(SCUBA:::.scuba_supported_classes("layer_names"), supported_classes)
})

test_that("an unwrapped Python AnnData object gets the Python stack hint", {
  requirement <- SCUBA:::.scuba_python_requirements$packages[1]
  for (python_class in c("anndata.AnnData", "anndata._core.anndata.AnnData")) {
    object <- structure(list(), class = c(python_class, "python.builtin.object"))
    for (generic in c("get_all_cells", "fetch_data")) {
      message <- unsupported_message(unsupported_calls[[generic]], object)
      expect_match(message, paste0(generic, "() does not know"), fixed = TRUE)
      expect_match(message, "did not wrap as AnnDataR6", fixed = TRUE)
      expect_match(message, requirement, fixed = TRUE)
    }
  }
})

test_that("other Python objects get the open_h5ad() / open_h5mu() hint", {
  object <- structure(list(), class = c("mudata._core.mudata.MuData", "python.builtin.object"))
  message <- unsupported_message(unsupported_calls$meta_varnames, object)
  expect_match(message, "reached through reticulate", fixed = TRUE)
  expect_no_match(message, "did not wrap", fixed = TRUE)
})

test_that("the error is the same when the Julia engine is selected", {
  object <- structure(list(), class = "scuba_not_a_class")
  expected <- c(
    fetch_data = unsupported_message(unsupported_calls$fetch_data, object),
    fetch_feature = unsupported_message(unsupported_calls$fetch_feature, object)
  )
  local_mocked_bindings(scuba_engine = function() "julia")
  expect_identical(unsupported_message(unsupported_calls$fetch_data, object), expected[["fetch_data"]])
  expect_identical(unsupported_message(unsupported_calls$fetch_feature, object), expected[["fetch_feature"]])
})
