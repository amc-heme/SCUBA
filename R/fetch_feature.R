#' Fetch feature expression data from single-cell objects
#'
#' @description
#' **This function is still in development.** Currently, it only has a full
#' implementation for Seurat objects. For all other object types, it acts as
#' a wrapper around [fetch_data()].
#'
#' This function retrieves feature expression data from single-cell objects.
#' For Seurat objects with BPCells-backed assay layers, this function uses
#' optimized direct matrix access for improved performance. For all other
#' cases, it delegates to [fetch_data()].
#'
#' Unlike [fetch_data()], this function is specifically designed for feature
#' retrieval and does not support fetching metadata or reduction coordinates.
#'
#' See our GitHub.io website for additional information and examples.
#'
#' @inheritParams object_param
#' @param features A character vector of feature names to retrieve from the
#'   object. Features should be specified with the assay key prefix if pulling
#'   from a non-default assay (e.g., "rna_FLT3"). To determine the key that
#'   corresponds to the assay to pull data from, run [all_keys()].
#' @param assay For Seurat objects, the assay to pull data from. If `NULL`,
#'   the default assay will be used.
#' @param layer The layer to pull data from. Layers are referred to as "slots"
#'   in Seurat objects v4 and earlier, and "assays" in SingleCellExperiment
#'   objects. If `NULL`, the default layer will be used.
#' @param cells A character vector of cell names to include, as they are named
#'   in the object (i.e. according to `colnames(object)`). If `NULL`, data will
#'   be returned for all cells in the object.
#'
#' @returns A data.frame with the requested `features` as columns and the cells
#'   as rows.
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' # Fetch feature expression data
#' fetch_feature(AML_Seurat, features = c("FLT3", "NPM1")) |> head()
#'
#' # Fetch from a specific layer
#' fetch_feature(AML_Seurat, features = "FLT3", layer = "data") |> head()
#'
fetch_feature <- function(
  object,
  features,
  assay = NULL,
  layer = NULL,
  cells = NULL
) {
  UseMethod("fetch_feature")
}

#' Default method for unsupported object types
#'
#' @noRd
#' @export
fetch_feature.default <- function(
  object,
  features,
  assay = NULL,
  layer = NULL,
  cells = NULL
) {
  warning(
    paste0(
      "fetch_feature does not know how to handle object of class ",
      paste(class(object), collapse = ", "),
      ". Currently supported classes: Seurat, SingleCellExperiment, and AnnDataR6."
    )
  )
}

#' @describeIn fetch_feature Seurat objects. Uses optimized BPCells access
#'   when available, otherwise delegates to [fetch_data()].
#'
#' @keywords internal
#'
#' @export
fetch_feature.Seurat <- function(
  object,
  features,
  assay = NULL,
  layer = NULL,
  cells = NULL
) {
  # Extract assay for each feature from key prefixes
  # Get all keys in the object (named vector: assay name -> key)
  object_keys <- SCUBA::all_keys(object)
  default_assay <- SeuratObject::DefaultAssay(object)

  # Set default layer if not provided
  layer <- layer %||% SCUBA::default_layer(object)

  # Set default cells if not provided
  cells <- cells %||% colnames(object)

  # Build a named list mapping each assay key to its features

  # Features without a matching key prefix are assigned to the default assay
  # Structure: list(assay_name = list(bare_feature = original_feature, ...))
  features_by_assay <- list()

  for (feature in features) {
    matched_assay <- NULL
    matched_feature <- feature

    # Check if feature starts with any known key prefix
    for (assay_name in names(object_keys)) {
      key <- object_keys[assay_name]
      # Keys in Seurat end with underscore, e.g., "rna_"
      # Match pattern: key prefix followed by feature name
      pattern <- paste0("^", key)
      if (grepl(pattern, feature, ignore.case = TRUE)) {
        matched_assay <- assay_name
        # Remove the key prefix to get the bare feature name
        matched_feature <- sub(pattern, "", feature, ignore.case = TRUE)
        break
      }
    }

    # If no key matched, assign to default assay
    if (is.null(matched_assay)) {
      matched_assay <- default_assay
    }

    # Add feature to the appropriate assay's list
    if (is.null(features_by_assay[[matched_assay]])) {
      features_by_assay[[matched_assay]] <- list()
    }
    # Store mapping: bare feature name -> original feature string
    features_by_assay[[matched_assay]][[matched_feature]] <- feature
  }

  # Fetch data for each assay and combine results
  fetched_data_list <-
    lapply(
      names(features_by_assay),
      function(assay_name) {
        # Get the keyless feature names for this assay
        keyless_features <- names(features_by_assay[[assay_name]])
        # Get the original keyed feature names for this assay
        original_features <-
          unlist(
            features_by_assay[[assay_name]],
            use.names = FALSE
          )

        # Check if this assay's layer is BPCells-backed
        if (is_bpcells(object, assay = assay_name, layer = layer)) {
          # Use optimized BPCells direct matrix access
          assay_data <-
            bp_cells_fetch_features(
              object = object,
              features = keyless_features,
              assay = assay_name,
              layer = layer,
              cells = cells
            )
        } else {
          # Standard Seurat assay - delegate to fetch_data
          # Pass keyed feature names so fetch_data can resolve them properly
          assay_data <-
            SCUBA::fetch_data(
              object = object,
              vars = original_features,
              layer = layer,
              cells = cells
            )
        }

        # Rename columns to use original keyed feature names
        # (bp_cells_fetch_features returns bare names, fetch_data may vary)
        if (ncol(assay_data) > 0) {
          # Map current column names back to original feature names
          col_mapping <- setNames(original_features, keyless_features)
          new_colnames <- sapply(
            colnames(assay_data),
            function(col) {
              if (col %in% names(col_mapping)) {
                col_mapping[[col]]
              } else {
                col
              }
            },
            USE.NAMES = FALSE
          )
          colnames(assay_data) <- new_colnames
        }

        assay_data
      }
    )

  # Combine all assay data frames by columns
  if (length(fetched_data_list) == 0) {
    # No data fetched - return empty data.frame with correct row names
    combined_data <- data.frame(row.names = cells)
  } else if (length(fetched_data_list) == 1) {
    combined_data <- fetched_data_list[[1]]
  } else {
    # Combine multiple data frames column-wise
    combined_data <- do.call(cbind, fetched_data_list)
  }

  # Reorder columns to match the original order of features as passed by user
  if (ncol(combined_data) > 0) {
    # Match original feature order
    col_order <- match(features, colnames(combined_data))
    col_order <- col_order[!is.na(col_order)]
    combined_data <- combined_data[, col_order, drop = FALSE]
  }

  combined_data
}

#' @describeIn fetch_feature SingleCellExperiment objects. Wrapper for
#'   [fetch_data()].
#'
#' @keywords internal
#'
#' @export
fetch_feature.SingleCellExperiment <- function(
  object,
  features,
  assay = NULL,
  layer = NULL,
  cells = NULL
) {
  # SingleCellExperiment does not use assay parameter in the same way
  # Delegate to fetch_data
  SCUBA::fetch_data(
    object = object,
    vars = features,
    layer = layer,
    cells = cells
  )
}

#' @describeIn fetch_feature AnnDataR6 objects. Wrapper for [fetch_data()].
#'
#' @keywords internal
#'
#' @export
fetch_feature.AnnDataR6 <- function(
  object,
  features,
  assay = NULL,
  layer = NULL,
  cells = NULL
) {
  # AnnDataR6 does not use assay parameter in the same way
  # Delegate to fetch_data
  SCUBA::fetch_data(
    object = object,
    vars = features,
    layer = layer,
    cells = cells
  )
}
