#' Access expression data from a SingleCellExperiment object
#'
#' This function is akin to the FetchData method for Seurat objects, and currently supports feature data only. Data from multiple assays (exeriments) can be selected by adding a prefix for the desired assay.
#'
#' @param object a SingleCellExperiment object.
#' @param features a character vector with one or more features for which expression data will be fetched. To select features from experiments other than the main experiment, use the name of the desired alternate experiment followed by an underscore. For example, to fetch a surface protein from an alternate experiment named ADT, use \code{ADT_{name of protein}}.
#' @param cells a character vector giving the cells for which to fetch data. If \code{NULL}, data will be returned for all cells in the object.
#' @param slot the slot (assay in SingleCellExperiment objects) for which to pull data. To view the list of available assays in your object, use \code{assayNames({your object})}.
#'
#' @return A data.frame object containing expression data for the feature(s) entered.
#' @importFrom SingleCellExperiment mainExpName altExps altExpNames
#' @importFrom SummarizedExperiment assays assayNames
#' @export
#'
feature_data <- function(object, features, cells = NULL, slot = NULL){
  # 1. Set default values
  # Slot (assay): defaults to the first assay stored in the object
  slot <- slot %||% assayNames(object)[1]
  # Cells: if NULL, use all cells in the object
  cells <- cells %||% colnames(object)

  # 2. For each feature, fetch data from object for the defined slot
  data <-
    lapply(
      features,
      function(feature, object, slot){
        # 2.1. Determine the experiment "key" of the feature passed
        # alphanumeric features to the left of an underscore will be treated as
        # the key, and will be compared against the main and alternate experiment
        # names.

        # Determine if the the entry is in "{experiment}_{feature}" format
        if (grepl("^[[:alnum:]]+_[[:alnum:]]+", feature)){
          # If so, extract experiment name from the assay
          exp <-
            gsub(
              pattern = '(^[[:alnum:]]+)_.*',
              replacement = "\\1",
              x = feature,
              fixed = FALSE
              )
          # Feature_name: the name of the feature with the experiment prefix
          # removed. Used when searching the indicated experiment for the feature.
          feature_name <-
            gsub(
              pattern = '^[[:alnum:]]+_(.*)',
              replacement = "\\1",
              x = feature,
              fixed = FALSE
              )
        } else {
          # If not, assume the intended experiment is the main experiment
          exp <- mainExpName(object)
          feature_name <- feature
        }

        # 2.2. Fetch data for feature from matrix
        # Determine if feature is in the main experiment or an alternate experiment
        if (exp == mainExpName(object)){
          # Determine if the feature name is in the expression matrix under
          # the defined slot
          if (feature_name %in% rownames(assays(object)[[slot]])){
            # Pull feature row and transpose to column
            assays(object)[[slot]][feature_name, , drop = FALSE] #|>
              # as.data.frame() |>
              # t()
          } else {
            # Display message to user if not found, and return nothing
            warning(
              glue::glue(
                "Feature {feature_name} not found in the experiment {mainExpName(object)}."
              )
            )

            NULL
          }
        } else if (exp %in% altExpNames(object)){
          # Load alternate experiment (also a SingleCellExperiment object)
          alt_exp_data <- altExps(object)[[exp]]

          # Determine if the feature name is in the expression matrix under
          # the defined slot
          if (feature_name %in% rownames(assays(alt_exp_data)[[slot]])){
            assays(alt_exp_data)[[slot]][feature_name, , drop = FALSE] #|>
              #as.data.frame() |>
              #t()
          } else {
            # Display message to user if not found, and return nothing
            warning(
              glue::glue(
                "Feature {feature_name} not found in the experiment {mainExpName(alt_exp_data)}."
              )
            )

            NULL
          }
        } else {
          warning(glue::glue("No experiment found matching the entered key {exp}_ (for {feature})."))

          NULL
        }
      },
      object,
      slot
    )

  #as.data.frame(data)
  data
}
