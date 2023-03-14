#' FetchData Equivalent for SingleCellExperiment Objects
#'
#' The SingleCellExperiment equivalent of the SeuratObject FetchData method.
#'
#' @param object A SingleCellExperiment object.
#' @param vars A character vector with desired features or metadata variables
#' to pull from the object. To include features from an experiment other than
#' the main experiment, use the name of the experiment as a prefix (i.e. AB_CD4
#' for a feature in the experiment "AB" named "CD4".)
#' @param slot The assay (equivalent of slot in Seruat objects) to pull data
#' from. To view the list of available assays in your object, use
#' \code{assayNames({your object})}.
#' @param cells A character vector of cells to include, as they are named in
#' the object (i.e. according to colNames(object)). If \code{NULL}, data will
#' be returned for all cells in the object.
#'
#' @return A data.frame object containing the requested expression
#' data or metadata.
#'
#' @importFrom SingleCellExperiment mainExpName altExps altExpNames
#' @importFrom SummarizedExperiment assays assayNames
#'
#' @export
#'
#' @method FetchData SingleCellExperiment
#'
FetchData.SingleCellExperiment <-
  function(
    object,
    vars,
    slot = NULL,
    cells = NULL
    ){
    # 1. Set default values
    # Slot (assay): defaults to the first assay stored in the object
    slot <- slot %||% assayNames(object)[1]
    # Cells: if NULL, use all cells in the object
    cells <- cells %||% colnames(object)

    # 2. Identify which experiments have keyed features requested for them
    # Experiments will be looped through, instead of looping through each var to
    # find an associated experiment

    # Get a list of all experiment "keys" in the object
    exp_names <- c(mainExpName(object), altExpNames(object))

    # Construct a list of experiments with the indices of vars that match
    # each experiment
    keyed_vars <-
      lapply(
        exp_names,
        function(exp){
          # grep returns the indices of matching vars
          grep(pattern = paste0('^', exp), x = vars)
        }
      )

    names(keyed_vars) <- exp_names

    # Subset list for experiments that have at least one matching var
    keyed_vars <-
      Filter(
        # Filter list for elements with any length
        f = length,
        x = keyed_vars
      )

    # 3. Loop through experiment and get data for the keyed vars in that experiment
    fetched_data <-
      lapply(
        names(keyed_vars),
        function(exp){
          print(exp)
          # Variables in current experiment
          exp_vars <- vars[keyed_vars[[exp]]]
          print(exp_vars)

          # Remove experiment key for feature retrieval
          keyless_vars <-
            gsub(
              pattern = paste0('^', exp, '_(.*)'),
              replacement = "\\1",
              x = exp_vars,
              fixed = FALSE
            )
          print(keyless_vars)

          # Retrieve data
          if (exp == mainExpName(object)){
            # For main experiment
            # Subset to variables that are included in the experiment, to avoid errors
            keyless_vars <- keyless_vars[keyless_vars %in% rownames(object)]
            data <-
              assays(object)[[slot]][keyless_vars, cells, drop = FALSE] |>
              # Must be a matrix for feature names to properly display as names in
              # the final list
              # (begins as a DelayedArray)
              as.matrix() |>
              t()

            # Add experiment key back in
            colnames(data) <- paste0(exp, "_", keyless_vars)
          } else {
            # For alternate experimnent(s)
            # Switch to SingleCellExperiment object for the alternate experiment
            alt_sce <- altExps(object)[[exp]]

            keyless_vars <- keyless_vars[keyless_vars %in% rownames(alt_sce)]

            data <-
              assays(alt_sce)[[slot]][keyless_vars, cells, drop = FALSE] |>
              as.matrix() |>
              t()

            print("as.list on data")
            print(str(data))
            print(str(as.list(as.data.frame(data))))

            # Add experiment key back in
            colnames(data) <- paste0(exp, "_", keyless_vars)
          }

          # Return as a list
          data <- as.list(as.data.frame(data))
          data
        }
      )

    # Nested list is returned, condense to a list (only unlist at the top level)
    fetched_data <- unlist(fetched_data, recursive = FALSE)

    str(fetched_data)

    # 4. Fetch metadata variables
    # Identify metadata variables
    remaining_vars <- vars[!vars %in% names(fetched_data)]
    metadata_vars <- remaining_vars[remaining_vars %in% names(colData(object))]

    # Fetch metadata vars and append to fetched_data
    fetched_data <-
      c(
        fetched_data,
        # SeuratObjects return metadata as a data.frame by default. This must be
        # done manually to properly append to fetched_data (data.frames can be
        # coerced to a list)
        as.data.frame(colData(object))[cells, metadata_vars, drop = FALSE]
      )

    # Handle ambiguous case of metadata variables also existing in the
    # main experiment (warn users that only metadata will be returned)
    ambiguous_meta_vars <- metadata_vars[metadata_vars %in% rownames(object)]
    if (length(ambiguous_meta_vars) > 0){
      warning(
        "The following variables were found in both object metadata and the main experiment: ",
        paste0(ambiguous_meta_vars, collapse = ", "),
        '\nOnly the metadata will be returned. To get feature data from the main experiment, please add the "key" of the main experiment to the feature (eg. ',
        paste0(mainExpName(object), "_", ambiguous_meta_vars[1]),
        ")",
        call. = FALSE
      )
    }

    # 5. Fetch data for vars in the main experiment that were not specified
    # with an assay key
    remaining_vars <- vars[!vars %in% names(fetched_data)]
    main_exp_vars <- remaining_vars[remaining_vars %in% rownames(object)]

    # Pull vars from main experiment
    main_exp_data <-
      assays(object)[[slot]][main_exp_vars, cells, drop = FALSE] |>
      as.matrix() |>
      t()

    # Explicitly specify feature names as column names
    colnames(main_exp_data) <- main_exp_vars

    main_exp_data <-
      main_exp_data |>
      as.data.frame() |>
      as.list()

    # Append data to the list of fetched data
    fetched_data <-
      c(
        fetched_data,
        main_exp_data
      )

    # 6. Handle variables that have not yet been fetched
    missing_vars <- vars[!vars %in% names(fetched_data)]

    print("Missing vars")
    print(missing_vars)

    if (length(missing_vars) > 0){
      # 6.1. Create a list to store the experiment(s) each missing feature is in
      # Empty list for storing data
      where_missing_vars <- vector(mode = 'list', length = length(missing_vars))
      names(where_missing_vars) <- missing_vars

      for (exp in altExpNames(object)){
        alt_sce <- altExps(object)[[exp]]

        # Determine which of the missing vars are in the current exp., if any
        missing_in_exp <-
          missing_vars[missing_vars %in% rownames(assays(alt_sce)[[slot]])]

        print(missing_in_exp)

        # For each variable found in this experiment, append the assay name to the
        # feature's entry in missing_in_exp (each entry is a vector)
        for (var in missing_in_exp){
          where_missing_vars[[var]] <-
            append(
              where_missing_vars[[var]],
              values = exp
            )
        }
      }

      # 6.2. Warn user if there are vars in multiple experiments. Do not pull data
      # in this case
      vars_multi_exp <-
        # Subset list from 6.1. for vars with more than one associated experiment
        Filter(
          f = function(x) {
            length(x) > 1
          },
          where_missing_vars
        ) |>
        names()

      if (length(vars_multi_exp) > 0){
        warning(
          "The following features were found in more than one alternate experiment. These features will not be included in the data returned: ",
          paste(vars_multi_exp, collapse = ', '),
          ". \n",
          "To include these features, please specify which experiment you would like to pull data from using the experiment name and an underscore (i.e. ",
          # Display an example with the experiment key added (using an exp that
          # the example is certain to be in)
          paste0(where_missing_vars[[vars_multi_exp[1]]][1], "_", vars_multi_exp[1]),
          ").",
          call. = FALSE,
          immediate. = TRUE
        )
      }

      # 6.3. Pull data for missing variables found in one alternate experiment
      # Update list of missing vars to exclude vars in one experiment
      missing_vars <-
        Filter(
          f = function(x){
            length(x) != 1
          },
          where_missing_vars
        ) |>
        names()

      # Subset missing vars for vars in one exp
      where_missing_vars <-
        Filter(
          f = function(x) {
            length(x) == 1
          },
          where_missing_vars
        )

      #
      for (var in names(where_missing_vars)){
        exp <- where_missing_vars[[var]]
        # Load alternate experiment
        alt_sce <- altExps(object)[[exp]]

        data <-
          assays(alt_sce)[[slot]][var, cells, drop = FALSE] |>
          # Only one var will be fetched at once in this case, so data can be added
          # as a vector to the list of fetched data
          as.vector()

        # Add experiment key to var
        keyed_var <- paste0(exp, "_", var)

        fetched_data[[keyed_var]] <- data

        # Edit vars to use the keyed var
        vars <-
          sub(
            pattern = paste0('^', var, '$'),
            replacement = keyed_var,
            x = vars
          )
      }
    }

    # 7. User warnings/errors
    ten_plus_message <-
      if (length(x = missing_vars) > 10) {
        paste0(' (10 out of ', length(x = missing_vars), ' shown)')
      } else {
        ''
      }

    if (length(x = missing_vars) == length(x = vars)) {
      stop(
        "None of the requested variables were found",
        ten_plus_message,
        ': ',
        paste(head(x = missing_vars, n = 10L), collapse = ', ')
      )
    } else if (length(x = missing_vars) > 0) {
      warning(
        "The following requested variables were not found",
        ten_plus_message,
        ': ',
        paste(head(x = missing_vars, n = 10L), collapse = ', ')
      )
    }

    # 8. Construct data.frame and return to user
    # Store names of vars fetched for downstream operations
    fetched_vars <- names(fetched_data)
    # Convert fetched_data to data.frame
    fetched_data <-
      as.data.frame(
        fetched_data,
        row.names = cells,
        stringsAsFactors = FALSE
      )

    # Re-order vars to reflect the order entered, instead of the order fetched
    data_order <-
      na.omit(
        object =
          pmatch(
            x = vars,
            table = fetched_vars
          )
      )

    print("data_order")
    print(data_order)

    if (length(x = data_order) > 1) {
      fetched_data <- fetched_data[, data_order]
    }

    # Change column names to reflect `vars` that were fetched
    colnames(x = fetched_data) <- vars[vars %in% fetched_vars]

    fetched_data
  }
