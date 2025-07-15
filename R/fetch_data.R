#' Fetch feature expression data, reduction coordinates, or metadata from single-cell objects
#' 
#' This function extends the behavior of SeuratObject's FetchData to other 
#' single-cell objects, allowing for expression data, metdadata, or reduction 
#' coordinates to be pulled using consistent syntax.
#' 
#' See our GitHub.io website for additional information and examples.
#'
#' @param object A single-cell object. Currently, Seurat, SingleCellExperiment, 
#' and anndata objects are supported. If a Seruat object is passed to this 
#' generic, the FetchData method from the SeruatObject method will be ran. 
#' For all other objects, methods that replicate the behavior of FetchData in 
#' that object will be ran.
#' @param vars A character vector with desired features, metadata variables, or
#' reduction dimensions to pull from the object. By default, features are 
#' returned from the default assay (or "experiment" in SingleCellExperiment 
#' objects, or "modality" in anndata objects). To pull feature expression data, 
#' the assay to pull data from should be defined using the "key" of the assay 
#' before the feature name. To determine the key that corresponds to the assay 
#' to pull data from, run `all_keys`. For more information, see 
#' [our user guide]().
#' @param layer For feature expression data, the layer to pull data from. 
#' Layers are referred to as "slots" in Seurat objects v4 and earlier, and 
#' "assays" in SingleCellExperiment objects.
#' @param cells A character vector of cells to include, as they are named in
#' the object (i.e. according to colNames(object)). If `NULL`, data will
#' be returned for all cells in the object.
#' @param ... Additional parameters, beyond the ones listed above, to be passed to S3 methods. This includes the following, all of which are documented in the parameter entries below:
#'  - fetch_data.Seurat: `slot` parameter
#' @param slot This parameter is added for backwards compatability with Seruat 
#' v4 and earlier. This was deprecated in Seurat version 5.0.0. *If you have Seruat v5.0.0 or later, this should not be used. This parameter should also not be used for SingleCellExperiment objects or anndata objects and will not work at all for these object classes. Use* `layer` *instead.*
#' 
#'
#' @returns A data.frame with the requested `vars` as columns and the cells as 
#' rows.
#'
#' @export
#'
#' @examples
#' # Feature expression data
#' fetch_data(AML_Seurat, vars = "rna_FLT3") |> str()
#' 
#' # Reduction coordinates
#' fetch_data(AML_Seurat, vars = c("UMAP_1", "UMAP_2")) |> str()
#' fetch_data(AML_Seurat, vars = c("PC_1", "PC_2", "PC_3")) |> str()
#' 
#' # Metadata
#' fetch_data(AML_Seurat, vars = c("condensed_cell_type", "Batch", "nCount_RNA")) |> str()
fetch_data <-
  function(
    object,
    vars = NULL,
    layer = NULL,
    cells = NULL,
    ...
  ){
    UseMethod("fetch_data")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
fetch_data.default <-
  function(
    object,
    vars = NULL,
    layer = NULL,
    cells = NULL,
    ...
  ){
    warning(
      paste0(
        "fetch_data does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn fetch_data Seurat objects. This will run FetchData from the 
#' SeuratObject package.
#' @export
fetch_data.Seurat <-
  function(
    object,
    vars = NULL,
    layer = NULL,
    cells = NULL,
    slot = NULL,
    ...
  ){
    # Check vars input
    # If more than 1000 features are requested, warn the user of potential 
    # performance issues
    if (!is.null(vars)){
      if (length(vars) >= 1000){
        warning(
          paste0(
            "A very large number of features was requested (", 
            length(vars), 
            " features). fetch_data is not intended to be used with feature ",
            "queries of this length. Data is returned in a dense format, so ",
            "the memory usage of the output may be very large. Also, this ",
            "query may take a while to complete."
            ),
          immediate. = TRUE
        )
      }
    }
    
    # This is a wrapper for SeuratObject::FetchData
    if (!is.null(slot)){
      # For backwards compatability with Seruat v4 and earlier, use `slot` 
      # parameter instead of layer.
      SeuratObject::FetchData(
        object = object,
        vars = vars,
        slot = slot,
        cells = cells
        )
    } else {
      # Seurat v5 and later
      SeuratObject::FetchData(
        object = object,
        vars = vars,
        layer = layer,
        cells = cells
        )
    }
  }

#' @describeIn fetch_data SingleCellExperiment objects
#'
#' @importFrom SingleCellExperiment mainExpName altExps altExpNames
#' @importFrom SummarizedExperiment assays assayNames
#'
#' @export
fetch_data.SingleCellExperiment <-
  function(
    object,
    vars,
    layer = NULL,
    cells = NULL,
    ...
  ){
    # 1. Check vars input ####
    # 1.1. Very large queries in vars ####
    # If more than 1000 features are requested, warn the user of potential 
    # performance issues
    if (length(vars) >= 1000){
      warning(
        paste0(
          "A very large number of features was requested (", 
          length(vars), 
          " features). fetch_data is not intended to be used with feature ",
          "queries of this length. Data is returned in a dense format, so ",
          "the memory usage of the output may be very large. Also, this ",
          "query may take a while to complete."
          ),
        call. = FALSE,
        immediate. = TRUE
      )
    }
    
    ## 1.2. Duplicate feature entries ####
    # Remove duplicates and warn the user that duplicates were entered
    if (any(duplicated(vars))){
      warning(
        paste0(
          "The following entries to `vars` are duplicates: '",
          paste(vars[duplicated(vars)], collapse = "', '"), 
          "'. Only one entry for each duplicated var will be returned."
          ),
        call. = FALSE,
        immediate. = TRUE
      )
      
      vars <- vars[!duplicated(vars)]
    }
    
    # 1. Set default values
    # Layer (assay): fill with default if null (via default_layer method)
    layer <- layer %||% default_layer(object)
    # Cells: if NULL, use all cells in the object
    cells <- cells %||% get_all_cells(object)
    
    # 2. Identify which experiments/reductions have keyed features requested
    # for them. Loop through experiments and reductions, instead of looping
    # through each var
    
    # Get a list of all "keys" in the object
    key_names <-
      c(mainExpName(object),
        altExpNames(object),
        reducedDimNames(object)
      )
    
    # Construct a list of experiments with the indices of vars that match
    # each experiment
    keyed_var_locations <-
      lapply(
        key_names,
        function(key){
          # grep returns the indices of matching vars
          grep(pattern = paste0('^', key), x = vars)
        }
      )
    
    names(keyed_var_locations) <- key_names
    
    # Subset list for experiments that have at least one matching var
    keyed_var_locations <-
      Filter(
        # Filter list for elements with any length
        f = length,
        x = keyed_var_locations
      )
    
    # 3. Loop through experiment and get data for the keyed vars in that experiment
    fetched_data <-
      lapply(
        names(keyed_var_locations),
        function(key){
          # Variables in current experiment/reduction
          key_vars <- vars[keyed_var_locations[[key]]]
          
          # Remove experiment key for feature retrieval
          keyless_vars <-
            gsub(
              pattern = paste0('^', key, '_(.*)'),
              replacement = "\\1",
              x = key_vars,
              fixed = FALSE
            )
          
          # Retrieve data
          if (key == mainExpName(object)){
            # For main experiment
            # Subset to variables that are included in the experiment,
            # to avoid errors
            keyless_vars <- keyless_vars[keyless_vars %in% rownames(object)]
            
            # Before pulling data, make sure the layer provided by the user
            # exists in the object. Throw an error if not
            if (!layer %in% names(assays(object))){
              stop(
                "Error for vars ",
                paste(keyless_vars, collapse = ", "),
                ": layer ",
                layer,
                " does not exist in the indicated experiment (",
                mainExpName(object),
                ")"
              )
            }
            
            data <-
              assays(object)[[layer]][keyless_vars, cells, drop = FALSE] |>
              # Must be a matrix for feature names to properly display as names in
              # the final list
              # (begins as a DelayedArray)
              as.matrix() |>
              t()
            
            # Add experiment key back in (assuming at least one variable was
            # found. An error will result if none are found at this point)
            if (length(keyless_vars) >= 1){
              colnames(data) <- paste0(key, "_", keyless_vars)
            }
          } else if (key %in% altExpNames(object)){
            # For alternate experimnent(s)
            # Switch to SingleCellExperiment object for the alternate experiment
            alt_sce <- altExps(object)[[key]]
            
            keyless_vars <- keyless_vars[keyless_vars %in% rownames(alt_sce)]
            
            # Before pulling data, make sure the layer provided by the user
            # exists in the object. Throw an error if not
            if (!layer %in% names(assays(alt_sce))){
              stop(
                "Error for vars",
                keyless_vars,
                ": layer ",
                layer,
                " does not exist in the indicated experiment (",
                mainExpName(alt_sce),
                ")"
              )
            }
            
            
            data <-
              assays(alt_sce)[[layer]][keyless_vars, cells, drop = FALSE] |>
              as.matrix() |>
              t()
            
            # Add experiment key back in (assuming at least one variable was
            # found. An error will result if none are found at this point)
            if (length(keyless_vars) >= 1){
              colnames(data) <- paste0(key, "_", keyless_vars)
            }
            
          } else if (key %in% reducedDimNames(object)){
            # For reductions
            # keyless_vars will be equal to the indices of the dimensions to
            # pull data from
            dims <- as.integer(keyless_vars)
            
            # Nonsensical dim inputs
            # It is possible to enter a dim that does not exist in the reduction
            # matrix, for example via a typo. These will cause an error 
            # To avoid this, dims not within the bounds of the matrix 
            # are filtered out
            # Upper bound defined using the second element of dims (number of
            # columns)
            dims <- dims[dims >= 1 & dims <= dim(reducedDims(object)[[key]])[2]]
            
            data <-
              reducedDims(object)[[key]][cells, dims]
          }
          
          # Return as a list
          data <- as.list(as.data.frame(data))
          data
        }
      )
    
    # Nested list is returned, condense to a list (only unlist at the top level)
    fetched_data <- unlist(fetched_data, recursive = FALSE)
    
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
        call. = FALSE,
        immediate. = TRUE
      )
    }
    
    # 5. Fetch data for vars in the main experiment that were not specified
    # with an assay key
    remaining_vars <- vars[!vars %in% names(fetched_data)]
    main_exp_vars <- remaining_vars[remaining_vars %in% rownames(object)]
    
    # Catch duplicate entries: it is possible to enter the same feature, 
    # with a key in one case and without a key in another case
    # for example, GAPDH and RNA_GAPDH. Only one feature should be returned in 
    # this case. 
    # To test, add the key of the main experiment to the variables to test 
    # if these variables were already fetched above
    keyed_main_exp_vars <- paste0(mainExpName(object), "_", main_exp_vars)
    # Construct relationship of keyed vars to the vars as entered for error
    # message reporting
    names(keyed_main_exp_vars) <- main_exp_vars
    
    duplicate_main_exp_vars <- 
      keyed_main_exp_vars[keyed_main_exp_vars %in% names(fetched_data)]
    
    if (length(duplicate_main_exp_vars) > 0){
      warning(
        paste0(
          "The entries to `vars` '", 
          paste(names(duplicate_main_exp_vars), collapse = "', '"),
          "' are the same as the entries '", 
          paste(duplicate_main_exp_vars, collapse = "', '"), 
          "'. Only one entry for each of these variables will be returned."
          ),
        call. = FALSE, 
        immediate. = TRUE
        )
      
      main_exp_vars <- 
        main_exp_vars[!main_exp_vars %in% names(duplicate_main_exp_vars)]
    }
    
    if (!layer %in% names(assays(object))){
      stop(
        "Error for vars",
        main_exp_vars,
        ": layer ",
        layer,
        " does not exist in the indicated experiment (",
        mainExpName(object),
        ")"
      )
    }
    
    # Pull vars from main experiment
    main_exp_data <-
      assays(object)[[layer]][main_exp_vars, cells, drop = FALSE] |>
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
    # If there were any entries found to be duplicated in 5, remove them from
    # the list of missing vars
    missing_vars <- 
      missing_vars[!missing_vars %in% names(duplicate_main_exp_vars)]
    
    
    if (length(missing_vars) > 0){
      # 6.1. Create a list to store the experiment(s) each missing feature is in
      # Empty list for storing data
      where_missing_vars <- vector(mode = 'list', length = length(missing_vars))
      names(where_missing_vars) <- missing_vars
      
      for (key in altExpNames(object)){
        alt_sce <- altExps(object)[[key]]
        
        # Determine which of the missing vars are in the current key., if any
        missing_in_exp <-
          missing_vars[missing_vars %in% rownames(assays(alt_sce)[[layer]])]
        
        # For each variable found in this experiment, append the assay name to the
        # feature's entry in missing_in_exp (each entry is a vector)
        for (var in missing_in_exp){
          where_missing_vars[[var]] <-
            append(
              where_missing_vars[[var]],
              values = key
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
          paste0(
            "The following features were found in more than one alternate ",
            "experiment. These features will not be included in the data ",
            "returned, since the query does not specify which experiment ",
            "to pull the features from: ",
            paste(vars_multi_exp, collapse = ', '),
            ". \n",
            "To include these features, please specify which experiment you ",
            "would like to pull data from using the experiment name and an ",
            "underscore (i.e. ",
            # Display an example with the experiment key added (using an key that
            # the example is certain to be in)
            paste0(where_missing_vars[[vars_multi_exp[1]]][1], "_", vars_multi_exp[1]),
            ").",
            call. = FALSE,
            immediate. = TRUE
          )
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
      
      # Subset missing vars for vars in one key
      where_missing_vars <-
        Filter(
          f = function(x) {
            length(x) == 1
          },
          where_missing_vars
        )
      
      for (var in names(where_missing_vars)){
        key <- where_missing_vars[[var]]
        # Load alternate experiment
        alt_sce <- altExps(object)[[key]]
        
        if (!layer %in% names(assays(alt_sce))){
          stop(
            "Error for var",
            var,
            ": layer ",
            layer,
            " does not exist in the indicated experiment (",
            mainExpName(alt_sce),
            ")"
            )
          }
        
        data <-
          assays(alt_sce)[[layer]][var, cells, drop = FALSE] |>
          # Only one var will be fetched at once in this case, so data can be added
          # as a vector to the list of fetched data
          as.vector()
        
        warning(
          paste0(
            var,
            " was passed to fetch_data without specifying the key of the ",
            "experiment to pull the feature from, and it is not present in ",
            "the main experiment. The feature was found in ", 
            sQuote(key), 
            " and successfully returned."
            ),
          call. = FALSE,
          immediate. = TRUE
        )
        
        # Add experiment key to var
        keyed_var <- paste0(key, "_", var)
        
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
    
    if (length(x = data_order) > 1) {
      fetched_data <- fetched_data[, data_order]
    }
    
    # Change column names to reflect `vars` that were fetched
    colnames(x = fetched_data) <- vars[vars %in% fetched_vars]
    
    fetched_data
  }

#' @describeIn fetch_data anndata Objects
#' @export
fetch_data.AnnDataR6 <-
  function(
    object,
    vars,
    layer = NULL,
    cells = NULL,
    ...
    ){
    library(reticulate)
    
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop(
        paste0(
          'Package "reticulate" must be installed to use this ',
          'function with anndata objects.'
          ),
        call. = FALSE
      )
    }
    
    # Establish Python package dependencies
    # Reticulate will automatically manage a Python environment with these 
    # packages, installing each if they are not already present
    py_require("anndata>=0.11.4")
    py_require("pandas>=2.0.0")
    py_require("numpy")
    py_require("scipy>=1.14.0")
    
    # Check vars input
    # If more than 1000 features are requested, warn the user of potential 
    # performance issues
    if (length(vars) >= 1000){
      warning(
        paste0(
          "A very large number of features was requested (", 
          length(vars), 
          " features). fetch_data is not intended to be used with feature ",
          "queries of this length. Data is returned in a dense format, so ",
          "the memory usage of the output may be very large. Also, this ",
          "query may take a while to complete."
          ),
        immediate. = TRUE
      )
    }
    
    # Source fetch_anndata python script
    python_path =
      system.file(
        "extdata",
        "Python",
        "fetch_data.py",
        package = "SCUBA"
        )
    
    py_objs <- reticulate::py_run_file(python_path)
    
    # Runs python fetch_anndata function
    py_objs$fetch_anndata(
      obj = object,
      fetch_vars = vars,
      cells = cells,
      layer = layer
    )
  }

#' @export
fetch_data.mudata._core.mudata.MuData <-
  function(
    object,
    vars,
    layer = NULL,
    cells = NULL,
    ...
  ){
    # mudata._core.mudata.MuData: possible class when loading 
    # Redirect to fetch_data.md._core.mudata.MuData method
    fetch_data.md._core.mudata.MuData(
      object = object,
      vars = vars,
      layer = layer,
      cells = cells
    )
  }

#' @describeIn fetch_data MuData Objects
#' @export
fetch_data.md._core.mudata.MuData <-
  function(
    object,
    vars,
    layer = NULL,
    cells = NULL,
    ...
  ){
    # MuData objects
    # Run Python fetch_data scripts
    library(reticulate)
    
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop(
        paste0(
          'Package "reticulate" must be installed to use this ',
          'function with anndata objects.'
        ),
        call. = FALSE
      )
    }
    
    # Establish Python package dependencies
    # Reticulate will automatically manage a Python environment with these 
    # packages, installing each if they are not already present
    py_require("anndata>=0.11.4")
    py_require("pandas>=2.0.0")
    py_require("numpy")
    py_require("scipy>=1.14.0")
    
    # Check vars input
    # If more than 1000 features are requested, warn the user of potential 
    # performance issues
    if (length(vars) >= 1000){
      warning(
        paste0(
          "A very large number of features was requested (", 
          length(vars), 
          " features). fetch_data is not intended to be used with feature ",
          "queries of this length. Data is returned in a dense format, so ",
          "the memory usage of the output may be very large. Also, this ",
          "query may take a while to complete."
        ),
        immediate. = TRUE
      )
    }
    
    # Source fetch_anndata python script 
    # (contains functions for anndata and MuData)
    python_path =
      system.file(
        "extdata",
        "Python",
        "fetch_data.py",
        package = "SCUBA"
      )
    
    py_objs <- reticulate::py_run_file(python_path)
    
    # Runs python fetch_anndata function
    py_objs$fetch_mudata(
      obj = object,
      fetch_vars = vars,
      cells = cells,
      layer = layer
    )
  }
