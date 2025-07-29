#' Fetch metadata from single-cell objects
#'
#' Returns object metadata for a specified set of cells.
#'
#' See our GitHub.io website for additional information and examples.
#'
#' @inheritParams object_param
#' @param vars metadata variables to pull from object. This must be defined,
#' unless "full_table" is set to `TRUE`.
#' @param cells cell IDs for which to pull metadata. If `NULL`, 
#' coordinates will be returned from all cells in the object. Cell IDs can be 
#' generated with [fetch_cells()].
#' @param full_table if `TRUE`, return the entire metadata table. This is `FALSE` by
#' default.
#' @param return_class class of data returned. Set to "dataframe" by default to
#' return a data.frame, and may also be set to "vector" to yield a vector of
#' values. This is ignored if "full_table" is set to `TRUE`.
#'
#' @rdname fetch_metadata
#'
#' @export
#' 
#' @examples
#' # Return several metadata variables as a data.frame
#' fetch_metadata(
#'   AML_Seurat, 
#'   vars = c("condensed_cell_type", "Batch", "nCount_RNA")
#'   ) |> str()
#'   
#' # Return data for a single metadata variable as a vector
#' fetch_metadata(
#'   AML_Seurat, 
#'   vars = "condensed_cell_type",
#'   return_class = "vector"
#'   ) |> str()
#' 
#' # Return all metadata 
#' fetch_metadata(
#'   AML_Seurat,
#'   full_table = TRUE
#'   ) |> str()
#' 
fetch_metadata <-
  function(
    object,
    vars = NULL,
    cells = NULL,
    full_table = FALSE,
    return_class = "dataframe"
  ){
    UseMethod("fetch_metadata")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
fetch_metadata.default <-
  function(
    object,
    vars = NULL,
    cells = NULL,
    full_table = FALSE,
    return_class = "dataframe"
  ){
    warning(
      paste0(
        "fetch_metadata does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat, SingleCellExperiment, anndata (AnnDataR6), and MuData (as loaded from SCUBA::load_data)."
      )
    )
  }

#' @describeIn fetch_metadata Seurat objects
#' @export
fetch_metadata.Seurat <-
  function(
    object,
    vars = NULL,
    cells = NULL,
    full_table = FALSE,
    return_class = "dataframe"
  ){
    # Check return_class for valid entries
    if (!return_class %in% c("dataframe", "vector")){
      stop('return_class must be either "dataframe" or "vector".')
    }

    # Check vars for valid entries (may be undefined only if full_table or FALSE)
    if (full_table == FALSE && is.null(vars)){
      stop("`vars` must be defined, unless `full_table` is TRUE.")
    }
    # Also warn if vars is defined when full_table is TRUE
    if (full_table == TRUE && !is.null(vars)){
      warning("`full_table` is TRUE, ignoring entries in `vars`.")
    }

    # Return full table if enabled
    if (full_table == TRUE){
      return(object@meta.data)
    }

    # Return an error if the variable name is not in the object metadata
    vars_present <- meta_varnames(object)
    # not_found: TRUE if a var exists in the object, FALSE if not
    not_found <- !vars %in% vars_present
    
    # Return an error if any or all variables entered are not found
    if (any(not_found)){
      if (all(not_found)){
        stop(
          paste0(
            "\nNone of the variables entered in `vars` were not found\n",
            "in the object. To view available entries for your object,\n",
            "run SCUBA::meta_varnames()."
            )
          )
      } else {
        # If some but not all variables entered are not found, report the 
        # variables that are not found
        stop(
          paste0(
            "\nThe following variables entered in `vars` were not found in \n",
            "the object: ",
            paste(vars[not_found]), ". \n\n",
            "To view available entries for your object, run SCUBA::meta_varnames()."
          )
        )
      }
    }
    
    # Cells: if undefined, pull data for all cells
    cells <- cells %||% get_all_cells(object)

    # Seurat objects: pull metadata with `[[`
    data <- object[[vars]][cells, ]

    # Use as.data.frame() if return_class is equal to "dataframe"
    # (if "vector", no additional operations necessary)
    if (return_class == "dataframe"){
      data <-
        as.data.frame(
          data
        )

      # Add cell names to rownames, and metadata names to colnames
      rownames(data) <- cells
      colnames(data) <- vars
    } else if (return_class == "vector"){
      # Return named vector with cell names
      names(data) <- cells
    }

    data
  }

#' @describeIn fetch_metadata SingleCellExperiment objects
#' @export
fetch_metadata.SingleCellExperiment <-
  function(
    object,
    vars = NULL,
    cells = NULL,
    full_table = FALSE,
    return_class = "dataframe"
  ){
    # Check return_class for valid entries
    if (!return_class %in% c("dataframe", "vector")){
      stop('return_class must be either "dataframe" or "vector".')
    }

    # Check vars for valid entries (may be undefined only if full_table or FALSE)
    if (full_table == FALSE && is.null(vars)){
      stop("`vars` must be defined, unless `full_table` is TRUE.")
    }
    # Also warn if vars is defined when full_table is TRUE
    if (full_table == TRUE && !is.null(vars)){
      warning("`full_table` is TRUE, ignoring entries in `vars`.")
    }

    # Return full table if enabled
    if (full_table == TRUE){
      return(colData(object))
    }
    
    # Return an error if the variable name is not in the object metadata
    vars_present <- meta_varnames(object)
    # not_found: TRUE if a var exists in the object, FALSE if not
    not_found <- !vars %in% vars_present
    
    # Return an error if any or all variables entered are not found
    if (any(not_found)){
      if (all(not_found)){
        stop(
          paste0(
            "\nNone of the variables entered in `vars` were not found\n",
            "in the object. To view available entries for your object,\n",
            "run SCUBA::meta_varnames()."
          )
        )
      } else {
        # If some but not all variables entered are not found, report the 
        # variables that are not found
        stop(
          paste0(
            "\nThe following variables entered in `vars` were not found in \n",
            "the object: ",
            paste(vars[not_found]), ". \n\n",
            "To view available entries for your object, run SCUBA::meta_varnames()."
          )
        )
      }
    }

    # Cells: if undefined, pull data for all cells
    cells <- cells %||% get_all_cells(object)

    # SingleCellExperiment objects: use colData
    data <-
      colData(object)[cells, vars]

    # Use as.data.frame() if return_class is equal to "dataframe"
    # (if "vector", no additional operations necessary)
    if (return_class == "dataframe"){
      data <-
        as.data.frame(
          data
        )

      # Add cell names to rownames, and metadata names to colnames
      rownames(data) <- cells
      colnames(data) <- vars
    } else if (return_class == "vector"){
      # Return named vector with cell names
      names(data) <- cells
    }

    data
  }



#' @describeIn fetch_metadata AnnDataR6 objects
#' @export
fetch_metadata.AnnDataR6 <-
  function(
    object,
    vars = NULL,
    cells = NULL,
    full_table = FALSE,
    return_class = "dataframe"
  ){
    # Check return_class for valid entries
    if (!return_class %in% c("dataframe", "vector")){
      stop('return_class must be either "dataframe" or "vector".')
    }

    # Check vars for valid entries (may be undefined only if full_table or FALSE)
    if (full_table == FALSE && is.null(vars)){
      stop("`vars` must be defined, unless `full_table` is TRUE.")
    }
    # Also warn if vars is defined when full_table is TRUE
    if (full_table == TRUE && !is.null(vars)){
      warning("`full_table` is TRUE, ignoring entries in `vars`.")
    }

    # Return full table if enabled
    if (full_table == TRUE){
      return(object$obs)
    }

    # Return an error if the variable name is not in the object metadata
    vars_present <- meta_varnames(object)
    # not_found: TRUE if a var exists in the object, FALSE if not
    not_found <- !vars %in% vars_present
    
    # Return an error if any or all variables entered were not found
    if (any(not_found)){
      if (all(not_found)){
        stop(
          paste0(
            "\nNone of the variables entered in `vars` were not found\n",
            "in the object. To view available entries for your object,\n",
            "run SCUBA::meta_varnames()."
          )
        )
      } else {
        # If some but not all variables entered were not found, report the 
        # variables that are not found
        stop(
          paste0(
            "\nThe following variables entered in `vars` were not found in \n",
            "the object: ",
            paste(vars[not_found]), ". \n\n",
            "To view available entries for your object, run SCUBA::meta_varnames()."
          )
        )
      }
    }
    
    # Cells: if undefined, pull data for all cells
    cells <- cells %||% get_all_cells(object)

    # Pull metadata
    # For anndata objects, use obs
    data <-
      object$obs[cells, vars]

    # Use as.data.frame() if return_class is equal to "dataframe"
    # (if "vector", no additional operations necessary)
    if (return_class == "dataframe"){
      data <-
        as.data.frame(
          data
        )

      # Add cell names to rownames, and metadata names to colnames
      rownames(data) <- cells
      colnames(data) <- vars
    } else if (return_class == "vector"){
      # Return named vector with cell names
      names(data) <- cells
    }

    data
  }


#' @describeIn fetch_metadata Mudata objects
#' @export
fetch_metadata.md._core.mudata.MuData <-
  function(
    object,
    vars = NULL,
    cells = NULL,
    full_table = FALSE,
    return_class = "dataframe"
  ){
    # MuData objects
    # Run fetch_metadata_mudata Python function from fetch_data Python script
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
    reticulate::py_require("anndata>=0.11.4")
    reticulate::py_require("pandas>=2.0.0")
    reticulate::py_require("numpy")
    reticulate::py_require("scipy>=1.14.0")
    reticulate::py_require("mudata>=0.3.1")
    
    # Source fetch_data python script 
    # (contains functions for anndata and MuData)
    python_path =
      system.file(
        "extdata",
        "Python",
        "fetch_data.py",
        package = "SCUBA"
      )
    
    py_objs <- reticulate::py_run_file(python_path)
    
    # Error handling: vector return class is not possible when more than 
    # one variable is requested
    if (return_class == "vector" & length(vars) > 1){
      stop(
        'If more than one variable is requested via `vars`, ',
        '`return_class` must be "dataframe".'
        )
    }
    
    # Check vars for valid entries (may be undefined only if full_table or FALSE)
    if (full_table == FALSE && is.null(vars)){
      stop("`vars` must be defined, unless `full_table` is TRUE.")
    }
    # Also warn if vars is defined when full_table is TRUE
    if (full_table == TRUE && !is.null(vars)){
      warning("`full_table` is TRUE, ignoring entries in `vars`.")
    }
    
    # When full_table is TRUE,
    # Pull_obs, return main obs table
    if (full_table == TRUE){
      object$pull_obs()
      
      table <- object$obs
      
      # Transform separator for modality-metadata combinations for column names
      # from Mudata format (":") to "_" for consistency with the format returned
      # by other SCUBA methods
      colnames(table) <-
        sapply(
          colnames(table),
          function(colname){
            key_match <- 
              SCUBA::match_key(
                colname, 
                keys = object$mod_names, 
                sep = ":"
                )
            
            if (!is.null(key_match$key) & !is.null(key_match$suffix)){
              paste0(key_match$key, "_", key_match$suffix)
            } else {
              colname
            }
          }
        )
      
      # Scrub "pandas.index" attribute from R data.frame
      table <- remove_pandas_index(table)
      
      return(table)
    }
    
    # Fetch metadata for requested vars
    data <- 
      py_objs$fetch_metadata_mudata(
        obj = object,
        meta_vars = vars,
        cells = cells
        )
    
    if (return_class == "dataframe"){
      # For data.frame returns, scrub "pandas.index" attribute from 
      # R data.frame, and replace with rownames if they do not already exist
      data <- remove_pandas_index(data)
      
      return(data)
    } else if (return_class == "vector"){
      # Return as a vector: extract column from one-column 
      # dataframe and use row names as vector names
      data_vec <- data[[1]]
      names(data_vec) <- rownames(data)
      
      return(data_vec)
    }
  }

#' @export
fetch_metadata.mudata._core.mudata.MuData <-
  function(
    object,
    vars = NULL,
    cells = NULL,
    full_table = FALSE,
    return_class = "dataframe"
  ){
    # mudata._core.mudata.MuData: possible class when loading 
    # Redirect to fetch_data.md._core.mudata.MuData method
    fetch_metadata.md._core.mudata.MuData(
      object = object,
      vars = vars,
      cells = cells,
      full_table = full_table,
      return_class = return_class
    )
  }