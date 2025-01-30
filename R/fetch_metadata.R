#' Fetch metadata
#'
#' Returns object metadata for a specified set of cells.
#'
#' @param object a single-cell object. Currently, Seurat,
#' SingleCellExperiment, and anndata objects are supported.
#' @param vars metadata variables to pull from object. This must be defined,
#' unless "full_table" is set to \code{TRUE}.
#' @param cells cells to include in the returned metadata. If unspecified,
#' metadata will be returned for all cells in the object.
#' @param full_table if \code{TRUE}, return the entire metadata table (\code{FALSE} by
#' default)
#' @param return_class class of data returned. Set to "dataframe" by default to
#' return a data.frame, and may also be set to "vector" to yield a vector of
#' values. This is ignored if "full_table" is set to \code{TRUE}.
#'
#' @rdname fetch_metadata
#'
#' @export
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
        ". Currently supported classes: Seurat and SingleCellExperiment."
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
