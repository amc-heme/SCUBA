#' @param vars List of all variables to fetch, use keyword \dQuote{ident} to
#' pull identity classes
#' @param cells Cells to collect data for (default is all cells)
#' @param slot Slot to pull feature data for
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @rdname FetchData
#' @method FetchData Seurat
#' @export
#'
#' @concept data-access
#'
#' @examples
#' pc1 <- FetchData(object = pbmc_small, vars = 'PC_1')
#' head(x = pc1)
#' head(x = FetchData(object = pbmc_small, vars = c('groups', 'ident')))
#'
FetchData.Seurat <- function(object, vars, cells = NULL, slot = 'data', ...) {
  # Seurat FetchData function
  object <- UpdateSlots(object = object)
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  if (is.null(x = vars)) {
    df <- EmptyDF(n = length(x = cells))
    rownames(x = df) <- cells
    return(df)
  }

  # Get a list of all objects to search through and their keys
  # object.keys includes keys for each assay and projection
  object.keys <- Key(object = object)
  # Find all vars that are keyed
  # at this stage, keyed.vars will be the indexes of requested features that map
  # to each *key* in the object (each *assay or reduction*). If a key does not
  # have any variables mapped to it it will be an integer of length zero.
  keyed.vars <- lapply(
    X = object.keys,
    FUN = function(key) {
      if (length(x = key) == 0 || nchar(x = key) == 0) {
        return(integer(length = 0L))
      }
      return(grep(pattern = paste0('^', key), x = vars))
    }
  )

  keyed.vars <- Filter(f = length, x = keyed.vars)

  # Determine if features requested are part of spatial transcriptomics assays
  ret.spatial2 <-
    vapply(
      X = names(x = keyed.vars),
      FUN = function(x) {
        return(inherits(x = object[[x]], what = 'FOV'))
      },
      FUN.VALUE = logical(length = 1L)
    )

  if (any(ret.spatial2) && !all(ret.spatial2)) {
    warning(
      "Data returned from spatial coordinates are incompatible with other data, returning only spatial coordinates",
      call. = FALSE,
      immediate. = TRUE
    )
    keyed.vars <- keyed.vars[ret.spatial2]
  }

  # Generate list of data from variables mapped to keys
  # Data structure returned is a list of lists
  data.fetched <- lapply(
    # This iterates through each assay/reuction
    X = names(x = keyed.vars),
    FUN = function(x) {
      vars.use <- vars[keyed.vars[[x]]]
      key.use <- object.keys[x]
      # Get data for each var in the current assay
      data.return <-
        # Procedure for fetching data depends on the class of the current assay/reduction
        if (inherits(x = object[[x]], what = 'DimReduc')) {
          tryCatch(
            expr = FetchData(object = object[[x]], vars = vars.use, cells = cells),
            error = function(e) {
              return(NULL)
            }
          )
        } else if (inherits(x = object[[x]], what = 'Assay')) {
          vars.use <- gsub(pattern = paste0('^', key.use), replacement = '', x = vars.use)
          data.assay <- GetAssayData(
            object = object,
            slot = slot,
            assay = x
          )
          # Check to see if vars in vars.use are in the rownames of the assay first
          # (this keeps errors from occurring when vars that are not features from
          # this assay start with what appears to be the assay key, for example
          # if there is a metadata variable named 'rna_mitochondrial_percentage"
          # and an RNA assay)
          vars.use <- vars.use[vars.use %in% rownames(x = data.assay)]
          data.vars <- t(x = as.matrix(data.assay[vars.use, cells, drop = FALSE]))
          if (ncol(data.vars) > 0) {
            colnames(x = data.vars) <- paste0(key.use, vars.use)
          }
          data.vars
        } else if (inherits(x = object[[x]], what = 'FOV')) {
          vars.use <- gsub(pattern = paste0('^', key.use), replacement = '', x = vars.use)
          FetchData(object = object[[x]], vars = vars.use, cells = cells)
        } else if (inherits(x = object[[x]], what = 'SpatialImage')) {
          vars.unkeyed <- gsub(pattern = paste0('^', key.use), replacement = '', x = vars.use)
          names(x = vars.use) <- vars.unkeyed
          coords <- GetTissueCoordinates(object = object[[x]])[cells, vars.unkeyed, drop = FALSE]
          colnames(x = coords) <- vars.use[colnames(x = coords)]
          coords
        }

      data.return <- as.list(x = as.data.frame(x = data.return))
      return(data.return)
    }
  )

  # Removes the main list structure, but retains the list of expression values
  # returned for each variable found (via recursive = FALSE)
  data.fetched <-
    unlist(x = data.fetched, recursive = FALSE)

  if (any(ret.spatial2)) {
    return(as.data.frame(x = data.fetched))
  }

  # Pull vars from object metadata
  # Identify metadata variables
  # Test if variables that have not already been fetched are metadata variables
  meta.vars <-
    vars[vars %in% colnames(x = object[[]]) & !(vars %in% names(x = data.fetched))]

  # Add metadata variables requested to the list of data fetched
  # Data.frame extracted will be coerced to a list
  data.fetched <- c(data.fetched, object[[meta.vars]][cells, , drop = FALSE])

  # Check if any metadata variables entered are also in the default assay
  # (possible source of ambiguity)
  meta.default <-
    meta.vars[meta.vars %in% rownames(x = GetAssayData(object = object, slot = slot))]
  if (length(x = meta.default)) {
    warning(
      "The following variables were found in both object metadata and the default assay: ",
      paste0(meta.default, collapse = ", "),
      "\nReturning metadata; if you want the feature, please use the assay's key (eg. ",
      paste0(Key(object = object[[DefaultAssay(object = object)]]), meta.default[1]),
      ")",
      call. = FALSE
    )
  }

  # Pull vars that are in the default assay, and are not already fetched
  default.vars <-
    vars[
      vars %in% rownames(x = GetAssayData(object = object, slot = slot)) &
        !(vars %in% names(x = data.fetched))
      ]

  # Attempt to subset the matrix for `slot` in the default assay for the variable
  data.fetched <- c(
    data.fetched,
    tryCatch(
      expr =
        as.data.frame(
          x =
            t(
              x = as.matrix(
                x = GetAssayData(
                  object = object,
                  slot = slot
                )[default.vars, cells, drop = FALSE]
              )
            )
        ),
      error = function(...) {
        # If the variable is not present in the matrix, an error will be returned.
        # Return NULL in this case, which will leave data.fetched unchanged.
        return(NULL)
      }
    )
  )

  # Pull identities if "ident" is passed to `vars`
  if ('ident' %in% vars && !'ident' %in% colnames(x = object[[]])) {
    data.fetched[['ident']] <- Idents(object = object)[cells]
  }

  # Try to find ambiguous vars (those not yet found with above methods)
  fetched <- names(x = data.fetched)
  vars.missing <- setdiff(x = vars, y = fetched)

  if (length(x = vars.missing) > 0) {
    # Search in alternative assays for variables not found in the main assay
    # or metadata
    # vars.alt: an empty list with a length equal to the number of variables not
    # yet found
    vars.alt <-
      vector(mode = 'list', length = length(x = vars.missing))
    names(x = vars.alt) <- vars.missing

    # Loop constructs a list of the assays each remaining feature was found in
    for (assay in FilterObjects(object = object, classes.keep = 'Assay')) {
      # Store which features from vars.missing are in the current assay
      vars.assay <- Filter(
        f = function(x) {
          features.assay <-
            rownames(
              x = GetAssayData(
                object = object,
                assay = assay,
                slot = slot
                )
              )
          return(x %in% features.assay)
        },
        x = vars.missing
      )

      # For each variable found in this assay, append the assay name to the
      # feature's entry in vars.alt.
      for (var in vars.assay) {
        vars.alt[[var]] <-
          append(
            x = vars.alt[[var]],
            values = assay
            )
      }
    }

    # Message for variables found in multiple alternative assays
    vars.many <- names(x = Filter(
      f = function(x) {
        # (Those with more than one assay entry in the list)
        return(length(x = x) > 1)
      },
      x = vars.alt
    ))
    if (length(x = vars.many) > 0) {
      warning(
        "Found the following features in more than one assay, excluding the default. We will not include these in the final data frame: ",
        paste(vars.many, collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }

    # Filter vars.missing for vars that are *not* in one alternate assay
    # (either none or multiple)
    vars.missing <- names(x = Filter(
      f = function(x) {
        return(length(x = x) != 1)
      },
      x = vars.alt
    ))

    # Pull vars found in only one alternative assay and return message
    vars.alt <- Filter(
      f = function(x) {
        return(length(x = x) == 1)
      },
      x = vars.alt
    )

    for (var in names(x = vars.alt)) {
      assay <- vars.alt[[var]]
      warning(
        'Could not find ',
        var,
        ' in the default search locations, found in ',
        assay,
        ' assay instead',
        immediate. = TRUE,
        call. = FALSE
      )
      keyed.var <- paste0(Key(object = object[[assay]]), var)
      data.fetched[[keyed.var]] <- as.vector(
        x = GetAssayData(object = object, assay = assay, slot = slot)[var, cells]
      )

      # Edit vars to include the missing var with its key
      vars <- sub(
        pattern = paste0('^', var, '$'),
        replacement = keyed.var,
        x = vars
      )
    }
    fetched <- names(x = data.fetched)
  }

  # Name the vars not found in a warning (or error if none of the vars
  # entered were found)
  m2 <- if (length(x = vars.missing) > 10) {
    paste0(' (10 out of ', length(x = vars.missing), ' shown)')
  } else {
    ''
  }

  if (length(x = vars.missing) == length(x = vars)) {
    stop(
      "None of the requested variables were found",
      m2,
      ': ',
      paste(head(x = vars.missing, n = 10L), collapse = ', ')
    )
  } else if (length(x = vars.missing) > 0) {
    warning(
      "The following requested variables were not found",
      m2,
      ': ',
      paste(head(x = vars.missing, n = 10L), collapse = ', ')
    )
  }

  # Assembled fetched vars in a data frame
  data.fetched <- as.data.frame(
    x = data.fetched,
    row.names = cells,
    stringsAsFactors = FALSE
  )

  data.order <-
    na.omit(
      object =
        pmatch(
          x = vars,
          # `fetched`: vector of variables for which data was found
          table = fetched
        )
    )
  if (length(x = data.order) > 1) {
    data.fetched <- data.fetched[, data.order]
  }

  # Add column names to data.frame
  colnames(x = data.fetched) <- vars[vars %in% fetched]
  return(data.fetched)
}
