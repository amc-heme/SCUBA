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
  object.keys <- Key(object = object)
  # Find all vars that are keyed
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
  ret.spatial2 <- vapply(
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

  # Generate list of data
  data.fetched <- lapply(
    X = names(x = keyed.vars),
    FUN = function(x) {
      vars.use <- vars[keyed.vars[[x]]]
      key.use <- object.keys[x]
      data.return <- if (inherits(x = object[[x]], what = 'DimReduc')) {
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
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)
  if (any(ret.spatial2)) {
    return(as.data.frame(x = data.fetched))
  }
  # Pull vars from object metadata
  meta.vars <- vars[vars %in% colnames(x = object[[]]) & !(vars %in% names(x = data.fetched))]
  data.fetched <- c(data.fetched, object[[meta.vars]][cells, , drop = FALSE])
  meta.default <- meta.vars[meta.vars %in% rownames(x = GetAssayData(object = object, slot = slot))]
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
  # Pull vars from the default assay
  default.vars <- vars[vars %in% rownames(x = GetAssayData(object = object, slot = slot)) & !(vars %in% names(x = data.fetched))]
  data.fetched <- c(
    data.fetched,
    tryCatch(
      expr = as.data.frame(x = t(x = as.matrix(x = GetAssayData(
        object = object,
        slot = slot
      )[default.vars, cells, drop = FALSE]))),
      error = function(...) {
        return(NULL)
      }
    )
  )
  # Pull identities
  if ('ident' %in% vars && !'ident' %in% colnames(x = object[[]])) {
    data.fetched[['ident']] <- Idents(object = object)[cells]
  }
  # Try to find ambiguous vars
  fetched <- names(x = data.fetched)
  vars.missing <- setdiff(x = vars, y = fetched)
  if (length(x = vars.missing) > 0) {
    # Search for vars in alternative assays
    vars.alt <- vector(mode = 'list', length = length(x = vars.missing))
    names(x = vars.alt) <- vars.missing
    for (assay in FilterObjects(object = object, classes.keep = 'Assay')) {
      vars.assay <- Filter(
        f = function(x) {
          features.assay <- rownames(x = GetAssayData(
            object = object,
            assay = assay,
            slot = slot
          ))
          return(x %in% features.assay)
        },
        x = vars.missing
      )
      for (var in vars.assay) {
        vars.alt[[var]] <- append(x = vars.alt[[var]], values = assay)
      }
    }
    # Vars found in multiple alternative assays are truly ambiguous, will not pull
    vars.many <- names(x = Filter(
      f = function(x) {
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
    vars.missing <- names(x = Filter(
      f = function(x) {
        return(length(x = x) != 1)
      },
      x = vars.alt
    ))
    # Pull vars found in only one alternative assay
    # Key this var to highlight that it was found in an alternate assay
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
      vars <- sub(
        pattern = paste0('^', var, '$'),
        replacement = keyed.var,
        x = vars
      )
    }
    fetched <- names(x = data.fetched)
  }
  # Name the vars not found in a warning (or error if no vars found)
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
  data.order <- na.omit(object = pmatch(
    x = vars,
    table = fetched
  ))
  if (length(x = data.order) > 1) {
    data.fetched <- data.fetched[, data.order]
  }
  colnames(x = data.fetched) <- vars[vars %in% fetched]
  return(data.fetched)
}
