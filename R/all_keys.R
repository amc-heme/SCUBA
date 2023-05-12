#' Get keys for all assays/reductions in an object
#'
#' Returns all keys in the object. For SingleCellExperiment objects, currently
#' returns only the names of all reductions in the object, as well as all
#' experiments.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param ... Currently unused.
#'
#' @rdname all_keys
#'
#' @export
all_keys <-
  function(
    object,
    ...
  ){
    UseMethod("all_keys")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
all_keys.default <-
  function(
    object
  ){
    warning(
      paste0(
        "all_keys does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
        )
      )
    }

#' @describeIn all_keys Seurat objects
#' @export
all_keys.Seurat <-
  function(
    object
  ){
    # Seurat objects: run SeuratObject.Key method
    Key(object)
  }

#' @describeIn all_keys SingleCellExperiment objects
#' @export
all_keys.SingleCellExperiment <-
  function(
    object
  ){
    # SingleCellExperiment objects: return names of experiments and reductions
    keys <-
      c(mainExpName(object),
      altExpNames(object),
      reducedDimNames(object)
      )

    # Add names to vector returned for consistency with SeuratObject.Key() method
    names(keys) <- keys

    keys
  }
