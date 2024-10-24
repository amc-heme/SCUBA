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

#' @describeIn all_keys SingleCellExperiment objects
#' @export
all_keys.AnnDataR6 <-
  function(
    object
  ){
    # Since anndata objects don't have defined locations for modalities and
    # reductions, downstream functions calling this method may not work
    # in the same way they do for Seurat and SingleCellExperiment objects.

    # Anndata objects: return names of experiments and reductions
    keys <-
      c("X",
        object$obsm_keys()
        )

    # Add names to vector returned for consistency with
    # SeuratObject.Key() method
    names(keys) <- keys

    keys
  }

#' @describeIn all_keys MuData objects
#' @export
all_keys.md._core.mudata.MuData <-
  function(
    object
  ) {
    
     keys <- 
       c("X",
         object$obsm_keys()
         )
     
     names(keys) <- keys
     
     keys
  }
#SCUBA::all_keys(SCUBA::AML_h5mu())