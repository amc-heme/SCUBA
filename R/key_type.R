#' Get keys for all assays/reductions in an object
#'
#' Given a key, key_type will determine whether the key corresponds to an
#' assay/experiment, or a reduction.
#'
#' @param object a single-cell object. Currently, Seurat and
#' SingleCellExperiment objects are supported.
#' @param key the key to assess.
#' @param ... Currently unused.
#'
#' @rdname key_type
#'
#' @export
key_type <-
  function(
    object,
    key,
    ...
  ){
    UseMethod("key_type")
  }

#' Function to display an error message when an unsupported object
#' type is detected
#'
#' @noRd
#' @export
key_type.default <-
  function(
    object,
    key
  ){
    warning(
      paste0(
        "key_type does not know how to handle object of class ",
        paste(class(object), collapse = ", "),
        ". Currently supported classes: Seurat and SingleCellExperiment."
      )
    )
  }

#' @describeIn key_type Seurat objects
#' @export
key_type.Seurat <-
  function(
    object,
    key
  ){
    # obj_ref: Object referred to by key. The name of the assay/reduction/etc
    # associated with a key (what you would use to pull the data via `[[`.
    # For "rna_", this would be "RNA" (the RNA assay). For "UMAP_", this would
    # be "umap".
    # (Code adapted from Seurat:::ExIPlot)
    obj_ref <- names(which(Key(object) == key))

    # Seurat objects: Test class of obj[[obj_ref]]
    if (inherits(x = object[[obj_ref]], what = 'Assay')){
      return("Assay")
    } else if (inherits(x = object[[obj_ref]], what = 'DimReduc')){
      return("Reduction")
    } else {
      return("Other")
    }
  }

#' @describeIn key_type SingleCellExperiment objects
#' @export
key_type.SingleCellExperiment <-
  function(
    object,
    key
  ){
    # SingleCellExperiment objects: search for key in experiments, reductions
    if (key %in% c(mainExpName(object), altExpNames(object))){
      # Return "Assay" for consistency with Seurat method
      return("Assay")
    } else if (key %in% reducedDimNames(object)){
      return("Reduction")
    } else {
      return("Other")
    }
  }
