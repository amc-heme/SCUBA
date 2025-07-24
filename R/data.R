#' Reference dataset used for testing and demonstration (Seurat)
#'
#' A reference dataset for accute myeloid leukemia was included in this package for demonstration and testing. The data was originally published in [Triana et al. 2021](https://doi.org/10.1038/s41590-021-01059-0). The SCUBA authors downsampled the original Seurat object to use in the package for automated testing, and converted it into other object formats. The cell types provided by Triana et al. were also condensed into 10 generalized cell types to facilitate demonstration of SCUBA visualization capabilities. Details on the operations performed from the original object are provided in [this script](https://github.com/amc-heme/SCUBA_Manuscript/blob/main/Demo_Object_Generation.Rmd) in the SCUBA manuscript repository.
#'
#' @source The dataset was obtained from the [Figshare](https://figshare.com/articles/dataset/Expression_of_197_surface_markers_and_462_mRNAs_in_15281_cells_from_blood_and_bone_marrow_from_a_young_healthy_donor/13398065/2). For more information on the operations performed on the original object, see the [SCUBA manuscript repository](https://github.com/amc-heme/SCUBA_Manuscript/blob/main/Demo_Object_Generation.Rmd). 
#' 
#' @examples
#' # Object summary
#' AML_Seurat
#' 
#' # Summary of metadata variables in object
#' meta_varnames(AML_Seurat)
"AML_Seurat"

#' Reference dataset used for testing and demonstration (SingleCellExperiment)
#'
#' @inherit AML_Seurat description source
#'
#' @details
#' Contrary to convention for loading data, `AML_SCE()` is called as a function, with parentheses. This is because the
#' dataset is loaded with `HDF5Array::loadHDF5SummarizedExperiment()` instead of the typical process of loading data from R packages. The dataset was saved using the `HDF5Array` package to test support of SingleCellExperiment objects supporting HDF5 storage saved using this package.
#'
#' @import HDF5Array
#'
#' @export
#' 
#' @usage AML_SCE()
#' 
#' @examples
#' # Object summary
#' AML_SCE()
#' 
#' # Summary of metadata variables in object
#' meta_varnames(AML_SCE())
AML_SCE <- function(){
  HDF5Array::loadHDF5SummarizedExperiment(
    dir = system.file("extdata", "AML_sce", package = "SCUBA")
  )
}


#' Reference dataset used for testing and demonstration (anndata)
#'
#' @inherit AML_Seurat description source
#'
#' @details
#' - `AML_h5ad()`: Loads the reference dataset in-memory.
#' - `AML_h5ad_backed()`: Loads the reference dataset on-disk (via `backed=r` in `anndata::read_h5ad`).
#' 
#' `AML_h5ad()` and `AML_h5ad_backed()` are called as functions, with parentheses. 
#' This is because the dataset is in Python, and can't be loaded using the 
#' typical process of loading data included with R packages.
#'
#' @rdname anndata_example
#'
#' @export
#'
#' @usage 
#' AML_h5ad()
#'
#' @examples
#' # These examples require a functional Python installation with prerequisite
#' # packages installed to work
#' # Please see our website for more details
#' # https://amc-heme.github.io/SCUBA/index.html#installation
#' # 
#' # The examples may take a while (about 10 seconds) to run the first 
#' # time they are executed in a session due to the time required to 
#' # initialize a Python environment. 
#' # R Studio does not display a spinner while these run, so the "run_examples"
#' # link may not appear to do anything until the examples are finished. 
#' AML_h5ad()
#' 
#' # Summary of metadata variables in object
#' meta_varnames(AML_h5ad())
AML_h5ad <- function(){
  anndata::read_h5ad(
    system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA")
  )
}

#' @rdname anndata_example
#' @usage AML_h5ad_backed()
#' @export
AML_h5ad_backed <- function(){
    anndata::read_h5ad(
        system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA"),
        backed = "r"
    )
}

#' @export
AML_h5mu <- function(){
  path <- system.file("extdata", "AML_mudata.h5mu", package = "SCUBA")
  
  py_require("mudata")
  
  md <- reticulate::import("mudata", as = "md", convert = TRUE)
  
  md$set_options(pull_on_update = FALSE)
  
  md$read_h5mu(path)
}
