#' AML Reference Dataset
#'
#' @description
#' An AML single cell reference dataset produced by \href{https://doi.org/10.1038/s41590-021-01059-0}{Triana et al. 2021}. This is a small subset of the original  data that is intended for use with tests in this package.
#'
#' @description
#' The dataset includes a pheresis and a bone marrow sample from a single young healthy donor. A panel of 462 genes is included in the dataset, along with 197 surface protein markers. A UMAP projection based on gene expression was computed from the original Seurat object downloaded from Figshare (see source below). The cell types provided by Triana et al. were condensed into 10 generalized cell types, and the object was randomly downsampled to include 25 cells from each cell type.
#'
#' @format
#' A Seurat object with 659 features and 250 cells in 2 assays. Included assays:
#' \describe{
#'   \item{RNA}{mRNA expression data for 462 genes used as markers for hematopoietic stem and progenitor (HSPC) cell types at varying stages of differentiation.}
#'   \item{AB}{Surface protein expression data for 197 surface protein markers associated with HSPCs.}
#' }
#'
#' @source <https://figshare.com/articles/dataset/Expression_of_197_surface_markers_and_462_mRNAs_in_15281_cells_from_blood_and_bone_marrow_from_a_young_healthy_donor/13398065>
"AML_Seurat"

#' AML Reference Dataset (SingleCellExperiment Format)
#'
#' @inherit AML_Seurat description source
#'
#' @usage AML_SCE()
#'
#' @details
#' The function described here is used to load the SingleCellExperiment object using the \code{HDF5Array} package. Unlike `AML_Seurat`, this dataset should be called as a function.
#'
#' @format
#' A SingleCellExperiment object with 659 features and 250 cells in 2 SummarizedExperiment objects. Included experiments:
#' \describe{
#'   \item{RNA (main experiment)}{mRNA expression data for 462 genes used as markers for hamatopoietic stem and progenitor (HSPC) cell types at varying stages of differentiation.}
#'   \item{AB}{Surface protein expression data for 197 surface protein markers associated with HSPCs.}
#' }
#'
#' @import HDF5Array
#'
#' @export
AML_SCE <- function(){
  HDF5Array::loadHDF5SummarizedExperiment(
    dir = system.file("extdata", "AML_sce", package = "SCUBA")
  )
}


#' AnnData object for testing
#'
#' @usage AML_h5ad()
#'
#' @details
#' The function described here is used to load the AnnData R6 object, it should be called as a function.
#'
#' @format
#' An AnnDataR6 object with 659 features and 250 cells. Included experiments:
#' \describe{
#'   \item{RNA (main experiment)}{mRNA expression data for 462 genes used as markers for hamatopoietic stem and progenitor (HSPC) cell types at varying stages of differentiation.}
#'   \item{AB}{Surface protein expression data for 197 surface protein markers associated with HSPCs.}
#' }
#'
#' @export
#'
AML_h5ad <- function(){
  anndata::read_h5ad(
    system.file("extdata", "AML_h5ad.h5ad", package = "SCUBA")
  )
}
