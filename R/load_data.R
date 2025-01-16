#' Load Data 
#'
#' @param file_path 
#' @param backed
#'
#' @return obj
#' @export
#'
#' @examples
#' load_data("~/files/test.h5mu", backed="r+")
load_data <- 
  function(
    file_path,
    backed = FALSE
  ){
    #Need to add backed=False
    file_ext <- tools::file_ext(file_path)
    
    # Handle MuData files
    if (file_ext == "h5mu") {
      md <- reticulate::import("mudata", as = "md", convert = TRUE)
      
      obj <- md$read_h5mu(file_path, backed)
    }
    # Handle AnnData Objects
    else if (file_ext == "h5ad") {
      obj <- anndata::read_h5ad(file_path)
    }
    else if (file_ext == "rds") {
      obj <-  readRDS(file_path)
    } else if (file_ext == "h5Seurat") {
      obj <- SeuratDisk::LoadH5Seurat(file_path)
    }
    else {
      stop("Unsupported file format or directory structure: ", file_ext)
    }
    return(obj)
    
    
    #TODO: Add 
    
    
  }