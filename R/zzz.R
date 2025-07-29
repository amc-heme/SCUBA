md <- NULL

.onLoad <- function(libname, pkgname){
  # Establish Python package dependencies
  # Reticulate will automatically manage a Python environment with these 
  # packages, installing each if they are not already present
  reticulate::py_require("anndata>=0.11.4")
  reticulate::py_require("pandas>=2.0.0")
  reticulate::py_require("numpy")
  reticulate::py_require("scipy>=1.14.0")
  reticulate::py_require("mudata>=0.3.1")
  
  # Import Python packages used by Reticulate
  md <<- 
    reticulate::import("mudata", as = "md", convert = TRUE, delay_load = TRUE) 
}