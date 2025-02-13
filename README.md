# SCUBA

SCUBA (*S*ingle *C*ell *U*nified *B*ack end *A*PI) is a unified data accession interface for single-cell object classes. The package streamlines R data analysis for Seurat, SingleCellExperiment, and anndata objects by providing a consistent interface for data access, exploration, and visualization.

## Pre-requisites

SCUBA relies on the [reticulate](https://rstudio.github.io/reticulate/) package to operate on AnnData objects. Please make sure you follow the instructions on the reticulate page before moving on.  

It is recommended to use conda environments to manage your python version that is being used by reticulate. A conda environment can be used on the current R session by using `reticulate::use_condaenv()`. If you want to make it automatically load on startup you can select your conda environment in RStudio > Options > Python. 

You must have the following packages installed and available in your Python environment:

- Numpy
- Scipy
- Pandas
- AnnData

## Installation

Run the command below to install SCUBA. BiocManager is used to automatically install Bioconductor dependencies (SCUBA is not a Bioconductor package).

If you plan to use SCUBA with anndata objects, use `dependencies = TRUE`. If you only plan to use SCUBA with Seurat and SingleCellExperiment objects, use `dependencies = FALSE`.

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Set dependencies to FALSE if you do not plan to use anndata objects
BiocManager::install("amc-heme/SCUBA", dependencies = TRUE)
```

## Citation
Showers,W.M., Desai,J., Engel,K.L., Smith,C., Jordan,C.T. and Gillen,A.E. (2024) SCUBA implements a storage format-agnostic API for single-cell data access in R. [10.12688/f1000research.154675.1](https://doi.org/10.12688/f1000research.154675.1).

## Problems
If any issues arise please file a Github issue on this repository. 

<!-- badges: start -->
  [![R-CMD-check](https://github.com/amc-heme/SCUBA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/amc-heme/SCUBA/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
