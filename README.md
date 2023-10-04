# SCUBA

###  *S*ingle *C*ell *U*nified *B*ack end *A*PI 

SCUBA is a unified data accession interface for single cell data formats. The package pulls data from Seurat, SingleCellExeriment, and Anndata objects for use in downstream plotting and analysis. Plotting functions that produce Seruat-style plots from all three object formats are also included.

### Pre-requisites

SCUBA relies on the [reticulate](https://rstudio.github.io/reticulate/) package to operate on AnnData objects. Please make sure you follow the instructions on the reticulate page before moving on.  

It is recommended to use conda environments to manage your python version that is being used by reticulate. A conda environment can be used on the current R session by using `reticulate::use_condaenv()`. If you want to make it automatically load on startup you can select your conda environment in RStudio > Options > Python. 

You must have the following packages installed and available in your python environment:
- Numpy
- Scipy
- Pandas
- AnnData


### Working with different objects in SCUBA
The following demonstrates how to use SCUBA to access data in various formats.

There are three primary data accession methods: `FetchData`, `fetch_metadata` and `fetch_reduction`. 

There are two object exploration methods: `meta_varnames` and `unique_values`


 
#### Priary Data Accession Methods

FetchData for Seurat objects 

```
FetchData(
      AML_Seurat,
      slot = "data",
      vars =
        c("ab_CD117-AB",
          "ab_CD123-AB",
          "ab_CD11c-AB",
          "rna_GAPDH",
          "rna_MEIS1",
          # Reductions
          "UMAP_1",
          "UMAP_2",
          # Nonexistent features
          "ab_CD900",
          # Metadata
          "nCount_RNA",
          "nFeature_RNA",
          # "Ambiguous" feature not in RNA assay
          "CD11a-AB"
        )
    )
```

For SingleCellExperiment

```
sce <- AML_SCE()
FetchData(
      sce,
      slot = "logcounts",
      vars =
        c("AB_CD117-AB",
          "AB_CD123-AB",
          "AB_CD11c-AB",
          "RNA_GAPDH",
          "RNA_MEIS1",
          # Reductions
          "UMAP_1",
          "UMAP_2",
          # Nonexistent features
          "AB_CD900",
          # Metadata
          "nCount_RNA",
          "nFeature_RNA",
          # "Ambiguous" feature not in RNA assay
          "CD11a-AB"
        )
    )
```

AnnData:

```
AML_anndata <- AML_h5ad()
FetchData(
      AML_anndata,
      vars =
        c("protein_CD123-AB",
          "protein_CD117-AB",
          "X_MEIS1",
          "X_GAPDH",
          # Reductions
          "X_umap_1",
          "X_umap_2",
          # Metadata
          "nCount_RNA",
          "nFeature_RNA"
        )
    )
```


### Citation
Insert paper title, and preprint link when available here. 

### Problems
If any issues arise please use Github issues on this repository. 
