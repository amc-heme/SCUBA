# Reference dataset used for testing and demonstration (anndata)

A reference dataset for accute myeloid leukemia was included in this
package for demonstration and testing. The data was originally published
in [Triana et al. 2021](https://doi.org/10.1038/s41590-021-01059-0). The
SCUBA authors downsampled the original Seurat object to use in the
package for automated testing, and converted it into other object
formats. The cell types provided by Triana et al. were also condensed
into 10 generalized cell types to facilitate demonstration of SCUBA
visualization capabilities. Details on the operations performed from the
original object are provided in [this
script](https://github.com/amc-heme/SCUBA_Manuscript/blob/main/Demo_Object_Generation.Rmd)
in the SCUBA manuscript repository.

## Usage

``` r
AML_h5ad()

AML_h5ad_backed()
```

## Source

The dataset was obtained from the
[Figshare](https://figshare.com/articles/dataset/Expression_of_197_surface_markers_and_462_mRNAs_in_15281_cells_from_blood_and_bone_marrow_from_a_young_healthy_donor/13398065/2).
For more information on the operations performed on the original object,
see the [SCUBA manuscript
repository](https://github.com/amc-heme/SCUBA_Manuscript/blob/main/Demo_Object_Generation.Rmd).

## Details

- `AML_h5ad()`: Loads the reference dataset in-memory.

- `AML_h5ad_backed()`: Loads the reference dataset on-disk (via
  `backed=r` in
  [`anndata::read_h5ad`](https://anndata.dynverse.org/reference/read_h5ad.html)).

`AML_h5ad()` and `AML_h5ad_backed()` are called as functions, with
parentheses. This is because the dataset is in Python, and can't be
loaded using the typical process of loading data included with R
packages.

## Examples

``` r
# These examples require a functional Python installation with prerequisite
# packages installed to work
# Please see our website for more details
# https://amc-heme.github.io/SCUBA/index.html#installation
# 
# The examples may take a while (about 10 seconds) to run the first 
# time they are executed in a session due to the time required to 
# initialize a Python environment. 
# R Studio does not display a spinner while these run, so the "run_examples"
# link may not appear to do anything until the examples are finished. 
AML_h5ad()
#> AnnData object with n_obs × n_vars = 250 × 462
#>     obs: 'nCount_RNA', 'nFeature_RNA', 'nCount_AB', 'nFeature_AB', 'nCount_BOTH', 'nFeature_BOTH', 'BOTH_snn_res.0.9', 'seurat_clusters', 'Prediction_Ind', 'BOTH_snn_res.1', 'ClusterID', 'Batch', 'x', 'y', 'x_mean', 'y_mean', 'cor', 'ct', 'prop', 'meandist', 'cDC', 'B.cells', 'Myelocytes', 'Erythroid', 'Megakaryocte', 'Ident', 'RNA_snn_res.0.4', 'condensed_cell_type'
#>     var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized'
#>     obsm: 'X_pca', 'X_umap', 'protein'
#>     layers: 'counts'

# Summary of metadata variables in object
meta_varnames(AML_h5ad())
#>  [1] "nCount_RNA"          "nFeature_RNA"        "nCount_AB"          
#>  [4] "nFeature_AB"         "nCount_BOTH"         "nFeature_BOTH"      
#>  [7] "BOTH_snn_res.0.9"    "seurat_clusters"     "Prediction_Ind"     
#> [10] "BOTH_snn_res.1"      "ClusterID"           "Batch"              
#> [13] "x"                   "y"                   "x_mean"             
#> [16] "y_mean"              "cor"                 "ct"                 
#> [19] "prop"                "meandist"            "cDC"                
#> [22] "B.cells"             "Myelocytes"          "Erythroid"          
#> [25] "Megakaryocte"        "Ident"               "RNA_snn_res.0.4"    
#> [28] "condensed_cell_type"
```
