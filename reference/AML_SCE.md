# Reference dataset used for testing and demonstration (SingleCellExperiment)

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
AML_SCE()
```

## Source

The dataset was obtained from the
[Figshare](https://figshare.com/articles/dataset/Expression_of_197_surface_markers_and_462_mRNAs_in_15281_cells_from_blood_and_bone_marrow_from_a_young_healthy_donor/13398065/2).
For more information on the operations performed on the original object,
see the [SCUBA manuscript
repository](https://github.com/amc-heme/SCUBA_Manuscript/blob/main/Demo_Object_Generation.Rmd).

## Details

Contrary to convention for loading data, `AML_SCE()` is called as a
function, with parentheses. This is because the dataset is loaded with
[`HDF5Array::loadHDF5SummarizedExperiment()`](https://rdrr.io/pkg/HDF5Array/man/saveHDF5SummarizedExperiment.html)
instead of the typical process of loading data from R packages. The
dataset was saved using the `HDF5Array` package to test support of
SingleCellExperiment objects supporting HDF5 storage saved using this
package.

## Examples

``` r
# Object summary
AML_SCE()
#> class: SingleCellExperiment 
#> dim: 462 250 
#> metadata(0):
#> assays(3): counts logcounts scaledata
#> rownames(462): ACTG1 ADGRG1 ... ZFAS1 ZFP36L2
#> rowData names(0):
#> colnames(250): 487013_1 39207_1 ... 545562_2 849364_2
#> colData names(30): orig.ident nCount_RNA ... condensed_cell_type ident
#> reducedDimNames(2): PCA UMAP
#> mainExpName: RNA
#> altExpNames(1): AB

# Summary of metadata variables in object
meta_varnames(AML_SCE())
#>  [1] "orig.ident"          "nCount_RNA"          "nFeature_RNA"       
#>  [4] "nCount_AB"           "nFeature_AB"         "nCount_BOTH"        
#>  [7] "nFeature_BOTH"       "BOTH_snn_res.0.9"    "seurat_clusters"    
#> [10] "Prediction_Ind"      "BOTH_snn_res.1"      "ClusterID"          
#> [13] "Batch"               "x"                   "y"                  
#> [16] "x_mean"              "y_mean"              "cor"                
#> [19] "ct"                  "prop"                "meandist"           
#> [22] "cDC"                 "B.cells"             "Myelocytes"         
#> [25] "Erythroid"           "Megakaryocte"        "Ident"              
#> [28] "RNA_snn_res.0.4"     "condensed_cell_type" "ident"              
```
