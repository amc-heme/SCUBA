# Reference dataset used for testing and demonstration (Seurat with BPCells-backed matrices)

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
AML_BPCells
```

## Format

An object of class `Seurat` with 462 rows and 250 columns.

## Source

The dataset was obtained from the
[Figshare](https://figshare.com/articles/dataset/Expression_of_197_surface_markers_and_462_mRNAs_in_15281_cells_from_blood_and_bone_marrow_from_a_young_healthy_donor/13398065/2).
For more information on the operations performed on the original object,
see the [SCUBA manuscript
repository](https://github.com/amc-heme/SCUBA_Manuscript/blob/main/Demo_Object_Generation.Rmd).

## Details

This object is identical to `AML_Seurat`, but the counts and data
matrices for the RNA and AB assays are stored using
[BPCells](https://bnprks.github.io/BPCells/) on-disk matrices. The
BPCells-backed matrices enable efficient storage and access of large
single-cell datasets while maintaining full compatibility with SCUBA's
data retrieval functions. The on-disk matrices are stored in the
`inst/extdata/AML_BPCells/` directory.

## Examples

``` r
# Object summary
AML_BPCells
#> An object of class Seurat 
#> 659 features across 250 samples within 2 assays 
#> Active assay: RNA (462 features, 462 variable features)
#>  3 layers present: data, counts, scale.data
#>  1 other assay present: AB
#>  2 dimensional reductions calculated: pca, umap

# Summary of metadata variables in object
meta_varnames(AML_BPCells)
#>  [1] "orig.ident"          "nCount_RNA"          "nFeature_RNA"       
#>  [4] "nCount_AB"           "nFeature_AB"         "nCount_BOTH"        
#>  [7] "nFeature_BOTH"       "BOTH_snn_res.0.9"    "seurat_clusters"    
#> [10] "Prediction_Ind"      "BOTH_snn_res.1"      "ClusterID"          
#> [13] "Batch"               "x"                   "y"                  
#> [16] "x_mean"              "y_mean"              "cor"                
#> [19] "ct"                  "prop"                "meandist"           
#> [22] "cDC"                 "B.cells"             "Myelocytes"         
#> [25] "Erythroid"           "Megakaryocte"        "Ident"              
#> [28] "RNA_snn_res.0.4"     "condensed_cell_type"
```
