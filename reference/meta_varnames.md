# Summarize object metadata

Returns the names of all metadata variables in an object.

## Usage

``` r
meta_varnames(object)

# S3 method for class 'Seurat'
meta_varnames(object)

# S3 method for class 'SingleCellExperiment'
meta_varnames(object)

# S3 method for class 'AnnDataR6'
meta_varnames(object)
```

## Arguments

- object:

  A single cell object. Currently, Seurat, SingleCellExpleriment, and
  anndata objects are supported.

## Methods (by class)

- `meta_varnames(Seurat)`: Seurat objects

- `meta_varnames(SingleCellExperiment)`: SingleCellExperiment objects

- `meta_varnames(AnnDataR6)`: Anndata objects

## Examples

``` r
# Seurat objects
meta_varnames(AML_Seurat)
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

# SingleCellExperiment objects
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

# anndata objects
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
