---
title: "Introduction to Single Cell RNA-Seq Part 4: Dimensionality reduction"
author: "UCD Bioinformatics Core"
date: "2023-12-06"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 4: Dimensionality reduction
Single cell (or nucleus) data are extremely high-dimensional. In order to reduce the complexity of analysis and remove sources of noise, dimensionality reduction is an important step in the analysis workflow. In this section, we will be using two dimension reduction methods: PCA and UMAP.


## Set up workspace

```r
library(Seurat)
library(ggplot2)
experiment.aggregate <- readRDS(file="scRNA_workshop-03.rds")
experiment.aggregate
```

```
## An object of class Seurat 
## 11292 features across 6313 samples within 1 assay 
## Active assay: RNA (11292 features, 7012 variable features)
```

```r
set.seed(12345)
```

## Perform dimensionality reduction with PCA
Principal Components Analysis (PCA) is a widely-used dimension reduction method. Each PC is a vector in the reduced-dimensional space that is orthogonal to all preceding PCs. The first of these explains the largest amount of variation and each subsequent PC explains slightly less than the preceding component. PCA is performed on the scaled data, and therefore uses only the variable features. 

```r
?RunPCA
```


```r
experiment.aggregate <- RunPCA(experiment.aggregate, npcs = 100)
```

While it is theoretically possible to calculate as many PCs as there are features in the data, typically 100 PCs is more than sufficient. In fact, many of these PCs may explain negligible amounts of variation. Seurat provides a number of ways to visualize the PCA results.

### Principal components plot
The PCA biplot is a scatter plot showing the placement of each cell on two selected components. By default, the first and second PC are used, but any two calculted PCs may be specified.

At this point in the analysis, since we are no longer performing QA and filtering, we can move to examining relationships between cells on a per-group rather than per-sample basis.

```r
DimPlot(experiment.aggregate,
        group.by = "group",
        reduction = "pca",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](04-dimensionality_reduction_files/figure-html/plot_pca-1.png)<!-- -->

The axes are unit-less; points (cells or nuclei) that are farther apart are more dissimilar on the displayed PC than points that are closer together.

### PCA loadings
Each PC can be imagined as a sort of meta-gene for which every cell has an expression value. The top genes associated with the reduction component (i.e. contributing to a cell's "expression level" of that meta-gene) can be plotted for a selected dimension(s) using the `VizDimLoadings` function.

```r
VizDimLoadings(experiment.aggregate,
               dims = 1,
               nfeatures = 25,
               reduction = "pca",
               ncol = 1) +
  theme_minimal(base_size = 8)
```

![](04-dimensionality_reduction_files/figure-html/viz_pca-1.png)<!-- -->

### Heat map
Heat maps display similar information. On the x-axis, cells are ordered by their embeddings ("expression" of the PC), while on the y-axis, genes are ordered by PC loading. When fewer than the total number of cells is selected, this results in selection of the cells with the largest absolute value embeddings, which emphasizes variation on the PC.

```r
DimHeatmap(experiment.aggregate,
           dims = 1,
           nfeatures = 25,
           cells = 500,
           reduction = "pca",
           balanced = TRUE,
           slot = "scale.data")
```

![](04-dimensionality_reduction_files/figure-html/heatmap_pca-1.png)<!-- -->

#### Explore
Re-import the original data and try modifying the ScaleData vars.to.regress argument. You could remove some variables, or add others. What happens? See how choices effect the plots.

```r
experiment.explore <- readRDS("scRNA_workshop-03.rds")
experiment.explore <- ScaleData(experiment.explore) # make changes here to explore the data
experiment.explore <- RunPCA(experiment.explore) # what happens if you adjust npcs?
VizDimLoadings(experiment.explore, dims = 1:2)
DimPlot(experiment.explore, reduction = "pca")
DimHeatmap(experiment.explore, dims = 1:6, cells = 500, balanced = TRUE) # adjust parameters
```

## Selecting PCs to use
To overcome the extensive technical noise in any single gene, Seurat clusters cells based on their PCA scores, with each PC essentially representing a meta-gene that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

### Elbow plot
An elbow plot displays the standard deviations (or approximate singular values if running PCAFast) of the principle components for easy identification of an elbow in the graph. This elbow often corresponds well with the significant PCs and is much faster to run.  This is the traditional approach to selecting principal components.

The appearance of elbow plots tends to be highly consistent across single cell / single nucleus experiments. Generally, the line approaches zero at around PC 50. This is a reasonable number of PCs to use for the downstream steps.

```r
ElbowPlot(experiment.aggregate, ndims = 100)
```

![](04-dimensionality_reduction_files/figure-html/elbow-1.png)<!-- -->

### JackStraw
The JackStraw function randomly permutes a subset of data, and calculates projected PCA scores for these genes. The PCA scores for these randomly permuted genes are then compared with the observed PCA scores to determine statistical significance. The end result is a p-value for each gene's association with each principal component.

PCs with a strong enrichment of low p-value genes are identified as significant components.

**The JackStraw permutation is computationally intensive and can be quite slow. Consider skipping this step and exploring the function when you have some extra time.**

```r
experiment.aggregate <- JackStraw(experiment.aggregate, dims = 100)
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:100)
JackStrawPlot(object = experiment.aggregate, dims = 1:100) +
  scale_color_viridis_d() +
  theme(legend.position="bottom")
```

![](04-dimensionality_reduction_files/figure-html/jackstraw-1.png)<!-- -->

Let's use the first 50 PCs.


```r
use.pcs <- 1:50
```

## UMAP
[Uniform Manifold Approximation and Projection](https://arxiv.org/pdf/1802.03426v3.pdf) (UMAP) is a dimensionality reduction method that is commonly used in single cell RNA-Seq analysis. Single cell data is extremely high-dimensional; UMAP calculates a nearest neighbor network describing the relationships between cells as measured by the PC loadings of variable genes and creates a low-dimensional space that preserves these relationships.

```r
# calculate UMAP
experiment.aggregate <- RunUMAP(experiment.aggregate,
                                dims = use.pcs)
```

While UMAP can be a general non-linear dimensionality reduction approach, it's most frequently used as a visualization technique. A UMAP biplot offers a very useful graphical representation of the relationships captured by the nearest neighbor graph.


```r
# UMAP colored by sample identity
DimPlot(experiment.aggregate,
        group.by = "group",
        reduction = "umap",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](04-dimensionality_reduction_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

## Prepare for the next section

#### Save object

```r
saveRDS(experiment.aggregate, file="scRNA_workshop-04.rds")
```

#### Download Rmd

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/05-integration.Rmd", "05-integration.Rmd")
```

#### Session information

```r
sessionInfo()
```

```
## R version 4.2.2 (2022-10-31)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.4.0      SeuratObject_4.1.3 Seurat_4.3.0      
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6          
##   [4] ellipsis_0.3.2         ggridges_0.5.4         rstudioapi_0.14       
##   [7] spatstat.data_3.0-0    farver_2.1.1           leiden_0.4.3          
##  [10] listenv_0.8.0          ggrepel_0.9.2          fansi_1.0.3           
##  [13] codetools_0.2-18       splines_4.2.2          cachem_1.0.6          
##  [16] knitr_1.41             polyclip_1.10-4        jsonlite_1.8.4        
##  [19] ica_1.0-3              cluster_2.1.4          png_0.1-8             
##  [22] uwot_0.1.14            shiny_1.7.3            sctransform_0.3.5     
##  [25] spatstat.sparse_3.0-0  compiler_4.2.2         httr_1.4.4            
##  [28] assertthat_0.2.1       Matrix_1.5-3           fastmap_1.1.0         
##  [31] lazyeval_0.2.2         cli_3.4.1              later_1.3.0           
##  [34] htmltools_0.5.3        tools_4.2.2            igraph_1.3.5          
##  [37] gtable_0.3.1           glue_1.6.2             RANN_2.6.1            
##  [40] reshape2_1.4.4         dplyr_1.0.10           Rcpp_1.0.9            
##  [43] scattermore_0.8        jquerylib_0.1.4        vctrs_0.5.1           
##  [46] nlme_3.1-160           spatstat.explore_3.0-5 progressr_0.11.0      
##  [49] lmtest_0.9-40          spatstat.random_3.0-1  xfun_0.35             
##  [52] stringr_1.4.1          globals_0.16.2         mime_0.12             
##  [55] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5.1         
##  [58] goftest_1.2-3          future_1.29.0          MASS_7.3-58.1         
##  [61] zoo_1.8-11             scales_1.2.1           promises_1.2.0.1      
##  [64] spatstat.utils_3.0-1   parallel_4.2.2         RColorBrewer_1.1-3    
##  [67] yaml_2.3.6             reticulate_1.28        pbapply_1.6-0         
##  [70] gridExtra_2.3          sass_0.4.4             stringi_1.7.8         
##  [73] highr_0.9              rlang_1.0.6            pkgconfig_2.0.3       
##  [76] matrixStats_0.63.0     evaluate_0.18          lattice_0.20-45       
##  [79] tensor_1.5             ROCR_1.0-11            purrr_0.3.5           
##  [82] labeling_0.4.2         patchwork_1.1.2        htmlwidgets_1.5.4     
##  [85] cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.32.1     
##  [88] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
##  [91] R6_2.5.1               generics_0.1.3         DBI_1.1.3             
##  [94] withr_2.5.0            pillar_1.8.1           fitdistrplus_1.1-8    
##  [97] survival_3.4-0         abind_1.4-5            sp_1.5-1              
## [100] tibble_3.1.8           future.apply_1.10.0    KernSmooth_2.23-20    
## [103] utf8_1.2.2             spatstat.geom_3.0-3    plotly_4.10.1         
## [106] rmarkdown_2.18         grid_4.2.2             data.table_1.14.6     
## [109] digest_0.6.30          xtable_1.8-4           tidyr_1.2.1           
## [112] httpuv_1.6.6           munsell_0.5.0          viridisLite_0.4.1     
## [115] bslib_0.4.1
```
