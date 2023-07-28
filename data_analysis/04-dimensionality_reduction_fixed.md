---
title: "Introduction to Single Cell RNA-Seq Part 4: Dimensionality reduction"
author: "UCD Bioinformatics Core"
date: "2023-07-28"
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

<div class='r_output'> An object of class Seurat 
 11292 features across 6312 samples within 1 assay 
 Active assay: RNA (11292 features, 7012 variable features)
</div>
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

<div class='r_output'> R version 4.3.1 (2023-06-16)
 Platform: aarch64-apple-darwin20 (64-bit)
 Running under: macOS Monterey 12.4
 
 Matrix products: default
 BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
 LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
 
 locale:
 [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
 
 time zone: America/Los_Angeles
 tzcode source: internal
 
 attached base packages:
 [1] stats     graphics  grDevices utils     datasets  methods   base     
 
 other attached packages:
 [1] ggplot2_3.4.2      SeuratObject_4.1.3 Seurat_4.3.0.1    
 
 loaded via a namespace (and not attached):
   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
   [4] rlang_1.1.1            magrittr_2.0.3         RcppAnnoy_0.0.21      
   [7] spatstat.geom_3.2-4    matrixStats_1.0.0      ggridges_0.5.4        
  [10] compiler_4.3.1         png_0.1-8              vctrs_0.6.3           
  [13] reshape2_1.4.4         stringr_1.5.0          pkgconfig_2.0.3       
  [16] fastmap_1.1.1          ellipsis_0.3.2         labeling_0.4.2        
  [19] utf8_1.2.3             promises_1.2.0.1       rmarkdown_2.23        
  [22] purrr_1.0.1            xfun_0.39              cachem_1.0.8          
  [25] jsonlite_1.8.7         goftest_1.2-3          highr_0.10            
  [28] later_1.3.1            spatstat.utils_3.0-3   irlba_2.3.5.1         
  [31] parallel_4.3.1         cluster_2.1.4          R6_2.5.1              
  [34] ica_1.0-3              spatstat.data_3.0-1    bslib_0.5.0           
  [37] stringi_1.7.12         RColorBrewer_1.1-3     reticulate_1.30       
  [40] parallelly_1.36.0      lmtest_0.9-40          jquerylib_0.1.4       
  [43] scattermore_1.2        Rcpp_1.0.11            knitr_1.43            
  [46] tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
  [49] sctransform_0.3.5      httpuv_1.6.11          Matrix_1.6-0          
  [52] splines_4.3.1          igraph_1.5.0           tidyselect_1.2.0      
  [55] abind_1.4-5            rstudioapi_0.15.0      yaml_2.3.7            
  [58] spatstat.random_3.1-5  codetools_0.2-19       miniUI_0.1.1.1        
  [61] spatstat.explore_3.2-1 listenv_0.9.0          lattice_0.21-8        
  [64] tibble_3.2.1           plyr_1.8.8             withr_2.5.0           
  [67] shiny_1.7.4.1          ROCR_1.0-11            evaluate_0.21         
  [70] Rtsne_0.16             future_1.33.0          survival_3.5-5        
  [73] polyclip_1.10-4        fitdistrplus_1.1-11    pillar_1.9.0          
  [76] KernSmooth_2.23-22     plotly_4.10.2          generics_0.1.3        
  [79] sp_2.0-0               munsell_0.5.0          scales_1.2.1          
  [82] globals_0.16.2         xtable_1.8-4           glue_1.6.2            
  [85] lazyeval_0.2.2         tools_4.3.1            data.table_1.14.8     
  [88] RANN_2.6.1             leiden_0.4.3           cowplot_1.1.1         
  [91] grid_4.3.1             tidyr_1.3.0            colorspace_2.1-0      
  [94] nlme_3.1-162           patchwork_1.1.2        cli_3.6.1             
  [97] spatstat.sparse_3.0-2  fansi_1.0.4            viridisLite_0.4.2     
 [100] dplyr_1.1.2            uwot_0.1.16            gtable_0.3.3          
 [103] sass_0.4.7             digest_0.6.33          progressr_0.13.0      
 [106] ggrepel_0.9.3          farver_2.1.1           htmlwidgets_1.6.2     
 [109] htmltools_0.5.5        lifecycle_1.0.3        httr_1.4.6            
 [112] mime_0.12              MASS_7.3-60
</div>