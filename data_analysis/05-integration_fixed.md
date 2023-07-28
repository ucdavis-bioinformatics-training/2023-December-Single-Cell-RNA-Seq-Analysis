---
title: "Introduction to Single Cell RNA-Seq Part 5: Integration"
author: "UCD Bioinformatics Core"
date: "2023-07-28"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---
# Introduction to Single Cell RNA-Seq Part 5: Integration

More and more experiments involve a large number of samples/datasets, that may have been prepared in separate batches. Or in the case where one would like to include or integrate publicly available datasets. It is important to properly integrate these datasets, and we will see the effect the integration has at the end of this documentation.

Most of the methods that were developed to integrate single cell datasets fall into two categories. The first is the "anchor" based approach. In this approach, the first step is to select a batch as the "anchor" and convert other batches to the "anchor" batch. Among these approaches are [MNN](https://github.com/MarioniLab/MNN2017), [iMAP](https://github.com/Svvord/iMAP), and [SCALEX](https://github.com/jsxlei/SCALEX). The advantage of the anchor-based approach is that different batches of cells can be studied under the same experimental conditions, and the disadvantage is that it is not possible to fully combine the features of each batch because the cell types contained in each batch are unknown.

The second approach is to transform all batches of data to a low-dimensional space to correct batch effects, such as implemented in [Scanorama](https://github.com/brianhie/scanorama), [Harmony](https://github.com/immunogenomics/harmony), [DESC](https://www.nature.com/articles/s41467-020-15851-3), [BBKNN](https://github.com/Teichlab/bbknn), [STACAS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8098019/) and [Seurat's integration](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8). This second approach has the advantage of extracting biologically relevant latent features and reducing the impact of noise, but it cannot be used for differential gene expression analysis. Many of these existing methods work well when the batches of datasets have the same cell types, however, they fail when there are different cell types involved in different datasets.

Very recently, a [new approach](https://www.mdpi.com/1422-0067/23/4/2082) has been developed that uses connected graphs and generative adversarial networks (GAN) to achieve the goal of eliminating nonbiological noise between batches of datasets. This new method has been demonstrated to work well both in the situation where datasets have the same cell types and in the situation where datasets may have different cell types.

In this workshop, we are going to look at Seurat's integration approach using reciprocal PCA, which is supurior to its first integration approach using canonical correlation analysis. The basic idea is to identify cross-dataset pairs cells that are in a matched biological state ("anchors"), and use them to correct technical differences between datasets. The integration method we use has been implemented in Seurat and you can find the details of the method in [its publication](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub).

## Set up workspace

```r
library(Seurat)
library(ggplot2)
set.seed(12345)
experiment.aggregate <- readRDS("scRNA_workshop-02.rds") # filtered object
```

## Prepare data for integration
Prior to integration samples should be processed independently. First, we split the filtered object by sample to create a list of Seurat objects.

```r
experiment.split <- SplitObject(experiment.aggregate, split.by = "orig.ident")
rm(experiment.aggregate)
```

Each object is then normalized.

```r
experiment.split <- lapply(experiment.split, function(sce){
  sce = NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce = CellCycleScoring(sce, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
})
```

## Select integration features
Integration features are genes that are repeatedly variable across the objects to integrate. These are used to scale the data and run the PCA.

```r
features <- SelectIntegrationFeatures(object.list = experiment.split)
```

## Scale data and run PCA
Once integration features have been identified, we can scale the data and run the PCA.

```r
experiment.split <- lapply(experiment.split,function(sce){
  sce = ScaleData(sce, features = features, vars.to.regress = c("S.Score", "G2M.Score", "percent_MT", "nFeature_RNA"))
  RunPCA(sce, features = features)
})
```

## Idenfity integration anchors
The integration anchors are pairs of cells that are mutual nearest neighbors on the . These may be calculated either for each Seurat object relative to a reference object, or pairwise between all objects if no reference is provided.

```r
anchors <- FindIntegrationAnchors(object.list = experiment.split, anchor.features = features, reduction = "rpca")
```

## Integrate

```r
experiment.integrated <- IntegrateData(anchorset = anchors)
```
The new experiment.integrated object has two assays: RNA and integrated. The RNA assay contains the normalized, scaled data from the individual experiment.split objects merged into a single table, while the data in the integrated assay has been scaled in such a way that it is no longer appropriate to use this assay for differential expression.

The authors recommend using the integrated assay for clustering and visualization (UMAP plots).

## Impact of integration
In the dimensionality reduction section we performed PCA on the complete experiment.aggregate object, where we used the vars.to.regress argument of the ScaleData function to adjust for cell cycle, nucleus integrity, and sequencing depth. The PCA biplot looked like this:

![Previous PCA plot](04-dimensionality_reduction_files/figure-html/plot_pca-1.png)

After integration, the appearance of the PCA biplot has changed; cells no longer separate by group.

```r
experiment.integrated <- ScaleData(experiment.integrated, assay="integrated")
experiment.integrated <- RunPCA(experiment.integrated, assay="integrated")
DimPlot(experiment.integrated,
        group.by = "group",
        reduction = "pca",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](05-integration_files/figure-html/PCA-1.png)<!-- -->

A similar effect can be seen in the UMAP. Previously, the un-integrated UMAP plot had this appearance:

![Previous UMAP plot](04-dimensionality_reduction_files/figure-html/unnamed-chunk-1-1.png)

After integration, the polyp and colorectal cancer cells are more co-localized on the biplot.

```r
experiment.integrated <- RunUMAP(experiment.integrated, dims = 1:50)
DimPlot(experiment.integrated,
        reduction = "umap",
        group.by = "group",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](05-integration_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

In the next section, we will use the integrated data to perform clustering.

### Visualize metadata


```r
lapply(c("nCount_RNA", "nFeature_RNA", "percent_MT"), function(feature){
  FeaturePlot(experiment.integrated,
              reduction = "umap",
              features = feature)
})
```

<div class='r_output'> [[1]]
</div>
![](05-integration_files/figure-html/meta-1.png)<!-- -->

<div class='r_output'> 
 [[2]]
</div>
![](05-integration_files/figure-html/meta-2.png)<!-- -->

<div class='r_output'> 
 [[3]]
</div>
![](05-integration_files/figure-html/meta-3.png)<!-- -->

### Visualize cell cycle phase


```r
DimPlot(experiment.integrated,
        reduction = "umap",
        group.by = "Phase",
        shuffle = TRUE) +
  scale_color_viridis_d()
```

![](05-integration_files/figure-html/phase-1.png)<!-- -->

## Prepare for the next section

#### Save the Seurat object


```r
saveRDS(experiment.integrated, file="scRNA_workshop-05.rds")
```

#### Download the Rmd document

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-June-Single-Cell-RNA-Seq-Analysis/main/data_analysis/06-clustering_celltype.Rmd", "06-clustering_celltype.Rmd")
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