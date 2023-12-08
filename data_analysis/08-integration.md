---
title: "Introduction to Single Cell RNA-Seq Part 8: Integration"
author: "UCD Bioinformatics Core"
date: "2023-12-07"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---
# Introduction to Single Cell RNA-Seq Part 8: Integration

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

![](08-integration_files/figure-html/PCA-1.png)<!-- -->

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

![](08-integration_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

In the next section, we will use the integrated data to perform clustering.

### Visualize metadata


```r
lapply(c("nCount_RNA", "nFeature_RNA", "percent_MT"), function(feature){
  FeaturePlot(experiment.integrated,
              reduction = "umap",
              features = feature)
})
```

```
## [[1]]
```

![](08-integration_files/figure-html/meta-1.png)<!-- -->

```
## 
## [[2]]
```

![](08-integration_files/figure-html/meta-2.png)<!-- -->

```
## 
## [[3]]
```

![](08-integration_files/figure-html/meta-3.png)<!-- -->

### Visualize cell cycle phase


```r
DimPlot(experiment.integrated,
        reduction = "umap",
        group.by = "Phase",
        shuffle = TRUE) +
  scale_color_viridis_d()
```

![](08-integration_files/figure-html/phase-1.png)<!-- -->

### Clusters using the integrated data


```r
experiment.integrated <- FindNeighbors(experiment.integrated, reduction = "pca", dims = 50)
experiment.integrated <- FindClusters(experiment.integrated,
                                     resolution = seq(0.1, 0.4, 0.1))
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6313
## Number of edges: 105606
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9894
## Number of communities: 17
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6313
## Number of edges: 105606
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9842
## Number of communities: 23
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6313
## Number of edges: 105606
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9800
## Number of communities: 28
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6313
## Number of edges: 105606
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9764
## Number of communities: 35
## Elapsed time: 0 seconds
```

```r
cluster.resolutions <- grep("res", colnames(experiment.integrated@meta.data), value = TRUE)
sapply(cluster.resolutions, function(res){
  length(levels(experiment.integrated@meta.data[,res]))
})
```

```
## integrated_snn_res.0.1 integrated_snn_res.0.2 integrated_snn_res.0.3 
##                     17                     23                     28 
## integrated_snn_res.0.4 
##                     35
```

```r
lapply(cluster.resolutions, function(res){
  DimPlot(experiment.integrated,
          group.by = res,
          reduction = "umap",
          shuffle = TRUE) +
    scale_color_viridis_d(option = "turbo")
})
```

```
## [[1]]
```

![](08-integration_files/figure-html/clusters.integrated-1.png)<!-- -->

```
## 
## [[2]]
```

![](08-integration_files/figure-html/clusters.integrated-2.png)<!-- -->

```
## 
## [[3]]
```

![](08-integration_files/figure-html/clusters.integrated-3.png)<!-- -->

```
## 
## [[4]]
```

![](08-integration_files/figure-html/clusters.integrated-4.png)<!-- -->

```r
lapply(cluster.resolutions, function(res){
         tmp = experiment.integrated@meta.data[,c(res, "group")]
         colnames(tmp) = c("cluster", "group")
         ggplot(tmp, aes(x = cluster, fill = group)) +
           geom_bar() +
           scale_fill_viridis_d(option = "mako") +
           theme_classic()
})
```

```
## [[1]]
```

![](08-integration_files/figure-html/clusters.integrated-5.png)<!-- -->

```
## 
## [[2]]
```

![](08-integration_files/figure-html/clusters.integrated-6.png)<!-- -->

```
## 
## [[3]]
```

![](08-integration_files/figure-html/clusters.integrated-7.png)<!-- -->

```
## 
## [[4]]
```

![](08-integration_files/figure-html/clusters.integrated-8.png)<!-- -->

```r
FeaturePlot(experiment.integrated,
            reduction = "umap",
            features = "KCNMA1")
```

![](08-integration_files/figure-html/clusters.integrated-9.png)<!-- -->

```r
VlnPlot(experiment.integrated,
        group.by = "integrated_snn_res.0.1",
        features = "KCNMA1",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](08-integration_files/figure-html/clusters.integrated-10.png)<!-- -->

```r
Idents(experiment.integrated) <- "integrated_snn_res.0.1"
experiment.integrated <- BuildClusterTree(experiment.integrated, dims = 1:50)
PlotClusterTree(experiment.integrated)
```

![](08-integration_files/figure-html/clusters.integrated-11.png)<!-- -->

```r
experiment.integrated <- RenameIdents(experiment.integrated,
                                     '4' = '1',
                                     '16' = '1',
                                     '11' = '1',
                                     '10' = '1')
experiment.integrated$res.0.1_merged <- Idents(experiment.integrated)
table(experiment.integrated$res.0.1_merged)
```

```
## 
##    1    0    2    3    5    6    7    8    9   12   13   14   15 
## 1770  527  490  472  460  449  406  379  344  272  265  257  222
```

```r
new.order <- as.character(sort(as.numeric(levels(experiment.integrated$res.0.1_merged))))
experiment.integrated$res.0.1_merged <- factor(experiment.integrated$res.0.1_merged, levels = new.order)
DimPlot(experiment.integrated,
        reduction = "umap",
        group.by = "res.0.1_merged",
        label = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](08-integration_files/figure-html/clusters.integrated-12.png)<!-- -->



```r
experiment.aggregate <- readRDS("scRNA_workshop-05.rds")
identical(rownames(experiment.aggregate@meta.data), rownames(experiment.integrated@meta.data))
```

```
## [1] TRUE
```

```r
experiment.aggregate <- AddMetaData(experiment.aggregate,
                                    metadata = experiment.integrated$integrated_snn_res.0.1,
                                    col.name = "integrated_snn_res.0.1")
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "integrated_snn_res.0.1",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](08-integration_files/figure-html/AddMetaData-1.png)<!-- -->

In this case, integration appears to have interfered with forming easily-interpreted clusters. There is very little relationship between location on the UMAP and cluster identity, which makes it harder to identify possible cell populations at a glance. We can add the integrated object's cluster identities to the un-integrated object if we choose, though projecting the integrated clustering onto the un-integrated UMAP is unlikely to be useful.

#### Save the Seurat object


```r
saveRDS(experiment.integrated, file="scRNA_workshop-08.rds")
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
##   [1] ggbeeswarm_0.6.0       Rtsne_0.16             colorspace_2.0-3      
##   [4] deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.4        
##   [7] rstudioapi_0.14        spatstat.data_3.0-0    farver_2.1.1          
##  [10] leiden_0.4.3           listenv_0.8.0          ggrepel_0.9.2         
##  [13] fansi_1.0.3            codetools_0.2-18       splines_4.2.2         
##  [16] cachem_1.0.6           knitr_1.41             polyclip_1.10-4       
##  [19] jsonlite_1.8.4         ica_1.0-3              cluster_2.1.4         
##  [22] png_0.1-8              uwot_0.1.14            shiny_1.7.3           
##  [25] sctransform_0.3.5      spatstat.sparse_3.0-0  compiler_4.2.2        
##  [28] httr_1.4.4             assertthat_0.2.1       Matrix_1.5-3          
##  [31] fastmap_1.1.0          lazyeval_0.2.2         cli_3.4.1             
##  [34] later_1.3.0            htmltools_0.5.3        tools_4.2.2           
##  [37] igraph_1.3.5           gtable_0.3.1           glue_1.6.2            
##  [40] RANN_2.6.1             reshape2_1.4.4         dplyr_1.0.10          
##  [43] Rcpp_1.0.9             scattermore_0.8        jquerylib_0.1.4       
##  [46] vctrs_0.5.1            ape_5.6-2              nlme_3.1-160          
##  [49] spatstat.explore_3.0-5 progressr_0.11.0       lmtest_0.9-40         
##  [52] spatstat.random_3.0-1  xfun_0.35              stringr_1.4.1         
##  [55] globals_0.16.2         mime_0.12              miniUI_0.1.1.1        
##  [58] lifecycle_1.0.3        irlba_2.3.5.1          goftest_1.2-3         
##  [61] future_1.29.0          MASS_7.3-58.1          zoo_1.8-11            
##  [64] scales_1.2.1           promises_1.2.0.1       spatstat.utils_3.0-1  
##  [67] parallel_4.2.2         RColorBrewer_1.1-3     yaml_2.3.6            
##  [70] reticulate_1.28        pbapply_1.6-0          gridExtra_2.3         
##  [73] ggrastr_1.0.1          sass_0.4.4             stringi_1.7.8         
##  [76] highr_0.9              rlang_1.0.6            pkgconfig_2.0.3       
##  [79] matrixStats_0.63.0     evaluate_0.18          lattice_0.20-45       
##  [82] tensor_1.5             ROCR_1.0-11            purrr_0.3.5           
##  [85] labeling_0.4.2         patchwork_1.1.2        htmlwidgets_1.5.4     
##  [88] cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.32.1     
##  [91] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
##  [94] R6_2.5.1               generics_0.1.3         DBI_1.1.3             
##  [97] withr_2.5.0            pillar_1.8.1           fitdistrplus_1.1-8    
## [100] survival_3.4-0         abind_1.4-5            sp_1.5-1              
## [103] tibble_3.1.8           future.apply_1.10.0    KernSmooth_2.23-20    
## [106] utf8_1.2.2             spatstat.geom_3.0-3    plotly_4.10.1         
## [109] rmarkdown_2.18         grid_4.2.2             data.table_1.14.6     
## [112] digest_0.6.30          xtable_1.8-4           tidyr_1.2.1           
## [115] httpuv_1.6.6           munsell_0.5.0          beeswarm_0.4.0        
## [118] viridisLite_0.4.1      vipor_0.4.5            bslib_0.4.1
```
