---
title: "Introduction to Single Cell RNA-Seq Part 6: Clustering and cell type assignment"
author: "UCD Bioinformatics Core"
date: "2023-07-28"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 6: Clustering and cell type assignment
Clustering and cell type assignment are critical steps in many single cell (or single nucleus) experiments. The power of single cell experiments is to capture the heterogeneity within samples as well as between them. Clustering permits the user to organize cells into clusters that correspond to biologically and experimentally relevant populations.


## Set up workspace

```r
library(Seurat)
library(kableExtra)
library(tidyr)
library(dplyr)
library(ggplot2)
library(HGNChelper)
set.seed(12345)
experiment.aggregate <- readRDS("scRNA_workshop-04.rds")
experiment.aggregate
```

<div class='r_output'> An object of class Seurat 
 11292 features across 6312 samples within 1 assay 
 Active assay: RNA (11292 features, 7012 variable features)
  2 dimensional reductions calculated: pca, umap
</div>
## Construct network
Seurat implements an graph-based clustering approach. Distances between the cells are calculated based on previously identified PCs.

The default method for identifying k-nearest neighbors is [annoy](https://github.com/spotify/annoy), an approximate nearest-neighbor approach that is widely used for high-dimensional analysis in many fields, including single-cell analysis. Extensive community benchmarking has shown that annoy substantially improves the speed and memory requirements of neighbor discovery, with negligible impact to downstream results.

```r
experiment.aggregate <- FindNeighbors(experiment.aggregate, reduction = "pca", dims = 1:50)
```

## Find clusters
The FindClusters function implements the neighbor based clustering procedure, and contains a resolution parameter that sets the granularity of the downstream clustering, with increased values leading to a greater number of clusters. This code produces a series of resolutions for us to investigate and choose from.

The clustering resolution parameter is unit-less and somewhat arbitrary. The resolutions used here were selected to produce a useable number of clusters in the example experiment. 

```r
experiment.aggregate <- FindClusters(experiment.aggregate,
                                     resolution = seq(0.1, 0.4, 0.1))
```

<div class='r_output'> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 6312
 Number of edges: 264108
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.9679
 Number of communities: 11
 Elapsed time: 0 seconds
 Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 6312
 Number of edges: 264108
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.9485
 Number of communities: 12
 Elapsed time: 0 seconds
 Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 6312
 Number of edges: 264108
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.9342
 Number of communities: 15
 Elapsed time: 0 seconds
 Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 6312
 Number of edges: 264108
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.9218
 Number of communities: 15
 Elapsed time: 0 seconds
</div>
Seurat adds the clustering information to the metadata table. Each FindClusters call generates a new column named with the assay, followed by "_snn_res.", and the resolution.

```r
cluster.resolutions <- grep("res", colnames(experiment.aggregate@meta.data), value = TRUE)
head(experiment.aggregate@meta.data[,cluster.resolutions]) %>%
  kable(caption = 'Cluster identities are added to the meta.data slot.') %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Cluster identities are added to the meta.data slot.</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> RNA_snn_res.0.1 </th>
   <th style="text-align:left;"> RNA_snn_res.0.2 </th>
   <th style="text-align:left;"> RNA_snn_res.0.3 </th>
   <th style="text-align:left;"> RNA_snn_res.0.4 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> AAACCCAAGTTATGGA_A001-C-007 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAACGCTTCTCTGCTG_A001-C-007 </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAAGAACGTGCTTATG_A001-C-007 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAAGAACGTTTCGCTC_A001-C-007 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAAGGATTCATTACCT_A001-C-007 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAAGTGACACGCTTAA_A001-C-007 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 5 </td>
  </tr>
</tbody>
</table>

## Explore clustering resolutions
The number of clusters produced increases with the clustering resolution.

```r
sapply(cluster.resolutions, function(res){
  length(levels(experiment.aggregate@meta.data[,res]))
})
```

<div class='r_output'> RNA_snn_res.0.1 RNA_snn_res.0.2 RNA_snn_res.0.3 RNA_snn_res.0.4 
              11              12              15              15
</div>
### Visualize clustering

Dimensionality reduction plots can be used to visualize the clustering results. On these plots, we can see how each clustering resolution aligns with patterns in the data revealed by dimensionality reductions.

#### UMAP


```r
# UMAP colored by cluster
lapply(cluster.resolutions, function(res){
  DimPlot(experiment.aggregate,
          group.by = res,
          reduction = "umap",
          shuffle = TRUE) +
    scale_color_viridis_d(option = "turbo")
})
```

<div class='r_output'> [[1]]
</div>
![](06-clustering_celltype_files/figure-html/UMAP-1.png)<!-- -->

<div class='r_output'> 
 [[2]]
</div>
![](06-clustering_celltype_files/figure-html/UMAP-2.png)<!-- -->

<div class='r_output'> 
 [[3]]
</div>
![](06-clustering_celltype_files/figure-html/UMAP-3.png)<!-- -->

<div class='r_output'> 
 [[4]]
</div>
![](06-clustering_celltype_files/figure-html/UMAP-4.png)<!-- -->

### Investigate the relationship between cluster identity and sample identity


```r
lapply(cluster.resolutions, function(res){
         tmp = experiment.aggregate@meta.data[,c(res, "group")]
         colnames(tmp) = c("cluster", "group")
         ggplot(tmp, aes(x = cluster, fill = group)) +
           geom_bar() +
           scale_fill_viridis_d(option = "mako") +
           theme_classic()
})
```

<div class='r_output'> [[1]]
</div>
![](06-clustering_celltype_files/figure-html/membership-1.png)<!-- -->

<div class='r_output'> 
 [[2]]
</div>
![](06-clustering_celltype_files/figure-html/membership-2.png)<!-- -->

<div class='r_output'> 
 [[3]]
</div>
![](06-clustering_celltype_files/figure-html/membership-3.png)<!-- -->

<div class='r_output'> 
 [[4]]
</div>
![](06-clustering_celltype_files/figure-html/membership-4.png)<!-- -->

### Investigate the relationship between cluster identity and metadata values
Here, example plots are displayed for the lowest resolution in order to save space. To see plots for each resolution, use `lapply()`.

```r
VlnPlot(experiment.aggregate,
        group.by = "RNA_snn_res.0.4",
        features = "nCount_RNA",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
VlnPlot(experiment.aggregate,
        group.by = "RNA_snn_res.0.4",
        features = "nFeature_RNA",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

```r
VlnPlot(experiment.aggregate,
        group.by = "RNA_snn_res.0.4",
        features = "percent_MT",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/unnamed-chunk-1-3.png)<!-- -->
 
### Visualize expression of genes of interest


```r
FeaturePlot(experiment.aggregate,
            reduction = "umap",
            features = "KCNMA1")
```

![](06-clustering_celltype_files/figure-html/feature-1.png)<!-- -->

```r
VlnPlot(experiment.aggregate,
        group.by = "RNA_snn_res.0.4",
        features = "KCNMA1",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/feature-2.png)<!-- -->

## Select a resolution
For now, let's use resolution 0.4. Over the remainder of this section, we will refine the clustering further.

```r
Idents(experiment.aggregate) <- "RNA_snn_res.0.4"
```

## Visualize cluster tree
Building a phylogenetic tree relating the 'average' cell from each group in default 'Ident' (currently "RNA_snn_res.0.1"). This tree is estimated based on a distance matrix constructed in either gene expression space or PCA space.

```r
experiment.aggregate <- BuildClusterTree(experiment.aggregate, dims = 1:50)
PlotClusterTree(experiment.aggregate)
```

![](06-clustering_celltype_files/figure-html/tree-1.png)<!-- -->


## Merge clusters
In many experiments, the clustering resolution does not need to be uniform across all of the cell types present. While for some cell types of interest fine detail may be desirable, for others, simply grouping them into a larger parent cluster is sufficient. Merging cluster is very straightforward.

```r
experiment.aggregate <- RenameIdents(experiment.aggregate,
                                     '11' = '10',
                                     '12' = '8')
experiment.aggregate$res.0.4_merged <- Idents(experiment.aggregate)
table(experiment.aggregate$res.0.4_merged)
```

<div class='r_output'> 
   10    8    0    1    2    3    4    5    6    7    9   13   14 
  155  238 1419  915  854  658  606  604  414  210  179   35   25
</div>
```r
experiment.aggregate@meta.data %>%
  ggplot(aes(x = res.0.4_merged, fill = group)) +
  geom_bar() +
  scale_fill_viridis_d(option = "mako") +
  theme_classic()
```

![](06-clustering_celltype_files/figure-html/merge-1.png)<!-- -->

```r
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "res.0.4_merged",
        label = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/merge-2.png)<!-- -->

```r
VlnPlot(experiment.aggregate,
        group.by = "res.0.4_merged",
        features = "KCNMA1") +
  scale_fill_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/merge-3.png)<!-- -->

## Reorder the clusters
Merging the clusters changed the order in which they appear on a plot. In order to reorder the clusters for plotting purposes take a look at the levels of the identity, then re-level as desired.

```r
levels(experiment.aggregate$res.0.4_merged)
```

<div class='r_output'>  [1] "10" "8"  "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "9"  "13" "14"
</div>
```r
# move one cluster to the first position
experiment.aggregate$res.0.4_merged <- relevel(experiment.aggregate$res.0.4_merged, "0")
levels(experiment.aggregate$res.0.4_merged)
```

<div class='r_output'>  [1] "0"  "10" "8"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "9"  "13" "14"
</div>
```r
# the color assigned to some clusters will change
VlnPlot(experiment.aggregate,
        group.by = "res.0.4_merged",
        features = "CAPN9",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/relevel-1.png)<!-- -->

```r
# re-level entire factor
new.order <- as.character(sort(as.numeric(levels(experiment.aggregate$res.0.4_merged))))
experiment.aggregate$res.0.4_merged <- factor(experiment.aggregate$res.0.4_merged, levels = new.order)
levels(experiment.aggregate$res.0.4_merged)
```

<div class='r_output'>  [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "13" "14"
</div>
```r
VlnPlot(experiment.aggregate,
        group.by = "res.0.4_merged",
        features = "MYRIP",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/relevel-2.png)<!-- -->

## Subcluster
While merging clusters reduces the resolution in some parts of the experiment, sub-clustering has the opposite effect. Let's produce sub-clusters for cluster 5.

```r
experiment.aggregate <- FindSubCluster(experiment.aggregate,
                                       graph.name = "RNA_snn",
                                       cluster = 5,
                                       subcluster.name = "subcluster")
```

<div class='r_output'> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 604
 Number of edges: 23551
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.7268
 Number of communities: 4
 Elapsed time: 0 seconds
</div>
```r
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "subcluster",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/subcluster-1.png)<!-- -->

## Subset experiment by cluster identity
After exploring and refining the cluster resolution, we may have identified some clusters that are composed of cells we aren't interested in. For example, if we have identified a cluster likely composed of contaminants, this cluster can be removed from the analysis. Alternatively, if a group of clusters have been identified as particularly of interest, these can be isolated and re-analyzed.

```r
# remove cluster 13
Idents(experiment.aggregate) <- experiment.aggregate$subcluster
experiment.tmp <- subset(experiment.aggregate, subcluster != "13")
DimPlot(experiment.tmp,
        reduction = "umap",
        group.by = "subcluster",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/subset-1.png)<!-- -->

```r
# retain cells belonging only to specified clusters
experiment.tmp <- subset(experiment.aggregate, subcluster %in% c("1", "2", "5_0", "5_1", "5_2", "5_3", "7"))
DimPlot(experiment.tmp,
        reduction = "umap",
        group.by = "subcluster",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/subset-2.png)<!-- -->

## Identify marker genes

Seurat provides several functions that can help you find markers that define clusters via differential expression:

* `FindMarkers` identifies markers for a cluster relative to all other clusters

* `FindAllMarkers` performs the find markers operation for all clusters

* `FindAllMarkersNode` defines all markers that split a node from the cluster tree

### FindMarkers


```r
markers <- FindMarkers(experiment.aggregate,
                       group.by = "subcluster",
                       ident.1 = "1")
length(which(markers$p_val_adj < 0.05)) # how many are significant?
```

<div class='r_output'> [1] 1477
</div>
```r
head(markers) %>%
  kable() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> p_val </th>
   <th style="text-align:right;"> avg_log2FC </th>
   <th style="text-align:right;"> pct.1 </th>
   <th style="text-align:right;"> pct.2 </th>
   <th style="text-align:right;"> p_val_adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CLCA1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 3.397040 </td>
   <td style="text-align:right;"> 0.484 </td>
   <td style="text-align:right;"> 0.056 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAPN9 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.310178 </td>
   <td style="text-align:right;"> 0.551 </td>
   <td style="text-align:right;"> 0.037 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANO7 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.254361 </td>
   <td style="text-align:right;"> 0.551 </td>
   <td style="text-align:right;"> 0.032 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MYRIP </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.589449 </td>
   <td style="text-align:right;"> 0.763 </td>
   <td style="text-align:right;"> 0.143 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> STXBP5-AS1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.246695 </td>
   <td style="text-align:right;"> 0.617 </td>
   <td style="text-align:right;"> 0.124 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HEPACAM2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.149583 </td>
   <td style="text-align:right;"> 0.495 </td>
   <td style="text-align:right;"> 0.052 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

The "pct.1" and "pct.2" columns record the proportion of cells with normalized expression above 0 in ident.1 and ident.2, respectively. The "p_val" is the raw p-value associated with the differential expression test, while the BH-adjusted value is found in "p_val_adj". Finally, "avg_logFC" is the average log fold change difference between the two groups.

Marker genes identified this way can be visualized in violin plots, feature plots, and heat maps.

```r
view.markers <- c(rownames(markers[markers$avg_log2FC > 0,])[1],
                  rownames(markers[markers$avg_log2FC < 0,])[1])
lapply(view.markers, function(marker){
  VlnPlot(experiment.aggregate,
          group.by = "subcluster",
          features = marker) +
    scale_fill_viridis_d(option = "turbo")
})
```

<div class='r_output'> [[1]]
</div>
![](06-clustering_celltype_files/figure-html/view.markers-1.png)<!-- -->

<div class='r_output'> 
 [[2]]
</div>
![](06-clustering_celltype_files/figure-html/view.markers-2.png)<!-- -->

```r
FeaturePlot(experiment.aggregate,
            features = view.markers,
            ncol = 2)
```

![](06-clustering_celltype_files/figure-html/view.markers-3.png)<!-- -->

```r
DoHeatmap(experiment.aggregate,
          group.by = "subcluster",
          features = view.markers,
          group.colors = viridis::turbo(length(unique(experiment.aggregate$subcluster))))
```

![](06-clustering_celltype_files/figure-html/view.markers-4.png)<!-- -->


### FindAllMarkers

FindAllMarkers can be used to automate this process across all clusters.


```r
markers <- FindAllMarkers(experiment.aggregate,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          thresh.use = 0.25)
tapply(markers$p_val_adj, markers$cluster, function(x){
  length(x < 0.05)
})
```

<div class='r_output'>  5_2   10  5_1    3    8    2    7    1    4    6    0   13   14    9  5_0  5_3 
  213  261  739 1180  342  707  641  436  750 1119  572  294  267  546  256  323
</div>
```r
head(markers) %>%
  kable() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> p_val </th>
   <th style="text-align:right;"> avg_log2FC </th>
   <th style="text-align:right;"> pct.1 </th>
   <th style="text-align:right;"> pct.2 </th>
   <th style="text-align:right;"> p_val_adj </th>
   <th style="text-align:left;"> cluster </th>
   <th style="text-align:left;"> gene </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MUC4 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.794845 </td>
   <td style="text-align:right;"> 0.845 </td>
   <td style="text-align:right;"> 0.465 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:left;"> MUC4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SRRM2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.068336 </td>
   <td style="text-align:right;"> 0.845 </td>
   <td style="text-align:right;"> 0.673 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:left;"> SRRM2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DUOX2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.299133 </td>
   <td style="text-align:right;"> 0.282 </td>
   <td style="text-align:right;"> 0.074 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:left;"> DUOX2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GCC2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.019849 </td>
   <td style="text-align:right;"> 0.800 </td>
   <td style="text-align:right;"> 0.559 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:left;"> GCC2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XAF1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.041879 </td>
   <td style="text-align:right;"> 0.373 </td>
   <td style="text-align:right;"> 0.124 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:left;"> XAF1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CEACAM5 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.077321 </td>
   <td style="text-align:right;"> 0.782 </td>
   <td style="text-align:right;"> 0.521 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:left;"> CEACAM5 </td>
  </tr>
</tbody>
</table>

```r
view.markers <- tapply(markers$gene, markers$cluster, function(x){head(x,1)})
# violin plots
lapply(view.markers, function(marker){
  VlnPlot(experiment.aggregate,
          group.by = "subcluster",
          features = marker) +
    scale_fill_viridis_d(option = "turbo")
})
```

<div class='r_output'> $`5_2`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-1.png)<!-- -->

<div class='r_output'> 
 $`10`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-2.png)<!-- -->

<div class='r_output'> 
 $`5_1`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-3.png)<!-- -->

<div class='r_output'> 
 $`3`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-4.png)<!-- -->

<div class='r_output'> 
 $`8`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-5.png)<!-- -->

<div class='r_output'> 
 $`2`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-6.png)<!-- -->

<div class='r_output'> 
 $`7`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-7.png)<!-- -->

<div class='r_output'> 
 $`1`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-8.png)<!-- -->

<div class='r_output'> 
 $`4`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-9.png)<!-- -->

<div class='r_output'> 
 $`6`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-10.png)<!-- -->

<div class='r_output'> 
 $`0`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-11.png)<!-- -->

<div class='r_output'> 
 $`13`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-12.png)<!-- -->

<div class='r_output'> 
 $`14`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-13.png)<!-- -->

<div class='r_output'> 
 $`9`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-14.png)<!-- -->

<div class='r_output'> 
 $`5_0`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-15.png)<!-- -->

<div class='r_output'> 
 $`5_3`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-16.png)<!-- -->

```r
# feature plots
lapply(view.markers, function(marker){
  FeaturePlot(experiment.aggregate,
              features = marker)
})
```

<div class='r_output'> $`5_2`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-17.png)<!-- -->

<div class='r_output'> 
 $`10`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-18.png)<!-- -->

<div class='r_output'> 
 $`5_1`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-19.png)<!-- -->

<div class='r_output'> 
 $`3`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-20.png)<!-- -->

<div class='r_output'> 
 $`8`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-21.png)<!-- -->

<div class='r_output'> 
 $`2`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-22.png)<!-- -->

<div class='r_output'> 
 $`7`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-23.png)<!-- -->

<div class='r_output'> 
 $`1`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-24.png)<!-- -->

<div class='r_output'> 
 $`4`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-25.png)<!-- -->

<div class='r_output'> 
 $`6`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-26.png)<!-- -->

<div class='r_output'> 
 $`0`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-27.png)<!-- -->

<div class='r_output'> 
 $`13`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-28.png)<!-- -->

<div class='r_output'> 
 $`14`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-29.png)<!-- -->

<div class='r_output'> 
 $`9`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-30.png)<!-- -->

<div class='r_output'> 
 $`5_0`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-31.png)<!-- -->

<div class='r_output'> 
 $`5_3`
</div>
![](06-clustering_celltype_files/figure-html/FindAllMarkers-32.png)<!-- -->

```r
# heat map
DoHeatmap(experiment.aggregate,
          group.by = "subcluster",
          features = view.markers,
          group.colors = viridis::turbo(length(unique(experiment.aggregate$subcluster))))
```

![](06-clustering_celltype_files/figure-html/FindAllMarkers-33.png)<!-- -->

#### Calculate mean marker expression within clusters
You may want to get an idea of the mean expression of markers in a cluster or group of clusters. The percent expressing is provided by FindMarkers and FindAllMarkers, along with the average log fold change, but not the expression value itself. The function below calculates a mean for the supplied marker in the named cluster(s) and all other groups. Please note that this function accesses the active identity.

```r
# ensure active identity is set to desired clustering resolution
Idents(experiment.aggregate) <- experiment.aggregate$subcluster
# define function
getGeneClusterMeans <- function(feature, idents){
  x = GetAssayData(experiment.aggregate)[feature,]
  m = tapply(x, Idents(experiment.aggregate) %in% idents, mean)
  names(m) = c("mean.out.of.idents", "mean.in.idents")
  return(m[c(2,1)])
}
# calculate means for a single marker
getGeneClusterMeans("SMOC2", c("1"))
```

<div class='r_output'>     mean.in.idents mean.out.of.idents 
         0.06728279         0.27537976
</div>
```r
# add means to marker table (example using subset)
markers.small <- markers[view.markers,]
means <- matrix(mapply(getGeneClusterMeans, view.markers, markers.small$cluster), ncol = 2, byrow = TRUE)
colnames(means) <- c("mean.in.cluster", "mean.out.of.cluster")
rownames(means) <- view.markers
markers.small <- cbind(markers.small, means)
markers.small[,c("cluster", "mean.in.cluster", "mean.out.of.cluster", "avg_log2FC", "p_val_adj")] %>%
  kable() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> cluster </th>
   <th style="text-align:right;"> mean.in.cluster </th>
   <th style="text-align:right;"> mean.out.of.cluster </th>
   <th style="text-align:right;"> avg_log2FC </th>
   <th style="text-align:right;"> p_val_adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MUC4 </td>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:right;"> 2.3939995 </td>
   <td style="text-align:right;"> 0.9712338 </td>
   <td style="text-align:right;"> 1.7948451 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ARHGAP15 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:right;"> 2.0999714 </td>
   <td style="text-align:right;"> 0.0420769 </td>
   <td style="text-align:right;"> 3.8168553 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC9A3 </td>
   <td style="text-align:left;"> 5_1 </td>
   <td style="text-align:right;"> 0.4786179 </td>
   <td style="text-align:right;"> 0.0157107 </td>
   <td style="text-align:right;"> 1.4194496 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AC012494.1 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 1.9829173 </td>
   <td style="text-align:right;"> 0.0198549 </td>
   <td style="text-align:right;"> 3.4055733 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CALD1 </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:right;"> 2.4889292 </td>
   <td style="text-align:right;"> 0.0321684 </td>
   <td style="text-align:right;"> 4.2588299 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CEMIP </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 1.7103302 </td>
   <td style="text-align:right;"> 0.1470160 </td>
   <td style="text-align:right;"> 2.9158868 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DIAPH3 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 0.7078035 </td>
   <td style="text-align:right;"> 0.1208302 </td>
   <td style="text-align:right;"> 1.3733030 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KCNMA1 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 2.9322926 </td>
   <td style="text-align:right;"> 0.1806122 </td>
   <td style="text-align:right;"> 3.9378447 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SCNN1B </td>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 2.7636801 </td>
   <td style="text-align:right;"> 0.8965142 </td>
   <td style="text-align:right;"> 1.9240480 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GDA </td>
   <td style="text-align:left;"> 5_1 </td>
   <td style="text-align:right;"> 1.4183410 </td>
   <td style="text-align:right;"> 0.6512665 </td>
   <td style="text-align:right;"> 0.8596434 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBFOX1 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:right;"> 2.4620327 </td>
   <td style="text-align:right;"> 0.8939034 </td>
   <td style="text-align:right;"> 1.8115188 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LRMP </td>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:right;"> 3.6264698 </td>
   <td style="text-align:right;"> 0.0292279 </td>
   <td style="text-align:right;"> 5.4904492 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KCNB2 </td>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:right;"> 1.5533307 </td>
   <td style="text-align:right;"> 0.0185579 </td>
   <td style="text-align:right;"> 3.6169997 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CA7 </td>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:right;"> 2.7325679 </td>
   <td style="text-align:right;"> 0.0402416 </td>
   <td style="text-align:right;"> 4.5910162 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LINC00278 </td>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:right;"> 1.0227013 </td>
   <td style="text-align:right;"> 0.3394740 </td>
   <td style="text-align:right;"> 1.3526996 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELENBP1 </td>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 2.4528770 </td>
   <td style="text-align:right;"> 0.9125702 </td>
   <td style="text-align:right;"> 1.6300180 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

## Advanced visualizations

Researchers may use the tree, markers, domain knowledge, and goals to identify useful clusters. This may mean adjusting PCA to use, choosing a new resolution, merging clusters together, sub-clustering, sub-setting, etc. You may also want to use automated cell type identification at this point, which will be discussed in the next section.

### Address overplotting
Single cell and single nucleus experiments may include so many cells that dimensionality reduction plots sometimes suffer from overplotting, where individual points are difficult to see. The following code addresses this by adjusting the size and opacity of the points.

```r
alpha.use <- 0.4
p <- DimPlot(experiment.aggregate,
             group.by="subcluster",
             pt.size=0.5,
             reduction = "umap",
             shuffle = TRUE) +
  scale_color_viridis_d(option = "turbo")
p$layers[[1]]$mapping$alpha <- alpha.use
p + scale_alpha_continuous(range = alpha.use, guide = FALSE)
```

![](06-clustering_celltype_files/figure-html/alpha-1.png)<!-- -->


### Highlight a subset


```r
DimPlot(experiment.aggregate,
        group.by = "subcluster",
        cells.highlight = CellsByIdentities(experiment.aggregate, idents = c("1", "2", "3")),
        cols.highlight = c(viridis::viridis(3))) +
  ggtitle("Selected clusters")
```

![](06-clustering_celltype_files/figure-html/highlight-1.png)<!-- -->

### Split dimensionality reduction plots


```r
DimPlot(experiment.aggregate,
        group.by = "subcluster",
        split.by = "group") +
  scale_color_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/split-1.png)<!-- -->

### Plot a subset of cells

Note that the object itself is unchanged by the subsetting operation.

```r
DimPlot(experiment.aggregate,
        group.by = "subcluster",
        reduction = "umap",
        cells = Cells(experiment.aggregate)[experiment.aggregate$orig.ident %in% "A001-C-007"]) +
  scale_color_viridis_d(option = "turbo") +
  ggtitle("A00-C-007 subcluster")
```

![](06-clustering_celltype_files/figure-html/plot.subset-1.png)<!-- -->

## Repeat with integrated data
But what about the integrated object? Let's import that one and repeat the clustering workflow.

```r
experiment.integrated <- readRDS("scRNA_workshop-05.rds")
experiment.integrated
```

<div class='r_output'> An object of class Seurat 
 13292 features across 6312 samples within 2 assays 
 Active assay: integrated (2000 features, 2000 variable features)
  1 other assay present: RNA
  2 dimensional reductions calculated: pca, umap
</div>
```r
experiment.integrated <- FindNeighbors(experiment.integrated, reduction = "pca", dims = 50)
experiment.integrated <- FindClusters(experiment.integrated,
                                     resolution = seq(0.1, 0.4, 0.1))
```

<div class='r_output'> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 6312
 Number of edges: 105574
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.9893
 Number of communities: 17
 Elapsed time: 0 seconds
 Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 6312
 Number of edges: 105574
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.9843
 Number of communities: 25
 Elapsed time: 0 seconds
 Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 6312
 Number of edges: 105574
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.9802
 Number of communities: 28
 Elapsed time: 0 seconds
 Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
 
 Number of nodes: 6312
 Number of edges: 105574
 
 Running Louvain algorithm...
 Maximum modularity in 10 random starts: 0.9767
 Number of communities: 30
 Elapsed time: 0 seconds
</div>
```r
cluster.resolutions <- grep("res", colnames(experiment.integrated@meta.data), value = TRUE)
sapply(cluster.resolutions, function(res){
  length(levels(experiment.integrated@meta.data[,res]))
})
```

<div class='r_output'> integrated_snn_res.0.1 integrated_snn_res.0.2 integrated_snn_res.0.3 
                     17                     25                     28 
 integrated_snn_res.0.4 
                     30
</div>
```r
lapply(cluster.resolutions, function(res){
  DimPlot(experiment.integrated,
          group.by = res,
          reduction = "umap",
          shuffle = TRUE) +
    scale_color_viridis_d(option = "turbo")
})
```

<div class='r_output'> [[1]]
</div>
![](06-clustering_celltype_files/figure-html/integrated-1.png)<!-- -->

<div class='r_output'> 
 [[2]]
</div>
![](06-clustering_celltype_files/figure-html/integrated-2.png)<!-- -->

<div class='r_output'> 
 [[3]]
</div>
![](06-clustering_celltype_files/figure-html/integrated-3.png)<!-- -->

<div class='r_output'> 
 [[4]]
</div>
![](06-clustering_celltype_files/figure-html/integrated-4.png)<!-- -->

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

<div class='r_output'> [[1]]
</div>
![](06-clustering_celltype_files/figure-html/integrated-5.png)<!-- -->

<div class='r_output'> 
 [[2]]
</div>
![](06-clustering_celltype_files/figure-html/integrated-6.png)<!-- -->

<div class='r_output'> 
 [[3]]
</div>
![](06-clustering_celltype_files/figure-html/integrated-7.png)<!-- -->

<div class='r_output'> 
 [[4]]
</div>
![](06-clustering_celltype_files/figure-html/integrated-8.png)<!-- -->

```r
FeaturePlot(experiment.integrated,
            reduction = "umap",
            features = "KCNMA1")
```

![](06-clustering_celltype_files/figure-html/integrated-9.png)<!-- -->

```r
VlnPlot(experiment.integrated,
        group.by = "integrated_snn_res.0.1",
        features = "KCNMA1",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/integrated-10.png)<!-- -->

```r
Idents(experiment.integrated) <- "integrated_snn_res.0.1"
experiment.integrated <- BuildClusterTree(experiment.integrated, dims = 1:50)
PlotClusterTree(experiment.integrated)
```

![](06-clustering_celltype_files/figure-html/integrated-11.png)<!-- -->

```r
experiment.integrated <- RenameIdents(experiment.integrated,
                                     '4' = '1',
                                     '16' = '1',
                                     '11' = '1',
                                     '10' = '1')
experiment.integrated$res.0.1_merged <- Idents(experiment.integrated)
table(experiment.integrated$res.0.1_merged)
```

<div class='r_output'> 
    1    0    2    3    5    6    7    8    9   12   13   14   15 
 1778  603  501  498  454  419  366  352  313  270  269  252  237
</div>
```r
new.order <- as.character(sort(as.numeric(levels(experiment.integrated$res.0.1_merged))))
experiment.integrated$res.0.1_merged <- factor(experiment.integrated$res.0.1_merged, levels = new.order)
DimPlot(experiment.integrated,
        reduction = "umap",
        group.by = "res.0.1_merged",
        label = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](06-clustering_celltype_files/figure-html/integrated-12.png)<!-- -->

In this case, integration appears to have interfered with forming easily-interpreted clusters. There is very little relationship between location on the UMAP and cluster identity, which makes it harder to identify possible cell populations at a glance. We can add the integrated object's cluster identities to the un-integrated object if we choose, though projecting the integrated clustering onto the un-integrated UMAP is unlikely to be useful.

```r
identical(rownames(experiment.aggregate@meta.data), rownames(experiment.integrated@meta.data))
```

<div class='r_output'> [1] TRUE
</div>
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

![](06-clustering_celltype_files/figure-html/AddMetaData-1.png)<!-- -->

## Automated cell type detection
[ScType](https://www.nature.com/articles/s41467-022-28803-w) is one of many available automated cell type detection algorithms. It has the advantage of being fast and flexible; it can be used with the large human and mouse cell type database provided by the authors, with a user-defined database, or with some combination of the two.

In this section, we will use ScType to assign our clusters from Seurat to a cell type based on a hierarchical external database. The database supplied with the package is human but users can supply their own data. More details are available on [Github](https://github.com/IanevskiAleksandr/sc-type).

### Source ScType functions from Github
The `source()` functions in the code box below run scripts stored in the ScType GitHub repository. These scripts add two additional functions, `gene_sets_prepare()` and `sctype_score()`, to the global environment.

```r
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

### Create marker database
The authors of ScType have provided an extensive database of human and mouse cell type markers in the following tissues:

* Immune system
* Pancreas
* Liver
* Eye
* Kidney
* Brain
* Lung
* Adrenal
* Heart
* Intestine
* Muscle
* Placenta
* Spleen
* Stomach
* Thymus

To set up the database for ScType, we run the `gene_sets_prepare()` function on a URL pointing to the full database, specifying the tissue subset to retrieve.

```r
db.URL <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- "Intestine"
gs.list <- gene_sets_prepare(db.URL, cell_type = "Intestine")
```


```r
View(gs.list)
```

The database is composed of two lists (negative and positive) of cell type markers for the specified tissue. In this case, no negative markers have been provided, so the vectors of gene names in that part of the database are empty (length 0).

In addition to using the database supplied with the package, or substituting another database, it is possible to augment the provided database by changing the provided markers for an existing cell type, or adding another cell type to the object.

Let's add a cell type and associated markers (markers from [PanglaoDB](https://panglaodb.se/markers.html?cell_type=%27Tuft%20cells%27#google_vignette)).

```r
gs.list$gs_positive$`Tuft cells` <- c("SUCNR1", "FABP1", "POU2F3", "SIGLECF", "CDHR2", "AVIL", "ESPN", "LRMP", "TRPM5", "DCLK1", "TAS1R3", "SOX9", "TUBB5", "CAMK2B", "GNAT3", "IL25", "PLCB2", "GFI1B", "ATOH1", "CD24A", "ASIC5", "KLF3", "KLF6", "DRD3", "NRADD", "GNG13", "NREP", "RGS2", "RAC2", "PTGS1", "IRF7", "FFAR3", "ALOX5", "TSLP", "IL4RA", "IL13RA1", "IL17RB", "PTPRC")
gs.list$gs_negative$`Tuft cells` <- NULL
```

### Score cells
ScType scores every cell, summarizing the transformed and weighted expression of markers for each cell type. This generates a matrix with the dimensions cell types x cells.

```r
es.max <- sctype_score(GetAssayData(experiment.aggregate, "scale"),
                       scaled = TRUE,
                       gs = gs.list$gs_positive,
                       gs2 = gs.list$gs_negative)
# cell type scores of first cell
es.max[,1]
```

<div class='r_output'>         Smooth muscle cells              Lymphoid cells 
                  -0.3240119                  -0.1071450 
 Intestinal epithelial cells                    ENS glia 
                  -0.4512850                  -0.1365531 
  Vascular endothelial cells               Stromal cells 
                  -0.1035601                  -0.1551877 
 Lymphatic endothelial cells                  Tuft cells 
                   2.1424191                   0.6259383
</div>
### Score clusters
The "ScType score" of each cluster is calculated by summing the cell-level scores within each cluster. The cell type with the largest (positive) score is the most highly enriched cell type for that cluster.

```r
clusters <- sort(unique(experiment.aggregate$subcluster))
tmp <- lapply(clusters, function(cluster){
  es.max.cluster = sort(rowSums(es.max[, experiment.aggregate$subcluster == cluster]), decreasing = TRUE)
  out = head(data.frame(subcluster = cluster,
                        ScType = names(es.max.cluster),
                        scores = es.max.cluster,
                        ncells = sum(experiment.aggregate$subcluster == cluster)))
  out$rank = 1:length(out$scores)
  return(out)
})
cluster.ScType <- do.call("rbind", tmp)
cluster.ScType %>%
  pivot_wider(id_cols = subcluster,
              names_from = ScType,
              values_from = scores) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped"), fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> subcluster </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Lymphatic endothelial cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ENS glia </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Lymphoid cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Vascular endothelial cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Stromal cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Smooth muscle cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Tuft cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Intestinal epithelial cells </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:right;"> 78.730449 </td>
   <td style="text-align:right;"> -5.668744 </td>
   <td style="text-align:right;"> -101.454919 </td>
   <td style="text-align:right;"> -110.558849 </td>
   <td style="text-align:right;"> -127.199691 </td>
   <td style="text-align:right;"> -128.045735 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -13.108991 </td>
   <td style="text-align:right;"> -93.884079 </td>
   <td style="text-align:right;"> -57.191352 </td>
   <td style="text-align:right;"> -92.148994 </td>
   <td style="text-align:right;"> 158.397337 </td>
   <td style="text-align:right;"> 231.59090 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:right;"> -13.871817 </td>
   <td style="text-align:right;"> 21.532493 </td>
   <td style="text-align:right;"> 332.447412 </td>
   <td style="text-align:right;"> 28.440689 </td>
   <td style="text-align:right;"> -20.085597 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 185.12259 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 17.856316 </td>
   <td style="text-align:right;"> 1.837349 </td>
   <td style="text-align:right;"> -3.655264 </td>
   <td style="text-align:right;"> -5.749935 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 331.44774 </td>
   <td style="text-align:right;"> -2.38527 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14 </td>
   <td style="text-align:right;"> -3.731304 </td>
   <td style="text-align:right;"> -3.296859 </td>
   <td style="text-align:right;"> -2.874269 </td>
   <td style="text-align:right;"> -2.554377 </td>
   <td style="text-align:right;"> -3.525360 </td>
   <td style="text-align:right;"> -6.504726 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> -21.761310 </td>
   <td style="text-align:right;"> 41.101170 </td>
   <td style="text-align:right;"> -32.580049 </td>
   <td style="text-align:right;"> -42.047865 </td>
   <td style="text-align:right;"> -92.753617 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -39.12461 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -13.451515 </td>
   <td style="text-align:right;"> -18.982844 </td>
   <td style="text-align:right;"> 29.052992 </td>
   <td style="text-align:right;"> -33.307831 </td>
   <td style="text-align:right;"> -41.194337 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 136.11698 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 124.674670 </td>
   <td style="text-align:right;"> -57.587776 </td>
   <td style="text-align:right;"> -33.265003 </td>
   <td style="text-align:right;"> -46.651246 </td>
   <td style="text-align:right;"> -37.033475 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -42.59110 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_0 </td>
   <td style="text-align:right;"> 5.900122 </td>
   <td style="text-align:right;"> 9.517507 </td>
   <td style="text-align:right;"> -28.617150 </td>
   <td style="text-align:right;"> -33.203823 </td>
   <td style="text-align:right;"> -41.819751 </td>
   <td style="text-align:right;"> -84.907883 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_1 </td>
   <td style="text-align:right;"> 11.163123 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 5.194384 </td>
   <td style="text-align:right;"> -8.194288 </td>
   <td style="text-align:right;"> -8.498858 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 26.64893 </td>
   <td style="text-align:right;"> 81.12562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:right;"> -17.674720 </td>
   <td style="text-align:right;"> -7.856129 </td>
   <td style="text-align:right;"> -7.735188 </td>
   <td style="text-align:right;"> -10.866366 </td>
   <td style="text-align:right;"> -3.486374 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -3.06727 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_3 </td>
   <td style="text-align:right;"> 23.029373 </td>
   <td style="text-align:right;"> -8.269773 </td>
   <td style="text-align:right;"> -6.931117 </td>
   <td style="text-align:right;"> -5.935493 </td>
   <td style="text-align:right;"> -8.491533 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -15.47474 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:right;"> 84.443721 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -15.401331 </td>
   <td style="text-align:right;"> -37.998667 </td>
   <td style="text-align:right;"> -30.351971 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 242.82685 </td>
   <td style="text-align:right;"> 261.76646 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:right;"> 42.524840 </td>
   <td style="text-align:right;"> 15.451665 </td>
   <td style="text-align:right;"> -28.274959 </td>
   <td style="text-align:right;"> -11.967535 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 5.926533 </td>
   <td style="text-align:right;"> 25.21208 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:right;"> 41.449012 </td>
   <td style="text-align:right;"> 38.372473 </td>
   <td style="text-align:right;"> -3.534333 </td>
   <td style="text-align:right;"> 248.335004 </td>
   <td style="text-align:right;"> 525.323465 </td>
   <td style="text-align:right;"> 506.841818 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 10.187099 </td>
   <td style="text-align:right;"> -11.246702 </td>
   <td style="text-align:right;"> -2.545120 </td>
   <td style="text-align:right;"> -17.273971 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -30.03886 </td>
   <td style="text-align:right;"> 541.82993 </td>
  </tr>
</tbody>
</table>

```r
cluster.ScType.top <- cluster.ScType %>%
  filter(rank == 1) %>%
  select(subcluster, ncells, ScType, scores)
rownames(cluster.ScType.top) <- cluster.ScType.top$subcluster
cluster.ScType.top <- cluster.ScType.top %>%
  rename("score" = scores)
```

Once the scores have been calculated for each cluster, they can be added to the Seurat object.

```r
subcluster.ScType <- cluster.ScType.top[experiment.aggregate$subcluster, "ScType"]
experiment.aggregate <- AddMetaData(experiment.aggregate,
                                    metadata = subcluster.ScType,
                                    col.name = "subcluster_ScType")
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "subcluster_ScType",
        shuffle = TRUE) +
  scale_color_viridis_d()
```

![](06-clustering_celltype_files/figure-html/add_ScType-1.png)<!-- -->

The ScType developer suggests that assignments with a score less than (number of cells in cluster)/4 are low confidence and should be set to unknown:

```r
cluster.ScType.top <- cluster.ScType.top %>%
  mutate(ScType_filtered = ifelse(score >= ncells / 4, ScType, "Unkown"))
subcluster.ScType.filtered <- cluster.ScType.top[experiment.aggregate$subcluster, "ScType_filtered"]
experiment.aggregate <- AddMetaData(experiment.aggregate,
                                    metadata = subcluster.ScType.filtered,
                                    col.name = "subcluster_ScType_filtered")
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "subcluster_ScType_filtered",
        shuffle = TRUE) +
  scale_color_viridis_d()
```

![](06-clustering_celltype_files/figure-html/adjust_confidence-1.png)<!-- -->

## Prepare for the next section

#### Save the Seurat object and download the next Rmd file

```r
# set the finalcluster to subcluster
experiment.aggregate$finalcluster <- experiment.aggregate$subcluster
# save object
saveRDS(experiment.aggregate, file="scRNA_workshop-06.rds")
```

#### Download the Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/07-doublet_detection.Rmd", "07-doublet_detection.Rmd")
```

#### Session Information

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
 [1] HGNChelper_0.8.1   ggplot2_3.4.2      dplyr_1.1.2        tidyr_1.3.0       
 [5] kableExtra_1.3.4   SeuratObject_4.1.3 Seurat_4.3.0.1    
 
 loaded via a namespace (and not attached):
   [1] RColorBrewer_1.1-3     rstudioapi_0.15.0      jsonlite_1.8.7        
   [4] magrittr_2.0.3         spatstat.utils_3.0-3   farver_2.1.1          
   [7] rmarkdown_2.23         vctrs_0.6.3            ROCR_1.0-11           
  [10] spatstat.explore_3.2-1 webshot_0.5.5          htmltools_0.5.5       
  [13] sass_0.4.7             sctransform_0.3.5      parallelly_1.36.0     
  [16] KernSmooth_2.23-22     bslib_0.5.0            htmlwidgets_1.6.2     
  [19] ica_1.0-3              plyr_1.8.8             plotly_4.10.2         
  [22] zoo_1.8-12             cachem_1.0.8           igraph_1.5.0          
  [25] mime_0.12              lifecycle_1.0.3        pkgconfig_2.0.3       
  [28] Matrix_1.6-0           R6_2.5.1               fastmap_1.1.1         
  [31] fitdistrplus_1.1-11    future_1.33.0          shiny_1.7.4.1         
  [34] digest_0.6.33          colorspace_2.1-0       patchwork_1.1.2       
  [37] tensor_1.5             irlba_2.3.5.1          labeling_0.4.2        
  [40] progressr_0.13.0       fansi_1.0.4            spatstat.sparse_3.0-2 
  [43] httr_1.4.6             polyclip_1.10-4        abind_1.4-5           
  [46] compiler_4.3.1         withr_2.5.0            viridis_0.6.3         
  [49] highr_0.10             MASS_7.3-60            tools_4.3.1           
  [52] lmtest_0.9-40          ape_5.7-1              zip_2.3.0             
  [55] httpuv_1.6.11          future.apply_1.11.0    goftest_1.2-3         
  [58] glue_1.6.2             nlme_3.1-162           promises_1.2.0.1      
  [61] grid_4.3.1             Rtsne_0.16             cluster_2.1.4         
  [64] reshape2_1.4.4         generics_0.1.3         gtable_0.3.3          
  [67] spatstat.data_3.0-1    data.table_1.14.8      sp_2.0-0              
  [70] xml2_1.3.5             utf8_1.2.3             spatstat.geom_3.2-4   
  [73] RcppAnnoy_0.0.21       ggrepel_0.9.3          RANN_2.6.1            
  [76] pillar_1.9.0           stringr_1.5.0          limma_3.56.2          
  [79] later_1.3.1            splines_4.3.1          lattice_0.21-8        
  [82] survival_3.5-5         deldir_1.0-9           tidyselect_1.2.0      
  [85] miniUI_0.1.1.1         pbapply_1.7-2          knitr_1.43            
  [88] gridExtra_2.3          svglite_2.1.1          scattermore_1.2       
  [91] xfun_0.39              matrixStats_1.0.0      stringi_1.7.12        
  [94] lazyeval_0.2.2         yaml_2.3.7             evaluate_0.21         
  [97] codetools_0.2-19       tibble_3.2.1           cli_3.6.1             
 [100] uwot_0.1.16            xtable_1.8-4           reticulate_1.30       
 [103] systemfonts_1.0.4      munsell_0.5.0          jquerylib_0.1.4       
 [106] Rcpp_1.0.11            globals_0.16.2         spatstat.random_3.1-5 
 [109] png_0.1-8              parallel_4.3.1         ellipsis_0.3.2        
 [112] listenv_0.9.0          viridisLite_0.4.2      scales_1.2.1          
 [115] ggridges_0.5.4         openxlsx_4.2.5.2       leiden_0.4.3          
 [118] purrr_1.0.1            rlang_1.1.1            cowplot_1.1.1         
 [121] rvest_1.0.3
</div>