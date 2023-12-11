---
title: "Introduction to Single Cell RNA-Seq Part 5: Clustering and cell type assignment"
author: "UCD Bioinformatics Core"
date: "2023-12-10"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 5: Clustering and cell type assignment
Clustering and cell type assignment are critical steps in many single cell (or single nucleus) experiments. The power of single cell experiments is to capture the heterogeneity within samples as well as between them. Clustering permits the user to organize cells into clusters that correspond to biologically and experimentally relevant populations.


## Set up workspace

```r
library(Seurat)
library(kableExtra)
library(tidyr)
library(dplyr)
library(ggplot2)
library(HGNChelper)
library(ComplexHeatmap)
set.seed(12345)
experiment.aggregate <- readRDS("scRNA_workshop-04.rds")
experiment.aggregate
```

```
## An object of class Seurat 
## 11292 features across 6313 samples within 1 assay 
## Active assay: RNA (11292 features, 7012 variable features)
##  2 dimensional reductions calculated: pca, umap
```

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

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6313
## Number of edges: 265916
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9680
## Number of communities: 10
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6313
## Number of edges: 265916
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9484
## Number of communities: 11
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6313
## Number of edges: 265916
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9339
## Number of communities: 13
## Elapsed time: 0 seconds
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 6313
## Number of edges: 265916
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9217
## Number of communities: 15
## Elapsed time: 0 seconds
```

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
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> 9 </td>
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
   <td style="text-align:left;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAAGGATTCATTACCT_A001-C-007 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 4 </td>
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

```
## RNA_snn_res.0.1 RNA_snn_res.0.2 RNA_snn_res.0.3 RNA_snn_res.0.4 
##              10              11              13              15
```

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

```
## [[1]]
```

![](05-clustering_celltype_files/figure-html/UMAP-1.png)<!-- -->

```
## 
## [[2]]
```

![](05-clustering_celltype_files/figure-html/UMAP-2.png)<!-- -->

```
## 
## [[3]]
```

![](05-clustering_celltype_files/figure-html/UMAP-3.png)<!-- -->

```
## 
## [[4]]
```

![](05-clustering_celltype_files/figure-html/UMAP-4.png)<!-- -->

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

```
## [[1]]
```

![](05-clustering_celltype_files/figure-html/membership-1.png)<!-- -->

```
## 
## [[2]]
```

![](05-clustering_celltype_files/figure-html/membership-2.png)<!-- -->

```
## 
## [[3]]
```

![](05-clustering_celltype_files/figure-html/membership-3.png)<!-- -->

```
## 
## [[4]]
```

![](05-clustering_celltype_files/figure-html/membership-4.png)<!-- -->

### Investigate the relationship between cluster identity and metadata values
Here, example plots are displayed for the lowest resolution in order to save space. To see plots for each resolution, use `lapply()`.

```r
VlnPlot(experiment.aggregate,
        group.by = "RNA_snn_res.0.4",
        features = "nCount_RNA",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
VlnPlot(experiment.aggregate,
        group.by = "RNA_snn_res.0.4",
        features = "nFeature_RNA",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

```r
VlnPlot(experiment.aggregate,
        group.by = "RNA_snn_res.0.4",
        features = "percent_MT",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/unnamed-chunk-1-3.png)<!-- -->
 
### Visualize expression of genes of interest


```r
FeaturePlot(experiment.aggregate,
            reduction = "umap",
            features = "KCNMA1")
```

![](05-clustering_celltype_files/figure-html/feature-1.png)<!-- -->

```r
VlnPlot(experiment.aggregate,
        group.by = "RNA_snn_res.0.4",
        features = "KCNMA1",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/feature-2.png)<!-- -->

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

![](05-clustering_celltype_files/figure-html/tree-1.png)<!-- -->


## Merge clusters
In many experiments, the clustering resolution does not need to be uniform across all of the cell types present. While for some cell types of interest fine detail may be desirable, for others, simply grouping them into a larger parent cluster is sufficient. Merging cluster is very straightforward.

```r
experiment.aggregate <- RenameIdents(experiment.aggregate,
                                     '11' = '10',
                                     '12' = '13')
experiment.aggregate$res.0.4_merged <- Idents(experiment.aggregate)
table(experiment.aggregate$res.0.4_merged)
```

```
## 
##   10   13    0    1    2    3    4    5    6    7    8    9   14 
##  107   58 1391 1009  916  657  655  613  413  189  180  102   23
```

```r
experiment.aggregate@meta.data %>%
  ggplot(aes(x = res.0.4_merged, fill = group)) +
  geom_bar() +
  scale_fill_viridis_d(option = "mako") +
  theme_classic()
```

![](05-clustering_celltype_files/figure-html/merge-1.png)<!-- -->

```r
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "res.0.4_merged",
        label = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/merge-2.png)<!-- -->

```r
VlnPlot(experiment.aggregate,
        group.by = "res.0.4_merged",
        features = "KCNMA1") +
  scale_fill_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/merge-3.png)<!-- -->

## Reorder the clusters
Merging the clusters changed the order in which they appear on a plot. In order to reorder the clusters for plotting purposes take a look at the levels of the identity, then re-level as desired.

```r
levels(experiment.aggregate$res.0.4_merged)
```

```
##  [1] "10" "13" "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "14"
```

```r
# move one cluster to the first position
experiment.aggregate$res.0.4_merged <- relevel(experiment.aggregate$res.0.4_merged, "0")
levels(experiment.aggregate$res.0.4_merged)
```

```
##  [1] "0"  "10" "13" "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "14"
```

```r
# the color assigned to some clusters will change
VlnPlot(experiment.aggregate,
        group.by = "res.0.4_merged",
        features = "CAPN9",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/relevel-1.png)<!-- -->

```r
# re-level entire factor
new.order <- as.character(sort(as.numeric(levels(experiment.aggregate$res.0.4_merged))))
experiment.aggregate$res.0.4_merged <- factor(experiment.aggregate$res.0.4_merged, levels = new.order)
levels(experiment.aggregate$res.0.4_merged)
```

```
##  [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "13" "14"
```

```r
VlnPlot(experiment.aggregate,
        group.by = "res.0.4_merged",
        features = "MYRIP",
        pt.size = 0.1) +
  scale_fill_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/relevel-2.png)<!-- -->

## Subcluster
While merging clusters reduces the resolution in some parts of the experiment, sub-clustering has the opposite effect. Let's produce sub-clusters for cluster 5.

```r
experiment.aggregate <- FindSubCluster(experiment.aggregate,
                                       graph.name = "RNA_snn",
                                       cluster = 5,
                                       subcluster.name = "subcluster")
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 613
## Number of edges: 23846
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.7230
## Number of communities: 4
## Elapsed time: 0 seconds
```

```r
experiment.aggregate$subcluster <- factor(experiment.aggregate$subcluster,
                                          levels = c(as.character(0:4),
                                                     "5_0", "5_1", "5_2", "5_3",
                                                     as.character(c(6:10, 13, 14))))
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "subcluster",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/subcluster-1.png)<!-- -->

```r
sort(unique(experiment.aggregate$subcluster))
```

```
##  [1] 0   1   2   3   4   5_0 5_1 5_2 5_3 6   7   8   9   10  13  14 
## Levels: 0 1 2 3 4 5_0 5_1 5_2 5_3 6 7 8 9 10 13 14
```

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

![](05-clustering_celltype_files/figure-html/subset-1.png)<!-- -->

```r
# retain cells belonging only to specified clusters
experiment.tmp <- subset(experiment.aggregate, subcluster %in% c("1", "2", "5_0", "5_1", "5_2", "5_3", "7"))
DimPlot(experiment.tmp,
        reduction = "umap",
        group.by = "subcluster",
        shuffle = TRUE) +
  scale_color_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/subset-2.png)<!-- -->

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

```
## [1] 1744
```

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
   <td style="text-align:left;"> AGBL4 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.173808 </td>
   <td style="text-align:right;"> 0.611 </td>
   <td style="text-align:right;"> 0.103 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EDAR </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.739874 </td>
   <td style="text-align:right;"> 0.383 </td>
   <td style="text-align:right;"> 0.023 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GRM8 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.186329 </td>
   <td style="text-align:right;"> 0.636 </td>
   <td style="text-align:right;"> 0.101 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CEMIP </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 3.035450 </td>
   <td style="text-align:right;"> 0.621 </td>
   <td style="text-align:right;"> 0.056 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CASC19 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.923663 </td>
   <td style="text-align:right;"> 0.532 </td>
   <td style="text-align:right;"> 0.081 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NKD1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.310743 </td>
   <td style="text-align:right;"> 0.529 </td>
   <td style="text-align:right;"> 0.085 </td>
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

```
## [[1]]
```

![](05-clustering_celltype_files/figure-html/view.markers-1.png)<!-- -->

```
## 
## [[2]]
```

![](05-clustering_celltype_files/figure-html/view.markers-2.png)<!-- -->

```r
FeaturePlot(experiment.aggregate,
            features = view.markers,
            ncol = 2)
```

![](05-clustering_celltype_files/figure-html/view.markers-3.png)<!-- -->

```r
DoHeatmap(experiment.aggregate,
          group.by = "subcluster",
          features = view.markers,
          group.colors = viridis::turbo(length(levels(experiment.aggregate$subcluster))))
```

![](05-clustering_celltype_files/figure-html/view.markers-4.png)<!-- -->
#### Improved heatmap
The Seurat `DoHeatmap` function provided by Seurat provides a convenient look at expression of selected genes. The ComplexHeatmap library generates heat maps with a much finer level of control.

```r
cluster.colors <- viridis::turbo(length(levels(experiment.aggregate$subcluster)))
names(cluster.colors) <- levels(experiment.aggregate$subcluster)
group.colors <- viridis::mako(length(levels(experiment.aggregate$group)))
names(group.colors) <- levels(experiment.aggregate$group)
top.annotation <- columnAnnotation(df = experiment.aggregate@meta.data[,c("group", "subcluster")],
                                   col = list(group = group.colors,
                                              subcluster = cluster.colors))
mat <- as.matrix(GetAssayData(experiment.aggregate[rownames(markers)[1:20],],
                              slot = "data"))
Heatmap(mat,
        name = "normalized\ncounts",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        top_annotation = top.annotation)
```

![](05-clustering_celltype_files/figure-html/ComplexHeatmap-1.png)<!-- -->


### FindAllMarkers
FindAllMarkers can be used to automate this process across all clusters.

```r
Idents(experiment.aggregate) <- "subcluster"
markers <- FindAllMarkers(experiment.aggregate,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          thresh.use = 0.25)
tapply(markers$p_val_adj, markers$cluster, function(x){
  length(x < 0.05)
})
```

```
##    0    1    2    3    4  5_0  5_1  5_2  5_3    6    7    8    9   10   13   14 
##  552  700  436  775 1187  250  749  223  333 1117  329  532  255  327  257  296
```

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
   <td style="text-align:left;"> RBFOX1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.801917 </td>
   <td style="text-align:right;"> 0.855 </td>
   <td style="text-align:right;"> 0.409 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> RBFOX1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NXPE1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.603500 </td>
   <td style="text-align:right;"> 0.958 </td>
   <td style="text-align:right;"> 0.597 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NXPE1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ADAMTSL1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.460447 </td>
   <td style="text-align:right;"> 0.911 </td>
   <td style="text-align:right;"> 0.512 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> ADAMTSL1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XIST </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.337624 </td>
   <td style="text-align:right;"> 0.850 </td>
   <td style="text-align:right;"> 0.359 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> XIST </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HNF1A-AS1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.220046 </td>
   <td style="text-align:right;"> 0.879 </td>
   <td style="text-align:right;"> 0.563 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> HNF1A-AS1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SATB2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.114086 </td>
   <td style="text-align:right;"> 0.937 </td>
   <td style="text-align:right;"> 0.665 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> SATB2 </td>
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

```
## $`0`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-1.png)<!-- -->

```
## 
## $`1`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-2.png)<!-- -->

```
## 
## $`2`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-3.png)<!-- -->

```
## 
## $`3`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-4.png)<!-- -->

```
## 
## $`4`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-5.png)<!-- -->

```
## 
## $`5_0`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-6.png)<!-- -->

```
## 
## $`5_1`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-7.png)<!-- -->

```
## 
## $`5_2`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-8.png)<!-- -->

```
## 
## $`5_3`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-9.png)<!-- -->

```
## 
## $`6`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-10.png)<!-- -->

```
## 
## $`7`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-11.png)<!-- -->

```
## 
## $`8`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-12.png)<!-- -->

```
## 
## $`9`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-13.png)<!-- -->

```
## 
## $`10`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-14.png)<!-- -->

```
## 
## $`13`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-15.png)<!-- -->

```
## 
## $`14`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-16.png)<!-- -->

```r
# feature plots
lapply(view.markers, function(marker){
  FeaturePlot(experiment.aggregate,
              features = marker)
})
```

```
## $`0`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-17.png)<!-- -->

```
## 
## $`1`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-18.png)<!-- -->

```
## 
## $`2`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-19.png)<!-- -->

```
## 
## $`3`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-20.png)<!-- -->

```
## 
## $`4`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-21.png)<!-- -->

```
## 
## $`5_0`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-22.png)<!-- -->

```
## 
## $`5_1`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-23.png)<!-- -->

```
## 
## $`5_2`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-24.png)<!-- -->

```
## 
## $`5_3`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-25.png)<!-- -->

```
## 
## $`6`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-26.png)<!-- -->

```
## 
## $`7`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-27.png)<!-- -->

```
## 
## $`8`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-28.png)<!-- -->

```
## 
## $`9`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-29.png)<!-- -->

```
## 
## $`10`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-30.png)<!-- -->

```
## 
## $`13`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-31.png)<!-- -->

```
## 
## $`14`
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers-32.png)<!-- -->

```r
# heat map
DoHeatmap(experiment.aggregate,
          group.by = "subcluster",
          features = view.markers,
          group.colors = viridis::turbo(length(unique(experiment.aggregate$subcluster))))
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers_heat-1.png)<!-- -->

```r
# ComplexHeatmap
mat <- as.matrix(GetAssayData(experiment.aggregate[view.markers,],
                              slot = "data"))
Heatmap(mat,
        name = "normalized\ncounts",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        column_split = experiment.aggregate$subcluster,
        top_annotation = top.annotation)
```

![](05-clustering_celltype_files/figure-html/FindAllMarkers_heat-2.png)<!-- -->

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

```
##     mean.in.idents mean.out.of.idents 
##          0.9897597          0.1035296
```

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
   <td style="text-align:left;"> RBFOX1 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:right;"> 2.4572798 </td>
   <td style="text-align:right;"> 0.9039857 </td>
   <td style="text-align:right;"> 1.8019166 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CEMIP </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 1.6549190 </td>
   <td style="text-align:right;"> 0.1118442 </td>
   <td style="text-align:right;"> 3.0354498 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KCNMA1 </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 2.9307725 </td>
   <td style="text-align:right;"> 0.1803269 </td>
   <td style="text-align:right;"> 3.9394788 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SCNN1B </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:right;"> 1.5462679 </td>
   <td style="text-align:right;"> 0.9425927 </td>
   <td style="text-align:right;"> 0.2633176 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AC012494.1 </td>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 1.9860483 </td>
   <td style="text-align:right;"> 0.0205298 </td>
   <td style="text-align:right;"> 3.4060180 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LINC00278 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 0.7379974 </td>
   <td style="text-align:right;"> 0.2781939 </td>
   <td style="text-align:right;"> 0.6933742 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC9A3 </td>
   <td style="text-align:left;"> 5_1 </td>
   <td style="text-align:right;"> 0.4735995 </td>
   <td style="text-align:right;"> 0.0159515 </td>
   <td style="text-align:right;"> 1.4183212 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MUC4 </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 1.7035789 </td>
   <td style="text-align:right;"> 0.8761755 </td>
   <td style="text-align:right;"> 1.0070648 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDH13 </td>
   <td style="text-align:left;"> 5_1 </td>
   <td style="text-align:right;"> 1.0662439 </td>
   <td style="text-align:right;"> 0.1748105 </td>
   <td style="text-align:right;"> 1.9017779 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GDA </td>
   <td style="text-align:left;"> 5_1 </td>
   <td style="text-align:right;"> 1.4427951 </td>
   <td style="text-align:right;"> 0.6509512 </td>
   <td style="text-align:right;"> 0.8816666 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CALD1 </td>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:right;"> 2.7875189 </td>
   <td style="text-align:right;"> 0.0426053 </td>
   <td style="text-align:right;"> 4.4588420 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CA7 </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:right;"> 2.7054238 </td>
   <td style="text-align:right;"> 0.0405927 </td>
   <td style="text-align:right;"> 4.5800381 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ARHGAP15 </td>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:right;"> 2.4875120 </td>
   <td style="text-align:right;"> 0.0532663 </td>
   <td style="text-align:right;"> 4.1066037 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PECAM1 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:right;"> 1.1234271 </td>
   <td style="text-align:right;"> 0.0086313 </td>
   <td style="text-align:right;"> 2.7698349 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KCNB2 </td>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:right;"> 2.0897079 </td>
   <td style="text-align:right;"> 0.0079379 </td>
   <td style="text-align:right;"> 4.4642037 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TIMP1 </td>
   <td style="text-align:left;"> 14 </td>
   <td style="text-align:right;"> 0.6126504 </td>
   <td style="text-align:right;"> 0.0439210 </td>
   <td style="text-align:right;"> 1.3887280 </td>
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

![](05-clustering_celltype_files/figure-html/alpha-1.png)<!-- -->


### Highlight a subset

```r
DimPlot(experiment.aggregate,
        group.by = "subcluster",
        cells.highlight = CellsByIdentities(experiment.aggregate, idents = c("1", "2", "3")),
        cols.highlight = c(viridis::viridis(3))) +
  ggtitle("Selected clusters")
```

![](05-clustering_celltype_files/figure-html/highlight-1.png)<!-- -->

### Split dimensionality reduction plots

```r
DimPlot(experiment.aggregate,
        group.by = "subcluster",
        split.by = "group") +
  scale_color_viridis_d(option = "turbo")
```

![](05-clustering_celltype_files/figure-html/split-1.png)<!-- -->

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

![](05-clustering_celltype_files/figure-html/plot.subset-1.png)<!-- -->


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

```
##         Smooth muscle cells              Lymphoid cells 
##                  -0.3242388                  -0.1071300 
## Intestinal epithelial cells                    ENS glia 
##                  -0.4535236                  -0.1363118 
##  Vascular endothelial cells               Stromal cells 
##                  -0.1040835                  -0.1545696 
## Lymphatic endothelial cells                  Tuft cells 
##                   2.1436243                   0.6228897
```

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
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Smooth muscle cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Stromal cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Intestinal epithelial cells </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Tuft cells </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:right;"> 63.3898475 </td>
   <td style="text-align:right;"> -10.989658 </td>
   <td style="text-align:right;"> -102.3898472 </td>
   <td style="text-align:right;"> -107.653784 </td>
   <td style="text-align:right;"> -121.581411 </td>
   <td style="text-align:right;"> -132.402777 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 16.2827734 </td>
   <td style="text-align:right;"> 59.078725 </td>
   <td style="text-align:right;"> -49.9188110 </td>
   <td style="text-align:right;"> -48.799502 </td>
   <td style="text-align:right;"> -86.591059 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -73.8491565 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -13.276428 </td>
   <td style="text-align:right;"> -93.9934571 </td>
   <td style="text-align:right;"> -57.218099 </td>
   <td style="text-align:right;"> 157.702670 </td>
   <td style="text-align:right;"> -92.469406 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 234.68964 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 136.8070556 </td>
   <td style="text-align:right;"> -55.682447 </td>
   <td style="text-align:right;"> -34.7475138 </td>
   <td style="text-align:right;"> -50.998283 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -33.896438 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -50.15324 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -12.641656 </td>
   <td style="text-align:right;"> -18.4670491 </td>
   <td style="text-align:right;"> 28.967538 </td>
   <td style="text-align:right;"> -40.548555 </td>
   <td style="text-align:right;"> -32.596855 </td>
   <td style="text-align:right;"> 136.5868558 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_0 </td>
   <td style="text-align:right;"> 5.6125848 </td>
   <td style="text-align:right;"> 7.635946 </td>
   <td style="text-align:right;"> -30.7812363 </td>
   <td style="text-align:right;"> -34.664302 </td>
   <td style="text-align:right;"> -85.843664 </td>
   <td style="text-align:right;"> -44.329471 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_1 </td>
   <td style="text-align:right;"> 10.0200260 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 5.4386675 </td>
   <td style="text-align:right;"> -7.988753 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -8.264155 </td>
   <td style="text-align:right;"> 82.0453510 </td>
   <td style="text-align:right;"> 25.48921 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_2 </td>
   <td style="text-align:right;"> -17.0017682 </td>
   <td style="text-align:right;"> -8.432776 </td>
   <td style="text-align:right;"> -9.1732053 </td>
   <td style="text-align:right;"> -11.226220 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -9.981352 </td>
   <td style="text-align:right;"> -1.6643326 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_3 </td>
   <td style="text-align:right;"> 25.4537227 </td>
   <td style="text-align:right;"> -1.156749 </td>
   <td style="text-align:right;"> -5.6889118 </td>
   <td style="text-align:right;"> -5.063386 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -7.022599 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -10.69587 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:right;"> 84.9356495 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -15.1909103 </td>
   <td style="text-align:right;"> -37.889131 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -30.310748 </td>
   <td style="text-align:right;"> 262.0270174 </td>
   <td style="text-align:right;"> 240.41732 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:right;"> 33.5274002 </td>
   <td style="text-align:right;"> 19.229131 </td>
   <td style="text-align:right;"> 2.4685436 </td>
   <td style="text-align:right;"> 24.056031 </td>
   <td style="text-align:right;"> 513.372529 </td>
   <td style="text-align:right;"> 538.392998 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 10.178259 </td>
   <td style="text-align:right;"> -16.7005497 </td>
   <td style="text-align:right;"> -2.699479 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -17.729215 </td>
   <td style="text-align:right;"> 538.1786695 </td>
   <td style="text-align:right;"> -30.89682 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:right;"> 0.6974037 </td>
   <td style="text-align:right;"> 28.218189 </td>
   <td style="text-align:right;"> 330.4094902 </td>
   <td style="text-align:right;"> 24.560271 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -17.921440 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 107.87937 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:right;"> -8.9357191 </td>
   <td style="text-align:right;"> 11.774279 </td>
   <td style="text-align:right;"> 0.5720622 </td>
   <td style="text-align:right;"> 227.698970 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -10.345356 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 59.46269 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:right;"> -10.8861167 </td>
   <td style="text-align:right;"> 14.900876 </td>
   <td style="text-align:right;"> -0.5725128 </td>
   <td style="text-align:right;"> -5.977856 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -8.817908 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 304.47051 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14 </td>
   <td style="text-align:right;"> 10.5955721 </td>
   <td style="text-align:right;"> 1.745013 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -2.670358 </td>
   <td style="text-align:right;"> -6.344036 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.8258375 </td>
   <td style="text-align:right;"> 10.12950 </td>
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

![](05-clustering_celltype_files/figure-html/add_ScType-1.png)<!-- -->

The ScType developer suggests that assignments with a score less than (number of cells in cluster)/4 are low confidence and should be set to unknown:

```r
cluster.ScType.top <- cluster.ScType.top %>%
  mutate(ScType_filtered = ifelse(score >= ncells / 4, ScType, "Unknown"))
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

![](05-clustering_celltype_files/figure-html/adjust_confidence-1.png)<!-- -->


```r
types <- grep("Unknown", unique(experiment.aggregate$subcluster_ScType_filtered), value = TRUE, invert = TRUE)
type.markers <- unlist(lapply(gs.list$gs_positive[types], function(m){
  rownames(markers[which(m %in% rownames(markers)),])
}))
type.markers <- type.markers[!duplicated(type.markers)]
type.markers.df <- data.frame(marker = type.markers,
                              type = gsub('[0-9]', '', names(type.markers)))
type.colors <- viridis::viridis(length(unique(type.markers.df$type)))
names(type.colors) <- unique(type.markers.df$type)
left.annotation <- rowAnnotation("marker" = type.markers.df$type, col = list(marker = type.colors))
mat <- as.matrix(GetAssayData(experiment.aggregate,
                              slot = "data")[type.markers.df$marker,])
Heatmap(mat,
        name = "normalized\ncounts",
        top_annotation = top.annotation,
        left_annotation = left.annotation,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE)
```

![](05-clustering_celltype_files/figure-html/cell_type_heat-1.png)<!-- -->

## Prepare for the next section

#### Save the Seurat object and download the next Rmd file

```r
# set the finalcluster to subcluster
experiment.aggregate$finalcluster <- experiment.aggregate$subcluster
# save object
saveRDS(experiment.aggregate, file="scRNA_workshop-05.rds")
```

#### Download the Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/06-de_enrichment.Rmd", "06-de_enrichment.Rmd")
```

#### Session Information

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
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] ComplexHeatmap_2.14.0 HGNChelper_0.8.1      ggplot2_3.4.0        
## [4] dplyr_1.0.10          tidyr_1.2.1           kableExtra_1.3.4     
## [7] SeuratObject_4.1.3    Seurat_4.3.0         
## 
## loaded via a namespace (and not attached):
##   [1] circlize_0.4.15        systemfonts_1.0.4      plyr_1.8.8            
##   [4] igraph_1.3.5           lazyeval_0.2.2         sp_1.5-1              
##   [7] splines_4.2.2          listenv_0.8.0          scattermore_0.8       
##  [10] digest_0.6.30          foreach_1.5.2          htmltools_0.5.3       
##  [13] viridis_0.6.2          fansi_1.0.3            magrittr_2.0.3        
##  [16] tensor_1.5             cluster_2.1.4          doParallel_1.0.17     
##  [19] ROCR_1.0-11            openxlsx_4.2.5.1       limma_3.54.0          
##  [22] globals_0.16.2         matrixStats_0.63.0     svglite_2.1.0         
##  [25] spatstat.sparse_3.0-0  colorspace_2.0-3       rvest_1.0.3           
##  [28] ggrepel_0.9.2          xfun_0.35              crayon_1.5.2          
##  [31] jsonlite_1.8.4         progressr_0.11.0       spatstat.data_3.0-0   
##  [34] ape_5.6-2              survival_3.4-0         zoo_1.8-11            
##  [37] iterators_1.0.14       glue_1.6.2             polyclip_1.10-4       
##  [40] gtable_0.3.1           webshot_0.5.4          leiden_0.4.3          
##  [43] GetoptLong_1.0.5       shape_1.4.6            future.apply_1.10.0   
##  [46] BiocGenerics_0.44.0    abind_1.4-5            scales_1.2.1          
##  [49] DBI_1.1.3              spatstat.random_3.0-1  miniUI_0.1.1.1        
##  [52] Rcpp_1.0.9             viridisLite_0.4.1      xtable_1.8-4          
##  [55] clue_0.3-63            reticulate_1.28        stats4_4.2.2          
##  [58] htmlwidgets_1.5.4      httr_1.4.4             RColorBrewer_1.1-3    
##  [61] ellipsis_0.3.2         ica_1.0-3              farver_2.1.1          
##  [64] pkgconfig_2.0.3        sass_0.4.4             uwot_0.1.14           
##  [67] deldir_1.0-6           utf8_1.2.2             labeling_0.4.2        
##  [70] tidyselect_1.2.0       rlang_1.0.6            reshape2_1.4.4        
##  [73] later_1.3.0            munsell_0.5.0          tools_4.2.2           
##  [76] cachem_1.0.6           cli_3.4.1              generics_0.1.3        
##  [79] ggridges_0.5.4         evaluate_0.18          stringr_1.4.1         
##  [82] fastmap_1.1.0          yaml_2.3.6             goftest_1.2-3         
##  [85] knitr_1.41             fitdistrplus_1.1-8     zip_2.2.2             
##  [88] purrr_0.3.5            RANN_2.6.1             pbapply_1.6-0         
##  [91] future_1.29.0          nlme_3.1-160           mime_0.12             
##  [94] ggrastr_1.0.1          xml2_1.3.3             compiler_4.2.2        
##  [97] rstudioapi_0.14        beeswarm_0.4.0         plotly_4.10.1         
## [100] png_0.1-8              spatstat.utils_3.0-1   tibble_3.1.8          
## [103] bslib_0.4.1            stringi_1.7.8          highr_0.9             
## [106] lattice_0.20-45        Matrix_1.5-3           vctrs_0.5.1           
## [109] pillar_1.8.1           lifecycle_1.0.3        spatstat.geom_3.0-3   
## [112] lmtest_0.9-40          jquerylib_0.1.4        GlobalOptions_0.1.2   
## [115] RcppAnnoy_0.0.20       data.table_1.14.6      cowplot_1.1.1         
## [118] irlba_2.3.5.1          httpuv_1.6.6           patchwork_1.1.2       
## [121] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20    
## [124] gridExtra_2.3          vipor_0.4.5            IRanges_2.32.0        
## [127] parallelly_1.32.1      codetools_0.2-18       MASS_7.3-58.1         
## [130] assertthat_0.2.1       rjson_0.2.21           withr_2.5.0           
## [133] sctransform_0.3.5      S4Vectors_0.36.0       parallel_4.2.2        
## [136] rmarkdown_2.18         Cairo_1.6-0            Rtsne_0.16            
## [139] spatstat.explore_3.0-5 shiny_1.7.3            ggbeeswarm_0.6.0
```
