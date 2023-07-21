---
title: "Introduction to Single Cell RNA-Seq Part 2: QA and filtering"
author: "UCD Bioinformatics Core"
date: "2023-07-21"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 2: QA and filtering


## Set up workspace
First, we need to load the required libraries.

```r
library(Seurat)
library(ggplot2)
library(tidyr)
library(kableExtra)
```

If you are continuing directly from part 1, the experiment.aggregate object is likely already in your workspace. In case you cleared your workspace at the end of the previous section, or are working on this project at a later date after re-starting R, you can use the `readRDS` function to read your saved Seurat object from part 1.

```r
experiment.aggregate <- readRDS("scRNA_workshop-01.rds")
experiment.aggregate
```

<div class='r_output'> An object of class Seurat 
 36601 features across 9452 samples within 1 assay 
 Active assay: RNA (36601 features, 0 variable features)
</div>
The seed is used to initialize pseudo-random functions. Some of the functions we will be using have pseudo-random elements. Setting a common seed ensures that all of us will get the same results, and that the results will remain stable when re-run.

```r
set.seed(12345)
```

## Mitochondrial gene expression
Filtering on the expression of genes from the mitochondrial genome is not appropriate in all cell types. However, in many tissues, low-quality / dying cells may exhibit extensive mitochondrial contamination. Even when not filtering on mitochondrial expression, the data can be interesting or informative.

The `PercentageFeatureSet` function calculates the proportion of counts originating from a set of features. Genes in the human mitochondrial genome begin with 'MT', while those in the mouse mitochondrial genome begin with 'mt'. These naming conventions make calculating percent mitochondrial very straightforward.

```r
experiment.aggregate$percent_MT <- PercentageFeatureSet(experiment.aggregate, pattern = "^MT-")
summary(experiment.aggregate$percent_MT)
```

<div class='r_output'>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.0000  0.2802  0.5601  0.8283  1.0834 14.2037
</div>
In this workshop, we are using the filtered feature barcode matrix. While this lowers the likelihood of encountering barcodes that are not cell-associated within our expression matrix, it is still good practice to perform quality assurance on the experiment.

## Display metadata by quantile
Using a few nested functions, we can produce prettier, more detailed, versions of the simple exploratory summary statistics we generated for the available metadata in the last section. In the code below, 10% quantile tables are produced for each metadata value, separated by sample identity.

```r
kable(do.call("cbind", tapply(experiment.aggregate$nFeature_RNA,
                              Idents(experiment.aggregate),
                              quantile, probs = seq(0,1,0.1))),
      caption = "10% Quantiles of Genes/Cell by Sample") %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>10% Quantiles of Genes/Cell by Sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> A001-C-007 </th>
   <th style="text-align:right;"> A001-C-104 </th>
   <th style="text-align:right;"> B001-A-301 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0% </td>
   <td style="text-align:right;"> 404 </td>
   <td style="text-align:right;"> 397.0 </td>
   <td style="text-align:right;"> 416.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10% </td>
   <td style="text-align:right;"> 495 </td>
   <td style="text-align:right;"> 505.0 </td>
   <td style="text-align:right;"> 577.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20% </td>
   <td style="text-align:right;"> 574 </td>
   <td style="text-align:right;"> 581.0 </td>
   <td style="text-align:right;"> 733.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 30% </td>
   <td style="text-align:right;"> 677 </td>
   <td style="text-align:right;"> 676.3 </td>
   <td style="text-align:right;"> 907.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 40% </td>
   <td style="text-align:right;"> 783 </td>
   <td style="text-align:right;"> 791.0 </td>
   <td style="text-align:right;"> 1105.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 50% </td>
   <td style="text-align:right;"> 927 </td>
   <td style="text-align:right;"> 959.0 </td>
   <td style="text-align:right;"> 1331.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 60% </td>
   <td style="text-align:right;"> 1129 </td>
   <td style="text-align:right;"> 1127.0 </td>
   <td style="text-align:right;"> 1564.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70% </td>
   <td style="text-align:right;"> 1382 </td>
   <td style="text-align:right;"> 1358.0 </td>
   <td style="text-align:right;"> 1810.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80% </td>
   <td style="text-align:right;"> 1797 </td>
   <td style="text-align:right;"> 1718.8 </td>
   <td style="text-align:right;"> 2119.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90% </td>
   <td style="text-align:right;"> 2592 </td>
   <td style="text-align:right;"> 2344.7 </td>
   <td style="text-align:right;"> 2567.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 100% </td>
   <td style="text-align:right;"> 12063 </td>
   <td style="text-align:right;"> 12064.0 </td>
   <td style="text-align:right;"> 8812.0 </td>
  </tr>
</tbody>
</table>

```r
kable(do.call("cbind", tapply(experiment.aggregate$nCount_RNA,
                              Idents(experiment.aggregate),
                              quantile, probs = seq(0,1,0.1))),
      caption = "10% Quantiles of UMI/Cell by Sample") %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>10% Quantiles of UMI/Cell by Sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> A001-C-007 </th>
   <th style="text-align:right;"> A001-C-104 </th>
   <th style="text-align:right;"> B001-A-301 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0% </td>
   <td style="text-align:right;"> 500.0 </td>
   <td style="text-align:right;"> 500.0 </td>
   <td style="text-align:right;"> 500.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10% </td>
   <td style="text-align:right;"> 576.5 </td>
   <td style="text-align:right;"> 588.0 </td>
   <td style="text-align:right;"> 698.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20% </td>
   <td style="text-align:right;"> 686.0 </td>
   <td style="text-align:right;"> 692.2 </td>
   <td style="text-align:right;"> 909.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 30% </td>
   <td style="text-align:right;"> 826.5 </td>
   <td style="text-align:right;"> 821.0 </td>
   <td style="text-align:right;"> 1179.9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 40% </td>
   <td style="text-align:right;"> 981.0 </td>
   <td style="text-align:right;"> 995.8 </td>
   <td style="text-align:right;"> 1514.6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 50% </td>
   <td style="text-align:right;"> 1203.5 </td>
   <td style="text-align:right;"> 1231.0 </td>
   <td style="text-align:right;"> 1913.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 60% </td>
   <td style="text-align:right;"> 1493.0 </td>
   <td style="text-align:right;"> 1514.6 </td>
   <td style="text-align:right;"> 2399.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70% </td>
   <td style="text-align:right;"> 1920.0 </td>
   <td style="text-align:right;"> 1909.2 </td>
   <td style="text-align:right;"> 2960.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80% </td>
   <td style="text-align:right;"> 2690.0 </td>
   <td style="text-align:right;"> 2606.6 </td>
   <td style="text-align:right;"> 3704.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90% </td>
   <td style="text-align:right;"> 4419.0 </td>
   <td style="text-align:right;"> 4052.5 </td>
   <td style="text-align:right;"> 5087.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 100% </td>
   <td style="text-align:right;"> 150805.0 </td>
   <td style="text-align:right;"> 149096.0 </td>
   <td style="text-align:right;"> 89743.0 </td>
  </tr>
</tbody>
</table>

```r
kable(round(do.call("cbind", tapply(experiment.aggregate$percent_MT,
                                    Idents(experiment.aggregate),
                                    quantile, probs = seq(0,1,0.1))),
            digits = 3),
      caption = "10% Quantiles of Percent Mitochondrial by Sample") %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>10% Quantiles of Percent Mitochondrial by Sample</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> A001-C-007 </th>
   <th style="text-align:right;"> A001-C-104 </th>
   <th style="text-align:right;"> B001-A-301 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0% </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10% </td>
   <td style="text-align:right;"> 0.266 </td>
   <td style="text-align:right;"> 0.350 </td>
   <td style="text-align:right;"> 0.118 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20% </td>
   <td style="text-align:right;"> 0.368 </td>
   <td style="text-align:right;"> 0.510 </td>
   <td style="text-align:right;"> 0.164 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 30% </td>
   <td style="text-align:right;"> 0.473 </td>
   <td style="text-align:right;"> 0.677 </td>
   <td style="text-align:right;"> 0.211 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 40% </td>
   <td style="text-align:right;"> 0.581 </td>
   <td style="text-align:right;"> 0.851 </td>
   <td style="text-align:right;"> 0.263 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 50% </td>
   <td style="text-align:right;"> 0.705 </td>
   <td style="text-align:right;"> 1.043 </td>
   <td style="text-align:right;"> 0.323 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 60% </td>
   <td style="text-align:right;"> 0.864 </td>
   <td style="text-align:right;"> 1.281 </td>
   <td style="text-align:right;"> 0.404 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70% </td>
   <td style="text-align:right;"> 1.061 </td>
   <td style="text-align:right;"> 1.537 </td>
   <td style="text-align:right;"> 0.512 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80% </td>
   <td style="text-align:right;"> 1.330 </td>
   <td style="text-align:right;"> 1.912 </td>
   <td style="text-align:right;"> 0.681 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90% </td>
   <td style="text-align:right;"> 1.717 </td>
   <td style="text-align:right;"> 2.735 </td>
   <td style="text-align:right;"> 0.947 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 100% </td>
   <td style="text-align:right;"> 14.204 </td>
   <td style="text-align:right;"> 13.172 </td>
   <td style="text-align:right;"> 3.481 </td>
  </tr>
</tbody>
</table>

## Visualize distribution of metadata values
Seurat has a number of convenient built-in functions for visualizing metadata. These functions produce ggplot objects, which can easily be modified using ggplot2. Of course, all of these visualizations can be reproduced with custom code as well, and we will include some examples of both modifying Seurat plots and generating plots from scratch as the analysis continues.

### Violin plots
The `VlnPlot` function produces a composite plot with one panel for each element of the "features" vector. The data are grouped by the provided identity; by default, this is the active identity of the object, which can be accessed using the `Idents()` function, or in the "active.ident" slot.

```r
VlnPlot(experiment.aggregate,
        features = c("nFeature_RNA", "nCount_RNA","percent_MT"),
        ncol = 1,
        pt.size = 0.3)
```

![](02-filtering_files/figure-html/violins-1.png)<!-- -->

#### Modifying Seurat plots
Modifying the ggplot objects produced by a Seurat plotting function works best on individual panels. Therefore, to recreate the function above with modifications, we can use `lapply` to create a list of plots. In some cases it may be more appropriate to create the plots individually so that different modifications can be applied to each plot.

```r
lapply(c("nFeature_RNA", "nCount_RNA","percent_MT"), function(feature){
  VlnPlot(experiment.aggregate,
          features = feature,
          pt.size = 0.01) +
    scale_fill_viridis_d(option = "mako") # default colors are not colorblind-friendly
})
```

<div class='r_output'> [[1]]
</div>
![](02-filtering_files/figure-html/violins_list-1.png)<!-- -->

<div class='r_output'> 
 [[2]]
</div>
![](02-filtering_files/figure-html/violins_list-2.png)<!-- -->

<div class='r_output'> 
 [[3]]
</div>
![](02-filtering_files/figure-html/violins_list-3.png)<!-- -->

```r
VlnPlot(experiment.aggregate, features = "nCount_RNA", pt.size = 0.01) + 
  scale_y_continuous(trans = "log10") +
  scale_fill_viridis_d(option = "mako") +
  ggtitle("log10(nCount_RNA)")
```

![](02-filtering_files/figure-html/violins_list-4.png)<!-- -->

These can later be stitched together with another library, like patchwork, or cowplot.

### Ridge plots
Ridge plots are very similar in appearance to violin plots turned on their sides. In some cases it may be more appropriate to create the plots individually so that appropriate transformations can be applied to each plot.

```r
RidgePlot(experiment.aggregate, features="nFeature_RNA") +
  scale_fill_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/ridgeplot-1.png)<!-- -->

```r
RidgePlot(experiment.aggregate, features="nCount_RNA") +
  scale_x_continuous(trans = "log10") + # "un-squish" the distribution
  scale_fill_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/ridgeplot-2.png)<!-- -->

```r
RidgePlot(experiment.aggregate, features="percent_MT") +
  scale_fill_viridis_d(option = "mako") +
  coord_cartesian(xlim = c(0, 10)) # zoom in on the lower end of the distribution
```

![](02-filtering_files/figure-html/ridgeplot-3.png)<!-- -->

### Custom plots
The Seurat built-in functions are useful and easy to interact with, but sometimes you may wish to visualize something for which a plotting function does not already exist. For example, we might want to see how many cells are expressing each gene over some UMI threshold.

The code below produces a ranked plot similar to the barcode inflection plots from the last section. On the x-axis are the genes arranged from most ubiquitously expressed to rarest. In a single cell dataset, many genes are expessed in a relatively small number of cells, or not at all. The y-axis displays the number of cells in which each gene is expressed.

**Note: This function is SLOW. You may want to skip this code block or run it while you take a break for a few minutes.**

```r
# retrieve count data
counts <- GetAssayData(experiment.aggregate)
# order genes from most to least ubiquitous
ranked.genes <- names(sort(Matrix::rowSums(counts >= 3), decreasing = TRUE))
# drop genes not expressed in any cell
ranked.genes <- ranked.genes[ranked.genes %in% names(which(Matrix::rowSums(counts >= 3) >= 1))]
# get number of cells in which gene is expressed for each sample
cell.counts <- sapply(ranked.genes, function(gene){
  tapply(counts[gene,], experiment.aggregate$orig.ident, function(x){
    sum(x >= 3)
  })
})
cell.counts <- as.data.frame(t(cell.counts))
cell.counts$gene <- rownames(cell.counts)
cell.counts$rank <- seq(1:dim(cell.counts)[1])
cell.counts %>%
  pivot_longer(cols = 1:3, names_to = "sample", values_to = "count") %>%
  ggplot(mapping = aes(x = rank, y = count, color = sample)) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_smooth() +
  scale_color_viridis_d(option = "mako") +
  theme_classic() +
  theme(legend.title = element_blank())
```

![](02-filtering_files/figure-html/gene_range-1.png)<!-- -->

### Scatter plots
Scatter plots allow us to visualize the relationships between the metadata variables.

```r
# mitochondrial vs UMI
FeatureScatter(experiment.aggregate,
               feature1 = "nCount_RNA",
               feature2 = "percent_MT",
               shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/relationships-1.png)<!-- -->

```r
# mitochondrial vs genes
FeatureScatter(experiment.aggregate,
               feature1 = "nFeature_RNA",
               feature2 = "percent_MT",
               shuffle = TRUE) +
  scale_color_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/relationships-2.png)<!-- -->

```r
# genes vs UMI
FeatureScatter(experiment.aggregate,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA",
               shuffle = TRUE)  +
  scale_color_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/relationships-3.png)<!-- -->

## Cell filtering
The goal of cell filtering is to remove cells with anomolous expression profiles, typically low UMI cells, which may correspond to low-quality cells or background barcodes that made it through the Cell Ranger filtration algorithm. It may also be appropriate to remove outlier cells with extremely high UMI counts.
In this case, the proposed cut-offs on the high end of the distributions are quite conservative, in part to reduce the size of the object and speed up analysis during the workshop.

The plots below display proposed filtering cut-offs.

```r
FeatureScatter(experiment.aggregate,
               feature1 = "nCount_RNA",
               feature2 = "percent_MT",
               shuffle = TRUE) +
  geom_vline(xintercept = c(1000, 25000)) +
  geom_hline(yintercept = 5) +
  scale_color_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/thresholds-1.png)<!-- -->

```r
FeatureScatter(experiment.aggregate,
               feature1 = "nFeature_RNA",
               feature2 = "percent_MT",
               shuffle = TRUE) +
  geom_vline(xintercept = c(500, 7500)) +
  geom_hline(yintercept = 5) +
  scale_color_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/thresholds-2.png)<!-- -->

```r
FeatureScatter(experiment.aggregate,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA",
               pt.size = 0.5,
               shuffle = TRUE)  +
  geom_vline(xintercept = c(1000, 25000)) +
  geom_hline(yintercept = c(500, 7500)) +
  scale_color_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/thresholds-3.png)<!-- -->

These filters can be put in place with the `subset` function.

```r
table(experiment.aggregate$orig.ident)
```

<div class='r_output'> 
 A001-C-007 A001-C-104 B001-A-301 
       1796       3142       4514
</div>
```r
# mitochondrial filter
experiment.filter <- subset(experiment.aggregate, percent_MT <= 5)
# UMI filter
experiment.filter <- subset(experiment.filter, nCount_RNA >= 1000 & nCount_RNA <= 25000)
# gene filter
experiment.filter <- subset(experiment.filter, nFeature_RNA >= 500 & nFeature_RNA <= 7500)
# filtering results
experiment.filter
```

<div class='r_output'> An object of class Seurat 
 36601 features across 6312 samples within 1 assay 
 Active assay: RNA (36601 features, 0 variable features)
</div>
```r
table(experiment.filter$orig.ident)
```

<div class='r_output'> 
 A001-C-007 A001-C-104 B001-A-301 
       1023       1858       3431
</div>
```r
# ridge plots
RidgePlot(experiment.filter, features="nFeature_RNA") +
  scale_fill_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/cell_filtering-1.png)<!-- -->

```r
RidgePlot(experiment.filter, features="nCount_RNA") +
  scale_x_continuous(trans = "log10") + 
  scale_fill_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/cell_filtering-2.png)<!-- -->

```r
RidgePlot(experiment.filter, features="percent_MT") +
  scale_fill_viridis_d(option = "mako")
```

![](02-filtering_files/figure-html/cell_filtering-3.png)<!-- -->

```r
# use filtered results from now on
experiment.aggregate <- experiment.filter
```

**Play with the filtering parameters, and see how the results change. Is there a set of parameters you feel is more appropriate? Why?**

## Feature filtering
When creating the base Seurat object, we had the opportunity filter out some genes using the "min.cells" argument. At the time, we set that to 0. Since we didn't filter our features then, we can apply a filter at this point. If we had filtered when the object was created, this would be an opportunity to be more aggressive. The custom code below provides a function that filters genes requiring a min.umi in at least min.cells, or takes a user-provided list of genes.

```r
# define function
FilterGenes <- function(object, min.umi = NA, min.cells = NA, genes = NULL) {
  genes.use = NA
  if (!is.null(genes)) {
    genes.use = intersect(rownames(object), genes)
    } else if (min.cells & min.umi) {
      num.cells = Matrix::rowSums(GetAssayData(object) >= min.umi)
      genes.use = names(num.cells[which(num.cells >= min.cells)])
    }
  object = object[genes.use,]
  object = LogSeuratCommand(object = object)
  return(object)
}
# apply filter
experiment.filter <- FilterGenes(object = experiment.aggregate, min.umi = 2, min.cells = 10)
# filtering results
experiment.filter
```

<div class='r_output'> An object of class Seurat 
 11292 features across 6312 samples within 1 assay 
 Active assay: RNA (11292 features, 0 variable features)
</div>
```r
experiment.aggregate <- experiment.filter
```

## Prepare for the next section

#### Save object

```r
saveRDS(experiment.aggregate, file="scRNA_workshop-02.rds")
```

#### Download Rmd

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/03-normalize_scale.Rmd", "03-normalize_scale.Rmd")
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
 [1] kableExtra_1.3.4   tidyr_1.3.0        ggplot2_3.4.2      SeuratObject_4.1.3
 [5] Seurat_4.3.0.1    
 
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
  [43] mgcv_1.9-0             httr_1.4.6             polyclip_1.10-4       
  [46] abind_1.4-5            compiler_4.3.1         withr_2.5.0           
  [49] highr_0.10             MASS_7.3-60            tools_4.3.1           
  [52] lmtest_0.9-40          httpuv_1.6.11          future.apply_1.11.0   
  [55] goftest_1.2-3          glue_1.6.2             nlme_3.1-162          
  [58] promises_1.2.0.1       grid_4.3.1             Rtsne_0.16            
  [61] cluster_2.1.4          reshape2_1.4.4         generics_0.1.3        
  [64] gtable_0.3.3           spatstat.data_3.0-1    data.table_1.14.8     
  [67] sp_2.0-0               xml2_1.3.5             utf8_1.2.3            
  [70] spatstat.geom_3.2-4    RcppAnnoy_0.0.21       ggrepel_0.9.3         
  [73] RANN_2.6.1             pillar_1.9.0           stringr_1.5.0         
  [76] later_1.3.1            splines_4.3.1          dplyr_1.1.2           
  [79] lattice_0.21-8         survival_3.5-5         deldir_1.0-9          
  [82] tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2         
  [85] knitr_1.43             gridExtra_2.3          svglite_2.1.1         
  [88] scattermore_1.2        xfun_0.39              matrixStats_1.0.0     
  [91] stringi_1.7.12         lazyeval_0.2.2         yaml_2.3.7            
  [94] evaluate_0.21          codetools_0.2-19       tibble_3.2.1          
  [97] cli_3.6.1              uwot_0.1.16            xtable_1.8-4          
 [100] reticulate_1.30        systemfonts_1.0.4      munsell_0.5.0         
 [103] jquerylib_0.1.4        Rcpp_1.0.11            globals_0.16.2        
 [106] spatstat.random_3.1-5  png_0.1-8              parallel_4.3.1        
 [109] ellipsis_0.3.2         listenv_0.9.0          viridisLite_0.4.2     
 [112] scales_1.2.1           ggridges_0.5.4         leiden_0.4.3          
 [115] purrr_1.0.1            rlang_1.1.1            cowplot_1.1.1         
 [118] rvest_1.0.3
</div>