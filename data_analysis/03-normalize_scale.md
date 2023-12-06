---
title: "Introduction to Single Cell RNA-Seq Part 3: Normalize and scale"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

# Introduction to Single Cell RNA-Seq Part 3: Normalize and scale


## Set up workspace
First, we need to load the required libraries.

```r
library(Seurat)
library(kableExtra)
experiment.aggregate <- readRDS("scRNA_workshop-02.rds")
experiment.aggregate
```

```
## An object of class Seurat 
## 11292 features across 6313 samples within 1 assay 
## Active assay: RNA (11292 features, 0 variable features)
```

```r
set.seed(12345)
```

## Normalize the data
After filtering, the next step is to normalize the data. We employ a global-scaling normalization method, LogNormalize, that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and then log-transforms the data.

```r
?NormalizeData
```


```r
experiment.aggregate <- NormalizeData(
  object = experiment.aggregate,
  normalization.method = "LogNormalize",
  scale.factor = 10000)
```

## Cell cycle assignment
Cell cycle phase can be a significant source of variation in single cell and single nucleus experiments. There are a number of automated cell cycle stage detection methods available for single cell data. For this workshop, we will be using the built-in Seurat cell cycle function, `CellCycleScoring`. This tool compares gene expression in each cell to a list of cell cycle marker genes and scores each barcode based on marker expression. The phase with the highest score is selected for each barcode. Seurat includes a list of cell cycle genes in human single cell data.

```r
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
```

For other species, a user-provided gene list may be substituted, or the orthologs of the human gene list used instead.

**Do not run the code below for human experiments!**

```r
# mouse code DO NOT RUN for human data
library(biomaRt)
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     mirror = "uswest")
  mouse = useEnsembl("ensembl",
                     dataset = "mmusculus_gene_ensembl",
                     mirror = "uswest")
  genes = getLDS(attributes = c("hgnc_symbol"),
                 filters = "hgnc_symbol",
                 values = x ,
                 mart = human,
                 attributesL = c("mgi_symbol"),
                 martL = mouse,
                 uniqueRows=T)
  humanx = unique(genes[, 2])
  print(head(humanx)) # print first 6 genes found to the screen
  return(humanx)
}
# convert lists to mouse orthologs
s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
```

Once an appropriate gene list has been identified, the `CellCycleScoring` function can be run.

```r
experiment.aggregate <- CellCycleScoring(experiment.aggregate,
                                         s.features = s.genes,
                                         g2m.features = g2m.genes,
                                         set.ident = TRUE)
table(experiment.aggregate@meta.data$Phase) %>%
  kable(caption = "Number of Cells in each Cell Cycle Stage",
        col.names = c("Stage", "Count"),
        align = "c") %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Number of Cells in each Cell Cycle Stage</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> Stage </th>
   <th style="text-align:center;"> Count </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> G1 </td>
   <td style="text-align:center;"> 3777 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> G2M </td>
   <td style="text-align:center;"> 1131 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> S </td>
   <td style="text-align:center;"> 1405 </td>
  </tr>
</tbody>
</table>

Because the "set.ident" argument was set to TRUE (this is also the default behavior), the active identity of the Seurat object was changed to the phase. To return the active identity to the sample identity, use the `Idents` function.

```r
table(Idents(experiment.aggregate))
```

```
## 
##    S  G2M   G1 
## 1405 1131 3777
```

```r
Idents(experiment.aggregate) <- "orig.ident"
table(Idents(experiment.aggregate))
```

```
## 
## A001-C-007 A001-C-104 B001-A-301 
##       1023       1859       3431
```

## Identify variable genes
The function FindVariableFeatures identifies the most highly variable genes (default 2000 genes) by fitting a line to the relationship of log(variance) and log(mean) using loess smoothing, uses this information to standardize the data, then calculates the variance of the standardized data.  This helps avoid selecting genes that only appear variable due to their expression level.

```r
?FindVariableFeatures
```


```r
experiment.aggregate <- FindVariableFeatures(
  object = experiment.aggregate,
  selection.method = "vst")
length(VariableFeatures(experiment.aggregate))
```

```
## [1] 2000
```

```r
top10 <- head(VariableFeatures(experiment.aggregate), 10)
top10
```

```
##  [1] "BEST4"      "CLCA4"      "SMOC2"      "NRG1"       "TPH1"      
##  [6] "LRMP"       "TRPM3"      "PTPRR"      "AC007493.1" "CACNA1A"
```

```r
var.feat.plot <- VariableFeaturePlot(experiment.aggregate)
var.feat.plot <- LabelPoints(plot = var.feat.plot, points = top10, repel = TRUE)
var.feat.plot
```

![](03-normalize_scale_files/figure-html/find_variable_genes-1.png)<!-- -->

**How do the results change if you use selection.method = "dispersion" or selection.method = "mean.var.plot"?**

FindVariableFeatures isn't the only way to set the "variable features" of a Seurat object. Another reasonable approach is to select a set of "minimally expressed" genes.

```r
min.value <- 2
min.cells <- 10

num.cells <- Matrix::rowSums(GetAssayData(experiment.aggregate, slot = "count") > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])
length(genes.use)
```

```
## [1] 7012
```

```r
VariableFeatures(experiment.aggregate) <- genes.use
```

## Scale the data
The `ScaleData` function scales and centers genes in the dataset. If variables are provided with the "vars.to.regress" argument, they are individually regressed against each gene, and the resulting residuals are then scaled and centered unless otherwise specified. We regress out cell cycle results S.Score and G2M.Score, mitochondrial RNA level (percent_MT), and the number of features (nFeature_RNA) as a proxy for sequencing depth.

```r
experiment.aggregate <- ScaleData(experiment.aggregate,
                                  vars.to.regress = c("S.Score", "G2M.Score", "percent_MT", "nFeature_RNA"))
```

## Prepare for the next section

#### Save object

```r
saveRDS(experiment.aggregate, file = "scRNA_workshop-03.rds")
```

#### Download Rmd

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/04-dimensionality_reduction.Rmd", "04-dimensionality_reduction.Rmd")
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] kableExtra_1.3.4   SeuratObject_4.1.3 Seurat_4.3.0      
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6          
##   [4] ellipsis_0.3.2         ggridges_0.5.4         rstudioapi_0.14       
##   [7] spatstat.data_3.0-0    farver_2.1.1           leiden_0.4.3          
##  [10] listenv_0.8.0          ggrepel_0.9.2          fansi_1.0.3           
##  [13] xml2_1.3.3             codetools_0.2-18       splines_4.2.2         
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
##  [46] vctrs_0.5.1            svglite_2.1.0          nlme_3.1-160          
##  [49] spatstat.explore_3.0-5 progressr_0.11.0       lmtest_0.9-40         
##  [52] spatstat.random_3.0-1  xfun_0.35              stringr_1.4.1         
##  [55] globals_0.16.2         rvest_1.0.3            mime_0.12             
##  [58] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5.1         
##  [61] goftest_1.2-3          future_1.29.0          MASS_7.3-58.1         
##  [64] zoo_1.8-11             scales_1.2.1           promises_1.2.0.1      
##  [67] spatstat.utils_3.0-1   parallel_4.2.2         RColorBrewer_1.1-3    
##  [70] yaml_2.3.6             reticulate_1.28        pbapply_1.6-0         
##  [73] gridExtra_2.3          ggplot2_3.4.0          sass_0.4.4            
##  [76] stringi_1.7.8          highr_0.9              systemfonts_1.0.4     
##  [79] rlang_1.0.6            pkgconfig_2.0.3        matrixStats_0.63.0    
##  [82] evaluate_0.18          lattice_0.20-45        tensor_1.5            
##  [85] ROCR_1.0-11            purrr_0.3.5            labeling_0.4.2        
##  [88] patchwork_1.1.2        htmlwidgets_1.5.4      cowplot_1.1.1         
##  [91] tidyselect_1.2.0       parallelly_1.32.1      RcppAnnoy_0.0.20      
##  [94] plyr_1.8.8             magrittr_2.0.3         R6_2.5.1              
##  [97] generics_0.1.3         DBI_1.1.3              withr_2.5.0           
## [100] pillar_1.8.1           fitdistrplus_1.1-8     survival_3.4-0        
## [103] abind_1.4-5            sp_1.5-1               tibble_3.1.8          
## [106] future.apply_1.10.0    KernSmooth_2.23-20     utf8_1.2.2            
## [109] spatstat.geom_3.0-3    plotly_4.10.1          rmarkdown_2.18        
## [112] grid_4.2.2             data.table_1.14.6      webshot_0.5.4         
## [115] digest_0.6.30          xtable_1.8-4           tidyr_1.2.1           
## [118] httpuv_1.6.6           munsell_0.5.0          viridisLite_0.4.1     
## [121] bslib_0.4.1
```
