---
title: "Introduction to Single Cell RNA-Seq Part 7: Doublet detection"
author: "UCD Bioinformatics Core"
date: "2023-12-07"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---


# Introduction to Single Cell RNA-Seq Part 7: Doublet detection
Doublets are cells that appear to be, but are not, real cells. There are two major types of doublets: heterotypic and homotypic. *Heterotypic doublets* are formed by cells with distinct transcriptional profiles. *Homotypic doublets* are formed by cells with similar transcriptional profiles. Heterotypic doublets are relatively easier to detect compared with homotypic doublets.

Depending on the protocols used to barcode single cells/nuclei, doublet rates vary significantly and it can reach as high as 40%. Experimental strategies have been developed to reduce the doublet rate, such as [cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1), and [MULTI-Seq](https://www.nature.com/articles/s41592-019-0433-8). However, these techniques require extra steps in sample preparation which leads to extra costs and time, and they do not guarantee to remove all doublets.

Naturally, removing doublets _in silico_ is very appealing and there have been many tools/methods developed to achieve this: [DoubletFinder](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30073-0), DoubletDetection(https://github.com/JonathanShor/DoubletDetection), [DoubletDecon](https://www.sciencedirect.com/science/article/pii/S2211124719312860), [demuxlet](https://www.nature.com/articles/nbt.4042), among others.

<p align = "center">
<img src="figures/doublets.jpg" alt="micribial" width="85%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
Xi, etc., Cell Systems, 2021, https://www.sciencedirect.com/science/article/pii/S2405471220304592
</p>

## Set up workspace


```r
library(Seurat)
library(DoubletFinder)
library(ggplot2)
set.seed(12345)
```

## Prepare data for DoubletFinder
[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) takes fully pre-processed data from Seurat (NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunUMAP) as input and the process should be done for each sample individually. The input data should be processed to remove low-quality cell clusters first.

```r
experiment.aggregate <- readRDS("scRNA_workshop-02.rds") # filtered object
experiment.split <- SplitObject(experiment.aggregate, split.by = "orig.ident")
rm(experiment.aggregate)
experiment.split <- lapply(experiment.split, function(sce){
  sce = NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce = CellCycleScoring(sce, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score", "percent_MT", "nFeature_RNA"))
  RunPCA(sce, npcs = 100)
})
```

## Parameter selection
In addition to the Seurat object, DoubletFinder takes a number of arguments. These are

* PCs: a vector of statistically significant PCs to use
* pN: the number of artificially generated doublets (default = 0.25)
* pK: PC neighborhood size used to compute network
* nExp: threshold used to make doublet/singlet call

We will use PCs 1-50 (based on the dimensionality reduction section) and the default value of 0.25 for pN.

### pK
The optimal value of pK varies between samples, and is impacted by the number of cell states and magnitude of transcriptional heterogeneity present in the data. The mean-variance normalized bimodality coefficient (BCmvn) is used as a measure of pK optimization. In experiments with known doublet frequencies, BCmvn is maximized by values of pK that produce the most accurate doublet-calling. In the code below, the pK value corresponding to the maxiumum BCmvn is selected for each sample.

```r
pK <- lapply(experiment.split, function(sce){
  sweep.res = paramSweep_v3(sce, PCs = 1:50, sct = FALSE)
  sweep.stats = summarizeSweep(sweep.res, GT = FALSE)
  BCmvn = find.pK(sweep.stats)
  as.numeric(BCmvn$pK[which(BCmvn$BCmetric == max(BCmvn$BCmetric))])
})
```

![](07-doublet_detection_files/figure-html/pK_select-1.png)<!-- -->![](07-doublet_detection_files/figure-html/pK_select-2.png)<!-- -->![](07-doublet_detection_files/figure-html/pK_select-3.png)<!-- -->

```r
pK
```

```
## $`A001-C-007`
## [1] 30
## 
## $`A001-C-104`
## [1] 30
## 
## $`B001-A-301`
## [1] 10
```

### nExp
In single cell data generated using microfluidics, the frequency of all multiplets (both homotypic and heterotypic) can be modeled with a Poisson distribution. The probability of capturing multiple cells in a droplet is a function of the loading density of the microfluidic device. However, because homotypic doublets are far more difficult to detect informatically, the Poisson distribution overestimates the number of *detectable* doublets.

Here we will use the 10x expected doublet rate of 8% for chips that were loaded with 16,500 cells (targeting 10,000 cells for recovery).

```r
nExp <- lapply(experiment.split, function(sce){
  round(0.08*dim(sce)[2])
})
nExp
```

```
## $`A001-C-007`
## [1] 82
## 
## $`A001-C-104`
## [1] 149
## 
## $`B001-A-301`
## [1] 274
```
#### Homotypic doublet adjustment
If cell type annotations are available, the following code can be used to adjust the estimated detectable doublet frequency downwards by adjusting for the expected frequency of homotypic doublets.

Let's read in the cell type assignments from the previous part and attach them to the split object.

```r
experiment.aggregate <- readRDS("scRNA_workshop-05.rds")
sctype.split <- lapply(SplitObject(experiment.aggregate, split.by = "orig.ident"),
                       function(sce){
  sce$subcluster_ScType
})
experiment.split <- lapply(seq_along(experiment.split), function(i){
  AddMetaData(experiment.split[[i]],
              metadata = sctype.split[[i]],
              col.name = "subcluster_ScType_filtered")
})
homotypic.prop <- lapply(experiment.split, function(sce){
  modelHomotypic(annotations = sce@meta.data$subcluster_ScType_filtered)
})
nExp.adj <- lapply(seq_along(nExp), function(i){
  round(nExp[[i]] * (1 - homotypic.prop[[i]]))
})
nExp.adj
```

```
## [[1]]
## [1] 30
## 
## [[2]]
## [1] 75
## 
## [[3]]
## [1] 164
```

## Doublet detection

```r
experiment.split <- lapply(seq_along(experiment.split), function(i){
  doubletFinder_v3(experiment.split[[i]],
                   PCs = 1:50,
                   pN = 0.25,
                   pK = as.numeric(as.character(pK[[i]])),
                   nExp = nExp.adj[[i]],
                   reuse.pANN = FALSE,
                   sct = FALSE)
})
```

## Aggreate doublet calls
The doublet calls are currently split over three objects. Let's aggregate them into a named character vector.

```r
calls.list <- lapply(experiment.split, function(sce){
  calls = sce@meta.data[,grep("DF.classifications", colnames(sce@meta.data))]
  names(calls) = rownames(sce@meta.data)
  calls
})
doublet.calls <- unlist(calls.list)
head(doublet.calls)
```

```
## AAACCCAAGTTATGGA_A001-C-007 AAACGCTTCTCTGCTG_A001-C-007 
##                   "Doublet"                   "Doublet" 
## AAAGAACGTGCTTATG_A001-C-007 AAAGAACGTTTCGCTC_A001-C-007 
##                   "Doublet"                   "Doublet" 
## AAAGGATTCATTACCT_A001-C-007 AAAGTGACACGCTTAA_A001-C-007 
##                   "Doublet"                   "Doublet"
```

## Add doublet calls to aggregated object
Now that the doublet calls are in a single named vector, they can be added as metadata to the aggregate Seurat object.

```r
experiment.aggregate <- AddMetaData(experiment.aggregate,
                                    metadata = doublet.calls,
                                    col.name = "doublet_call")
DimPlot(experiment.aggregate,
        reduction = "umap",
        group.by = "doublet_call",
        pt.size = 0.1,
        shuffle = TRUE) +
  scale_color_manual(values = c("red", "gray")) +
  theme(plot.title = element_blank())
```

![](07-doublet_detection_files/figure-html/add_calls-1.png)<!-- -->

## Remove doublets

```r
experiment.aggregate <- subset(experiment.aggregate, doublet_call == "Singlet")
experiment.aggregate$doublet_call <- NULL # remove redundant column from metadata
```

## Prepare for the next section

#### Save object

```r
saveRDS(experiment.aggregate, file = "scRNA_workshop-07.rds")
```

#### Download Rmd document

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/08-integration.Rmd", "08-integration.Rmd")
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
##  [1] ROCR_1.0-11         KernSmooth_2.23-20  fields_14.1        
##  [4] viridis_0.6.2       viridisLite_0.4.1   spam_2.9-1         
##  [7] ggplot2_3.4.0       DoubletFinder_2.0.3 SeuratObject_4.1.3 
## [10] Seurat_4.3.0       
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
##  [34] htmltools_0.5.3        tools_4.2.2            dotCall64_1.0-2       
##  [37] igraph_1.3.5           gtable_0.3.1           glue_1.6.2            
##  [40] RANN_2.6.1             reshape2_1.4.4         dplyr_1.0.10          
##  [43] maps_3.4.1             Rcpp_1.0.9             scattermore_0.8       
##  [46] jquerylib_0.1.4        vctrs_0.5.1            nlme_3.1-160          
##  [49] spatstat.explore_3.0-5 progressr_0.11.0       lmtest_0.9-40         
##  [52] spatstat.random_3.0-1  xfun_0.35              stringr_1.4.1         
##  [55] globals_0.16.2         mime_0.12              miniUI_0.1.1.1        
##  [58] lifecycle_1.0.3        irlba_2.3.5.1          goftest_1.2-3         
##  [61] future_1.29.0          MASS_7.3-58.1          zoo_1.8-11            
##  [64] scales_1.2.1           promises_1.2.0.1       spatstat.utils_3.0-1  
##  [67] parallel_4.2.2         RColorBrewer_1.1-3     yaml_2.3.6            
##  [70] reticulate_1.28        pbapply_1.6-0          gridExtra_2.3         
##  [73] sass_0.4.4             stringi_1.7.8          highr_0.9             
##  [76] rlang_1.0.6            pkgconfig_2.0.3        matrixStats_0.63.0    
##  [79] evaluate_0.18          lattice_0.20-45        tensor_1.5            
##  [82] purrr_0.3.5            labeling_0.4.2         patchwork_1.1.2       
##  [85] htmlwidgets_1.5.4      cowplot_1.1.1          tidyselect_1.2.0      
##  [88] parallelly_1.32.1      RcppAnnoy_0.0.20       plyr_1.8.8            
##  [91] magrittr_2.0.3         R6_2.5.1               generics_0.1.3        
##  [94] DBI_1.1.3              withr_2.5.0            pillar_1.8.1          
##  [97] fitdistrplus_1.1-8     survival_3.4-0         abind_1.4-5           
## [100] sp_1.5-1               tibble_3.1.8           future.apply_1.10.0   
## [103] crayon_1.5.2           utf8_1.2.2             spatstat.geom_3.0-3   
## [106] plotly_4.10.1          rmarkdown_2.18         grid_4.2.2            
## [109] data.table_1.14.6      digest_0.6.30          xtable_1.8-4          
## [112] tidyr_1.2.1            httpuv_1.6.6           munsell_0.5.0         
## [115] bslib_0.4.1
```
