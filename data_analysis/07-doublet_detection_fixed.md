---
title: "Introduction to Single Cell RNA-Seq Part 7: Doublet detection"
author: "UCD Bioinformatics Core"
date: "2023-07-28"
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

<div class='r_output'> $`A001-C-007`
 [1] 30
 
 $`A001-C-104`
 [1] 1
 
 $`B001-A-301`
 [1] 10
</div>
### nExp
In single cell data generated using microfluidics, the frequency of all multiplets (both homotypic and heterotypic) can be modeled with a Poisson distribution. The probability of capturing multiple cells in a droplet is a function of the loading density of the microfluidic device. However, because homotypic doublets are far more difficult to detect informatically, the Poisson distribution overestimates the number of *detectable* doublets.

Here we will use the 10x expected doublet rate of 8% for chips that were loaded with 16,500 cells (targeting 10,000 cells for recovery).

```r
nExp <- lapply(experiment.split, function(sce){
  round(0.08*dim(sce)[2])
})
nExp
```

<div class='r_output'> $`A001-C-007`
 [1] 82
 
 $`A001-C-104`
 [1] 149
 
 $`B001-A-301`
 [1] 274
</div>#### Homotypic doublet adjustment
If cell type annotations are available, the following code can be used to adjust the estimated detectable doublet frequency downwards by adjusting for the expected frequency of homotypic doublets.

Let's read in the cell type assignments from the previous part and attach them to the split object.

```r
experiment.aggregate <- readRDS("scRNA_workshop-06.rds")
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

<div class='r_output'> [[1]]
 [1] 30
 
 [[2]]
 [1] 88
 
 [[3]]
 [1] 163
</div>
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

<div class='r_output'> AAACCCAAGTTATGGA_A001-C-007 AAACGCTTCTCTGCTG_A001-C-007 
                   "Doublet"                   "Doublet" 
 AAAGAACGTGCTTATG_A001-C-007 AAAGAACGTTTCGCTC_A001-C-007 
                   "Doublet"                   "Doublet" 
 AAAGGATTCATTACCT_A001-C-007 AAAGTGACACGCTTAA_A001-C-007 
                   "Doublet"                   "Doublet"
</div>
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
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/08-de_enrichment.Rmd", "08-de_enrichment.Rmd")
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
  [1] ROCR_1.0-11         KernSmooth_2.23-22  fields_14.1        
  [4] viridis_0.6.3       viridisLite_0.4.2   spam_2.9-1         
  [7] ggplot2_3.4.2       DoubletFinder_2.0.3 SeuratObject_4.1.3 
 [10] Seurat_4.3.0.1     
 
 loaded via a namespace (and not attached):
   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
   [4] rlang_1.1.1            magrittr_2.0.3         RcppAnnoy_0.0.21      
   [7] spatstat.geom_3.2-4    matrixStats_1.0.0      ggridges_0.5.4        
  [10] compiler_4.3.1         maps_3.4.1             png_0.1-8             
  [13] vctrs_0.6.3            reshape2_1.4.4         stringr_1.5.0         
  [16] pkgconfig_2.0.3        fastmap_1.1.1          ellipsis_0.3.2        
  [19] labeling_0.4.2         utf8_1.2.3             promises_1.2.0.1      
  [22] rmarkdown_2.23         purrr_1.0.1            xfun_0.39             
  [25] cachem_1.0.8           jsonlite_1.8.7         goftest_1.2-3         
  [28] highr_0.10             later_1.3.1            spatstat.utils_3.0-3  
  [31] irlba_2.3.5.1          parallel_4.3.1         cluster_2.1.4         
  [34] R6_2.5.1               ica_1.0-3              spatstat.data_3.0-1   
  [37] bslib_0.5.0            stringi_1.7.12         RColorBrewer_1.1-3    
  [40] reticulate_1.30        parallelly_1.36.0      lmtest_0.9-40         
  [43] jquerylib_0.1.4        scattermore_1.2        Rcpp_1.0.11           
  [46] knitr_1.43             tensor_1.5             future.apply_1.11.0   
  [49] zoo_1.8-12             sctransform_0.3.5      httpuv_1.6.11         
  [52] Matrix_1.6-0           splines_4.3.1          igraph_1.5.0          
  [55] tidyselect_1.2.0       abind_1.4-5            rstudioapi_0.15.0     
  [58] yaml_2.3.7             spatstat.random_3.1-5  codetools_0.2-19      
  [61] miniUI_0.1.1.1         spatstat.explore_3.2-1 listenv_0.9.0         
  [64] lattice_0.21-8         tibble_3.2.1           plyr_1.8.8            
  [67] withr_2.5.0            shiny_1.7.4.1          evaluate_0.21         
  [70] Rtsne_0.16             future_1.33.0          survival_3.5-5        
  [73] polyclip_1.10-4        fitdistrplus_1.1-11    pillar_1.9.0          
  [76] plotly_4.10.2          generics_0.1.3         sp_2.0-0              
  [79] munsell_0.5.0          scales_1.2.1           globals_0.16.2        
  [82] xtable_1.8-4           glue_1.6.2             lazyeval_0.2.2        
  [85] tools_4.3.1            data.table_1.14.8      RANN_2.6.1            
  [88] dotCall64_1.0-2        leiden_0.4.3           cowplot_1.1.1         
  [91] grid_4.3.1             tidyr_1.3.0            colorspace_2.1-0      
  [94] nlme_3.1-162           patchwork_1.1.2        cli_3.6.1             
  [97] spatstat.sparse_3.0-2  fansi_1.0.4            dplyr_1.1.2           
 [100] uwot_0.1.16            gtable_0.3.3           sass_0.4.7            
 [103] digest_0.6.33          progressr_0.13.0       ggrepel_0.9.3         
 [106] farver_2.1.1           htmlwidgets_1.6.2      htmltools_0.5.5       
 [109] lifecycle_1.0.3        httr_1.4.6             mime_0.12             
 [112] MASS_7.3-60
</div>