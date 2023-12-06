---
title: "Introduction to Single Cell RNA-Seq Part 1: Create Seurat object"
author: "UCD Bioinformatics Core"
date: "2023-12-06"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 1: Create Seurat object


Our first Markdown document concentrates on getting data into R and setting up our initial object. We will also replicate some of the tables and figures found in the Cellranger web summary.

## Load packages
We will start each section by loading the libraries necessary for that portion of the analysis.

```r
library(Seurat)     # single cell RNA-Seq analysis
library(kableExtra) # format tables
library(ggplot2)   # create graphics
library(viridis)   # accessible color palettes
```

## Experiment metadata
The metadata we have available for this subset of the [Becker experiment](https://www.nature.com/articles/s41588-022-01088-x) during this workshop is very basic; we don't have a record of patient identifiers, biopsy dates, treatment course, or prognosis. Instead, for each sample, we know the group (healthy, polyp, or cancerous tissue) and the sequencing run, which we can derive from the read header.
Let's create a data table containing this information.

```r
experiment.metadata <- data.frame(id = c("A001-C-007",
                                         "A001-C-104",
                                         "B001-A-301"),
                                  group = c("Colorectal Cancer",
                                            "Polyp",
                                            "Normal"),
                                  run = c("A00509:126:HTLFWDMXX:1",
                                          "A00509:116:HTLNJDMXX:1",
                                          "A00509:113:HTNCWDMXX:1"))
experiment.metadata %>%
  kable() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> id </th>
   <th style="text-align:left;"> group </th>
   <th style="text-align:left;"> run </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> A001-C-007 </td>
   <td style="text-align:left;"> Colorectal Cancer </td>
   <td style="text-align:left;"> A00509:126:HTLFWDMXX:1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A001-C-104 </td>
   <td style="text-align:left;"> Polyp </td>
   <td style="text-align:left;"> A00509:116:HTLNJDMXX:1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> B001-A-301 </td>
   <td style="text-align:left;"> Normal </td>
   <td style="text-align:left;"> A00509:113:HTNCWDMXX:1 </td>
  </tr>
</tbody>
</table>

## Create metrics tables
The **expression_data_cellranger.zip** file that we have downloaded in previous step contains the single cell matrix files and HDF5 files for three single nuclei RNASeq samples from [Becker et al., 2022](https://www.nature.com/articles/s41588-022-01088-x). After un-compressing the file, please make sure that you see three folders (A001-C-007, A001-C-104, and B001-A-301) in the same folder as this R markdown file. If the three folders are located elsewhere, please change the assignment of "dataset.loc" in the code box below to reflect the location of your data.

```r
experiment.name <- "Becker 2022 colorectal cancer continuum"
dataset.loc <- "./"
```

In this section, the metrics_summary.csv files produced by Cellranger are used to create a single table summarizing the sequencing metrics for each sample.

```r
sample.metrics <- lapply(experiment.metadata$id, function(id){
  metrics = read.csv(file.path(dataset.loc,
                               paste0(id,"/outs"),
                               "metrics_summary.csv"),
                     colClasses = "character")
})
experiment.metrics <- do.call("rbind", sample.metrics)
rownames(experiment.metrics) <- experiment.metadata$id

sequencing.metrics <- data.frame(t(experiment.metrics[,c(1:19)]))

rownames(sequencing.metrics) <- gsub("\\."," ", rownames(sequencing.metrics))

sequencing.metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Overview", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Sequencing Characteristics", 4, 9, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 10, 19, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Cell Ranger Results</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> A001.C.007 </th>
   <th style="text-align:left;"> A001.C.104 </th>
   <th style="text-align:left;"> B001.A.301 </th>
  </tr>
 </thead>
<tbody>
  <tr grouplength="3"><td colspan="4" style="background-color: #666; color: #fff;"><strong>Overview</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Estimated Number of Cells </td>
   <td style="text-align:left;"> 1,798 </td>
   <td style="text-align:left;"> 3,144 </td>
   <td style="text-align:left;"> 4,514 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Mean Reads per Cell </td>
   <td style="text-align:left;"> 77,438 </td>
   <td style="text-align:left;"> 151,221 </td>
   <td style="text-align:left;"> 38,935 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median Genes per Cell </td>
   <td style="text-align:left;"> 927 </td>
   <td style="text-align:left;"> 959 </td>
   <td style="text-align:left;"> 1,331 </td>
  </tr>
  <tr grouplength="6"><td colspan="4" style="background-color: #666; color: #fff;"><strong>Sequencing Characteristics</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Number of Reads </td>
   <td style="text-align:left;"> 139,233,487 </td>
   <td style="text-align:left;"> 475,437,350 </td>
   <td style="text-align:left;"> 175,752,014 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Valid Barcodes </td>
   <td style="text-align:left;"> 97.4% </td>
   <td style="text-align:left;"> 95.4% </td>
   <td style="text-align:left;"> 98.5% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Sequencing Saturation </td>
   <td style="text-align:left;"> 75.5% </td>
   <td style="text-align:left;"> 84.3% </td>
   <td style="text-align:left;"> 69.0% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Q30 Bases in Barcode </td>
   <td style="text-align:left;"> 97.4% </td>
   <td style="text-align:left;"> 96.4% </td>
   <td style="text-align:left;"> 96.6% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Q30 Bases in RNA Read </td>
   <td style="text-align:left;"> 95.6% </td>
   <td style="text-align:left;"> 94.4% </td>
   <td style="text-align:left;"> 94.1% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Q30 Bases in UMI </td>
   <td style="text-align:left;"> 97.5% </td>
   <td style="text-align:left;"> 96.4% </td>
   <td style="text-align:left;"> 96.5% </td>
  </tr>
  <tr grouplength="10"><td colspan="4" style="background-color: #666; color: #fff;"><strong>Mapping Characteristics</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped to Genome </td>
   <td style="text-align:left;"> 94.3% </td>
   <td style="text-align:left;"> 92.1% </td>
   <td style="text-align:left;"> 89.0% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Genome </td>
   <td style="text-align:left;"> 72.3% </td>
   <td style="text-align:left;"> 49.1% </td>
   <td style="text-align:left;"> 83.8% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Intergenic Regions </td>
   <td style="text-align:left;"> 8.2% </td>
   <td style="text-align:left;"> 7.4% </td>
   <td style="text-align:left;"> 5.0% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Intronic Regions </td>
   <td style="text-align:left;"> 26.9% </td>
   <td style="text-align:left;"> 20.6% </td>
   <td style="text-align:left;"> 40.2% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Exonic Regions </td>
   <td style="text-align:left;"> 37.2% </td>
   <td style="text-align:left;"> 21.2% </td>
   <td style="text-align:left;"> 38.6% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Transcriptome </td>
   <td style="text-align:left;"> 61.1% </td>
   <td style="text-align:left;"> 39.4% </td>
   <td style="text-align:left;"> 73.4% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Antisense to Gene </td>
   <td style="text-align:left;"> 2.4% </td>
   <td style="text-align:left;"> 1.9% </td>
   <td style="text-align:left;"> 4.8% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Fraction Reads in Cells </td>
   <td style="text-align:left;"> 29.6% </td>
   <td style="text-align:left;"> 37.7% </td>
   <td style="text-align:left;"> 36.8% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Total Genes Detected </td>
   <td style="text-align:left;"> 23,930 </td>
   <td style="text-align:left;"> 25,327 </td>
   <td style="text-align:left;"> 25,462 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median UMI Counts per Cell </td>
   <td style="text-align:left;"> 1,201 </td>
   <td style="text-align:left;"> 1,231 </td>
   <td style="text-align:left;"> 1,913 </td>
  </tr>
</tbody>
</table>

This roughly replicates the table that appears in the Cellranger web summary file.

## Create Seurat object
We will be using [Seurat](http://satijalab.org/seurat/) (currently version 4) as the basis of our single cell (or nucleus) RNA-Seq analysis. Seurat is a popular R package that is designed for QC, analysis, and exploration of single cell data, which aims to enable users to identify and interpret sources of heterogeneity from single cell transcriptomic measurements, and to integrate diverse types of single cell data.
In addition to the standard Seurat workflow, this documentation makes use of some custom code, and brings in functions from other packages. For additional information on Seurat standard workflows, see the authors' [tutorials](https://satijalab.org/seurat/vignettes.html).

#### Read in expression matrix
First, we read in data from each individual sample folder.

```r
expression.data <- lapply(experiment.metadata$id, function(id){
  sample.matrix = Read10X_h5(file.path(dataset.loc, id, "/outs","filtered_feature_bc_matrix.h5"))
  colnames(sample.matrix) = paste(sapply(strsplit(colnames(sample.matrix),split="-"), '[[', 1L), id, sep="_")
  sample.matrix
})
names(expression.data) <- experiment.metadata$id
```


```r
View(expression.data)
```

#### Merge matrices

```r
aggregate.data <- do.call("cbind", expression.data)
```

#### Create object
The `CreateSeuratObject` function allows feature (gene) and cell filtering by minimum cell and feature counts. We will set these to 0 for now in order to explore manual filtering more fully in part 2.

```r
experiment.aggregate <- CreateSeuratObject(
  aggregate.data,
  project = experiment.name,
  min.cells = 0,
  min.features = 0,
  names.field = 2, # tells Seurat which part of the cell identifier contains the sample name
  names.delim = "\\_")
```

## Add metadata
We can now attach the metadata in our table to the Seurat object.

#### Match metadata to expression matrix
The columns of the expression matrix correspond to the cells in the experiment. When we created the Seurat object, the "names.field" and "names.delim" arguments allowed Seurat to infer sample identity from the cell names. This information is stored in a variable called "orig.ident."

```r
levels(experiment.aggregate$orig.ident)
```

```
## [1] "A001-C-007" "A001-C-104" "B001-A-301"
```

These sample identifiers are stored in the experiment.metadata object as well, which allows us to match the other metadata contained within that table to the correct cells within the Seurat object.

```r
sample.index <- match(experiment.aggregate$orig.ident, experiment.metadata$id)
```

#### Attach metadata
The AddMetaData function returns a new Seurat object with an additional column in the metadata table containing the new information.

```r
experiment.aggregate <- AddMetaData(experiment.aggregate,
                                    metadata = experiment.metadata$group[sample.index],
            col.name = "group")
experiment.aggregate$group <- factor(experiment.aggregate$group,
                                     levels = c("Normal", "Polyp", "Colorectal Cancer"))
experiment.aggregate <- AddMetaData(experiment.aggregate,
                                    metadata = experiment.metadata$run[sample.index],
                                    col.name = "run")
experiment.aggregate$run <- factor(experiment.aggregate$run,
                                   levels = c("A00509:113:HTNCWDMXX:1",
                                              "A00509:116:HTLNJDMXX:1",
                                              "A00509:126:HTLFWDMXX:1"))
```

## Explore the Seurat object
A Seurat object is a complex data structure containing the data from a single cell or single nucleus assay and **all** of the information associated with the experiment, including annotations, analysis, and more. This data structure was developed by the authors of the Seurat analysis package, for use with their pipeline.

```r
View(experiment.aggregate)
```

Most Seurat functions take the object as an argument, and return either a new Seurat object or a ggplot object (a visualization). As the analysis continues, more and more data will be added to the object.


```r
slotNames(experiment.aggregate)
```

```
##  [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"      
##  [6] "neighbors"    "reductions"   "images"       "project.name" "misc"        
## [11] "version"      "commands"     "tools"
```

```r
experiment.aggregate@assays # a slot is accessed with the @ symbol
```

```
## $RNA
## Assay data with 36601 features for 9456 cells
## First 10 features:
##  MIR1302-2HG, FAM138A, OR4F5, AL627309.1, AL627309.3, AL627309.2,
## AL627309.5, AL627309.4, AP006222.2, AL732372.1
```

- Which slots are empty, and which contain data?
- What type of object is the content of the meta.data slot?
- What metadata is available?

There is often more than one way to interact with the information stored in each of a Seurat objects many slots. The default behaviors of different access functions are described in the help documentation.


```r
# which slot is being accessed here? find another way to produce the result
head(experiment.aggregate[[]])
```

```
##                             orig.ident nCount_RNA nFeature_RNA
## AAACCCAAGTTATGGA_A001-C-007 A001-C-007       2078         1549
## AAACCCACAACGCCCA_A001-C-007 A001-C-007        854          687
## AAACCCACAGAAGTTA_A001-C-007 A001-C-007        541          467
## AAACCCAGTCAGTCCG_A001-C-007 A001-C-007        605          538
## AAACGAAGTTGGTGTT_A001-C-007 A001-C-007        954          772
## AAACGCTAGGAGCAAA_A001-C-007 A001-C-007        734          629
##                                         group                    run
## AAACCCAAGTTATGGA_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
## AAACCCACAACGCCCA_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
## AAACCCACAGAAGTTA_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
## AAACCCAGTCAGTCCG_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
## AAACGAAGTTGGTGTT_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
## AAACGCTAGGAGCAAA_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
```

The use of syntax is often a matter of personal preference. In the interest of clarity, this documentation will generally use the more explicit syntax, with a few exceptions.

## Barcode inflection plots

Imagine the barcode rank plot from the Cell Ranger web summary. That graphic plots the number of UMIs against the barcode rank, and typically has a sharp inflection point where the number of UMIs drops dramatically. These points can represent a transition between cell types from a higher RNA content population to a lower RNA content population, or from cell-associated barcodes to background.

The Seurat `BarcodeInflectionsPlot` provides a similar graphic. In this case, because we are using the filtered barcode matrix, rather than all barcodes, much of the background is absent from the plot.


```r
experiment.aggregate <- CalculateBarcodeInflections(experiment.aggregate)
BarcodeInflectionsPlot(experiment.aggregate) +
  scale_color_viridis_d()
```

![](01-create_object_files/figure-html/barcode_inflection_plot-1.png)<!-- -->

Adding a log-scale transformation to the x-axis increases the resemblance to the Cell Ranger plot. Values on the y-axis are already log-transformed.


```r
BarcodeInflectionsPlot(experiment.aggregate) +
  scale_x_continuous(trans = "log10") +
  scale_color_viridis_d()
```

![](01-create_object_files/figure-html/barcode_inflection_plot_log-1.png)<!-- -->

## Prepare for the next section

#### Save object

```r
saveRDS(experiment.aggregate, file="scRNA_workshop-01.rds")
```

#### Download Rmd

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/02-filtering.Rmd", "02-filtering.Rmd")
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
## [1] viridis_0.6.2      viridisLite_0.4.1  ggplot2_3.4.0      kableExtra_1.3.4  
## [5] SeuratObject_4.1.3 Seurat_4.3.0      
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6          
##   [4] ellipsis_0.3.2         ggridges_0.5.4         rstudioapi_0.14       
##   [7] spatstat.data_3.0-0    farver_2.1.1           leiden_0.4.3          
##  [10] listenv_0.8.0          bit64_4.0.5            ggrepel_0.9.2         
##  [13] fansi_1.0.3            xml2_1.3.3             codetools_0.2-18      
##  [16] splines_4.2.2          cachem_1.0.6           knitr_1.41            
##  [19] polyclip_1.10-4        jsonlite_1.8.4         ica_1.0-3             
##  [22] cluster_2.1.4          png_0.1-8              uwot_0.1.14           
##  [25] shiny_1.7.3            sctransform_0.3.5      spatstat.sparse_3.0-0 
##  [28] compiler_4.2.2         httr_1.4.4             assertthat_0.2.1      
##  [31] Matrix_1.5-3           fastmap_1.1.0          lazyeval_0.2.2        
##  [34] cli_3.4.1              later_1.3.0            htmltools_0.5.3       
##  [37] tools_4.2.2            igraph_1.3.5           gtable_0.3.1          
##  [40] glue_1.6.2             RANN_2.6.1             reshape2_1.4.4        
##  [43] dplyr_1.0.10           Rcpp_1.0.9             scattermore_0.8       
##  [46] jquerylib_0.1.4        vctrs_0.5.1            svglite_2.1.0         
##  [49] nlme_3.1-160           spatstat.explore_3.0-5 progressr_0.11.0      
##  [52] lmtest_0.9-40          spatstat.random_3.0-1  xfun_0.35             
##  [55] stringr_1.4.1          globals_0.16.2         rvest_1.0.3           
##  [58] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
##  [61] irlba_2.3.5.1          goftest_1.2-3          future_1.29.0         
##  [64] MASS_7.3-58.1          zoo_1.8-11             scales_1.2.1          
##  [67] promises_1.2.0.1       spatstat.utils_3.0-1   parallel_4.2.2        
##  [70] RColorBrewer_1.1-3     yaml_2.3.6             reticulate_1.28       
##  [73] pbapply_1.6-0          gridExtra_2.3          sass_0.4.4            
##  [76] stringi_1.7.8          highr_0.9              systemfonts_1.0.4     
##  [79] rlang_1.0.6            pkgconfig_2.0.3        matrixStats_0.63.0    
##  [82] evaluate_0.18          lattice_0.20-45        tensor_1.5            
##  [85] ROCR_1.0-11            purrr_0.3.5            labeling_0.4.2        
##  [88] patchwork_1.1.2        htmlwidgets_1.5.4      bit_4.0.5             
##  [91] cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.32.1     
##  [94] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
##  [97] R6_2.5.1               generics_0.1.3         DBI_1.1.3             
## [100] withr_2.5.0            pillar_1.8.1           fitdistrplus_1.1-8    
## [103] survival_3.4-0         abind_1.4-5            sp_1.5-1              
## [106] tibble_3.1.8           future.apply_1.10.0    hdf5r_1.3.7           
## [109] KernSmooth_2.23-20     utf8_1.2.2             spatstat.geom_3.0-3   
## [112] plotly_4.10.1          rmarkdown_2.18         grid_4.2.2            
## [115] data.table_1.14.6      webshot_0.5.4          digest_0.6.30         
## [118] xtable_1.8-4           tidyr_1.2.1            httpuv_1.6.6          
## [121] munsell_0.5.0          bslib_0.4.1
```
