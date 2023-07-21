---
title: "Introduction to Single Cell RNA-Seq Part 1: Create Seurat object"
author: "UCD Bioinformatics Core"
date: "2023-07-21"
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
   <td style="text-align:left;"> 1,796 </td>
   <td style="text-align:left;"> 3,142 </td>
   <td style="text-align:left;"> 4,514 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Mean Reads per Cell </td>
   <td style="text-align:left;"> 77,524 </td>
   <td style="text-align:left;"> 151,317 </td>
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
   <td style="text-align:left;"> 37.6% </td>
   <td style="text-align:left;"> 36.8% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Total Genes Detected </td>
   <td style="text-align:left;"> 23,930 </td>
   <td style="text-align:left;"> 25,326 </td>
   <td style="text-align:left;"> 25,462 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median UMI Counts per Cell </td>
   <td style="text-align:left;"> 1,204 </td>
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

<div class='r_output'> [1] "A001-C-007" "A001-C-104" "B001-A-301"
</div>
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
experiment.aggregate <- AddMetaData(experiment.aggregate,
                                    metadata = experiment.metadata$run[sample.index],
                                    col.name = "run")
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

<div class='r_output'>  [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"      
  [6] "neighbors"    "reductions"   "images"       "project.name" "misc"        
 [11] "version"      "commands"     "tools"
</div>
```r
experiment.aggregate@assays # a slot is accessed with the @ symbol
```

<div class='r_output'> $RNA
 Assay data with 36601 features for 9452 cells
 First 10 features:
  MIR1302-2HG, FAM138A, OR4F5, AL627309.1, AL627309.3, AL627309.2,
 AL627309.5, AL627309.4, AP006222.2, AL732372.1
</div>
- Which slots are empty, and which contain data?
- What type of object is the content of the meta.data slot?
- What metadata is available?

There is often more than one way to interact with the information stored in each of a Seurat objects many slots. The default behaviors of different access functions are described in the help documentation.


```r
# which slot is being accessed here? find another way to produce the result
head(experiment.aggregate[[]])
```

<div class='r_output'>                             orig.ident nCount_RNA nFeature_RNA
 AAACCCAAGTTATGGA_A001-C-007 A001-C-007       2078         1549
 AAACCCACAACGCCCA_A001-C-007 A001-C-007        854          687
 AAACCCACAGAAGTTA_A001-C-007 A001-C-007        541          467
 AAACCCAGTCAGTCCG_A001-C-007 A001-C-007        605          538
 AAACGAAGTTGGTGTT_A001-C-007 A001-C-007        954          772
 AAACGCTAGGAGCAAA_A001-C-007 A001-C-007        734          629
                                         group                    run
 AAACCCAAGTTATGGA_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
 AAACCCACAACGCCCA_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
 AAACCCACAGAAGTTA_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
 AAACCCAGTCAGTCCG_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
 AAACGAAGTTGGTGTT_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
 AAACGCTAGGAGCAAA_A001-C-007 Colorectal Cancer A00509:126:HTLFWDMXX:1
</div>
The use of syntax is often a matter of personal preference. In the interest of clarity, this documentation will generally use the more explicit syntax, with a few exceptions.

## Barcode inflection plots

Imagine the barcode rank plot from the Cell Ranger web summary. That graphic plots the number of UMIs against the barcode rank, and typically has a sharp inflection point where the number of UMIs drops dramatically. These points can represent a transition between cell types from a higher RNA content population to a lower RNA content population, or from cell-associated barcodes to background.

The Seurat `BarcodeInflectionsPlot` provides a similar graphic. In this case, because we are using the filtered barcode matrix, rather than all barcodes, much of the background is absent from the plot.


```r
experiment.aggregate <- CalculateBarcodeInflections(experiment.aggregate)
BarcodeInflectionsPlot(experiment.aggregate) +
  scale_color_viridis_d(option = "mako")
```

![](01-create_object_files/figure-html/barcode_inflection_plot-1.png)<!-- -->

Adding a log-scale transformation to the x-axis increases the resemblance to the Cell Ranger plot. Values on the y-axis are already log-transformed.


```r
BarcodeInflectionsPlot(experiment.aggregate) +
  scale_x_continuous(trans = "log10") +
  scale_color_viridis_d(option = "mako")
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
 [1] viridis_0.6.3      viridisLite_0.4.2  ggplot2_3.4.2      kableExtra_1.3.4  
 [5] SeuratObject_4.1.3 Seurat_4.3.0.1    
 
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
  [46] compiler_4.3.1         bit64_4.0.5            withr_2.5.0           
  [49] highr_0.10             MASS_7.3-60            tools_4.3.1           
  [52] lmtest_0.9-40          httpuv_1.6.11          future.apply_1.11.0   
  [55] goftest_1.2-3          glue_1.6.2             nlme_3.1-162          
  [58] promises_1.2.0.1       grid_4.3.1             Rtsne_0.16            
  [61] cluster_2.1.4          reshape2_1.4.4         generics_0.1.3        
  [64] hdf5r_1.3.8            gtable_0.3.3           spatstat.data_3.0-1   
  [67] tidyr_1.3.0            data.table_1.14.8      sp_2.0-0              
  [70] xml2_1.3.5             utf8_1.2.3             spatstat.geom_3.2-4   
  [73] RcppAnnoy_0.0.21       ggrepel_0.9.3          RANN_2.6.1            
  [76] pillar_1.9.0           stringr_1.5.0          later_1.3.1           
  [79] splines_4.3.1          dplyr_1.1.2            lattice_0.21-8        
  [82] bit_4.0.5              survival_3.5-5         deldir_1.0-9          
  [85] tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2         
  [88] knitr_1.43             gridExtra_2.3          svglite_2.1.1         
  [91] scattermore_1.2        xfun_0.39              matrixStats_1.0.0     
  [94] stringi_1.7.12         lazyeval_0.2.2         yaml_2.3.7            
  [97] evaluate_0.21          codetools_0.2-19       tibble_3.2.1          
 [100] cli_3.6.1              uwot_0.1.16            xtable_1.8-4          
 [103] reticulate_1.30        systemfonts_1.0.4      munsell_0.5.0         
 [106] jquerylib_0.1.4        Rcpp_1.0.11            globals_0.16.2        
 [109] spatstat.random_3.1-5  png_0.1-8              parallel_4.3.1        
 [112] ellipsis_0.3.2         listenv_0.9.0          scales_1.2.1          
 [115] ggridges_0.5.4         leiden_0.4.3           purrr_1.0.1           
 [118] rlang_1.1.1            cowplot_1.1.1          rvest_1.0.3
</div>