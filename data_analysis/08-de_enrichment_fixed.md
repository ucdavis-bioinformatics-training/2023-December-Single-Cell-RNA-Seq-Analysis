---
title: "Introduction to Single Cell RNA-Seq Part 8: Enrichment and model-based differential expression"
author: "Bioinformatics Core"
date: "2023-07-28"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 8: Enrichment and model-based differential expression


## Set up workspace

```r
library(Seurat)
library(limma)
library(topGO)
library(dplyr)
library(kableExtra)
set.seed(12345)
experiment.aggregate <- readRDS("scRNA_workshop-07.rds")
Idents(experiment.aggregate) <- "finalcluster"
```

## 1. Gene Ontology (GO) Enrichment of Genes Expressed in a Cluster
[Gene Ontology](http://geneontology.org/docs/ontology-documentation/) provides a controlled vocabulary for describing gene products.  Here we use enrichment analysis to identify GO terms that are over-represented among the gene expressed in cells in a given cluster. 

```r
cluster10 <- subset(experiment.aggregate, idents = '10')
expr <- as.matrix(GetAssayData(cluster10))

# Select genes that are expressed > 0 in at least half of cells
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.5)]
all.genes <- rownames(expr)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
	GOdata <- new("topGOdata",
		ontology = "BP", # use biological process ontology
		allGenes = geneList,
		geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
# Test for enrichment using Fisher's Exact Test
	resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
	GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> GO.ID </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Term </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Annotated </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Significant </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Expected </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Fisher </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0006325 </td>
   <td style="text-align:left;"> chromatin organization </td>
   <td style="text-align:right;"> 465 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 6.02 </td>
   <td style="text-align:left;"> 2.3e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0048024 </td>
   <td style="text-align:left;"> regulation of mRNA splicing, via spliceosome </td>
   <td style="text-align:right;"> 81 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1.05 </td>
   <td style="text-align:left;"> 0.00061 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0061061 </td>
   <td style="text-align:left;"> muscle structure development </td>
   <td style="text-align:right;"> 361 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 4.67 </td>
   <td style="text-align:left;"> 0.00077 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0031580 </td>
   <td style="text-align:left;"> membrane raft distribution </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:left;"> 0.00098 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006357 </td>
   <td style="text-align:left;"> regulation of transcription by RNA polymerase II </td>
   <td style="text-align:right;"> 1456 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 18.85 </td>
   <td style="text-align:left;"> 0.00138 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0031047 </td>
   <td style="text-align:left;"> RNA-mediated gene silencing </td>
   <td style="text-align:right;"> 95 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1.23 </td>
   <td style="text-align:left;"> 0.00141 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0009791 </td>
   <td style="text-align:left;"> post-embryonic development </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.84 </td>
   <td style="text-align:left;"> 0.00149 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1903613 </td>
   <td style="text-align:left;"> regulation of protein tyrosine phosphatase activity </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:left;"> 0.00162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1903978 </td>
   <td style="text-align:left;"> regulation of microglial cell activation </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:left;"> 0.00162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045060 </td>
   <td style="text-align:left;"> negative thymic T cell selection </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:left;"> 0.00162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0070229 </td>
   <td style="text-align:left;"> negative regulation of lymphocyte apoptotic process </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.30 </td>
   <td style="text-align:left;"> 0.00310 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0072711 </td>
   <td style="text-align:left;"> cellular response to hydroxyurea </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:left;"> 0.00335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0034968 </td>
   <td style="text-align:left;"> histone lysine methylation </td>
   <td style="text-align:right;"> 81 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1.05 </td>
   <td style="text-align:left;"> 0.00392 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0009653 </td>
   <td style="text-align:left;"> anatomical structure morphogenesis </td>
   <td style="text-align:right;"> 1495 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 19.36 </td>
   <td style="text-align:left;"> 0.00420 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0060009 </td>
   <td style="text-align:left;"> Sertoli cell development </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:left;"> 0.00442 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0030889 </td>
   <td style="text-align:left;"> negative regulation of B cell proliferation </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:left;"> 0.00442 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045059 </td>
   <td style="text-align:left;"> positive thymic T cell selection </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:left;"> 0.00442 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0001771 </td>
   <td style="text-align:left;"> immunological synapse formation </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:left;"> 0.00442 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006376 </td>
   <td style="text-align:left;"> mRNA splice site selection </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.35 </td>
   <td style="text-align:left;"> 0.00494 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0010720 </td>
   <td style="text-align:left;"> positive regulation of cell development </td>
   <td style="text-align:right;"> 251 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 3.25 </td>
   <td style="text-align:left;"> 0.00521 </td>
  </tr>
</tbody>
</table>
* Annotated: number of genes (out of all.genes) that are annotated with that GO term
* Significant: number of genes that are annotated with that GO term and meet our criteria for "expressed"
* Expected: Under random chance, number of genes that would be expected to be annotated with that GO term and meeting our criteria for "expressed"
* Fisher: (Raw) p-value from Fisher's Exact Test

## 2. Model-based DE analysis in limma
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) is an R package for differential expression analysis of bulk RNASeq and microarray data.  We apply it here to single cell data.

Limma can be used to fit any linear model to expression data and is useful for analyses that go beyond two-group comparisons.  A detailed tutorial of model specification in limma is available [here](https://ucdavis-bioinformatics-training.github.io/2021-June-RNA-Seq-Analysis/data_analysis/DE_Analysis_mm_with_quizzes) and in the [limma User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf).


```r
# filter genes to those expressed in at least 10% of cells
keep <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.1)]
expr2 <- expr[keep,]

# Set up "design matrix" with statistical model
cluster10$proper.group <- make.names(cluster10$group)
mm <- model.matrix(~0 + proper.group + S.Score + G2M.Score + percent_MT + nFeature_RNA, data = cluster10[[]])
head(mm) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ACGGTTAGTCTCACAA_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0333657 </td>
   <td style="text-align:right;"> 0.7864107 </td>
   <td style="text-align:right;"> 0.5656109 </td>
   <td style="text-align:right;"> 1187 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACGTACAAGGAACTCG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.4162480 </td>
   <td style="text-align:right;"> 0.5902857 </td>
   <td style="text-align:right;"> 1.1111111 </td>
   <td style="text-align:right;"> 824 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACTGTCCAGGCGTTAG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0600477 </td>
   <td style="text-align:right;"> -0.0284044 </td>
   <td style="text-align:right;"> 0.5043712 </td>
   <td style="text-align:right;"> 1634 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACTGTCCGTATAGGAT_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0098200 </td>
   <td style="text-align:right;"> 0.0345753 </td>
   <td style="text-align:right;"> 1.0645376 </td>
   <td style="text-align:right;"> 1097 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AGAACCTTCGTTGTTT_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.5558142 </td>
   <td style="text-align:right;"> 0.5858386 </td>
   <td style="text-align:right;"> 0.6190581 </td>
   <td style="text-align:right;"> 2437 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AGCGTATAGCCTCAGC_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0548458 </td>
   <td style="text-align:right;"> -0.0674695 </td>
   <td style="text-align:right;"> 3.0966767 </td>
   <td style="text-align:right;"> 880 </td>
  </tr>
</tbody>
</table>

```r
tail(mm) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> TGGGTTACAAGAATGT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.1117535 </td>
   <td style="text-align:right;"> -0.1049327 </td>
   <td style="text-align:right;"> 0.2403846 </td>
   <td style="text-align:right;"> 883 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TGTTTGTTCACTACTT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.1121379 </td>
   <td style="text-align:right;"> -0.1127140 </td>
   <td style="text-align:right;"> 0.4361371 </td>
   <td style="text-align:right;"> 1064 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTCCGTGTCCGCTGTT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.1004205 </td>
   <td style="text-align:right;"> -0.0142675 </td>
   <td style="text-align:right;"> 0.6378455 </td>
   <td style="text-align:right;"> 994 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTCTTCCAGTCCCAAT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0498378 </td>
   <td style="text-align:right;"> -0.1072860 </td>
   <td style="text-align:right;"> 0.3696858 </td>
   <td style="text-align:right;"> 818 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTCTTCCCAGCGTAGA_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0536123 </td>
   <td style="text-align:right;"> 0.0009001 </td>
   <td style="text-align:right;"> 0.5504587 </td>
   <td style="text-align:right;"> 843 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTTACGTGTGTCTTAG_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0939035 </td>
   <td style="text-align:right;"> -0.1239464 </td>
   <td style="text-align:right;"> 0.5756579 </td>
   <td style="text-align:right;"> 882 </td>
  </tr>
</tbody>
</table>

```r
# Fit model in limma
fit <- lmFit(expr2, mm)
head(coef(fit)) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CCNL2 </td>
   <td style="text-align:right;"> 0.2623488 </td>
   <td style="text-align:right;"> 0.2428768 </td>
   <td style="text-align:right;"> 0.2945962 </td>
   <td style="text-align:right;"> -0.2984829 </td>
   <td style="text-align:right;"> 0.4853440 </td>
   <td style="text-align:right;"> 0.0505625 </td>
   <td style="text-align:right;"> 0.0002827 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK11B </td>
   <td style="text-align:right;"> 0.6027014 </td>
   <td style="text-align:right;"> 0.1003271 </td>
   <td style="text-align:right;"> 0.5807651 </td>
   <td style="text-align:right;"> -0.4460424 </td>
   <td style="text-align:right;"> 0.7724371 </td>
   <td style="text-align:right;"> -0.2130050 </td>
   <td style="text-align:right;"> 0.0000367 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK11A </td>
   <td style="text-align:right;"> 0.6918667 </td>
   <td style="text-align:right;"> -0.0418740 </td>
   <td style="text-align:right;"> 0.1530063 </td>
   <td style="text-align:right;"> -0.9900331 </td>
   <td style="text-align:right;"> 0.7322619 </td>
   <td style="text-align:right;"> -0.0752472 </td>
   <td style="text-align:right;"> 0.0002515 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NADK </td>
   <td style="text-align:right;"> 0.9304212 </td>
   <td style="text-align:right;"> 0.5667858 </td>
   <td style="text-align:right;"> 0.6113672 </td>
   <td style="text-align:right;"> -0.2316442 </td>
   <td style="text-align:right;"> 0.0363555 </td>
   <td style="text-align:right;"> -0.2690529 </td>
   <td style="text-align:right;"> -0.0001733 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GNB1 </td>
   <td style="text-align:right;"> 0.3846527 </td>
   <td style="text-align:right;"> 0.7198024 </td>
   <td style="text-align:right;"> 0.5191101 </td>
   <td style="text-align:right;"> -0.4141845 </td>
   <td style="text-align:right;"> 1.2323367 </td>
   <td style="text-align:right;"> 0.3057632 </td>
   <td style="text-align:right;"> 0.0002222 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SKI </td>
   <td style="text-align:right;"> 0.5390208 </td>
   <td style="text-align:right;"> 0.2955985 </td>
   <td style="text-align:right;"> 0.4723097 </td>
   <td style="text-align:right;"> -0.2405039 </td>
   <td style="text-align:right;"> -0.1385236 </td>
   <td style="text-align:right;"> -0.2108464 </td>
   <td style="text-align:right;"> 0.0000390 </td>
  </tr>
</tbody>
</table>

```r
# Test 'Normal' - 'Colorectal.Cancer'
contr <- makeContrasts(proper.groupNormal - proper.groupColorectal.Cancer, levels = colnames(coef(fit)))
contr %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal - proper.groupColorectal.Cancer </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> proper.groupColorectal.Cancer </td>
   <td style="text-align:right;"> -1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupNormal </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupPolyp </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S.Score </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G2M.Score </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> percent_MT </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nFeature_RNA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

```r
fit2 <- contrasts.fit(fit, contrasts = contr)
fit2 <- eBayes(fit2)
out <- topTable(fit2, n = Inf, sort.by = "P")
head(out, 30) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> logFC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AveExpr </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> t </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> P.Value </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> adj.P.Val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> B </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> XIST </td>
   <td style="text-align:right;"> 2.1117828 </td>
   <td style="text-align:right;"> 0.8617130 </td>
   <td style="text-align:right;"> 12.386400 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 45.063615 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC26A2 </td>
   <td style="text-align:right;"> 2.2815414 </td>
   <td style="text-align:right;"> 1.2586768 </td>
   <td style="text-align:right;"> 10.464502 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 33.716835 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PIGR </td>
   <td style="text-align:right;"> 2.0629994 </td>
   <td style="text-align:right;"> 1.3021222 </td>
   <td style="text-align:right;"> 8.735935 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 23.734924 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MMP12 </td>
   <td style="text-align:right;"> -1.6891241 </td>
   <td style="text-align:right;"> 0.4547259 </td>
   <td style="text-align:right;"> -8.639084 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 23.189319 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCND3 </td>
   <td style="text-align:right;"> 2.2031489 </td>
   <td style="text-align:right;"> 1.2335298 </td>
   <td style="text-align:right;"> 8.270839 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 21.133169 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC26A3 </td>
   <td style="text-align:right;"> 1.6256208 </td>
   <td style="text-align:right;"> 0.7653037 </td>
   <td style="text-align:right;"> 8.186660 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 20.667525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TYMP </td>
   <td style="text-align:right;"> -1.6659571 </td>
   <td style="text-align:right;"> 0.6145426 </td>
   <td style="text-align:right;"> -7.998999 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 19.635834 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PHGR1 </td>
   <td style="text-align:right;"> 1.6337462 </td>
   <td style="text-align:right;"> 0.8138926 </td>
   <td style="text-align:right;"> 7.780492 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 18.446331 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDE3A </td>
   <td style="text-align:right;"> 1.3916984 </td>
   <td style="text-align:right;"> 0.6038827 </td>
   <td style="text-align:right;"> 7.726625 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 18.155144 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CKB </td>
   <td style="text-align:right;"> 1.8743347 </td>
   <td style="text-align:right;"> 1.2263496 </td>
   <td style="text-align:right;"> 7.667410 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 17.836017 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GUCA2A </td>
   <td style="text-align:right;"> 1.3650407 </td>
   <td style="text-align:right;"> 0.5487307 </td>
   <td style="text-align:right;"> 7.360595 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 16.199629 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UTY </td>
   <td style="text-align:right;"> -1.3039318 </td>
   <td style="text-align:right;"> 0.6730543 </td>
   <td style="text-align:right;"> -6.327526 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 8.00e-07 </td>
   <td style="text-align:right;"> 10.934230 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FKBP5 </td>
   <td style="text-align:right;"> 1.8112118 </td>
   <td style="text-align:right;"> 1.3595449 </td>
   <td style="text-align:right;"> 6.315513 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 8.00e-07 </td>
   <td style="text-align:right;"> 10.875522 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MUC12 </td>
   <td style="text-align:right;"> 1.2426659 </td>
   <td style="text-align:right;"> 0.7475248 </td>
   <td style="text-align:right;"> 6.224798 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 1.10e-06 </td>
   <td style="text-align:right;"> 10.434283 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PARP14 </td>
   <td style="text-align:right;"> -1.5799500 </td>
   <td style="text-align:right;"> 1.0758847 </td>
   <td style="text-align:right;"> -6.217622 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 1.10e-06 </td>
   <td style="text-align:right;"> 10.399535 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZG16 </td>
   <td style="text-align:right;"> 1.0408636 </td>
   <td style="text-align:right;"> 0.4193840 </td>
   <td style="text-align:right;"> 6.109961 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 1.70e-06 </td>
   <td style="text-align:right;"> 9.881052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TNFAIP2 </td>
   <td style="text-align:right;"> -1.3112046 </td>
   <td style="text-align:right;"> 0.6137222 </td>
   <td style="text-align:right;"> -6.103989 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 1.70e-06 </td>
   <td style="text-align:right;"> 9.852445 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNF213 </td>
   <td style="text-align:right;"> -1.5300930 </td>
   <td style="text-align:right;"> 1.5047272 </td>
   <td style="text-align:right;"> -6.101504 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 1.70e-06 </td>
   <td style="text-align:right;"> 9.840552 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLA2G7 </td>
   <td style="text-align:right;"> -1.0776176 </td>
   <td style="text-align:right;"> 0.3860114 </td>
   <td style="text-align:right;"> -5.958340 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 3.20e-06 </td>
   <td style="text-align:right;"> 9.160027 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MBNL1 </td>
   <td style="text-align:right;"> 1.6166508 </td>
   <td style="text-align:right;"> 2.4884494 </td>
   <td style="text-align:right;"> 5.802791 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 6.60e-06 </td>
   <td style="text-align:right;"> 8.431854 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MT-CO2 </td>
   <td style="text-align:right;"> -1.1954459 </td>
   <td style="text-align:right;"> 2.3583166 </td>
   <td style="text-align:right;"> -5.745928 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 8.30e-06 </td>
   <td style="text-align:right;"> 8.168660 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SON </td>
   <td style="text-align:right;"> -1.4467357 </td>
   <td style="text-align:right;"> 1.3057203 </td>
   <td style="text-align:right;"> -5.625887 </td>
   <td style="text-align:right;"> 1e-07 </td>
   <td style="text-align:right;"> 1.41e-05 </td>
   <td style="text-align:right;"> 7.618453 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ARHGAP15 </td>
   <td style="text-align:right;"> 1.5448224 </td>
   <td style="text-align:right;"> 2.0602791 </td>
   <td style="text-align:right;"> 5.590987 </td>
   <td style="text-align:right;"> 1e-07 </td>
   <td style="text-align:right;"> 1.59e-05 </td>
   <td style="text-align:right;"> 7.459889 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC20A1 </td>
   <td style="text-align:right;"> -0.9744211 </td>
   <td style="text-align:right;"> 0.3808700 </td>
   <td style="text-align:right;"> -5.502835 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 2.31e-05 </td>
   <td style="text-align:right;"> 7.062245 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAML2 </td>
   <td style="text-align:right;"> 1.5276346 </td>
   <td style="text-align:right;"> 1.5398544 </td>
   <td style="text-align:right;"> 5.431854 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 3.05e-05 </td>
   <td style="text-align:right;"> 6.745075 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCDC88A </td>
   <td style="text-align:right;"> -1.2557448 </td>
   <td style="text-align:right;"> 0.8337547 </td>
   <td style="text-align:right;"> -5.422734 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 3.05e-05 </td>
   <td style="text-align:right;"> 6.704523 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LINC00996 </td>
   <td style="text-align:right;"> -1.1693561 </td>
   <td style="text-align:right;"> 0.5466776 </td>
   <td style="text-align:right;"> -5.419324 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 3.05e-05 </td>
   <td style="text-align:right;"> 6.689372 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAMD9L </td>
   <td style="text-align:right;"> -1.0672631 </td>
   <td style="text-align:right;"> 0.5365975 </td>
   <td style="text-align:right;"> -5.410692 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 3.06e-05 </td>
   <td style="text-align:right;"> 6.651046 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP1B1 </td>
   <td style="text-align:right;"> -1.1246590 </td>
   <td style="text-align:right;"> 0.5904788 </td>
   <td style="text-align:right;"> -5.357671 </td>
   <td style="text-align:right;"> 3e-07 </td>
   <td style="text-align:right;"> 3.78e-05 </td>
   <td style="text-align:right;"> 6.416535 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAT1 </td>
   <td style="text-align:right;"> -1.7454575 </td>
   <td style="text-align:right;"> 1.5955237 </td>
   <td style="text-align:right;"> -5.334894 </td>
   <td style="text-align:right;"> 3e-07 </td>
   <td style="text-align:right;"> 4.06e-05 </td>
   <td style="text-align:right;"> 6.316266 </td>
  </tr>
</tbody>
</table>

**Output columns:**

* logFC: log fold change (since we are working with Seurat's natural log transformed data, will be natural log fold change)
* AveExpr: Average expression across all cells in expr2
* t: logFC divided by its standard error
* P.Value: Raw p-value (based on t) from test that logFC differs from 0
* adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
* B: log-odds that gene is DE 

## Save files

```r
write.csv(GenTable(GOdata, Fisher = resultFisher), file = "cluster10_GOdata.csv")
write.csv(out, file = "cluster10_Normal-Colorectal.Cancer_topTable.csv")
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
 [1] stats4    stats     graphics  grDevices utils     datasets  methods  
 [8] base     
 
 other attached packages:
  [1] org.Hs.eg.db_3.17.0  kableExtra_1.3.4     dplyr_1.1.2         
  [4] topGO_2.52.0         SparseM_1.81         GO.db_3.17.0        
  [7] AnnotationDbi_1.62.2 IRanges_2.34.1       S4Vectors_0.38.1    
 [10] Biobase_2.60.0       graph_1.78.0         BiocGenerics_0.46.0 
 [13] limma_3.56.2         SeuratObject_4.1.3   Seurat_4.3.0.1      
 
 loaded via a namespace (and not attached):
   [1] RColorBrewer_1.1-3      rstudioapi_0.15.0       jsonlite_1.8.7         
   [4] magrittr_2.0.3          spatstat.utils_3.0-3    rmarkdown_2.23         
   [7] zlibbioc_1.46.0         vctrs_0.6.3             ROCR_1.0-11            
  [10] memoise_2.0.1           spatstat.explore_3.2-1  RCurl_1.98-1.12        
  [13] webshot_0.5.5           htmltools_0.5.5         sass_0.4.7             
  [16] sctransform_0.3.5       parallelly_1.36.0       KernSmooth_2.23-22     
  [19] bslib_0.5.0             htmlwidgets_1.6.2       ica_1.0-3              
  [22] plyr_1.8.8              plotly_4.10.2           zoo_1.8-12             
  [25] cachem_1.0.8            igraph_1.5.0            mime_0.12              
  [28] lifecycle_1.0.3         pkgconfig_2.0.3         Matrix_1.6-0           
  [31] R6_2.5.1                fastmap_1.1.1           GenomeInfoDbData_1.2.10
  [34] fitdistrplus_1.1-11     future_1.33.0           shiny_1.7.4.1          
  [37] digest_0.6.33           colorspace_2.1-0        patchwork_1.1.2        
  [40] tensor_1.5              irlba_2.3.5.1           RSQLite_2.3.1          
  [43] progressr_0.13.0        fansi_1.0.4             spatstat.sparse_3.0-2  
  [46] httr_1.4.6              polyclip_1.10-4         abind_1.4-5            
  [49] compiler_4.3.1          bit64_4.0.5             DBI_1.1.3              
  [52] highr_0.10              MASS_7.3-60             tools_4.3.1            
  [55] lmtest_0.9-40           httpuv_1.6.11           future.apply_1.11.0    
  [58] goftest_1.2-3           glue_1.6.2              nlme_3.1-162           
  [61] promises_1.2.0.1        grid_4.3.1              Rtsne_0.16             
  [64] cluster_2.1.4           reshape2_1.4.4          generics_0.1.3         
  [67] gtable_0.3.3            spatstat.data_3.0-1     tidyr_1.3.0            
  [70] data.table_1.14.8       xml2_1.3.5              XVector_0.40.0         
  [73] sp_2.0-0                utf8_1.2.3              spatstat.geom_3.2-4    
  [76] RcppAnnoy_0.0.21        ggrepel_0.9.3           RANN_2.6.1             
  [79] pillar_1.9.0            stringr_1.5.0           later_1.3.1            
  [82] splines_4.3.1           lattice_0.21-8          bit_4.0.5              
  [85] survival_3.5-5          deldir_1.0-9            tidyselect_1.2.0       
  [88] Biostrings_2.68.1       miniUI_0.1.1.1          pbapply_1.7-2          
  [91] knitr_1.43              gridExtra_2.3           svglite_2.1.1          
  [94] scattermore_1.2         xfun_0.39               matrixStats_1.0.0      
  [97] stringi_1.7.12          lazyeval_0.2.2          yaml_2.3.7             
 [100] evaluate_0.21           codetools_0.2-19        tibble_3.2.1           
 [103] cli_3.6.1               uwot_0.1.16             systemfonts_1.0.4      
 [106] xtable_1.8-4            reticulate_1.30         munsell_0.5.0          
 [109] jquerylib_0.1.4         GenomeInfoDb_1.36.1     Rcpp_1.0.11            
 [112] globals_0.16.2          spatstat.random_3.1-5   png_0.1-8              
 [115] parallel_4.3.1          ellipsis_0.3.2          blob_1.2.4             
 [118] ggplot2_3.4.2           bitops_1.0-7            listenv_0.9.0          
 [121] viridisLite_0.4.2       scales_1.2.1            ggridges_0.5.4         
 [124] crayon_1.5.2            leiden_0.4.3            purrr_1.0.1            
 [127] rlang_1.1.1             rvest_1.0.3             KEGGREST_1.40.0        
 [130] cowplot_1.1.1
</div>