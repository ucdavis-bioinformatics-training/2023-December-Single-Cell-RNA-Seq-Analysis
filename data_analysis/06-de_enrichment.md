---
title: "Introduction to Single Cell RNA-Seq Part 6: Enrichment and model-based differential expression"
author: "Bioinformatics Core"
date: "2023-12-11"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 6: Enrichment and model-based differential expression


## Set up workspace

```r
library(Seurat)
library(limma)
library(topGO)
library(dplyr)
library(kableExtra)
set.seed(12345)
experiment.aggregate <- readRDS("scRNA_workshop-05.rds")
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
   <td style="text-align:left;"> GO:0010613 </td>
   <td style="text-align:left;"> positive regulation of cardiac muscle hypertrophy </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.34 </td>
   <td style="text-align:left;"> 0.00031 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006325 </td>
   <td style="text-align:left;"> chromatin organization </td>
   <td style="text-align:right;"> 465 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 8.31 </td>
   <td style="text-align:left;"> 0.00047 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0030036 </td>
   <td style="text-align:left;"> actin cytoskeleton organization </td>
   <td style="text-align:right;"> 462 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 8.26 </td>
   <td style="text-align:left;"> 0.00057 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0034968 </td>
   <td style="text-align:left;"> histone lysine methylation </td>
   <td style="text-align:right;"> 81 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1.45 </td>
   <td style="text-align:left;"> 0.00059 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:2001014 </td>
   <td style="text-align:left;"> regulation of skeletal muscle cell differentiation </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.20 </td>
   <td style="text-align:left;"> 0.00083 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0000381 </td>
   <td style="text-align:left;"> regulation of alternative mRNA splicing, via spliceosome </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.75 </td>
   <td style="text-align:left;"> 0.00086 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006338 </td>
   <td style="text-align:left;"> chromatin remodeling </td>
   <td style="text-align:right;"> 270 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 4.82 </td>
   <td style="text-align:left;"> 0.00108 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0048813 </td>
   <td style="text-align:left;"> dendrite morphogenesis </td>
   <td style="text-align:right;"> 95 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1.70 </td>
   <td style="text-align:left;"> 0.00152 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045445 </td>
   <td style="text-align:left;"> myoblast differentiation </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1.25 </td>
   <td style="text-align:left;"> 0.00152 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1901888 </td>
   <td style="text-align:left;"> regulation of cell junction assembly </td>
   <td style="text-align:right;"> 124 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 2.22 </td>
   <td style="text-align:left;"> 0.00167 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1901897 </td>
   <td style="text-align:left;"> regulation of relaxation of cardiac muscle </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:left;"> 0.00186 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0031666 </td>
   <td style="text-align:left;"> positive regulation of lipopolysaccharide-mediated signaling... </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:left;"> 0.00186 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0071320 </td>
   <td style="text-align:left;"> cellular response to cAMP </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.54 </td>
   <td style="text-align:left;"> 0.00188 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0061049 </td>
   <td style="text-align:left;"> cell growth involved in cardiac muscle cell development </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.27 </td>
   <td style="text-align:left;"> 0.00218 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0010595 </td>
   <td style="text-align:left;"> positive regulation of endothelial cell migration </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1.36 </td>
   <td style="text-align:left;"> 0.00232 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0002063 </td>
   <td style="text-align:left;"> chondrocyte development </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.29 </td>
   <td style="text-align:left;"> 0.00264 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045906 </td>
   <td style="text-align:left;"> negative regulation of vasoconstriction </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:left;"> 0.00306 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0098909 </td>
   <td style="text-align:left;"> regulation of cardiac muscle cell action potential involved ... </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:left;"> 0.00306 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0048842 </td>
   <td style="text-align:left;"> positive regulation of axon extension involved in axon guida... </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:left;"> 0.00306 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0021785 </td>
   <td style="text-align:left;"> branchiomotor neuron axon guidance </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:left;"> 0.00306 </td>
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
   <td style="text-align:right;"> 0.0347548 </td>
   <td style="text-align:right;"> 0.7909218 </td>
   <td style="text-align:right;"> 0.5656109 </td>
   <td style="text-align:right;"> 1187 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACTGTCCAGGCGTTAG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0513069 </td>
   <td style="text-align:right;"> -0.0204757 </td>
   <td style="text-align:right;"> 0.5043712 </td>
   <td style="text-align:right;"> 1634 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AGCGTATAGCCTCAGC_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0576611 </td>
   <td style="text-align:right;"> -0.0549822 </td>
   <td style="text-align:right;"> 3.0966767 </td>
   <td style="text-align:right;"> 880 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AGCGTATAGGTATCTC_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.1080488 </td>
   <td style="text-align:right;"> -0.0416714 </td>
   <td style="text-align:right;"> 0.4595588 </td>
   <td style="text-align:right;"> 783 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AGGACTTGTACCTATG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0181025 </td>
   <td style="text-align:right;"> -0.0642517 </td>
   <td style="text-align:right;"> 0.4804393 </td>
   <td style="text-align:right;"> 986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AGGGTTTGTTCTTCAT_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0057069 </td>
   <td style="text-align:right;"> -0.0475972 </td>
   <td style="text-align:right;"> 1.3483146 </td>
   <td style="text-align:right;"> 937 </td>
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
   <td style="text-align:left;"> TCCACCAAGTAACGTA_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0872095 </td>
   <td style="text-align:right;"> -0.1505728 </td>
   <td style="text-align:right;"> 0.1785183 </td>
   <td style="text-align:right;"> 2096 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCCACCAGTGTCCACG_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0697658 </td>
   <td style="text-align:right;"> -0.1145196 </td>
   <td style="text-align:right;"> 0.3147954 </td>
   <td style="text-align:right;"> 1147 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGCTTGAGACTTGTC_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0989612 </td>
   <td style="text-align:right;"> -0.1146185 </td>
   <td style="text-align:right;"> 0.1570475 </td>
   <td style="text-align:right;"> 1471 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCTCCGACACAGTACT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0437678 </td>
   <td style="text-align:right;"> -0.1455716 </td>
   <td style="text-align:right;"> 0.1921230 </td>
   <td style="text-align:right;"> 2383 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TGATGGTCATAACTCG_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.1108661 </td>
   <td style="text-align:right;"> -0.1465694 </td>
   <td style="text-align:right;"> 0.5749668 </td>
   <td style="text-align:right;"> 1413 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTCTTCCCAGCGTAGA_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0361969 </td>
   <td style="text-align:right;"> 0.0067655 </td>
   <td style="text-align:right;"> 0.5504587 </td>
   <td style="text-align:right;"> 843 </td>
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
   <td style="text-align:right;"> 0.0579580 </td>
   <td style="text-align:right;"> 0.4644867 </td>
   <td style="text-align:right;"> 0.5707252 </td>
   <td style="text-align:right;"> -0.5229024 </td>
   <td style="text-align:right;"> 0.7195925 </td>
   <td style="text-align:right;"> 0.2619560 </td>
   <td style="text-align:right;"> 0.0001080 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK11B </td>
   <td style="text-align:right;"> 0.5915809 </td>
   <td style="text-align:right;"> 0.1722522 </td>
   <td style="text-align:right;"> 0.5923884 </td>
   <td style="text-align:right;"> -0.6623504 </td>
   <td style="text-align:right;"> 0.7849839 </td>
   <td style="text-align:right;"> -0.3260968 </td>
   <td style="text-align:right;"> 0.0000319 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK11A </td>
   <td style="text-align:right;"> 0.9175843 </td>
   <td style="text-align:right;"> 0.2172913 </td>
   <td style="text-align:right;"> 0.4343150 </td>
   <td style="text-align:right;"> -1.0275000 </td>
   <td style="text-align:right;"> 0.9939193 </td>
   <td style="text-align:right;"> -0.0600846 </td>
   <td style="text-align:right;"> 0.0000147 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NADK </td>
   <td style="text-align:right;"> 0.3736396 </td>
   <td style="text-align:right;"> 0.2646946 </td>
   <td style="text-align:right;"> 0.0662475 </td>
   <td style="text-align:right;"> -0.2493852 </td>
   <td style="text-align:right;"> 0.1639109 </td>
   <td style="text-align:right;"> 0.0637775 </td>
   <td style="text-align:right;"> -0.0000001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GNB1 </td>
   <td style="text-align:right;"> 0.7772654 </td>
   <td style="text-align:right;"> 0.5788333 </td>
   <td style="text-align:right;"> 0.8958150 </td>
   <td style="text-align:right;"> -0.5519191 </td>
   <td style="text-align:right;"> 0.9004079 </td>
   <td style="text-align:right;"> -0.0616478 </td>
   <td style="text-align:right;"> 0.0002935 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SKI </td>
   <td style="text-align:right;"> 0.7171411 </td>
   <td style="text-align:right;"> 0.1954765 </td>
   <td style="text-align:right;"> 0.4870630 </td>
   <td style="text-align:right;"> -0.7409301 </td>
   <td style="text-align:right;"> -0.1685434 </td>
   <td style="text-align:right;"> -0.2098351 </td>
   <td style="text-align:right;"> 0.0000133 </td>
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
   <td style="text-align:right;"> 1.9827318 </td>
   <td style="text-align:right;"> 0.6682736 </td>
   <td style="text-align:right;"> 11.870704 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 37.251244 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GUCA2A </td>
   <td style="text-align:right;"> 1.6727802 </td>
   <td style="text-align:right;"> 0.4946151 </td>
   <td style="text-align:right;"> 10.607512 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 31.021088 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MMP12 </td>
   <td style="text-align:right;"> -2.0976536 </td>
   <td style="text-align:right;"> 0.6819202 </td>
   <td style="text-align:right;"> -8.045790 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 18.347651 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PHGR1 </td>
   <td style="text-align:right;"> 1.7128625 </td>
   <td style="text-align:right;"> 0.6230695 </td>
   <td style="text-align:right;"> 7.990444 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 18.079041 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CKB </td>
   <td style="text-align:right;"> 2.0697765 </td>
   <td style="text-align:right;"> 1.1605407 </td>
   <td style="text-align:right;"> 7.813288 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 17.222231 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC26A2 </td>
   <td style="text-align:right;"> 1.8826649 </td>
   <td style="text-align:right;"> 0.9641306 </td>
   <td style="text-align:right;"> 6.787927 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000004 </td>
   <td style="text-align:right;"> 12.377261 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC26A3 </td>
   <td style="text-align:right;"> 1.5570563 </td>
   <td style="text-align:right;"> 0.6022747 </td>
   <td style="text-align:right;"> 6.661061 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000006 </td>
   <td style="text-align:right;"> 11.794657 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PIGR </td>
   <td style="text-align:right;"> 1.7503636 </td>
   <td style="text-align:right;"> 1.1102687 </td>
   <td style="text-align:right;"> 6.460438 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000015 </td>
   <td style="text-align:right;"> 10.882534 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HIVEP3 </td>
   <td style="text-align:right;"> -1.4376532 </td>
   <td style="text-align:right;"> 0.7340200 </td>
   <td style="text-align:right;"> -5.994125 </td>
   <td style="text-align:right;"> 0.0e+00 </td>
   <td style="text-align:right;"> 0.0000121 </td>
   <td style="text-align:right;"> 8.810834 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NR3C2 </td>
   <td style="text-align:right;"> 1.6376297 </td>
   <td style="text-align:right;"> 1.0291416 </td>
   <td style="text-align:right;"> 5.792722 </td>
   <td style="text-align:right;"> 1.0e-07 </td>
   <td style="text-align:right;"> 0.0000276 </td>
   <td style="text-align:right;"> 7.939339 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UTY </td>
   <td style="text-align:right;"> -1.3631601 </td>
   <td style="text-align:right;"> 0.7485786 </td>
   <td style="text-align:right;"> -5.484691 </td>
   <td style="text-align:right;"> 3.0e-07 </td>
   <td style="text-align:right;"> 0.0000928 </td>
   <td style="text-align:right;"> 6.636837 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZG16 </td>
   <td style="text-align:right;"> 0.9778633 </td>
   <td style="text-align:right;"> 0.3534387 </td>
   <td style="text-align:right;"> 5.484091 </td>
   <td style="text-align:right;"> 3.0e-07 </td>
   <td style="text-align:right;"> 0.0000928 </td>
   <td style="text-align:right;"> 6.634339 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MUC13 </td>
   <td style="text-align:right;"> 1.2139028 </td>
   <td style="text-align:right;"> 0.5911549 </td>
   <td style="text-align:right;"> 5.163231 </td>
   <td style="text-align:right;"> 1.1e-06 </td>
   <td style="text-align:right;"> 0.0003496 </td>
   <td style="text-align:right;"> 5.320392 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDE4D </td>
   <td style="text-align:right;"> 1.6637089 </td>
   <td style="text-align:right;"> 1.5119900 </td>
   <td style="text-align:right;"> 5.120095 </td>
   <td style="text-align:right;"> 1.3e-06 </td>
   <td style="text-align:right;"> 0.0003901 </td>
   <td style="text-align:right;"> 5.147303 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TNFAIP2 </td>
   <td style="text-align:right;"> -1.4572991 </td>
   <td style="text-align:right;"> 1.0361014 </td>
   <td style="text-align:right;"> -5.104445 </td>
   <td style="text-align:right;"> 1.4e-06 </td>
   <td style="text-align:right;"> 0.0003901 </td>
   <td style="text-align:right;"> 5.084722 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HSP90AA1 </td>
   <td style="text-align:right;"> -1.3268169 </td>
   <td style="text-align:right;"> 0.7386294 </td>
   <td style="text-align:right;"> -5.041036 </td>
   <td style="text-align:right;"> 1.8e-06 </td>
   <td style="text-align:right;"> 0.0004794 </td>
   <td style="text-align:right;"> 4.832355 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MACF1 </td>
   <td style="text-align:right;"> -1.2728166 </td>
   <td style="text-align:right;"> 2.5152218 </td>
   <td style="text-align:right;"> -5.024663 </td>
   <td style="text-align:right;"> 1.9e-06 </td>
   <td style="text-align:right;"> 0.0004838 </td>
   <td style="text-align:right;"> 4.767506 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SATB2 </td>
   <td style="text-align:right;"> 1.2335513 </td>
   <td style="text-align:right;"> 0.6184811 </td>
   <td style="text-align:right;"> 4.931185 </td>
   <td style="text-align:right;"> 2.9e-06 </td>
   <td style="text-align:right;"> 0.0006781 </td>
   <td style="text-align:right;"> 4.399765 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNF213 </td>
   <td style="text-align:right;"> -1.4619869 </td>
   <td style="text-align:right;"> 1.4301922 </td>
   <td style="text-align:right;"> -4.899682 </td>
   <td style="text-align:right;"> 3.3e-06 </td>
   <td style="text-align:right;"> 0.0007080 </td>
   <td style="text-align:right;"> 4.276803 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PARP14 </td>
   <td style="text-align:right;"> -1.5705769 </td>
   <td style="text-align:right;"> 1.4065849 </td>
   <td style="text-align:right;"> -4.895788 </td>
   <td style="text-align:right;"> 3.3e-06 </td>
   <td style="text-align:right;"> 0.0007080 </td>
   <td style="text-align:right;"> 4.261639 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IL7R </td>
   <td style="text-align:right;"> -1.0858376 </td>
   <td style="text-align:right;"> 0.3242462 </td>
   <td style="text-align:right;"> -4.865053 </td>
   <td style="text-align:right;"> 3.8e-06 </td>
   <td style="text-align:right;"> 0.0007591 </td>
   <td style="text-align:right;"> 4.142219 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FOXP1 </td>
   <td style="text-align:right;"> 1.5111842 </td>
   <td style="text-align:right;"> 1.6488493 </td>
   <td style="text-align:right;"> 4.856249 </td>
   <td style="text-align:right;"> 3.9e-06 </td>
   <td style="text-align:right;"> 0.0007591 </td>
   <td style="text-align:right;"> 4.108095 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMN </td>
   <td style="text-align:right;"> 0.8531736 </td>
   <td style="text-align:right;"> 0.2741440 </td>
   <td style="text-align:right;"> 4.731821 </td>
   <td style="text-align:right;"> 6.6e-06 </td>
   <td style="text-align:right;"> 0.0011476 </td>
   <td style="text-align:right;"> 3.630123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EPCAM </td>
   <td style="text-align:right;"> 1.1109426 </td>
   <td style="text-align:right;"> 0.5789830 </td>
   <td style="text-align:right;"> 4.715942 </td>
   <td style="text-align:right;"> 7.0e-06 </td>
   <td style="text-align:right;"> 0.0011476 </td>
   <td style="text-align:right;"> 3.569703 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLA2G7 </td>
   <td style="text-align:right;"> -1.1812964 </td>
   <td style="text-align:right;"> 0.5448965 </td>
   <td style="text-align:right;"> -4.705319 </td>
   <td style="text-align:right;"> 7.3e-06 </td>
   <td style="text-align:right;"> 0.0011476 </td>
   <td style="text-align:right;"> 3.529356 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZBTB20 </td>
   <td style="text-align:right;"> 1.5396424 </td>
   <td style="text-align:right;"> 1.7559376 </td>
   <td style="text-align:right;"> 4.700057 </td>
   <td style="text-align:right;"> 7.5e-06 </td>
   <td style="text-align:right;"> 0.0011476 </td>
   <td style="text-align:right;"> 3.509395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RAPGEF1 </td>
   <td style="text-align:right;"> -1.5180128 </td>
   <td style="text-align:right;"> 1.5356339 </td>
   <td style="text-align:right;"> -4.698717 </td>
   <td style="text-align:right;"> 7.5e-06 </td>
   <td style="text-align:right;"> 0.0011476 </td>
   <td style="text-align:right;"> 3.504313 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAT1 </td>
   <td style="text-align:right;"> -1.8924435 </td>
   <td style="text-align:right;"> 2.0416164 </td>
   <td style="text-align:right;"> -4.696409 </td>
   <td style="text-align:right;"> 7.6e-06 </td>
   <td style="text-align:right;"> 0.0011476 </td>
   <td style="text-align:right;"> 3.495562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TYMP </td>
   <td style="text-align:right;"> -1.3860201 </td>
   <td style="text-align:right;"> 0.9203430 </td>
   <td style="text-align:right;"> -4.688993 </td>
   <td style="text-align:right;"> 7.8e-06 </td>
   <td style="text-align:right;"> 0.0011476 </td>
   <td style="text-align:right;"> 3.467469 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SCNN1B </td>
   <td style="text-align:right;"> 0.9079763 </td>
   <td style="text-align:right;"> 0.3004799 </td>
   <td style="text-align:right;"> 4.668610 </td>
   <td style="text-align:right;"> 8.5e-06 </td>
   <td style="text-align:right;"> 0.0011916 </td>
   <td style="text-align:right;"> 3.390400 </td>
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

## A note on pseudobulk DE

Pseudobulk differential expression uses count data summed across all cells in each sample (typically within each cell type or cluster).  Unlike cell-level DE, pseudobulk DE *requires biological replicates* so we won't perform it on this dataset.

Once counts are summed, pseudobulk data are analyzed like bulk RNASeq data.

Pseudobulk DE may result in better false discovery rate control than cell-level DE, as shown [here](https://www.nature.com/articles/s41467-021-25960-2).

The Seurat function `AggregateExpression()` can be used to sum counts as described [here](https://satijalab.org/seurat/articles/de_vignette).

A tutorial on using limma for bulk RNASeq is available [here](https://ucdavis-bioinformatics-training.github.io/2023-June-RNA-Seq-Analysis/data_analysis/DE_Analysis_mm_with_quizzes).

## Prepare for the next section

#### Download Rmd document

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/07-doublet_detection.Rmd", "07-doublet_detection.Rmd")
```

#### Session Information

```r
sessionInfo()
```

<div class='r_output'> R version 4.3.2 (2023-10-31 ucrt)
 Platform: x86_64-w64-mingw32/x64 (64-bit)
 Running under: Windows 10 x64 (build 19045)
 
 Matrix products: default
 
 
 locale:
 [1] LC_COLLATE=English_United States.utf8 
 [2] LC_CTYPE=English_United States.utf8   
 [3] LC_MONETARY=English_United States.utf8
 [4] LC_NUMERIC=C                          
 [5] LC_TIME=English_United States.utf8    
 
 time zone: America/Los_Angeles
 tzcode source: internal
 
 attached base packages:
 [1] stats4    stats     graphics  grDevices utils     datasets  methods  
 [8] base     
 
 other attached packages:
  [1] org.Hs.eg.db_3.17.0  kableExtra_1.3.4     dplyr_1.1.4         
  [4] topGO_2.52.0         SparseM_1.81         GO.db_3.17.0        
  [7] AnnotationDbi_1.62.2 IRanges_2.34.1       S4Vectors_0.38.2    
 [10] Biobase_2.60.0       graph_1.78.0         BiocGenerics_0.46.0 
 [13] limma_3.56.2         Seurat_5.0.1         SeuratObject_5.0.1  
 [16] sp_2.1-2            
 
 loaded via a namespace (and not attached):
   [1] RColorBrewer_1.1-3      rstudioapi_0.15.0       jsonlite_1.8.8         
   [4] magrittr_2.0.3          spatstat.utils_3.0-4    rmarkdown_2.25         
   [7] zlibbioc_1.46.0         vctrs_0.6.5             ROCR_1.0-11            
  [10] memoise_2.0.1           spatstat.explore_3.2-5  RCurl_1.98-1.13        
  [13] webshot_0.5.5           htmltools_0.5.7         sass_0.4.8             
  [16] sctransform_0.4.1       parallelly_1.36.0       KernSmooth_2.23-22     
  [19] bslib_0.6.1             htmlwidgets_1.6.4       ica_1.0-3              
  [22] plyr_1.8.9              plotly_4.10.3           zoo_1.8-12             
  [25] cachem_1.0.8            igraph_1.5.1            mime_0.12              
  [28] lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.6-4           
  [31] R6_2.5.1                fastmap_1.1.1           GenomeInfoDbData_1.2.10
  [34] fitdistrplus_1.1-11     future_1.33.0           shiny_1.8.0            
  [37] digest_0.6.33           colorspace_2.1-0        patchwork_1.1.3        
  [40] tensor_1.5              RSpectra_0.16-1         irlba_2.3.5.1          
  [43] RSQLite_2.3.4           progressr_0.14.0        fansi_1.0.6            
  [46] spatstat.sparse_3.0-3   httr_1.4.7              polyclip_1.10-6        
  [49] abind_1.4-5             compiler_4.3.2          bit64_4.0.5            
  [52] DBI_1.1.3               fastDummies_1.7.3       highr_0.10             
  [55] MASS_7.3-60             tools_4.3.2             lmtest_0.9-40          
  [58] httpuv_1.6.13           future.apply_1.11.0     goftest_1.2-3          
  [61] glue_1.6.2              nlme_3.1-164            promises_1.2.1         
  [64] grid_4.3.2              Rtsne_0.17              cluster_2.1.6          
  [67] reshape2_1.4.4          generics_0.1.3          gtable_0.3.4           
  [70] spatstat.data_3.0-3     tidyr_1.3.0             data.table_1.14.10     
  [73] xml2_1.3.6              XVector_0.40.0          utf8_1.2.4             
  [76] spatstat.geom_3.2-7     RcppAnnoy_0.0.21        ggrepel_0.9.4          
  [79] RANN_2.6.1              pillar_1.9.0            stringr_1.5.1          
  [82] spam_2.10-0             RcppHNSW_0.5.0          later_1.3.2            
  [85] splines_4.3.2           lattice_0.22-5          bit_4.0.5              
  [88] survival_3.5-7          deldir_2.0-2            tidyselect_1.2.0       
  [91] Biostrings_2.68.1       miniUI_0.1.1.1          pbapply_1.7-2          
  [94] knitr_1.45              gridExtra_2.3           svglite_2.1.3          
  [97] scattermore_1.2         xfun_0.41               matrixStats_1.1.0      
 [100] stringi_1.8.2           lazyeval_0.2.2          yaml_2.3.7             
 [103] evaluate_0.23           codetools_0.2-19        tibble_3.2.1           
 [106] cli_3.6.1               uwot_0.1.16             systemfonts_1.0.5      
 [109] xtable_1.8-4            reticulate_1.34.0       munsell_0.5.0          
 [112] jquerylib_0.1.4         GenomeInfoDb_1.36.4     Rcpp_1.0.11            
 [115] globals_0.16.2          spatstat.random_3.2-2   png_0.1-8              
 [118] parallel_4.3.2          ellipsis_0.3.2          blob_1.2.4             
 [121] ggplot2_3.4.4           dotCall64_1.1-1         bitops_1.0-7           
 [124] listenv_0.9.0           viridisLite_0.4.2       scales_1.3.0           
 [127] ggridges_0.5.4          crayon_1.5.2            leiden_0.4.3.1         
 [130] purrr_1.0.2             rlang_1.1.2             rvest_1.0.3            
 [133] KEGGREST_1.40.1         cowplot_1.1.1
</div>