# Single   CEll   Variational   Aneuploidy aNalysis  (SCEVAN)

## Introduction

SCEVAN is an R tool that starting from the raw count matrix of scRNA data automatically classify malignant cells and non-malignant cells of tumor microenviroment present in the biopsy and characterize the clonal structure of tumour cells, specifically analysing the presence of subpopulations and characterising the specific and shared alterations of each subpopulation. The aim of the tool is to automate the entire analysis by allowing it to be performed in a very simple and completely unsupervised way. Analyses carried out on more than 100 samples and 93,000 cells show better classification with an F1 score for all samples of 0.90 compared to 0.63 obtained with the state-of-the-art copyKAT, with a better classification in 63% of the samples and worse in only 23% of the samples.

## Installation
```
library(devtools)
install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
```

## Usage

A single call (pipelineCNA) allows the execution of the entire analysis.

- ***count_mtx*** : Count matrix with genes on rows (both Gene Symbol or Ensembl ID are allowed) and cells on columns.
- ***sample*** : Sample name to save results (optional)
- ***par_cores*** : Number of cores to run the pipeline  (optional)
- ***norm_cells*** : vectors of normal cells if the classification is already known and you are only interested in the clonal structure (optional)
- ***SUBCLONES*** : Boolean value TRUE if you are interested in analysing the clonal structure and FALSE if you are only interested in the classification of malignant and non-malignant cells (optional)

```
## GET TUMOR CELLS AND SUBCLONES
results <- pipelineCNA(count_mtx)
```

## Usage examples (vignettes)

- [Intratumoral heterogeneity in Glioblastoma](http://htmlpreview.github.io/?https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/IntratumoralHeterogeneityInGlioblastoma.html)
- [](http://htmlpreview.github.io/)

## Output

The pipeline returns in ***results*** a data frame containing for each cell its classification and which clonal subpopulation it belongs to, and in the output folder the following plots:

### Heatmap classification of tumor cells

Heatmap of the Copy Number Alteration matrix with classification of non-malignant and malignant cells.

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106heatmap.png)

### Heatmap of tumour cell subclones

Heatmap of the Copy Number Alteration matrix with the clonal subpopulations found.

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106heatmap_subclones.png)

### Consensus plot

Compact plot of the alterations present in each subpopulation.

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106consensus.png)

### OncoPrint-like plot

OncoPrint-like plot that highlighting specific alteration, shared alteration between subclones, or clonal alteration.

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106OncoHeat.png)

### DE analysis in specific alterations

Vulcano plot obtained from differential expression analysis of the genes belonging to the specific alterations found.

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106-DEchr3_subclones.png)

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106-DEchr12_subclones.png)

### Pathway Analysis of subclones

REACTOME pathways activity obtained with GSEA for each subclone in contrast to the others for each subclone.

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106pathwayAnalysis_subclones1.png)
![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106pathwayAnalysis_subclones2.png)
![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106pathwayAnalysis_subclones3.png)
