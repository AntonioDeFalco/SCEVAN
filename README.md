# Single   CEll   Variational   Aneuploidy aNalysis  (SCEVAN)

## Introduction

## Installation
```
library(devtools)
install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
```

## Usage
- count_mtx : count matrix with genes on rows (both Gene Symbol or Ensembl ID are allowed) and cells on columns.

```
## GET TUMOR CELLS AND SUBCLONES
results <- pipelineCNA(count_mtx, sample = sample, par_cores = 20, gr_truth = gr_truth, SUBCLONES = TRUE)
```

Heatmap classification of tumor cells

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106heatmap.png)

Heatmap of tumour cell subclones

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106heatmap_subclones.png)

Consensus plot

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106consensus.png)

OncoPrint-like plot

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106OncoHeat.png)

DE analysis in specific alterations

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106-DEchr3_subclones.png)

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106-DEchr12_subclones.png)

Pathway Analysis of subclones

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106pathwayAnalysis_subclones1.png)
![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106pathwayAnalysis_subclones2.png)
![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/output/MGH106pathwayAnalysis_subclones3.png)
