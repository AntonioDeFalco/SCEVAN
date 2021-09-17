# Single   CEll   Variational   Aneuploidy aNalysis  (SCEVAN)

## Introduction

## Installation
```
library(devtools)
install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
```

## Usage

```
## GET TUMOR CELLS AND SUBCLONES
results <- pipelineCNA(count_mtx, sample = sample, par_cores = 20, gr_truth = gr_truth, SUBCLONES = TRUE)
```

Heatmap classification of tumor cells

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125heatmap.jpeg)

Heatmap of tumour cell subclones

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125heatmap_subclones.jpeg)

Fishplot subclones of tumor cells

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125fishplot_subclones.jpeg)

Consensus plot

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125consensus.png)

CN plot

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125plotCNline.png)

Tsne

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125tsne.jpeg)

DE analysis in specific alterations

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125-DE1.png)

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125-DE2.png)

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125-DE3.png)

Pathway Analysis of subclones

![image](https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/images/MGH125pathwayAnalysis_subclones.png)

## Use Cases

* [Example (MGH125)](https://htmlpreview.github.io/?https://raw.githubusercontent.com/AntonioDeFalco/SCEVAN/main/example.html?token=ACYAPINNYS2IVGAIMIUPRZDA5VWVU)