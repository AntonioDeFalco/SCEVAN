# CNAvegaMC

## Introduction

## Installation
```
library(devtools)
install_github("AntonioDeFalco/CNAvegaMC")
library(CNAvegaMC)
```

## Usage

```
## GET TUMOR CELLS AND SUBCLONES
results <- pipelineCNA(count_mtx, sample = sample, par_cores = 20, gr_truth = gr_truth, SUBCLONES = TRUE)
```
Heatmap classification of tumor cells

![image](https://github.com/AntonioDeFalco/CNAvegaMC/blob/main/vignettes/images/MGH125heatmap.jpeg)

Heatmap of tumour cell subclones

![image](https://github.com/AntonioDeFalco/CNAvegaMC/blob/main/vignettes/images/MGH125heatmap_subclones.jpeg)

## Use Cases

* [Example (TO DO)](https://github.com/AntonioDeFalco/CNAvegaMC)