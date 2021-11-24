# Single   CEll   Variational   Aneuploidy aNalysis  (SCEVAN)

<p align="center">
  <img width="560" height="560" src="https://github.com/AntonioDeFalco/SCEVAN/blob/main/SCEVAN.png">
</p>

Preprint Link: [A fast variational algorithm to detect the clonal copy number substructure of tumors from single-cell data](https://www.biorxiv.org/content/10.1101/2021.11.20.469390v1)

## Introduction

SCEVAN is an R package that starting from the raw count matrix of scRNA data automatically classifies the cells present in the biopsy by segregating non-malignant cells of tumor microenviroment from the malignant cells and also characterizes the clonal structure of these malignant cells. It identfies cell subpopulations with different copy number architecture and reports g the specific and shared alterations of each subpopulation. The aim of the tool is to automate the entire analysis by allowing it to be performed in a very simple and completely unsupervised way. Analyses carried out on 106 samples and 93332 cells show better classification with an F1 score for all samples of 0.90 compared to 0.63 obtained with the state-of-the-art tools. It also explits a greedy multichannel segmentation algorithms making it particularly fast even for large datasets. 

## Installation

```
library(devtools)
install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
```

## Usage

### Single-sample analysis
A single call (pipelineCNA) allows the execution of the entire analysis of classification and characterization of clonal structure.

- ***count_mtx*** : Count matrix with genes on rows (both Gene Symbol or Ensembl ID are allowed) and cells on columns.
- ***sample*** : Sample name to save results (optional)
- ***par_cores*** : Number of cores to run the pipeline  (optional)
- ***norm_cells*** : vectors of normal cells if the classification is already known and you are only interested in the clonal structure (optional)
- ***SUBCLONES*** : Boolean value TRUE if you are interested in analysing the clonal structure and FALSE if you are only interested in the classification of malignant and non-malignant cells (optional)

```
results <- pipelineCNA(count_mtx)
```

### Multi-sample analysis
A single call (compareClonalStructure) allows the comparison of clonal profiles of the different samples.

- ***count_mtx1*** : Count matrix of sample 1.
- ***count_mtx2*** : Count matrix of sample 2.
- ***samp_1*** : Name of sample 1.
- ***samp_2*** : Name of sample 2.
- ***par_cores*** : Number of cores to run the pipeline  (optional)

```
compareClonalStructure(count_mtx1, count_mtx2, samp_1, samp_2)
```

## Usage examples (vignettes)

- [Intratumoral heterogeneity](http://htmlpreview.github.io/?https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/IntratumoralHeterogeneityInGlioblastoma.html)
- [Comparison of clonal profiles](http://htmlpreview.github.io/?https://github.com/AntonioDeFalco/SCEVAN/blob/main/vignettes/ComparisonOfClonalProfiles.html)

## Sample Datasets

We provide some pre-processed samples used in the examples (vignettes):

- ***MGH106.RData*** : scRNA data of MGH106 sample from the public dataset of Gliobastoma (GSE131928), you can download the pre-processed data from [here](https://www.dropbox.com/s/b9udpvhnc2ez9pc/MGH106_data.RData?dl=0)

- ***HNSCC26.RData*** : scRNA data of HNSCC26 Primary and HNSCC26 Lymph Node sample from the public dataset of Head&Neck cancer (GSE10332), you can download the pre-processed data from [here](https://www.dropbox.com/s/6zns12amobs39g8/HNSCC26_data.RData?dl=0)

## Citation

> 
>@article {De Falco2021.11.20.469390,\
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	author = {De Falco, Antonio and Caruso, Francesca P and Su, Xiao Dong and Iavarone, Antonio and Ceccarelli, Michele},\
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	title = {A fast variational algorithm to detect the clonal copy number substructure of tumors from single-cell data},\
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	elocation-id = {2021.11.20.469390},\
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	year = {2021},\
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	doi = {10.1101/2021.11.20.469390},\
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	publisher = {Cold Spring Harbor Laboratory},\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	URL = { https://www.biorxiv.org/content/early/2021/11/22/2021.11.20.469390 },  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	eprint = { https://www.biorxiv.org/content/early/2021/11/22/2021.11.20.469390.full.pdf }, \
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	journal = {bioRxiv}\
>}
