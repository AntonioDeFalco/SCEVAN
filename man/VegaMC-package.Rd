\name{VegaMC-package}

\alias{VegaMC-package}

\alias{VegaMC}

\docType{package}

\title{
    VegaMC: A Package Implementing a Variational Piecewise Smooth Model for 
    Identification of Driver Chromosomal Imbalances in Cancer.
}

\description{
    VegaMC enables the detection of driver chromosomal imbalances 
    (deletions, amplifications and loss of heterozygosities (LOHs)) from 
    array comparative genomic hybridization (aCGH) data. VegaMC performs 
    a joint segmentation of aCGH data. Segmented regions are then used 
    into a statistical framework to distinguish between driver and 
    passenger mutations. In this way, significant imbalances can be 
    detected by the associated p-value. VegaMC has been implemented to be 
    easily integrated with the output produced by PennCNV. VegaMC produces 
    in output two web pages allowing a rapid navigation between both 
    detected regions and altered genes. In the web page summarizing the 
    altered genes, the user finds the link to the respective Ensembl gene
    web page.
}

\details{
    \tabular{ll}{
        Package: \tab VegaMC\cr
        Type: \tab Package\cr
        Version: \tab 3.9.3\cr
        License: \tab GPL-2\cr
        LazyLoad: \tab yes\cr
    }
}

\examples{
    ## Copy the example dataset in current folder
    file.copy(system.file("example/breast_Affy500K.txt", package="VegaMC"), 
                            ".")

    ## Analyse data and save results in sorted.txt file
    results <- vegaMC("breast_Affy500K.txt", "results", html=FALSE, 
                        getGenes=FALSE)
}
