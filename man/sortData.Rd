\name{sortData}

\alias{sortData}

\title{
    This function sorts a dataset file by the genomic position of the probes.
}

\description{
    This function sorts a dataset file by the genomic position of the probes. 
    This function makes very easy the integration of VegaMC with the output of 
    PennCNV tool.
}

\usage{
    sortData(dataset, output_file_name = "")
}

\arguments{
    \item{dataset}{Dataset file.}
    \item{output_file_name}{Name of the file in which sorted data are stored.}
}

\value{
    This function returns the input matrix ordered by the genomic
    position of the probes.
}

\references{
    Morganella S., and Ceccarelli M. VegaMC: a R/bioconductor package for fast downstream analysis of large array comparative genomic hybridization datasets. Bioinformatics, 28(19):2512-4 (2012).
}

\author{
    Sandro Morganella
}

\note{
    This function allows to sort a dataset by the genomic position. The input file
    must have the chromosome and the position in column two and three respectively.
    This format follows the standard output of PennCNV. An example of file can be 
    found in inst/example folder.
}
\examples{

    ## Copy the example dataset in current folder
    file.copy(system.file("example/breast_Affy500K.txt", package="VegaMC"), 
                            ".")


    ## Sort data and save results in sorted.txt file
    sortData("breast_Affy500K.txt", "sorted.txt");

}
